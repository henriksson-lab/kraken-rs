use std::collections::HashSet;
use std::io;

use crate::hash::murmurhash3;
use crate::minimizer::MinimizerScanner;
use crate::seq::BatchSequenceReader;
use crate::types::*;

const RANGE_SECTIONS: usize = 1024;
const RANGE_MASK: u64 = (RANGE_SECTIONS - 1) as u64;
const MAX_N: usize = RANGE_SECTIONS;
const DEFAULT_N: usize = 4;
const DEFAULT_BLOCK_SIZE: usize = 30 * 1024 * 1024;

/// Options for capacity estimation.
pub struct EstimateOptions {
    pub k: usize,
    pub l: usize,
    pub n: usize,
    pub input_is_protein: bool,
    pub threads: usize,
    pub block_size: usize,
    pub spaced_seed_mask: u64,
    pub toggle_mask: u64,
}

impl Default for EstimateOptions {
    fn default() -> Self {
        EstimateOptions {
            k: 0,
            l: 0,
            n: DEFAULT_N,
            input_is_protein: false,
            threads: 1,
            block_size: DEFAULT_BLOCK_SIZE,
            spaced_seed_mask: DEFAULT_SPACED_SEED_MASK,
            toggle_mask: DEFAULT_TOGGLE_MASK,
        }
    }
}

pub fn parse_command_line(args: &[String], opts: &mut EstimateOptions) -> io::Result<()> {
    let mut i = 1usize;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "-?" => return Err(io::Error::new(io::ErrorKind::InvalidInput, usage(0))),
            "-p" => {
                i += 1;
                let sig = args
                    .get(i)
                    .and_then(|s| s.parse::<isize>().ok())
                    .unwrap_or(0);
                if sig < 1 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "must have at least 1 thread",
                    ));
                }
                opts.threads = sig as usize;
            }
            "-B" => {
                i += 1;
                let sig = args
                    .get(i)
                    .and_then(|s| s.parse::<isize>().ok())
                    .unwrap_or(0);
                if sig < 1 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "block size must be positive",
                    ));
                }
                opts.block_size = sig as usize;
            }
            "-n" => {
                i += 1;
                let sig = args
                    .get(i)
                    .and_then(|s| s.parse::<isize>().ok())
                    .unwrap_or(0);
                if sig < 1 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "n must be positive integer",
                    ));
                }
                if sig as usize > MAX_N {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("n must be no more than {}", MAX_N),
                    ));
                }
                opts.n = sig as usize;
            }
            "-k" => {
                i += 1;
                let sig = args
                    .get(i)
                    .and_then(|s| s.parse::<isize>().ok())
                    .unwrap_or(0);
                if sig < 1 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "k must be positive integer",
                    ));
                }
                opts.k = sig as usize;
            }
            "-l" => {
                i += 1;
                let sig = args
                    .get(i)
                    .and_then(|s| s.parse::<isize>().ok())
                    .unwrap_or(0);
                if sig < 1 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "l must be positive integer",
                    ));
                }
                if sig > 31 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "l must be no more than 31",
                    ));
                }
                opts.l = sig as usize;
            }
            "-X" => opts.input_is_protein = true,
            "-S" => {
                i += 1;
                opts.spaced_seed_mask = args
                    .get(i)
                    .and_then(|s| u64::from_str_radix(s, 2).ok())
                    .unwrap_or(DEFAULT_SPACED_SEED_MASK);
            }
            "-T" => {
                i += 1;
                opts.toggle_mask = args
                    .get(i)
                    .and_then(|s| u64::from_str_radix(s, 2).ok())
                    .unwrap_or(DEFAULT_TOGGLE_MASK);
            }
            _ => {}
        }
        i += 1;
    }

    if opts.spaced_seed_mask != DEFAULT_SPACED_SEED_MASK {
        crate::utilities::expand_spaced_seed_mask(
            &mut opts.spaced_seed_mask,
            if opts.input_is_protein {
                BITS_PER_CHAR_PRO as i32
            } else {
                BITS_PER_CHAR_DNA as i32
            },
        );
    }
    if opts.k == 0 || opts.l == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "missing mandatory integer parameter",
        ));
    }
    if opts.k < opts.l {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "k cannot be less than l",
        ));
    }
    Ok(())
}

pub fn usage(_exit_code: i32) -> String {
    [
        "Usage: estimate_capacity <options>",
        "",
        "Options (*mandatory):",
        "* -k INT        Set length of k-mers",
        "* -l INT        Set length of minimizers",
        "  -n INT        Set maximum qualifying hash code",
        "  -X            Input sequences are proteins",
        "  -S BITSTRING  Spaced seed mask",
        "  -T BITSTRING  Minimizer ordering toggle mask",
        "  -B INT        Read block size",
        "  -p INT        Number of threads",
    ]
    .join("\n")
}

/// Process a single sequence for capacity estimation.
pub fn process_sequence(seq: &str, opts: &EstimateOptions, sets: &mut [HashSet<u64>]) {
    let mut seq_str = seq.to_string();

    // Add terminator for protein sequences if not already there
    if opts.input_is_protein && !seq_str.ends_with('*') {
        seq_str.push('*');
    }

    let mut scanner = MinimizerScanner::new(
        opts.k as isize,
        opts.l as isize,
        opts.spaced_seed_mask,
        !opts.input_is_protein,
        opts.toggle_mask,
        CURRENT_REVCOM_VERSION,
    );
    scanner.load_sequence(&seq_str, 0, usize::MAX);

    while let Some(minimizer) = scanner.next_minimizer() {
        if scanner.is_ambiguous() {
            continue;
        }
        let hash_code = murmurhash3(minimizer);
        let idx = (hash_code & RANGE_MASK) as usize;
        if idx < opts.n {
            sets[idx].insert(minimizer);
        }
    }
}

pub fn process_sequences(opts: &EstimateOptions) -> io::Result<usize> {
    let mut vector_of_sets: Vec<Vec<HashSet<u64>>> = (0..opts.threads)
        .map(|_| (0..opts.n).map(|_| HashSet::new()).collect())
        .collect();
    let mut reader = BatchSequenceReader::new(None)?;

    let mut thread_idx = 0usize;
    while reader.load_block(opts.block_size) {
        while let Some(seq) = reader.next_sequence() {
            process_sequence(&seq.seq, opts, &mut vector_of_sets[thread_idx]);
            thread_idx = (thread_idx + 1) % opts.threads.max(1);
        }
    }

    let mut sets = vector_of_sets.swap_remove(0);
    for worker_sets in vector_of_sets {
        for j in 0..opts.n {
            sets[j].extend(worker_sets[j].iter().copied());
        }
    }

    let mut sum_set_sizes = 0usize;
    for s in &sets {
        sum_set_sizes += s.len();
    }
    sum_set_sizes += 1;
    Ok((sum_set_sizes * RANGE_SECTIONS) / opts.n)
}

pub fn estimate_capacity_main(args: &[String]) -> io::Result<usize> {
    let mut opts = EstimateOptions::default();
    parse_command_line(args, &mut opts)?;
    process_sequences(&opts)
}

/// Estimate the number of distinct minimizers in the input (from stdin).
/// Returns the estimated capacity.
pub fn estimate_capacity(opts: &EstimateOptions) -> io::Result<usize> {
    let mut reader = BatchSequenceReader::new(None)?;
    let mut sets: Vec<HashSet<u64>> = (0..opts.n).map(|_| HashSet::new()).collect();

    while reader.load_block(opts.block_size) {
        while let Some(seq) = reader.next_sequence() {
            process_sequence(&seq.seq, opts, &mut sets);
        }
    }

    let sum_set_sizes: usize = sets.iter().map(|s| s.len()).sum::<usize>() + 1;
    Ok(sum_set_sizes * RANGE_SECTIONS / opts.n)
}

/// Estimate capacity from a list of files (instead of stdin).
pub fn estimate_capacity_from_files(
    filenames: &[String],
    opts: &EstimateOptions,
) -> io::Result<usize> {
    let mut sets: Vec<HashSet<u64>> = (0..opts.n).map(|_| HashSet::new()).collect();

    for filename in filenames {
        let mut reader = BatchSequenceReader::new(Some(filename))?;
        while reader.load_block(opts.block_size) {
            while let Some(seq) = reader.next_sequence() {
                process_sequence(&seq.seq, opts, &mut sets);
            }
        }
    }

    let sum_set_sizes: usize = sets.iter().map(|s| s.len()).sum::<usize>() + 1;
    Ok(sum_set_sizes * RANGE_SECTIONS / opts.n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_command_line() {
        let args = vec![
            "estimate_capacity".to_string(),
            "-k".to_string(),
            "35".to_string(),
            "-l".to_string(),
            "31".to_string(),
            "-n".to_string(),
            "4".to_string(),
            "-B".to_string(),
            "1024".to_string(),
            "-p".to_string(),
            "2".to_string(),
        ];
        let mut opts = EstimateOptions::default();
        parse_command_line(&args, &mut opts).unwrap();
        assert_eq!(opts.k, 35);
        assert_eq!(opts.l, 31);
        assert_eq!(opts.n, 4);
        assert_eq!(opts.block_size, 1024);
        assert_eq!(opts.threads, 2);
    }

    #[test]
    fn test_usage() {
        let text = usage(0);
        assert!(text.contains("Usage: estimate_capacity <options>"));
        assert!(text.contains("-k INT"));
    }

    #[test]
    fn test_process_sequence() {
        let opts = EstimateOptions {
            k: 5,
            l: 3,
            n: RANGE_SECTIONS,
            ..Default::default()
        };
        let mut sets: Vec<HashSet<u64>> = (0..opts.n).map(|_| HashSet::new()).collect();
        process_sequence("ACGTACGTACGT", &opts, &mut sets);
        assert!(sets.iter().map(|s| s.len()).sum::<usize>() > 0);
    }
}

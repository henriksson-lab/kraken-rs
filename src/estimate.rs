use std::collections::HashSet;
use std::io;

use crate::hash::murmurhash3;
use crate::minimizer::MinimizerScanner;
use crate::seq::BatchSequenceReader;
use crate::types::*;


const RANGE_SECTIONS: usize = 1024;
const RANGE_MASK: u64 = (RANGE_SECTIONS - 1) as u64;
const DEFAULT_N: usize = 4;

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
            block_size: 30 * 1024 * 1024,
            spaced_seed_mask: DEFAULT_SPACED_SEED_MASK,
            toggle_mask: DEFAULT_TOGGLE_MASK,
        }
    }
}

/// Process a single sequence for capacity estimation.
fn process_sequence(seq: &str, opts: &EstimateOptions, sets: &mut [HashSet<u64>]) {
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

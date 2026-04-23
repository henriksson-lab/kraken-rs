use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use rayon::prelude::*;

use crate::seq::BatchSequenceReader;
use crate::types::Sequence;

const DEFAULT_WINDOW_SIZE: i32 = 64;
const DEFAULT_THRESHOLD: i32 = 20;
const DEFAULT_LINE_WIDTH: usize = 72;
const DEFAULT_BATCH_SIZE: usize = 64;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PerfectInterval {
    pub start: i32,
    pub finish: i32,
    pub left: i32,
    pub right: i32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MaskRange {
    pub start: i32,
    pub finish: i32,
}

#[derive(Debug, Clone)]
pub struct SDust {
    pub kmers: VecDeque<i32>,
    pub perfect_intervals: Vec<PerfectInterval>,
    pub ranges: Vec<MaskRange>,
    pub seq: Sequence,
    pub cw: [i32; 64],
    pub cv: [i32; 64],
    pub rw: i32,
    pub rv: i32,
    pub l: i32,
    pub window_size: i32,
    pub threshold: i32,
    pub replace_masked_with: Option<u8>,
}

impl Default for SDust {
    fn default() -> Self {
        Self {
            kmers: VecDeque::new(),
            perfect_intervals: Vec::new(),
            ranges: Vec::new(),
            seq: Sequence::default(),
            cw: [0; 64],
            cv: [0; 64],
            rw: 0,
            rv: 0,
            l: 0,
            window_size: DEFAULT_WINDOW_SIZE,
            threshold: DEFAULT_THRESHOLD,
            replace_masked_with: None,
        }
    }
}

impl SDust {
    pub fn default_new() -> Self {
        Self::default()
    }

    pub fn new(window_size: i32, threshold: i32, replace_masked_with: Option<u8>) -> Self {
        Self {
            window_size,
            threshold,
            replace_masked_with,
            ..Default::default()
        }
    }

    pub fn reset(&mut self) {
        self.kmers.clear();
        self.perfect_intervals.clear();
        self.ranges.clear();
        self.cw.fill(0);
        self.cv.fill(0);
        self.l = 0;
        self.rw = 0;
        self.rv = 0;
    }

    fn process_masked_nucleotide(&self, c: u8) -> u8 {
        self.replace_masked_with
            .unwrap_or_else(|| c.to_ascii_lowercase())
    }
}

fn asc2dna(c: u8) -> u8 {
    match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' | b'U' | b'u' => 3,
        _ => 4,
    }
}

pub fn usage(prog: &str) -> String {
    format!(
        "usage: {prog} [-T | -level threshold] [-W | -window window size] [-i | -in input file]\n [-w | -width sequence width] [-o | -out output file] [-r | -replace-masked-with char]\n [-f | -outfmt output format] [-t | -threads threads]\n"
    )
}

pub fn file_exists(filename: &str) -> bool {
    Path::new(filename).exists()
}

pub fn stricasecmp(s1: &str, s2: &str) -> bool {
    s1.eq_ignore_ascii_case(s2)
}

pub fn shift_window(sd: &mut SDust, t: i32) {
    let mut s;
    if sd.kmers.len() as i32 >= sd.window_size - 2 {
        s = sd.kmers.pop_front().unwrap();
        sd.cw[s as usize] -= 1;
        sd.rw -= sd.cw[s as usize];
        if sd.l > sd.kmers.len() as i32 {
            sd.l -= 1;
            sd.cv[s as usize] -= 1;
            sd.rv -= sd.cv[s as usize];
        }
    }
    sd.kmers.push_back(t);
    sd.l += 1;
    sd.rw += sd.cw[t as usize];
    sd.cw[t as usize] += 1;
    sd.rv += sd.cv[t as usize];
    sd.cv[t as usize] += 1;
    if sd.cv[t as usize] * 10 > sd.threshold * 2 {
        loop {
            let idx = sd.kmers.len() - sd.l as usize;
            s = sd.kmers[idx];
            sd.cv[s as usize] -= 1;
            sd.rv -= sd.cv[s as usize];
            sd.l -= 1;
            if s == t {
                break;
            }
        }
    }
}

pub fn save_masked_regions(sd: &mut SDust, window_start: i32) {
    let mut saved = false;
    if sd.perfect_intervals.is_empty() || sd.perfect_intervals.last().unwrap().start >= window_start
    {
        return;
    }
    let p = *sd.perfect_intervals.last().unwrap();
    if let Some(last_range) = sd.ranges.last_mut() {
        let start = last_range.start;
        let finish = last_range.finish;
        if p.start <= finish {
            *last_range = MaskRange {
                start,
                finish: p.finish.max(finish),
            };
            saved = true;
        }
    }
    if !saved {
        sd.ranges.push(MaskRange {
            start: p.start,
            finish: p.finish,
        });
    }
    while let Some(last) = sd.perfect_intervals.last() {
        if last.start < window_start {
            sd.perfect_intervals.pop();
        } else {
            break;
        }
    }
}

pub fn find_perfect(sd: &mut SDust, window_start: i32) {
    let mut cv = sd.cv;
    let mut max_left = 0;
    let mut max_right = 0;
    let mut new_left;
    let mut new_right = sd.rv;

    for i in (0..=(sd.kmers.len() as i32 - sd.l - 1)).rev() {
        let kmer = sd.kmers[i as usize];
        new_right += cv[kmer as usize];
        cv[kmer as usize] += 1;
        new_left = sd.kmers.len() as i32 - i - 1;
        if new_right * 10 > sd.threshold * new_left {
            let mut j = 0usize;
            while j < sd.perfect_intervals.len()
                && sd.perfect_intervals[j].start >= i + window_start
            {
                let p = sd.perfect_intervals[j];
                if max_right == 0 || p.right * max_left > max_right * p.left {
                    max_left = p.left;
                    max_right = p.right;
                }
                j += 1;
            }
            if max_right == 0 || new_right * max_left >= max_right * new_left {
                max_left = new_left;
                max_right = new_right;
                let p = PerfectInterval {
                    start: i + window_start,
                    finish: sd.kmers.len() as i32 + 2 + window_start,
                    left: new_left,
                    right: new_right,
                };
                sd.perfect_intervals.insert(j, p);
            }
        }
    }
}

pub fn run_symmetric_dust(sd: &mut SDust, seq: &mut [u8], size: usize, _offset: i32) {
    let mut triplet = 0i32;
    let mut window_start = 0i32;
    let mut l = 0i32;
    for &base_char in seq.iter().take(size) {
        let base = asc2dna(base_char);
        if base < 4 {
            l += 1;
            triplet = ((triplet << 2) | i32::from(base)) & 63;
            if l >= 3 {
                window_start = (l - sd.window_size).max(0);
                save_masked_regions(sd, window_start);
                shift_window(sd, triplet);
                if sd.rw * 10 > sd.l * sd.threshold {
                    find_perfect(sd, window_start);
                }
            }
        }
    }
    while !sd.perfect_intervals.is_empty() {
        save_masked_regions(sd, window_start);
        window_start += 1;
    }
}

pub fn print_fasta(seq: &Sequence, out: &mut dyn Write, width: usize) -> io::Result<()> {
    write!(out, ">{}", seq.header)?;
    if !seq.comment.is_empty() {
        write!(out, " {}", seq.comment)?;
    }
    writeln!(out)?;
    let bytes = seq.seq.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() {
        let end = (i + width).min(bytes.len());
        out.write_all(&bytes[i..end])?;
        writeln!(out)?;
        i = end;
    }
    out.flush()
}

pub fn mask(sd: &mut SDust) -> &mut SDust {
    let mut seq_bytes = std::mem::take(&mut sd.seq.seq).into_bytes();
    let mut i = 0usize;
    while i < seq_bytes.len() {
        if asc2dna(seq_bytes[i]) != 4 {
            let start = i;
            loop {
                seq_bytes[i] = seq_bytes[i].to_ascii_uppercase();
                if i + 1 == seq_bytes.len() || asc2dna(seq_bytes[i + 1]) == 4 {
                    break;
                }
                i += 1;
            }
            run_symmetric_dust(sd, &mut seq_bytes[start..=i], i - start + 1, start as i32);
            for range in &sd.ranges {
                for pos in range.start..range.finish {
                    seq_bytes[start + pos as usize] =
                        sd.process_masked_nucleotide(seq_bytes[start + pos as usize]);
                }
            }
            sd.reset();
        }
        i += 1;
    }
    sd.seq.seq = String::from_utf8(seq_bytes).expect("masked sequence should stay ASCII");
    sd
}

fn mask_owned(
    mut seq: Sequence,
    window_size: i32,
    threshold: i32,
    replace_masked_with: Option<u8>,
) -> Sequence {
    let mut sd = SDust::new(window_size, threshold, replace_masked_with);
    sd.seq = std::mem::take(&mut seq);
    mask(&mut sd);
    sd.seq
}

pub fn k2mask_main(args: &[String]) -> io::Result<()> {
    let prog = args.first().map(String::as_str).unwrap_or("k2mask");
    let mut line_width = DEFAULT_LINE_WIDTH;
    let mut threads = 1usize;
    let mut infile = "/dev/stdin".to_string();
    let mut outfile = "/dev/stdout".to_string();
    let mut window_size = DEFAULT_WINDOW_SIZE;
    let mut threshold = DEFAULT_THRESHOLD;
    let mut replace_masked_with: Option<u8> = None;

    let mut i = 1usize;
    while i < args.len() {
        match args[i].as_str() {
            "-W" | "--window" => {
                i += 1;
                window_size = args
                    .get(i)
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(window_size);
            }
            "-T" | "--level" => {
                i += 1;
                threshold = args
                    .get(i)
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(threshold);
            }
            "-i" | "--in" => {
                i += 1;
                infile = args.get(i).cloned().unwrap_or_else(|| infile.clone());
            }
            "-o" | "--out" => {
                i += 1;
                outfile = args.get(i).cloned().unwrap_or_else(|| outfile.clone());
            }
            "-w" | "--width" => {
                i += 1;
                line_width = args
                    .get(i)
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(DEFAULT_LINE_WIDTH);
            }
            "-f" | "--outfmt" => {
                i += 1;
                let arg = args.get(i).map(String::as_str).unwrap_or("fasta");
                if !stricasecmp(arg, "fasta") {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("{prog}:  currently only supports outputting FASTA."),
                    ));
                }
            }
            "-t" | "--threads" => {
                i += 1;
                threads = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(1);
            }
            "-r" | "--replace-masked-with" => {
                i += 1;
                let arg = args.get(i).map(String::as_str).unwrap_or_default();
                if arg.len() != 1 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("{prog}: -r expects a single character, {arg} given."),
                    ));
                }
                replace_masked_with = Some(arg.as_bytes()[0]);
            }
            "-h" | "--help" => {
                return Err(io::Error::new(io::ErrorKind::InvalidInput, usage(prog)));
            }
            _ => {
                return Err(io::Error::new(io::ErrorKind::InvalidInput, usage(prog)));
            }
        }
        i += 1;
    }

    let mut out: Box<dyn Write> = if outfile == "/dev/stdout" {
        Box::new(io::stdout())
    } else {
        Box::new(File::create(&outfile)?)
    };
    let mut reader = if infile == "/dev/stdin" {
        BatchSequenceReader::new(None)?
    } else {
        BatchSequenceReader::new(Some(&infile))?
    };

    let worker_threads = threads.saturating_sub(1).max(1);
    while reader.load_batch(DEFAULT_BATCH_SIZE.max(worker_threads)) {
        let mut batch = Vec::new();
        while let Some(seq) = reader.next_sequence() {
            batch.push(seq.clone());
        }
        let processed: Vec<Sequence> = if threads > 1 && batch.len() > 1 {
            rayon::ThreadPoolBuilder::new()
                .num_threads(worker_threads)
                .build()
                .map_err(|e| io::Error::other(e.to_string()))?
                .install(|| {
                    batch
                        .into_par_iter()
                        .map(|seq| mask_owned(seq, window_size, threshold, replace_masked_with))
                        .collect()
                })
        } else {
            batch
                .into_iter()
                .map(|seq| mask_owned(seq, window_size, threshold, replace_masked_with))
                .collect()
        };
        for seq in &processed {
            print_fasta(seq, &mut out, line_width)?;
        }
    }
    out.flush()
}

/// Compatibility wrapper retained for existing callers.
pub fn mask_low_complexity(seq: &mut [u8], threshold: usize) {
    let original = std::str::from_utf8(seq).expect("mask_low_complexity expects ASCII DNA input");
    let mut sd = SDust::new(DEFAULT_WINDOW_SIZE, threshold as i32, None);
    sd.seq.seq = original.to_string();
    mask(&mut sd);
    seq.copy_from_slice(sd.seq.seq.as_bytes());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_file_exists_and_stricasecmp() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        assert!(file_exists(tmp.path().to_str().unwrap()));
        assert!(stricasecmp("FaStA", "fasta"));
        assert!(!stricasecmp("fasta", "fastq"));
    }

    #[test]
    fn test_shift_window_updates_counts() {
        let mut sd = SDust::default();
        shift_window(&mut sd, 3);
        shift_window(&mut sd, 3);
        assert_eq!(sd.kmers.len(), 2);
        assert_eq!(sd.cw[3], 2);
        assert!(sd.rw >= 1);
    }

    #[test]
    fn test_mask_low_complexity_masks_repeats() {
        let mut seq = b"AAAAAAAAAAAAAAAAAAAA".to_vec();
        mask_low_complexity(&mut seq, DEFAULT_THRESHOLD as usize);
        assert!(seq.iter().any(u8::is_ascii_lowercase));
    }

    #[test]
    fn test_print_fasta_wraps_lines() {
        let seq = Sequence {
            header: "seq1".to_string(),
            comment: "comment".to_string(),
            seq: "ACGTACGT".to_string(),
            ..Default::default()
        };
        let mut out = Vec::new();
        print_fasta(&seq, &mut out, 4).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            ">seq1 comment\nACGT\nACGT\n"
        );
    }

    #[test]
    fn test_k2mask_main_writes_masked_fasta() {
        let input = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(input.path(), b">r1\nAAAAAAAAAAAAAAA\n>r2\nACGTNACGT\n").unwrap();
        let output = tempfile::NamedTempFile::new().unwrap();
        let args = vec![
            "k2mask".to_string(),
            "-i".to_string(),
            input.path().to_str().unwrap().to_string(),
            "-o".to_string(),
            output.path().to_str().unwrap().to_string(),
            "-w".to_string(),
            "6".to_string(),
            "-t".to_string(),
            "2".to_string(),
        ];
        k2mask_main(&args).unwrap();
        let text = std::fs::read_to_string(output.path()).unwrap();
        assert!(text.contains(">r1"));
        assert!(text.contains(">r2"));
        assert!(text.contains("aaaaaa") || text.contains("AAAAAA"));
    }

    #[test]
    fn test_mask_preserves_non_dna_gaps() {
        let mut sd = SDust::default();
        sd.seq.seq = "AAAAANNNNAAAAA".to_string();
        mask(&mut sd);
        assert_eq!(&sd.seq.seq[5..9], "NNNN");
    }

    #[test]
    fn test_usage_contains_expected_flags() {
        let text = usage("k2mask");
        assert!(text.contains("usage: k2mask"));
        assert!(text.contains("-replace-masked-with"));
    }
}

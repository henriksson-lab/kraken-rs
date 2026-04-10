use std::collections::VecDeque;
use crate::types::*;

struct MinimizerData {
    candidate: u64,
    pos: isize,
}

/// Minimizer scanner — extracts minimizers from DNA or protein sequences.
/// Exact port of C++ `MinimizerScanner` from `mmscanner.h/cc`.
pub struct MinimizerScanner {
    seq: Vec<u8>,
    k: isize,
    l: isize,
    str_pos: usize,
    start: usize,
    finish: usize,
    spaced_seed_mask: u64,
    dna: bool,
    toggle_mask: u64,
    lmer: u64,
    lmer_mask: u64,
    last_minimizer: u64,
    loaded_ch: isize,
    queue: VecDeque<MinimizerData>,
    queue_pos: isize,
    last_ambig: u64,
    lookup_table: [u8; 256],
    revcom_version: i32,
}

impl MinimizerScanner {
    pub fn new(
        k: isize,
        l: isize,
        spaced_seed_mask: u64,
        dna: bool,
        toggle_mask: u64,
        revcom_version: i32,
    ) -> Self {
        let bits_per_char = if dna { BITS_PER_CHAR_DNA } else { BITS_PER_CHAR_PRO };

        let mut lmer_mask: u64 = 1;
        lmer_mask <<= l as u64 * bits_per_char as u64;
        lmer_mask -= 1;

        let toggle_mask = toggle_mask & lmer_mask;

        let mut lookup_table = [0xFFu8; 256];

        if dna {
            Self::set_lookup(&mut lookup_table, b'A', 0x00);
            Self::set_lookup(&mut lookup_table, b'C', 0x01);
            Self::set_lookup(&mut lookup_table, b'G', 0x02);
            Self::set_lookup(&mut lookup_table, b'T', 0x03);
        } else {
            // Reduced 15-letter protein alphabet
            Self::set_lookup(&mut lookup_table, b'*', 0x00);
            Self::set_lookup(&mut lookup_table, b'U', 0x00);
            Self::set_lookup(&mut lookup_table, b'O', 0x00);
            Self::set_lookup(&mut lookup_table, b'A', 0x01);
            Self::set_lookup(&mut lookup_table, b'N', 0x02);
            Self::set_lookup(&mut lookup_table, b'Q', 0x02);
            Self::set_lookup(&mut lookup_table, b'S', 0x02);
            Self::set_lookup(&mut lookup_table, b'C', 0x03);
            Self::set_lookup(&mut lookup_table, b'D', 0x04);
            Self::set_lookup(&mut lookup_table, b'E', 0x04);
            Self::set_lookup(&mut lookup_table, b'F', 0x05);
            Self::set_lookup(&mut lookup_table, b'G', 0x06);
            Self::set_lookup(&mut lookup_table, b'H', 0x07);
            Self::set_lookup(&mut lookup_table, b'I', 0x08);
            Self::set_lookup(&mut lookup_table, b'L', 0x08);
            Self::set_lookup(&mut lookup_table, b'K', 0x09);
            Self::set_lookup(&mut lookup_table, b'P', 0x0a);
            Self::set_lookup(&mut lookup_table, b'R', 0x0b);
            Self::set_lookup(&mut lookup_table, b'M', 0x0c);
            Self::set_lookup(&mut lookup_table, b'V', 0x0c);
            Self::set_lookup(&mut lookup_table, b'T', 0x0d);
            Self::set_lookup(&mut lookup_table, b'W', 0x0e);
            Self::set_lookup(&mut lookup_table, b'Y', 0x0f);
        }

        MinimizerScanner {
            seq: Vec::new(),
            k,
            l,
            str_pos: 0,
            start: 0,
            finish: 0,
            spaced_seed_mask,
            dna,
            toggle_mask,
            lmer: 0,
            lmer_mask,
            last_minimizer: !0u64,
            loaded_ch: 0,
            queue: VecDeque::new(),
            queue_pos: 0,
            last_ambig: 0,
            lookup_table,
            revcom_version,
        }
    }

    fn set_lookup(table: &mut [u8; 256], ch: u8, val: u8) {
        table[ch as usize] = val;
        table[ch.to_ascii_lowercase() as usize] = val;
    }

    pub fn load_sequence(&mut self, seq: &str, start: usize, finish: usize) {
        let bytes = seq.as_bytes();
        self.seq.clear();
        self.seq.extend_from_slice(bytes);
        self.start = start;
        let mut finish = finish;
        if finish > self.seq.len() {
            finish = self.seq.len();
        }
        self.finish = finish;
        self.str_pos = self.start;

        if (self.finish as isize - self.start as isize) + 1 < self.l {
            self.str_pos = self.finish;
        }

        self.queue.clear();
        self.queue_pos = 0;
        self.loaded_ch = 0;
        self.last_minimizer = !0u64;
        self.last_ambig = 0;
    }

    /// Returns the next minimizer, or None if exhausted.
    pub fn next_minimizer(&mut self) -> Option<u64> {
        if self.str_pos >= self.finish {
            return None;
        }

        let bits_per_char = if self.dna { BITS_PER_CHAR_DNA } else { BITS_PER_CHAR_PRO };
        let ambig_code: u64 = (1u64 << bits_per_char) - 1;
        let mut changed_minimizer = false;

        while !changed_minimizer {
            // Incorporate next character(s)
            if self.loaded_ch == self.l {
                self.loaded_ch -= 1;
            }
            while self.loaded_ch < self.l && self.str_pos < self.finish {
                self.loaded_ch += 1;
                self.lmer <<= bits_per_char;
                self.last_ambig <<= bits_per_char;

                let lookup_code = self.lookup_table[self.seq[self.str_pos] as usize];
                self.str_pos += 1;

                if lookup_code == 0xFF {
                    self.queue.clear();
                    self.queue_pos = 0;
                    self.lmer = 0;
                    self.loaded_ch = 0;
                    self.last_ambig |= ambig_code;
                } else {
                    self.lmer |= lookup_code as u64;
                }
                self.lmer &= self.lmer_mask;
                self.last_ambig &= self.lmer_mask;

                if (self.str_pos - self.start) >= self.k as usize && self.loaded_ch < self.l {
                    return Some(self.last_minimizer);
                }
            }

            if self.loaded_ch < self.l {
                return None;
            }

            let canonical_lmer = if self.dna {
                self.canonical_representation(self.lmer, self.l as u8)
            } else {
                self.lmer
            };

            let canonical_lmer = if self.spaced_seed_mask != 0 {
                canonical_lmer & self.spaced_seed_mask
            } else {
                canonical_lmer
            };

            let candidate_lmer = canonical_lmer ^ self.toggle_mask;

            if self.k == self.l {
                self.last_minimizer = candidate_lmer ^ self.toggle_mask;
                return Some(self.last_minimizer);
            }

            // Sliding window minimum
            while !self.queue.is_empty() && self.queue.back().unwrap().candidate > candidate_lmer {
                self.queue.pop_back();
            }

            let data = MinimizerData {
                candidate: candidate_lmer,
                pos: self.queue_pos,
            };

            if self.queue.is_empty() && self.queue_pos >= self.k - self.l {
                changed_minimizer = true;
            }

            self.queue.push_back(data);

            if self.queue.front().unwrap().pos < self.queue_pos - self.k + self.l {
                self.queue.pop_front();
                changed_minimizer = true;
            }

            if self.queue_pos == self.k - self.l {
                changed_minimizer = true;
            }
            self.queue_pos += 1;

            if self.str_pos >= self.k as usize {
                break;
            }
        }

        assert!(!self.queue.is_empty());
        self.last_minimizer = self.queue.front().unwrap().candidate ^ self.toggle_mask;
        Some(self.last_minimizer)
    }

    /// Check if the current minimizer position is ambiguous.
    pub fn is_ambiguous(&self) -> bool {
        self.queue_pos < self.k - self.l || self.last_ambig != 0
    }

    /// Reverse complement of a k-mer (DNA only).
    /// Exact port of C++ `MinimizerScanner::reverse_complement()`.
    fn reverse_complement(&self, mut kmer: u64, n: u8) -> u64 {
        // Reverse bit pairs
        kmer = ((kmer & 0xCCCCCCCCCCCCCCCC) >> 2)
             | ((kmer & 0x3333333333333333) << 2);
        kmer = ((kmer & 0xF0F0F0F0F0F0F0F0) >> 4)
             | ((kmer & 0x0F0F0F0F0F0F0F0F) << 4);
        kmer = ((kmer & 0xFF00FF00FF00FF00) >> 8)
             | ((kmer & 0x00FF00FF00FF00FF) << 8);
        kmer = ((kmer & 0xFFFF0000FFFF0000) >> 16)
             | ((kmer & 0x0000FFFF0000FFFF) << 16);
        kmer = kmer.rotate_left(32);

        // Complement
        if self.revcom_version == 0 {
            (!kmer) & ((1u64 << (n as u64 * 2)) - 1)
        } else {
            ((!kmer) >> (64 - n as u64 * 2)) & ((1u64 << (n as u64 * 2)) - 1)
        }
    }

    /// Canonical representation = min(kmer, reverse_complement(kmer)).
    fn canonical_representation(&self, kmer: u64, n: u8) -> u64 {
        let revcom = self.reverse_complement(kmer, n);
        if kmer < revcom { kmer } else { revcom }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn collect_minimizers(seq: &str, k: isize, l: isize) -> (Vec<u64>, Vec<bool>) {
        let mut scanner = MinimizerScanner::new(
            k, l, DEFAULT_SPACED_SEED_MASK, true,
            DEFAULT_TOGGLE_MASK, CURRENT_REVCOM_VERSION,
        );
        scanner.load_sequence(seq, 0, usize::MAX);
        let mut mins = Vec::new();
        let mut ambs = Vec::new();
        while let Some(m) = scanner.next_minimizer() {
            mins.push(m);
            ambs.push(scanner.is_ambiguous());
        }
        (mins, ambs)
    }

    #[test]
    fn test_minimizer_produces_output() {
        let seq = "ACGATCGACGACGATCGATCGATCGATCGATCGATCGATCG";
        let (mins, _) = collect_minimizers(seq, 10, 5);
        assert!(!mins.is_empty(), "Should produce minimizers for a 41-char sequence with k=10 l=5");
    }

    #[test]
    fn test_minimizer_consistency() {
        // Same input should always give same output
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let (m1, a1) = collect_minimizers(seq, 15, 8);
        let (m2, a2) = collect_minimizers(seq, 15, 8);
        assert_eq!(m1, m2);
        assert_eq!(a1, a2);
    }

    #[test]
    fn test_minimizer_ambiguous_bases() {
        // Sequence with N's should report ambiguous positions
        let seq = "NNNNNACGATCGACGATCGATCGATCGATCGATCGATCG";
        let (mins, ambs) = collect_minimizers(seq, 10, 5);
        assert!(!mins.is_empty());
        // First minimizers should be ambiguous (N's at start)
        assert!(ambs[0], "First minimizer should be ambiguous due to leading N's");
        // Later minimizers should not be ambiguous
        let first_clear = ambs.iter().position(|&a| !a);
        assert!(first_clear.is_some(), "Should have some non-ambiguous minimizers");
    }

    #[test]
    fn test_minimizer_all_same_bases() {
        // Homopolymer: all A's. Should produce minimizers.
        let seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let (mins, _) = collect_minimizers(seq, 10, 5);
        assert!(!mins.is_empty());
        // All minimizers from a homopolymer should be the same
        let first = mins[0];
        for &m in &mins[1..] {
            assert_eq!(m, first, "Homopolymer should yield identical minimizers");
        }
    }

    #[test]
    fn test_minimizer_k_equals_l() {
        // When k==l, should short-circuit to direct l-mer output
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let (mins, _) = collect_minimizers(seq, 10, 10);
        assert!(!mins.is_empty());
        // Should have (len - k + 1) minimizers
        assert_eq!(mins.len(), seq.len() - 10 + 1);
    }

    #[test]
    fn test_minimizer_count() {
        // For a clean sequence of length L with k-mer length k,
        // we should get exactly L - k + 1 minimizers (one per k-mer position)
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 40 bp
        let k = 10isize;
        let l = 5isize;
        let (mins, _) = collect_minimizers(seq, k, l);
        assert_eq!(mins.len(), seq.len() - k as usize + 1,
            "Expected {} minimizers, got {}", seq.len() - k as usize + 1, mins.len());
    }

    #[test]
    fn test_minimizer_different_params() {
        let seq = "ACGATCGACGACGATCGATCGATCGATCGATCGATCGATCG";
        // Different k/l should generally give different minimizer sequences
        let (m1, _) = collect_minimizers(seq, 10, 5);
        let (m2, _) = collect_minimizers(seq, 15, 8);
        // They should have different lengths at minimum
        assert_ne!(m1.len(), m2.len());
    }
}

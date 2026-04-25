use crate::types::*;

#[derive(Clone, Copy)]
struct MinimizerData {
    candidate: u64,
    pos: isize,
}

/// Sliding-window deque for minimizer candidates.
///
/// The window can hold at most `k - l + 1` entries; for any realistic
/// `(k, l)` pair this is well below 64. We use a fixed-size inline buffer
/// with a power-of-two capacity so the wrap-around indexing compiles to a
/// bitwise AND. Replaces `std::collections::VecDeque` from earlier revisions
/// to remove a level of indirection and a bounds-check from the inner loop.
const Q_CAP: usize = 64;
const Q_MASK: usize = Q_CAP - 1;

struct MinQueue {
    data: [MinimizerData; Q_CAP],
    head: usize,
    len: usize,
}

impl MinQueue {
    #[inline]
    fn new() -> Self {
        Self {
            data: [MinimizerData {
                candidate: 0,
                pos: 0,
            }; Q_CAP],
            head: 0,
            len: 0,
        }
    }

    #[inline]
    fn clear(&mut self) {
        self.head = 0;
        self.len = 0;
    }

    #[inline]
    fn is_empty(&self) -> bool {
        self.len == 0
    }

    #[inline]
    fn front(&self) -> &MinimizerData {
        debug_assert!(self.len > 0);
        unsafe { self.data.get_unchecked(self.head & Q_MASK) }
    }

    #[inline]
    fn back(&self) -> &MinimizerData {
        debug_assert!(self.len > 0);
        unsafe { self.data.get_unchecked((self.head + self.len - 1) & Q_MASK) }
    }

    #[inline]
    fn push_back(&mut self, x: MinimizerData) {
        debug_assert!(self.len < Q_CAP);
        let idx = (self.head + self.len) & Q_MASK;
        unsafe {
            *self.data.get_unchecked_mut(idx) = x;
        }
        self.len += 1;
    }

    #[inline]
    fn pop_back(&mut self) {
        debug_assert!(self.len > 0);
        self.len -= 1;
    }

    #[inline]
    fn pop_front(&mut self) {
        debug_assert!(self.len > 0);
        self.head = (self.head + 1) & Q_MASK;
        self.len -= 1;
    }
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
    queue: MinQueue,
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
        let bits_per_char = if dna {
            BITS_PER_CHAR_DNA
        } else {
            BITS_PER_CHAR_PRO
        };

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
            queue: MinQueue::new(),
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
    #[inline]
    pub fn next_minimizer(&mut self) -> Option<u64> {
        if self.dna {
            self.next_minimizer_inner::<{ BITS_PER_CHAR_DNA as u32 }, true>()
        } else {
            self.next_minimizer_inner::<{ BITS_PER_CHAR_PRO as u32 }, false>()
        }
    }

    /// Inner body of `next_minimizer` parameterized on the alphabet so LLVM
    /// can constant-fold `bits_per_char` (shift amounts) and the DNA branch.
    /// Also inlined into both monomorphizations.
    #[inline]
    fn next_minimizer_inner<const BITS: u32, const DNA: bool>(&mut self) -> Option<u64> {
        if self.str_pos >= self.finish {
            return None;
        }

        const fn ambig(bits: u32) -> u64 {
            (1u64 << bits) - 1
        }
        let ambig_code: u64 = ambig(BITS);
        let mut changed_minimizer = false;

        while !changed_minimizer {
            // Incorporate next character(s)
            if self.loaded_ch == self.l {
                self.loaded_ch -= 1;
            }
            while self.loaded_ch < self.l && self.str_pos < self.finish {
                self.loaded_ch += 1;
                self.lmer <<= BITS;
                self.last_ambig <<= BITS;

                // Safety: str_pos < finish <= seq.len(), and lookup_table has 256 entries covering all u8
                let lookup_code = unsafe {
                    *self
                        .lookup_table
                        .get_unchecked(*self.seq.get_unchecked(self.str_pos) as usize)
                };
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

            let canonical_lmer = if DNA {
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
            while !self.queue.is_empty() && self.queue.back().candidate > candidate_lmer {
                self.queue.pop_back();
            }

            if self.queue.is_empty() && self.queue_pos >= self.k - self.l {
                changed_minimizer = true;
            }

            self.queue.push_back(MinimizerData {
                candidate: candidate_lmer,
                pos: self.queue_pos,
            });

            if self.queue.front().pos < self.queue_pos - self.k + self.l {
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

        debug_assert!(!self.queue.is_empty());
        self.last_minimizer = self.queue.front().candidate ^ self.toggle_mask;
        Some(self.last_minimizer)
    }

    /// Check if the current minimizer position is ambiguous.
    #[inline]
    pub fn is_ambiguous(&self) -> bool {
        self.queue_pos < self.k - self.l || self.last_ambig != 0
    }

    /// Reverse complement of a k-mer (DNA only).
    /// Exact port of C++ `MinimizerScanner::reverse_complement()`.
    #[inline]
    fn reverse_complement(&self, mut kmer: u64, n: u8) -> u64 {
        // Reverse bit pairs
        kmer = ((kmer & 0xCCCCCCCCCCCCCCCC) >> 2) | ((kmer & 0x3333333333333333) << 2);
        kmer = ((kmer & 0xF0F0F0F0F0F0F0F0) >> 4) | ((kmer & 0x0F0F0F0F0F0F0F0F) << 4);
        kmer = ((kmer & 0xFF00FF00FF00FF00) >> 8) | ((kmer & 0x00FF00FF00FF00FF) << 8);
        kmer = ((kmer & 0xFFFF0000FFFF0000) >> 16) | ((kmer & 0x0000FFFF0000FFFF) << 16);
        kmer = kmer.rotate_left(32);

        // Complement
        if self.revcom_version == 0 {
            (!kmer) & ((1u64 << (n as u64 * 2)) - 1)
        } else {
            ((!kmer) >> (64 - n as u64 * 2)) & ((1u64 << (n as u64 * 2)) - 1)
        }
    }

    /// Canonical representation = min(kmer, reverse_complement(kmer)).
    #[inline]
    fn canonical_representation(&self, kmer: u64, n: u8) -> u64 {
        let revcom = self.reverse_complement(kmer, n);
        if kmer < revcom {
            kmer
        } else {
            revcom
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn collect_minimizers(seq: &str, k: isize, l: isize) -> (Vec<u64>, Vec<bool>) {
        let mut scanner = MinimizerScanner::new(
            k,
            l,
            DEFAULT_SPACED_SEED_MASK,
            true,
            DEFAULT_TOGGLE_MASK,
            CURRENT_REVCOM_VERSION,
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
        assert!(
            !mins.is_empty(),
            "Should produce minimizers for a 41-char sequence with k=10 l=5"
        );
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
        assert!(
            ambs[0],
            "First minimizer should be ambiguous due to leading N's"
        );
        // Later minimizers should not be ambiguous
        let first_clear = ambs.iter().position(|&a| !a);
        assert!(
            first_clear.is_some(),
            "Should have some non-ambiguous minimizers"
        );
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
        assert_eq!(
            mins.len(),
            seq.len() - k as usize + 1,
            "Expected {} minimizers, got {}",
            seq.len() - k as usize + 1,
            mins.len()
        );
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

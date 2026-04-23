use crate::hyperloglog::HyperLogLogPlusMinus;
use std::collections::HashSet;

/// ReadCounts with HyperLogLog for distinct k-mer estimation (default mode).
#[derive(Clone, Default)]
pub struct ReadCounts {
    pub n_reads: u64,
    pub n_kmers: u64,
    pub kmers: HyperLogLogPlusMinus,
}

impl ReadCounts {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn read_count(&self) -> u64 {
        self.n_reads
    }

    pub fn increment_read_count(&mut self) {
        self.n_reads += 1;
    }

    pub fn kmer_count(&self) -> u64 {
        self.n_kmers
    }

    pub fn distinct_kmer_count(&self) -> u64 {
        self.kmers.size()
    }

    pub fn add_kmer(&mut self, kmer: u64) {
        self.n_kmers += 1;
        self.kmers.insert(kmer);
    }

    pub fn merge(&mut self, other: &ReadCounts) {
        self.n_reads += other.n_reads;
        self.n_kmers += other.n_kmers;
        self.kmers.merge(&other.kmers);
    }
}

impl std::ops::AddAssign for ReadCounts {
    fn add_assign(&mut self, other: Self) {
        self.merge(&other);
    }
}

/// ReadCounts with exact counting via HashSet (when EXACT_COUNTING is needed).
#[derive(Clone, Default)]
pub struct ReadCountsExact {
    pub n_reads: u64,
    pub n_kmers: u64,
    pub kmers: HashSet<u64>,
}

impl ReadCountsExact {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn read_count(&self) -> u64 {
        self.n_reads
    }

    pub fn increment_read_count(&mut self) {
        self.n_reads += 1;
    }

    pub fn kmer_count(&self) -> u64 {
        self.n_kmers
    }

    pub fn distinct_kmer_count(&self) -> u64 {
        self.kmers.len() as u64
    }

    pub fn add_kmer(&mut self, kmer: u64) {
        self.n_kmers += 1;
        self.kmers.insert(kmer);
    }

    pub fn merge(&mut self, other: &ReadCountsExact) {
        self.n_reads += other.n_reads;
        self.n_kmers += other.n_kmers;
        self.kmers.extend(&other.kmers);
    }
}

/// Type alias matching C++ READCOUNTER (default = HLL-based).
pub type ReadCounter = ReadCounts;

/// Map from taxid to read counters.
pub type TaxonCounters = std::collections::BTreeMap<u64, ReadCounter>;

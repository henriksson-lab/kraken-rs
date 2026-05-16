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
    /// Construct an empty `ReadCounts` (zero reads, zero k-mers, empty HLL).
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of reads counted for this taxon.
    pub fn read_count(&self) -> u64 {
        self.n_reads
    }

    /// Increment the read counter by one. Mirrors C++ `incrementReadCount()`.
    pub fn increment_read_count(&mut self) {
        self.n_reads += 1;
    }

    /// Total number of k-mers observed (with multiplicity).
    pub fn kmer_count(&self) -> u64 {
        self.n_kmers
    }

    /// Estimated number of distinct k-mers (HyperLogLog cardinality).
    pub fn distinct_kmer_count(&self) -> u64 {
        self.kmers.size()
    }

    /// Record a k-mer: bumps the total count and inserts it into the HLL
    /// sketch used for the distinct-k-mer estimate.
    pub fn add_kmer(&mut self, kmer: u64) {
        self.n_kmers += 1;
        self.kmers.insert(kmer);
    }

    /// Combine another counter into this one — used to aggregate per-thread
    /// counters at the end of classification. Mirrors C++ `operator+=`.
    pub fn merge(&mut self, other: &ReadCounts) {
        self.n_reads += other.n_reads;
        self.n_kmers += other.n_kmers;
        self.kmers.merge(&other.kmers);
    }
}

impl std::ops::AddAssign for ReadCounts {
    /// Forwards to [`merge`](Self::merge) so `a += b` matches the C++
    /// `operator+=` semantics.
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
    /// Construct an empty exact counter.
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of reads counted for this taxon.
    pub fn read_count(&self) -> u64 {
        self.n_reads
    }

    /// Increment the read counter by one.
    pub fn increment_read_count(&mut self) {
        self.n_reads += 1;
    }

    /// Total number of k-mers observed (with multiplicity).
    pub fn kmer_count(&self) -> u64 {
        self.n_kmers
    }

    /// Exact distinct-k-mer count (set cardinality).
    pub fn distinct_kmer_count(&self) -> u64 {
        self.kmers.len() as u64
    }

    /// Record a k-mer: bumps the total count and inserts it into the
    /// exact-set of distinct k-mers.
    pub fn add_kmer(&mut self, kmer: u64) {
        self.n_kmers += 1;
        self.kmers.insert(kmer);
    }

    /// Combine another exact counter into this one. Mirrors the C++
    /// `operator+=` plus the set-union overload at the bottom of
    /// `readcounts.h`.
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

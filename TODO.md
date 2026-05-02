# TODO

## High priority

- [x] Deterministic build mode — safe-prefix parallel insertion algorithm
- [x] Multi-threaded classification — rayon within-batch parallelism with ordered output
- [x] Remove C++ FFI entirely — deleted `ffi.rs`, removed `cc` build dep, pure Rust

## Medium priority

- [x] Full genome download pipeline — assembly_summary.txt parsing, genome download with taxid-prefixed headers, tar extraction for taxdump
- [x] Gzip/bzip2 compressed input — auto-detect and decompress `.gz`/`.bz2` input files
- [x] FASTQ output format — classified/unclassified output preserves input format
- [x] Second pass verification — all core C++ functions covered (~105/120, remaining are daemon/BLAST/utilities)

## Lower priority

- [x] Daemon mode — stdin-based persistent server for cached database access
- [x] BLAST-to-FASTA converter — BLAST binary database to FASTA conversion (index parsing, 2-bit nucleotide decoding, protein decoding)

## Performance optimization ideas

- [x] Switch HashMap to ahash — replaced SipHash with ahash in hit_counts, TaxonCounts, taxonomy ID map
- [x] Reuse output buffer — classify_sequence writes to caller-owned String, avoiding per-sequence allocation
- [x] Reuse scanner Vec — load_sequence reuses existing buffer capacity via extend_from_slice
- [x] Use write! instead of format! — avoid intermediate String allocations in hitlist building
- [~] Borrow sequence in scanner — investigated; load_sequence copy is negligible (<1% of time). The real cost is per-character scanning (41% CPU). Buffer reuse already applied.
- [x] SIMD-accelerated minimizer scanning — tried SSSE3 precompute of lookup codes (16 bytes at a time). Did NOT help: the 256-byte lookup table is always hot in L1 cache, and the extra memory pass offset the gains. Instead, `get_unchecked` bounds elimination in the inner loop gave ~3% improvement.
- [ ] Batch hash lookups with prefetch — prefetch hash table entries for the next minimizer while processing the current one (cache-line optimization).

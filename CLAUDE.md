# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Rust port of Kraken 2 — a taxonomic sequence classification system for DNA/RNA sequences. The original C++ code lives under `kraken2/` as a reference. The Rust rewrite lives in `src/` as a single crate exposing both a library (`kraken2_pure_rs`, `src/lib.rs`) and a CLI binary (`kraken2`, `src/bin/kraken2.rs` — uses `clap` subcommands `classify`, `build`, `estimate`, `download`, `inspect`, `lookup`, `clean`, `dump`, `mask`, `blast`, `mmtest`).

Output must be byte-identical to the original C++ code. Database files (`hash.k2d`, `taxo.k2d`, `opts.k2d`) are wire-compatible in both directions — a database built by either implementation must work in the other.

The repo has a second registered worktree at `/home/mahogny/github/claude/kraken-rs`.

## Build Commands

```bash
# Build Rust project
cargo build

# Build in release mode (for benchmarks, use native CPU)
RUSTFLAGS="-C target-cpu=native" cargo build --release

# Run tests (compares against C++-generated reference data in tests/reference/)
cargo test

# Run a single test
cargo test test_murmurhash3_known_values

# Build original C++ (for reference/comparison only)
make -C kraken2/src all
```

Build requirements: Rust toolchain only. No C++ compiler needed for the Rust crate.

`tests/integration_test.rs` runs against pre-generated reference outputs in `tests/reference/` (`hash.k2d`, `taxo.k2d`, `opts.k2d`, `classify_all.txt`, `classify_covid.txt`, `estimate_capacity.txt`). When changing anything that touches database layout, hashing, minimizer extraction, or classify/report output, regenerate the C++ reference and re-run these — diverging silently here means breaking on-disk compatibility.

## Rust Module Structure (`src/`)

| Module | Status | C++ Equivalent |
|--------|--------|----------------|
| `types.rs` | Done | `kraken2_data.h`, `kv_store.h` — IndexOptions (#[repr(C)]), TaxonomyNode, CompactHashCell, type aliases |
| `hash.rs` | Done | `kv_store.h` — MurmurHash3 finalizer |
| `utilities.rs` | Done | `utilities.h/cc` — ExpandSpacedSeedMask, SplitString |
| `aa_translate.rs` | Done | `aa_translate.h/cc` — 6-frame DNA-to-protein translation |
| `minimizer.rs` | Done | `mmscanner.h/cc` — MinimizerScanner with sliding window, spaced seeds, revcom |
| `ffi.rs` | Legacy shim | Thin wrappers kept for older internal callers; no C ABI any more — prefer the native module APIs |
| `omp.rs` | Stub | Single-thread fill-ins for OpenMP calls still referenced in ported code (`omp_get_thread_num`, etc.); real parallelism uses `rayon`/`parking_lot` |
| `taxonomy.rs` | Done | `taxonomy.h/cc` — NCBITaxonomy, Taxonomy, LCA, BFS conversion, binary I/O |
| `compact_hash.rs` | Done | `compact_hash.h/cc` — CompactHashTable with linear probing, zone locks, binary I/O |
| `mmap_file.rs` | Done | `mmap_file.h/cc` — memmap2 wrapper |
| `seq.rs` | Done | `seqreader.h/cc` — BatchSequenceReader for FASTA/FASTQ |
| `reports.rs` | Done | `reports.h/cc` — Kraken-style and MPA-style report generation |
| `hyperloglog.rs` | Done | `hyperloglogplus.h/cc` — HyperLogLogPlusMinus with Ertl estimator |
| `readcounts.rs` | Done | `readcounts.h` — ReadCounts with HLL and exact variants |
| `classify.rs` | Done | `classify.cc` — ClassifySequence, ResolveTree, hitlist formatting, full pipeline |
| `build_db.rs` | Done | `build_db.cc` — database building with fast mode, SetMinimizerLCA |
| `estimate.rs` | Done | `estimate_capacity.cc` — range-partitioned minimizer counting |
| `lookup.rs` | Done | `lookup_accession_numbers.cc` — accession-to-taxid lookup |
| `dust.rs` | Done | `k2mask.cc` — SDust low-complexity masking |
| `download.rs` | Done | Perl/Bash download scripts — NCBI downloads via ureq |
| `dump_table.rs` | Done | `dump_table.cc` — dump hash table contents to text |
| `blast.rs` | Done | `blast_to_fasta.c`, `blast_defline.c`, `blast_utils.c` — BLAST DB → FASTA conversion with optional taxid annotation |
| `mmtest.rs` | Done | small minimizer self-test command |
| `lib.rs` | — | Re-exports every module above as `pub mod`; library users consume `kraken2_pure_rs::<module>` |

## Critical Binary Format Details

These formats must be exactly matched for database compatibility:

**IndexOptions** (64 bytes): `k(8) l(8) spaced_seed_mask(8) toggle_mask(8) dna_db(1) _pad(7) minimum_acceptable_hash_value(8) revcom_version(4) db_version(4) db_type(4) _pad(4)`

**CompactHashTable** (hash.k2d): `capacity(8B) size(8B) key_bits(8B) value_bits(8B) cells(capacity × 4B)` — all header fields are `size_t` (8 bytes). LINEAR_PROBING mode (step=1) is the default.

**Taxonomy** (taxo.k2d): `"K2TAXDAT\0"(9B) node_count(8B) name_data_len(8B) rank_data_len(8B) nodes(N × 56B) name_data rank_data`

**TaxonomyNode**: 7 × u64 = 56 bytes: parent_id, first_child, child_count, name_offset, rank_offset, external_id, godparent_id

## Original C++ Architecture (`kraken2/src/`)

| Module | Role |
|--------|------|
| `classify.cc` | Main classifier — reads sequences, queries hash table, resolves taxonomy, writes reports |
| `build_db.cc` | Builds the compact hash table database from sequence libraries |
| `compact_hash.h/cc` | Probabilistic key-value store: 32-bit cells with linear probing, 256 zone locks |
| `taxonomy.h/cc` | Taxonomy tree with NCBI-to-internal ID mapping (BFS ordering) and LCA computation |
| `mmscanner.h/cc` | Minimizer extraction: deque sliding window, DNA 2-bit/protein 4-bit encoding, revcom |
| `seqreader.h/cc` | Batched FASTA/FASTQ reader (wraps `kseq.h`) |
| `reports.h/cc` | Kraken-style and MetaPhlAn-style report output |
| `hyperloglogplus.h/cc` | HyperLogLog++ for distinct k-mer counting |
| `mmap_file.h/cc` | Memory-mapped file I/O for database loading |
| `aa_translate.h/cc` | DNA-to-protein translation (AGCT ordering, not ACGT) |

**Data flow:** `estimate_capacity` -> `build_db` (creates database) -> `classify` (queries database)

## Key Dependencies

- `rayon` — replaces OpenMP parallel for
- `parking_lot` — replaces omp_lock_t zone locks
- `memmap2` — replaces MMapFile
- `noodles` — FASTA/FASTQ parsing (replaces kseq.h)
- `ureq` — HTTP downloads (replaces wget/rsync)
- `flate2`/`bzip2` — compression (replaces zlib/subprocess calls)
- `clap` — CLI parsing (replaces Perl wrapper argument handling)

## Working notes

- The release profile sets `lto = true` and `codegen-units = 1` — do not silently disable these for benchmarks; they materially change measured numbers.
- Performance comparisons in `BENCHMARKS.md` are translation-evaluation, not marketing claims; do not promote them as general performance statements when editing the README.
- `kraken2/` (the original C++ tree) is intentionally vendored for diffing behaviour — when porting a fix or feature, treat the C++ file as the spec and reproduce its observable behaviour, *including* known bugs, unless explicitly told to deviate.

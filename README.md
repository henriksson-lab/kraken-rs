# kraken2

A Rust port of [Kraken 2](https://github.com/DerrickWood/kraken2) — a taxonomic sequence classification system for DNA/RNA sequences.

This crate provides both a command-line tool and a library for classifying biological sequences against a Kraken 2 database. 

This is a translation of the original code and not the authorative implementation. This code should generate bitwise
equal output to the original. Please report any deviations

The aim of this project is to increase performance, especially by providing this code through a type-safe library interface.
The code can also be compiled to be used for webassembly.

## Installation

```bash
cargo install --path .
```

Or build from source:

```bash
cargo build --release
```

For best performance, build with native CPU optimizations:

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

### Build Requirements

- Rust toolchain (1.70+)
- C++11 compiler (for FFI bridge during transition period)
- OpenMP runtime (`libgomp`)
- zlib (`libz`)

## CLI Usage

### Classify sequences

```bash
# Classify sequences against a database
kraken2 classify -d /path/to/db -O output.txt input.fa

# With confidence threshold
kraken2 classify -d /path/to/db -T 0.2 -O output.txt input.fa

# Paired-end reads
kraken2 classify -d /path/to/db -P -O output.txt reads_1.fq reads_2.fq

# Generate a report
kraken2 classify -d /path/to/db -R report.txt -O output.txt input.fa

# Memory-mapped database (lower RAM usage)
kraken2 classify -d /path/to/db --memory-mapping -O output.txt input.fa

# Print scientific names instead of taxids
kraken2 classify -d /path/to/db --use-names -O output.txt input.fa
```

### Build a database

```bash
# Estimate capacity
cat library/*.fna | kraken2 estimate -k 35 -l 31

# Build database
cat library/*.fna | kraken2 build \
  -H db/hash.k2d -t db/taxo.k2d -o db/opts.k2d \
  -m seqid2taxid.map -n taxonomy/ \
  -k 35 -l 31 -c 92526
```

### Inspect a database

```bash
kraken2 inspect -d /path/to/db
kraken2 inspect -d /path/to/db --skip-counts
```

## Library Usage

```rust
use kraken2::classify::{ClassifyOptions, classify_sequence, resolve_tree, run_classify};
use kraken2::compact_hash::CompactHashTable;
use kraken2::taxonomy::Taxonomy;
use kraken2::types::IndexOptions;
use std::collections::HashMap;

// Load a database
let idx_opts = IndexOptions::read_from_file("db/opts.k2d").unwrap();
let mut taxonomy = Taxonomy::from_file("db/taxo.k2d", false).unwrap();
taxonomy.generate_external_to_internal_id_map();
let hash = CompactHashTable::from_file("db/hash.k2d", false).unwrap();

// Or use the high-level pipeline
let opts = ClassifyOptions::default();
let stats = run_classify(
    &["input.fa".to_string()],
    "db/hash.k2d", "db/taxo.k2d", "db/opts.k2d",
    &opts,
    Some("output.txt"),  // kraken output
    None,                // classified output
    None,                // unclassified output
    None,                // report
).unwrap();

println!("{} sequences classified", stats.total_classified);
```

### Building a taxonomy programmatically

```rust
use kraken2::taxonomy::NCBITaxonomy;

let mut ncbi_tax = NCBITaxonomy::new("nodes.dmp", "names.dmp").unwrap();
ncbi_tax.mark_node(2697049); // SARS-CoV-2
ncbi_tax.mark_node(694009);  // SARS coronavirus
ncbi_tax.convert_to_kraken_taxonomy("taxo.k2d").unwrap();
```

### Working with minimizers

```rust
use kraken2::minimizer::MinimizerScanner;
use kraken2::types::*;

let mut scanner = MinimizerScanner::new(
    35, 31,
    DEFAULT_SPACED_SEED_MASK,
    true,  // DNA mode
    DEFAULT_TOGGLE_MASK,
    CURRENT_REVCOM_VERSION,
);

scanner.load_sequence("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", 0, usize::MAX);
while let Some(minimizer) = scanner.next_minimizer() {
    if !scanner.is_ambiguous() {
        println!("Minimizer: {}", minimizer);
    }
}
```

## Benchmarks

Classification of 3,800 viral sequences (21.9 MB) against a small test database (38 reference genomes), on x86_64 Linux (40 cores):

| Threads | C++ (avg 3 runs) | Rust (avg 3 runs) | Ratio |
|---------|------------------|-------------------|-------|
| 1 | 3,019 ms | 3,019 ms | 1.00x |
| 4 | 859 ms | 925 ms | 1.08x |
| 8 | — | 557 ms | — |

Rust single-threaded and multi-threaded performance are on par with the C++ original. Multi-threading uses rayon for within-batch parallelism with ordered output.

*Benchmarked with `RUSTFLAGS="-C target-cpu=native"` on the same hardware. Database: 38 viral genomes, k=35, l=31.*

## Output Format

The classification output follows the standard Kraken 2 format:

```
C/U  header  taxid  length  hitlist
```

- `C` = classified, `U` = unclassified
- `header` = sequence ID
- `taxid` = NCBI taxonomy ID (0 if unclassified)
- `length` = sequence length in bp (or `len1|len2` for paired reads)
- `hitlist` = space-separated `taxid:count` pairs

## Database Compatibility

This crate reads and writes database files (`hash.k2d`, `taxo.k2d`, `opts.k2d`) that are binary-compatible with the original Kraken 2 C++ implementation. Databases built with either version can be used interchangeably.

## Ecosystem Integration

This crate is designed to work well with the Rust bioinformatics ecosystem:

- **[noodles](https://crates.io/crates/noodles)** — Used internally for FASTA/FASTQ parsing. Compatible with noodles-based pipelines.
- **[rayon](https://crates.io/crates/rayon)** — Used for data parallelism (replaces OpenMP from the C++ version).
- **[memmap2](https://crates.io/crates/memmap2)** — Memory-mapped database access for low-RAM systems.
- **[flate2](https://crates.io/crates/flate2)** / **[bzip2](https://crates.io/crates/bzip2)** — In-process decompression (no external gzip/bzip2 needed).
- **[ureq](https://crates.io/crates/ureq)** — HTTP downloads for NCBI data (replaces wget/rsync).

The library API accepts `Read` trait objects, so it integrates naturally with any Rust I/O source (files, network streams, in-memory buffers).

## License

MIT

# kraken2

A Rust port of [Kraken 2](https://github.com/DerrickWood/kraken2) — a taxonomic sequence classification system for DNA/RNA sequences.

This crate provides both a command-line tool and a library for classifying biological sequences against a Kraken 2 database.

* 2026-04-23: Possibly ready for real world data. If you decide to use it, ensure to compare the output on our own data with the original Kraken 2 software, as bugs may still remain
* Possibly about 10% faster than the original (on par)

## This is an LLM-mediated faithful (hopefully) translation, not the original code! 

Most users should probably first see if the existing original code works for them, unless they have reason otherwise. The original source
may have newer features and it has had more love in terms of fixing bugs. In fact, we aim to replicate bugs if they are present, for the
sake of reproducibility! (but then we might have added a few more in the process)

There are however cases when you might prefer this Rust version. We generally agree with [this manifesto](https://rewrites.bio/) but more specifically:
* We have had many issues with ensuring that our software works using existing containers (Docker, PodMan, Singularity). One size does not fit all and it eats our resources trying to keep up with every way of delivering software
* Common package managers do not work well. It was great when we had a few Linux distributions with stable procedures, but now there are just too many ecosystems (Homebrew, Conda). Conda has an NP-complete resolver which does not scale. Homebrew is only so-stable. And our dependencies in Python still break. These can no longer be considered professional serious options. Meanwhile, Cargo enables multiple versions of packages to be available, even within the same program(!)
* The future is the web. We deploy software in the web browser, and until now that has meant Javascript. This is a language where even the == operator is broken. Typescript is one step up, but a game changer is the ability to compile Rust code into webassembly, enabling performance and sharing of code with the backend. Translating code to Rust enables new ways of deployment and running code in the browser has especial benefits for science - researchers do not have deep pockets to run servers, so pushing compute to the user enables deployment that otherwise would be impossible
* Old CLI-based utilities are bad for the environment(!). A large amount of compute resources are spent creating and communicating via small files, which we can bypass by using code as libraries. Even better, we can avoid frequent reloading of databases by hoisting this stage, with up to 100x speedups in some cases. Less compute means faster compute and less electricity wasted
* LLM-mediated translations may actually be safer to use than the original code. This article shows that [running the same code on different operating systems can give somewhat different answers](https://doi.org/10.1038/nbt.3820). This is a gap that Rust+Cargo can reduce. Typesafe interfaces also reduce coding mistakes and error handling, as opposed to typical command-line scripting

But:

* **This approach should still be considered experimental**. The LLM technology is immature and has sharp corners. But there are opportunities to reap, and the genie is not going back into the bottle. This translation is as much aimed to learn how to improve the technology and get feedback on the results.
* Translations are not endorsed by the original authors unless otherwise noted. **Do not send bug reports to the original developers**. Use our Github issues page instead.
* **Do not treat any benchmark here as a general performance claim**. They are used to evaluate the translation on specific workloads. Current measurements are kept in [BENCHMARKS.md](BENCHMARKS.md), and the current `classify` comparison is a fair optimized-build benchmark for this repository rather than a marketing number
* **Check the original Github pages for information about the package**. This README is kept sparse on purpose. It is not meant to be the primary source of information
* **If you are the author of the original code and wish to move to Rust, you can obtain ownership of this repository and crate**. Until then, our commitment is to offer an as-faithful-as-possible translation of a snapshot of your code. If we find serious bugs, we will report them to you. Otherwise we will just replicate them, to ensure comparability across studies that claim to use package XYZ v.666. Think of this like a fancy Ubuntu .deb-package of your software - that is how we treat it

This blurb might be out of date. Go to [this page](https://github.com/henriksson-lab/rustification) for the latest information and further information about how we approach translation


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

## CLI Usage

### Classify sequences

```bash
# Classify sequences against a database
kraken2 classify -d /path/to/db -O output.txt input.fa

# With confidence threshold
kraken2 classify -d /path/to/db -T 0.2 -O output.txt input.fa

# Paired-end reads (two files)
kraken2 classify -d /path/to/db -P -O output.txt reads_1.fq reads_2.fq

# Paired-end reads (interleaved single file)
kraken2 classify -d /path/to/db -S -O output.txt interleaved.fq

# Multi-threaded
kraken2 classify -d /path/to/db -p 8 -O output.txt input.fa

# Generate a report
kraken2 classify -d /path/to/db -R report.txt -O output.txt input.fa

# Gzip/bzip2 compressed input (auto-detected)
kraken2 classify -d /path/to/db -O output.txt input.fa.gz

# Memory-mapped database (lower RAM usage)
kraken2 classify -d /path/to/db --memory-mapping -O output.txt input.fa

# Print scientific names instead of taxids
kraken2 classify -d /path/to/db --use-names -O output.txt input.fa
```

### Build a database

```bash
# Download taxonomy from NCBI
kraken2 download -d mydb --taxonomy

# Download a genomic library (e.g. viral)
kraken2 download -d mydb --library viral

# Estimate hash table capacity
cat mydb/library/viral/library.fna | kraken2 estimate -k 35 -l 31

# Build the database
cat mydb/library/viral/library.fna | kraken2 build \
  -H mydb/hash.k2d -t mydb/taxo.k2d -o mydb/opts.k2d \
  -m mydb/library/viral/prelim_map.txt -n mydb/taxonomy \
  -k 35 -l 31 -c 92526

# Clean intermediate files
kraken2 clean -d mydb
```

### Inspect a database

```bash
kraken2 inspect -d /path/to/db
kraken2 inspect -d /path/to/db --skip-counts
```

## Library Usage

### Classify sequences (load once, classify many)

```rust
use kraken2_rs::classify::{ClassifyDb, ClassifyOptions};
use kraken2_rs::types::Sequence;

// Load database once
let db = ClassifyDb::from_directory("path/to/db").unwrap();
let opts = ClassifyOptions::default();

// Classify a batch of sequences
let sequences = vec![
    Sequence {
        header: "read1".to_string(),
        seq: "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_string(),
        ..Default::default()
    },
];

let results = db.classify(&sequences, &opts);
for r in &results {
    println!("{}: taxid={}", if r.classified { "C" } else { "U" }, r.external_id);
}

// Classify more batches with the same loaded database — no reload cost
let results2 = db.classify(&more_sequences, &opts);

// Or classify a single sequence
let result = db.classify_one(&single_seq, &opts);
```

### Read sequences from any source

```rust
use kraken2_rs::seq::BatchSequenceReader;

// From a file (auto-detects .gz and .bz2)
let mut reader = BatchSequenceReader::new(Some("input.fa.gz")).unwrap();

// From in-memory data
let data = b">seq1\nACGT\n>seq2\nTTTT\n";
let mut reader = BatchSequenceReader::from_reader(std::io::Cursor::new(data));

reader.load_block(1_000_000);
while let Some(seq) = reader.next_sequence() {
    println!("{}: {} bp", seq.header, seq.seq.len());
}
```

### Build a taxonomy programmatically

```rust
use kraken2_rs::taxonomy::NCBITaxonomy;

let mut ncbi_tax = NCBITaxonomy::new("nodes.dmp", "names.dmp").unwrap();
ncbi_tax.mark_node(2697049); // SARS-CoV-2
ncbi_tax.mark_node(694009);  // SARS coronavirus
ncbi_tax.convert_to_kraken_taxonomy("taxo.k2d").unwrap();
```

### Work with minimizers

```rust
use kraken2_rs::minimizer::MinimizerScanner;
use kraken2_rs::types::*;

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

## Output Format

## Benchmark Notes

The current stored benchmark is in [BENCHMARKS.md](BENCHMARKS.md).

At the moment, the most recent real-data `classify` benchmark on a `600,000` read-pair subset from `/husky/henriksson/atrandi/rawdata/241206_novaseq_wgs3/filtered` shows Rust slightly faster than the original C++ on the repository reference DB:

- `-p 1`: C++ `22.54 s`, Rust `19.68 s`
- `-p 4`: C++ `6.73 s`, Rust `5.82 s`
- peak RSS on a measured `-p 4` run: C++ `187,568 kB`, Rust `59,812 kB` (`68.1%` less for Rust)

Those numbers are specific to that workload and should be read as translation-evaluation data, not as a blanket speed claim. For this repository they are a fair benchmark: both sides used normal optimized builds, neither side used `-march=native`, and checking Rust with and without LTO showed only a small effect.

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

- **[noodles](https://crates.io/crates/noodles)** — Compatible with noodles-based sequence processing pipelines.
- **[rayon](https://crates.io/crates/rayon)** — Data parallelism (replaces OpenMP from the C++ version).
- **[memmap2](https://crates.io/crates/memmap2)** — Memory-mapped database access for low-RAM systems.
- **[flate2](https://crates.io/crates/flate2)** / **[bzip2](https://crates.io/crates/bzip2)** — In-process decompression (no external gzip/bzip2 needed).
- **[ureq](https://crates.io/crates/ureq)** — HTTP downloads for NCBI data (replaces wget/rsync).

The library API accepts `Read` trait objects, so it integrates naturally with any Rust I/O source (files, network streams, in-memory buffers).

## License

MIT (derived from the original source)

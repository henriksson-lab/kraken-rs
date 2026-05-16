//! `kraken2_pure_rs` — a pure-Rust port of [Kraken 2], a taxonomic sequence
//! classification system for DNA/RNA reads.
//!
//! This crate exposes every module of the port as a `pub mod` so downstream users can
//! depend on individual pieces (e.g. just [`taxonomy`], [`compact_hash`], or
//! [`minimizer`]) without going through the CLI. The `kraken2` binary in
//! `src/bin/kraken2.rs` is a thin wrapper around [`cli::run_cli_from`].
//!
//! On-disk database files (`hash.k2d`, `taxo.k2d`, `opts.k2d`) are wire-compatible with
//! the original C++ Kraken 2 in both directions — see each module's docs for the exact
//! binary layout.
//!
//! [Kraken 2]: https://github.com/DerrickWood/kraken2

// Allow C-style indexing loops in code ported directly from C++
#![allow(clippy::needless_range_loop)]
// Allow many-argument functions matching C++ API signatures
#![allow(clippy::too_many_arguments)]

pub mod aa_translate;
pub mod blast;
pub mod build_db;
pub mod classify;
pub mod cli;
pub mod compact_hash;
pub mod download;
pub mod dump_table;
pub mod dust;
pub mod estimate;
pub mod ffi;
pub mod hash;
pub mod hyperloglog;
pub mod lookup;
pub mod minimizer;
pub mod mmap_file;
pub mod mmtest;
pub mod omp;
pub mod readcounts;
pub mod reports;
pub mod seq;
pub mod taxonomy;
pub mod types;
pub mod utilities;

// Allow C-style indexing loops in code ported directly from C++
#![allow(clippy::needless_range_loop)]
// Allow many-argument functions matching C++ API signatures
#![allow(clippy::too_many_arguments)]

pub mod aa_translate;
pub mod blast;
pub mod build_db;
pub mod classify;
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

// Allow C-style indexing loops in code ported directly from C++
#![allow(clippy::needless_range_loop)]
// Allow many-argument functions matching C++ API signatures
#![allow(clippy::too_many_arguments)]

pub mod types;
pub mod hash;
pub mod utilities;
pub mod aa_translate;
pub mod minimizer;
pub mod taxonomy;
pub mod compact_hash;
pub mod mmap_file;
pub mod seq;
pub mod reports;
pub mod hyperloglog;
pub mod readcounts;
pub mod classify;
pub mod build_db;
pub mod estimate;
pub mod lookup;
pub mod dust;
pub mod blast;
pub mod download;

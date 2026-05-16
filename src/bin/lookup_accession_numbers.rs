//! CLI shim for the `lookup_accession_numbers` binary.
//!
//! Thin wrapper that forwards `argv` to
//! [`kraken2_pure_rs::lookup::lookup_accession_numbers_main`]. Resolved
//! `seqid\ttaxid` lines are written to stdout; an `unmapped.txt` file is
//! produced in the current directory for any unresolved accessions, matching
//! the original C++ `lookup_accession_numbers` executable.

use std::process;

/// Run the `lookup_accession_numbers` command. Prints an error message and
/// exits with a non-zero status if the underlying driver fails.
fn main() {
    let args: Vec<String> = std::env::args().collect();
    if let Err(e) = kraken2_pure_rs::lookup::lookup_accession_numbers_main(&args) {
        eprintln!("{e}");
        process::exit(1);
    }
}

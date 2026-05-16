//! CLI entry point for the `blast_to_fasta` tool.
//!
//! Thin wrapper around [`kraken2_pure_rs::blast::blast_to_fasta_main`].
//! See the `blast` module for full BLAST-database → FASTA conversion logic.

use std::process;

/// Forwards process arguments to `blast::blast_to_fasta_main` and exits
/// non-zero on error.
fn main() {
    let args: Vec<String> = std::env::args().collect();
    if let Err(e) = kraken2_pure_rs::blast::blast_to_fasta_main(&args) {
        eprintln!("{e}");
        process::exit(1);
    }
}

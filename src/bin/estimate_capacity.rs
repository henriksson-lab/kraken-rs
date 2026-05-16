//! CLI shim for the `estimate_capacity` binary.
//!
//! Thin wrapper that forwards `argv` to
//! [`kraken2_pure_rs::estimate::estimate_capacity_main`] and prints the
//! returned distinct-minimizer capacity estimate to stdout. Mirrors the
//! original C++ `estimate_capacity` executable.

use std::process;

/// Run the `estimate_capacity` command. Prints the estimated capacity on stdout,
/// or the error message on stderr with a non-zero exit code on failure.
fn main() {
    let args: Vec<String> = std::env::args().collect();
    match kraken2_pure_rs::estimate::estimate_capacity_main(&args) {
        Ok(capacity) => println!("{capacity}"),
        Err(e) => {
            eprintln!("{e}");
            process::exit(1);
        }
    }
}

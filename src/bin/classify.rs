//! Tiny CLI shim that forwards `classify` invocations to the library's
//! `classify_main`. Matches the original `classify` binary entry point in
//! `kraken2/src/classify.cc`.

use std::process;

/// Collect command-line arguments and dispatch to `classify_main`,
/// exiting with status 1 on error.
fn main() {
    let args: Vec<String> = std::env::args().collect();
    if let Err(e) = kraken2_pure_rs::classify::classify_main(&args) {
        eprintln!("{e}");
        process::exit(1);
    }
}

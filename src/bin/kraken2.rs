//! Thin `kraken2` binary entry point.
//!
//! All real work happens in the [`kraken2_pure_rs::cli`] module; this shim just
//! forwards the process argv into [`kraken2_pure_rs::cli::run_cli_from`] and
//! propagates the returned status code to the OS.

/// Forward `std::env::args()` into the library CLI dispatcher and exit with its
/// returned status code.
fn main() {
    std::process::exit(kraken2_pure_rs::cli::run_cli_from(std::env::args()));
}

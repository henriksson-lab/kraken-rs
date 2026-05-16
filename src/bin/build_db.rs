//! CLI shim for `build_db`. Forwards argv to
//! [`kraken2_pure_rs::build_db::build_db_main`] and exits non-zero on error.

use std::process;

/// Forward argv to the library entry point and translate any error into a
/// non-zero exit code (matching the original C++ binary's behaviour).
fn main() {
    let args: Vec<String> = std::env::args().collect();
    if let Err(e) = kraken2_pure_rs::build_db::build_db_main(&args) {
        eprintln!("{e}");
        process::exit(1);
    }
}

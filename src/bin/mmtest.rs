//! Standalone CLI shim around [`kraken2_pure_rs::mmtest::mmtest_main`]. Prints
//! the minimizer self-test output to stdout; mirrors the original C++
//! `mmtest` executable.

/// Run the minimizer self-test and print its output to stdout.
fn main() {
    print!("{}", kraken2_pure_rs::mmtest::mmtest_main());
}

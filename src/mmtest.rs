use crate::minimizer::MinimizerScanner;
use crate::types::{CURRENT_REVCOM_VERSION, DEFAULT_TOGGLE_MASK};

/// Minimizer self-test: scan a fixed 13-base sequence with k=10, l=5 and
/// return one hex-encoded minimizer per line. Port of the C++ `mmtest.cc`
/// program; useful as a smoke test that the scanner still produces the
/// same byte stream as the reference implementation.
pub fn mmtest_main() -> String {
    let seq = "ACGATCGACGACG";
    let mut scanner =
        MinimizerScanner::new(10, 5, 0, true, DEFAULT_TOGGLE_MASK, CURRENT_REVCOM_VERSION);
    scanner.load_sequence(seq, 0, usize::MAX);
    let mut out = String::new();
    while let Some(mmp) = scanner.next_minimizer() {
        out.push_str(&format!("{:016x}\n", mmp));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mmtest_main() {
        let out = mmtest_main();
        assert!(!out.is_empty());
        assert!(out.lines().all(|line| line.len() == 16));
    }
}

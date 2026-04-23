use crate::minimizer::MinimizerScanner;
use crate::types::{CURRENT_REVCOM_VERSION, DEFAULT_TOGGLE_MASK};

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

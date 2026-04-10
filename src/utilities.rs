/// Expand a spaced seed mask by a given bit expansion factor.
/// Each bit in the original mask is expanded to `factor` bits.
/// E.g., 010110 expanded by factor 2 becomes 001100111100.
///
/// Exact port of C++ `ExpandSpacedSeedMask()` from `utilities.cc`.
pub fn expand_spaced_seed_mask(mask: &mut u64, bit_expansion_factor: i32) {
    let mut new_mask: u64 = 0;
    let bits: u64 = (1u64 << bit_expansion_factor) - 1;

    for i in (0..64 / bit_expansion_factor).rev() {
        new_mask <<= bit_expansion_factor;
        if (*mask >> i) & 1 != 0 {
            new_mask |= bits;
        }
    }
    *mask = new_mask;
}

/// Split a string by a delimiter, up to max_fields pieces.
/// Exact port of C++ `SplitString()` from `utilities.cc`.
pub fn split_string(s: &str, delim: &str, max_fields: usize) -> Vec<String> {
    let mut output = Vec::new();
    let mut pos1 = 0;
    let mut field_ct = 0;

    while field_ct < max_fields {
        field_ct += 1;
        if let Some(pos2) = s[pos1..].find(delim) {
            output.push(s[pos1..pos1 + pos2].to_string());
            pos1 += pos2 + delim.len();
        } else {
            output.push(s[pos1..].to_string());
            break;
        }
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_expand_spaced_seed_mask() {
        // 010110 expanded by 2 -> 001100111100
        let mut mask = 0b010110u64;
        expand_spaced_seed_mask(&mut mask, 2);
        assert_eq!(mask, 0b001100111100);

        // All zeros stays zero
        let mut mask = 0u64;
        expand_spaced_seed_mask(&mut mask, 2);
        assert_eq!(mask, 0);

        // All ones: 1111 expanded by 3 -> 111111111111
        let mut mask = 0b1111u64;
        expand_spaced_seed_mask(&mut mask, 3);
        assert_eq!(mask, 0b111111111111);

        // 111000 expanded by 2 -> 110000000000 (wait, let me think)
        // Actually: bit 0=0, bit 1=0, bit 2=0, bit 3=1, bit 4=1, bit 5=1
        // Expansion iterates i from (64/factor - 1) down to 0
        // For factor=2, 64/2=32, i from 31..0
        // For mask=0b111000 (=56), bits 3,4,5 are set
        let mut mask = 0b111000u64;
        expand_spaced_seed_mask(&mut mask, 2);
        assert_eq!(mask, 0b111111000000);
    }

    #[test]
    fn test_split_string() {
        assert_eq!(split_string("a\tb\tc", "\t", usize::MAX), vec!["a", "b", "c"]);
        assert_eq!(split_string("a\tb\tc", "\t", 2), vec!["a", "b"]);
        assert_eq!(split_string("hello", "\t", usize::MAX), vec!["hello"]);
        assert_eq!(split_string("a::b::c", "::", usize::MAX), vec!["a", "b", "c"]);
    }
}

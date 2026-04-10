/// MurmurHash3 64-bit finalizer.
/// Exact port of the C++ `MurmurHash3()` from `kv_store.h`.
#[inline]
pub fn murmurhash3(mut key: u64) -> u64 {
    key ^= key >> 33;
    key = key.wrapping_mul(0xff51afd7ed558ccd);
    key ^= key >> 33;
    key = key.wrapping_mul(0xc4ceb9fe1a85ec53);
    key ^= key >> 33;
    key
}

#[cfg(test)]
mod tests {
    use super::*;

    // Reference values verified against C++ MurmurHash3 implementation
    const REFERENCE: &[(u64, u64)] = &[
        (0, 0),
        (1, 12994781566227106604),
        (2, 4233148493373801447),
        (42, 9297814886316923340),
        (255, 1297215527019907880),
        (256, 15401116602503760187),
        (1000, 12610127409379334721),
        (u64::MAX, 7256831767414464289),
        (u64::MAX - 1, 4216938840244723755),
        (0xdeadbeef, 15153440252345589164),
        (0xcafebabe, 987948222317082758),
        (0x123456789abcdef0, 1781385183907554200),
    ];

    #[test]
    fn test_murmurhash3_known_values() {
        for &(key, expected) in REFERENCE {
            assert_eq!(murmurhash3(key), expected,
                "MurmurHash3 mismatch for key {}", key);
        }
    }

    #[test]
    fn test_murmurhash3_properties() {
        // Hash of 0 should be 0 (identity for this finalizer)
        assert_eq!(murmurhash3(0), 0);
        // Different inputs should give different outputs
        assert_ne!(murmurhash3(1), murmurhash3(2));
        // Consistent across calls
        assert_eq!(murmurhash3(42), murmurhash3(42));
    }
}

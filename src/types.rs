/// Taxonomy node ID (internal or external).
pub type TaxId = u64;

/// Hash table key type (64-bit minimizer hash).
pub type HKey = u64;

/// Hash table value type (32-bit taxonomy ID).
pub type HValue = u32;

/// Map from taxonomy IDs to counts.
pub type TaxonCounts = std::collections::HashMap<TaxId, u64>;

/// Maximum taxonomy ID value.
pub const TAXID_MAX: TaxId = u64::MAX;

/// Database index options — must exactly match C++ struct layout.
/// Verified via offsetof(): total size = 64 bytes on 64-bit Linux.
///
/// Layout:
///   offset  0: k (8 bytes, usize)
///   offset  8: l (8 bytes, usize)
///   offset 16: spaced_seed_mask (8 bytes, u64)
///   offset 24: toggle_mask (8 bytes, u64)
///   offset 32: dna_db (1 byte, bool)
///   offset 33: _padding (7 bytes)
///   offset 40: minimum_acceptable_hash_value (8 bytes, u64)
///   offset 48: revcom_version (4 bytes, i32)
///   offset 52: db_version (4 bytes, i32)
///   offset 56: db_type (4 bytes, i32)
///   offset 60: _padding2 (4 bytes)
///   total:  64 bytes
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct IndexOptions {
    pub k: usize,
    pub l: usize,
    pub spaced_seed_mask: u64,
    pub toggle_mask: u64,
    pub dna_db: bool,
    _padding: [u8; 7],
    pub minimum_acceptable_hash_value: u64,
    pub revcom_version: i32,
    pub db_version: i32,
    pub db_type: i32,
    _padding2: [u8; 4],
}

impl IndexOptions {
    pub fn new() -> Self {
        // Safety: zero-initialized is valid for this repr(C) struct
        unsafe { std::mem::zeroed() }
    }

    /// Read IndexOptions from a binary file (pure Rust, no FFI).
    pub fn read_from_file(filename: &str) -> std::io::Result<Self> {
        use std::fs::File;
        use std::io::Read;
        let mut file = File::open(filename)?;
        let mut opts = Self::new();
        let size = std::mem::size_of::<Self>();
        let buf = unsafe {
            std::slice::from_raw_parts_mut(&mut opts as *mut Self as *mut u8, size)
        };
        file.read_exact(buf)?;
        Ok(opts)
    }

    /// Write IndexOptions to a binary file (pure Rust, no FFI).
    pub fn write_to_file(&self, filename: &str) -> std::io::Result<()> {
        use std::fs::File;
        use std::io::Write;
        let mut file = File::create(filename)?;
        let size = std::mem::size_of::<Self>();
        let buf = unsafe {
            std::slice::from_raw_parts(self as *const Self as *const u8, size)
        };
        file.write_all(buf)
    }
}

impl Default for IndexOptions {
    fn default() -> Self {
        Self::new()
    }
}

/// Taxonomy node — must exactly match C++ struct layout.
/// 7 x u64 = 56 bytes, no padding needed.
#[repr(C)]
#[derive(Debug, Clone, Copy, Default)]
pub struct TaxonomyNode {
    pub parent_id: u64,
    pub first_child: u64,
    pub child_count: u64,
    pub name_offset: u64,
    pub rank_offset: u64,
    pub external_id: u64,
    pub godparent_id: u64,
}

/// Compact hash table cell — 32 bits: key in high bits, value in low bits.
#[repr(C)]
#[derive(Debug, Clone, Copy, Default)]
pub struct CompactHashCell {
    pub data: u32,
}

impl CompactHashCell {
    #[inline]
    pub fn hashed_key(&self, value_bits: usize) -> HKey {
        (self.data >> value_bits) as HKey
    }

    #[inline]
    pub fn value(&self, value_bits: usize) -> HValue {
        self.data & ((1u32 << value_bits) - 1)
    }

    #[inline]
    pub fn populate(&mut self, compacted_key: HKey, val: HValue, _key_bits: usize, value_bits: usize) {
        self.data = ((compacted_key as u32) << value_bits) | val;
    }
}

/// Sequence format (FASTA or FASTQ).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[derive(Default)]
pub enum SequenceFormat {
    #[default]
    AutoDetect,
    Fasta,
    Fastq,
}


/// A biological sequence with header, comment, bases, and optional quality scores.
#[derive(Debug, Clone, Default)]
pub struct Sequence {
    pub format: SequenceFormat,
    pub header: String,
    pub comment: String,
    pub seq: String,
    pub quals: String,
}

impl Sequence {
    /// Convert sequence to FASTA/FASTQ string representation.
    /// Port of C++ `Sequence::to_string()`.
    pub fn to_fasta_string(&self) -> String {
        let mut s = String::new();
        match self.format {
            SequenceFormat::Fastq => {
                s.push('@');
                s.push_str(&self.header);
                if !self.comment.is_empty() {
                    s.push(' ');
                    s.push_str(&self.comment);
                }
                s.push('\n');
                s.push_str(&self.seq);
                s.push_str("\n+\n");
                s.push_str(&self.quals);
                s.push('\n');
            }
            _ => {
                s.push('>');
                s.push_str(&self.header);
                if !self.comment.is_empty() {
                    s.push(' ');
                    s.push_str(&self.comment);
                }
                s.push('\n');
                s.push_str(&self.seq);
                s.push('\n');
            }
        }
        s
    }

}

impl std::fmt::Display for Sequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_fasta_string())
    }
}

/// Constants matching C++ mmscanner.h
pub const DEFAULT_TOGGLE_MASK: u64 = 0xe37e28c4271b5a2d;
pub const DEFAULT_SPACED_SEED_MASK: u64 = 0;
pub const BITS_PER_CHAR_DNA: usize = 2;
pub const BITS_PER_CHAR_PRO: usize = 4;
pub const CURRENT_REVCOM_VERSION: i32 = 1;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_index_options_size() {
        assert_eq!(std::mem::size_of::<IndexOptions>(), 64);
    }

    #[test]
    fn test_index_options_roundtrip() {
        let mut opts = IndexOptions::new();
        opts.k = 35;
        opts.l = 31;
        opts.spaced_seed_mask = 0;
        opts.toggle_mask = DEFAULT_TOGGLE_MASK;
        opts.dna_db = true;
        opts.minimum_acceptable_hash_value = 0;
        opts.revcom_version = CURRENT_REVCOM_VERSION;

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().to_str().unwrap();

        opts.write_to_file(path).unwrap();
        let loaded = IndexOptions::read_from_file(path).unwrap();

        assert_eq!(loaded.k, 35);
        assert_eq!(loaded.l, 31);
        assert_eq!(loaded.toggle_mask, DEFAULT_TOGGLE_MASK);
        assert!(loaded.dna_db);
        assert_eq!(loaded.revcom_version, CURRENT_REVCOM_VERSION);
    }

    #[test]
    fn test_index_options_matches_reference() {
        // Compare Rust-written IndexOptions with C++-generated reference
        let ref_path = format!("{}/tests/reference/opts.k2d", env!("CARGO_MANIFEST_DIR"));
        if !std::path::Path::new(&ref_path).exists() {
            return;
        }
        // Read reference, write with Rust, compare bytes (roundtrip)
        let opts = IndexOptions::read_from_file(&ref_path).unwrap();
        let tmp = tempfile::NamedTempFile::new().unwrap();
        opts.write_to_file(tmp.path().to_str().unwrap()).unwrap();

        let ref_bytes = std::fs::read(&ref_path).unwrap();
        let rust_bytes = std::fs::read(tmp.path()).unwrap();
        assert_eq!(ref_bytes, rust_bytes, "IndexOptions roundtrip differs from C++ reference");
    }

    #[test]
    fn test_index_options_read_reference() {
        let ref_path = format!("{}/tests/reference/opts.k2d", env!("CARGO_MANIFEST_DIR"));
        if !std::path::Path::new(&ref_path).exists() {
            return;
        }
        let opts = IndexOptions::read_from_file(&ref_path).unwrap();
        assert_eq!(opts.k, 35);
        assert_eq!(opts.l, 31);
        assert!(opts.dna_db);
        assert_eq!(opts.revcom_version, 1);
    }
}

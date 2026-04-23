use std::fs::File;
use std::io::{self, Read, Write};
use std::sync::atomic::{AtomicUsize, Ordering};

use parking_lot::Mutex;

use crate::hash::murmurhash3;
use crate::mmap_file::MMapFile;
use crate::types::*;

const LOCK_ZONES: usize = 256;

/// Compact hash table — probabilistic key-value store with 32-bit cells.
/// Exact port of C++ `CompactHashTable` from `compact_hash.h/cc`.
///
/// Each cell packs a truncated key hash in the high bits and a value in the low bits.
/// Uses linear probing (LINEAR_PROBING mode, step=1).
/// Thread-safe writes via 256 zone locks.
pub struct CompactHashTable {
    capacity: usize,
    size: AtomicUsize,
    key_bits: usize,
    value_bits: usize,
    table: Vec<CompactHashCell>,
    zone_locks: Vec<Mutex<()>>,
    file_backed: bool,
    _mmap: Option<MMapFile>,
}

impl CompactHashTable {
    /// Create a new empty hash table with the given parameters.
    pub fn new(capacity: usize, key_bits: usize, value_bits: usize) -> Self {
        assert_eq!(
            key_bits + value_bits,
            32,
            "sum of key bits and value bits must equal 32"
        );
        assert!(key_bits > 0, "key bits cannot be zero");
        assert!(value_bits > 0, "value bits cannot be zero");

        let zone_locks: Vec<Mutex<()>> = (0..LOCK_ZONES).map(|_| Mutex::new(())).collect();

        CompactHashTable {
            capacity,
            size: AtomicUsize::new(0),
            key_bits,
            value_bits,
            table: vec![CompactHashCell { data: 0 }; capacity],
            zone_locks,
            file_backed: false,
            _mmap: None,
        }
    }

    /// Load a hash table from a binary file.
    /// Format: capacity(8B) size(8B) key_bits(8B) value_bits(8B) cells(capacity × 4B)
    pub fn from_file(filename: &str, memory_mapping: bool) -> io::Result<Self> {
        if memory_mapping {
            Self::from_file_mmap(filename)
        } else {
            Self::from_file_read(filename)
        }
    }

    pub fn from_cstr(filename: &str, memory_mapping: bool) -> io::Result<Self> {
        Self::from_file(filename, memory_mapping)
    }

    fn from_file_read(filename: &str) -> io::Result<Self> {
        let mut file = File::open(filename)?;
        let mut buf8 = [0u8; 8];

        file.read_exact(&mut buf8)?;
        let capacity = usize::from_le_bytes(buf8);
        file.read_exact(&mut buf8)?;
        let size = usize::from_le_bytes(buf8);
        file.read_exact(&mut buf8)?;
        let key_bits = usize::from_le_bytes(buf8);
        file.read_exact(&mut buf8)?;
        let value_bits = usize::from_le_bytes(buf8);

        let mut cell_bytes = vec![0u8; capacity * 4];
        file.read_exact(&mut cell_bytes)?;

        let table = unsafe {
            let mut table = vec![CompactHashCell { data: 0 }; capacity];
            std::ptr::copy_nonoverlapping(
                cell_bytes.as_ptr(),
                table.as_mut_ptr() as *mut u8,
                cell_bytes.len(),
            );
            table
        };

        Ok(CompactHashTable {
            capacity,
            size: AtomicUsize::new(size),
            key_bits,
            value_bits,
            table,
            zone_locks: (0..LOCK_ZONES).map(|_| Mutex::new(())).collect(),
            file_backed: false,
            _mmap: None,
        })
    }

    fn from_file_mmap(filename: &str) -> io::Result<Self> {
        let mmap = MMapFile::open_read_only(filename)?;
        let data = mmap.as_slice();

        let capacity = usize::from_le_bytes(data[0..8].try_into().unwrap());
        let size = usize::from_le_bytes(data[8..16].try_into().unwrap());
        let key_bits = usize::from_le_bytes(data[16..24].try_into().unwrap());
        let value_bits = usize::from_le_bytes(data[24..32].try_into().unwrap());

        let cell_data = &data[32..];
        let expected_len = capacity * 4;
        if cell_data.len() < expected_len {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Capacity mismatch in hash table file",
            ));
        }

        let table = unsafe {
            let mut table = vec![CompactHashCell { data: 0 }; capacity];
            std::ptr::copy_nonoverlapping(
                cell_data.as_ptr(),
                table.as_mut_ptr() as *mut u8,
                expected_len,
            );
            table
        };

        Ok(CompactHashTable {
            capacity,
            size: AtomicUsize::new(size),
            key_bits,
            value_bits,
            table,
            zone_locks: (0..LOCK_ZONES).map(|_| Mutex::new(())).collect(),
            file_backed: true,
            _mmap: Some(mmap),
        })
    }

    /// Write the hash table to a binary file.
    /// Format matches C++ `CompactHashTable::WriteTable()`.
    pub fn write_table(&self, filename: &str) -> io::Result<()> {
        let mut file = File::create(filename)?;

        let capacity = self.capacity as u64;
        let size = self.size.load(Ordering::Relaxed) as u64;
        let key_bits = self.key_bits as u64;
        let value_bits = self.value_bits as u64;

        file.write_all(&capacity.to_le_bytes())?;
        file.write_all(&size.to_le_bytes())?;
        file.write_all(&key_bits.to_le_bytes())?;
        file.write_all(&value_bits.to_le_bytes())?;

        let cell_bytes = unsafe {
            std::slice::from_raw_parts(
                self.table.as_ptr() as *const u8,
                self.table.len() * std::mem::size_of::<CompactHashCell>(),
            )
        };
        file.write_all(cell_bytes)?;

        Ok(())
    }

    /// Look up a value by key. Returns 0 if not found.
    /// Exact port of C++ `CompactHashTable::Get()`.
    pub fn get(&self, key: u64) -> HValue {
        let hc = murmurhash3(key);
        let compacted_key = (hc >> (32 + self.value_bits)) as HKey;
        let mut idx = (hc % self.capacity as u64) as usize;
        let first_idx = idx;
        let mut step: usize = 0;

        loop {
            if self.table[idx].value(self.value_bits) == 0 {
                break; // empty cell
            }
            if self.table[idx].hashed_key(self.value_bits) == compacted_key {
                return self.table[idx].value(self.value_bits);
            }
            if step == 0 {
                step = self.second_hash(hc);
            }
            idx = (idx + step) % self.capacity;
            if idx == first_idx {
                break; // exhausted table
            }
        }
        0
    }

    /// Find the index where a key would be stored. Returns Some(idx) if found.
    /// Exact port of C++ `CompactHashTable::FindIndex()`.
    pub fn find_index(&self, key: u64) -> Option<usize> {
        let hc = murmurhash3(key);
        let compacted_key = (hc >> (32 + self.value_bits)) as HKey;
        let mut idx = (hc % self.capacity as u64) as usize;
        let first_idx = idx;
        let mut step: usize = 0;

        loop {
            if self.table[idx].value(self.value_bits) == 0 {
                return None; // empty cell
            }
            if self.table[idx].hashed_key(self.value_bits) == compacted_key {
                return Some(idx);
            }
            if step == 0 {
                step = self.second_hash(hc);
            }
            idx = (idx + step) % self.capacity;
            if idx == first_idx {
                break;
            }
        }
        None
    }

    /// Compare-and-set: if the cell's value matches *old_value, update to new_value.
    /// Otherwise, set *old_value to the current cell value.
    /// Returns true if the set was successful.
    /// Exact port of C++ `CompactHashTable::CompareAndSet()`.
    pub fn compare_and_set(&self, key: u64, new_value: HValue, old_value: &mut HValue) -> bool {
        if self.file_backed || new_value == 0 {
            return false;
        }

        let hc = murmurhash3(key);
        let compacted_key = (hc >> (32 + self.value_bits)) as HKey;
        let mut idx = (hc % self.capacity as u64) as usize;
        let first_idx = idx;
        let mut step: usize = 0;
        let mut set_successful = false;

        loop {
            let zone = idx % LOCK_ZONES;
            let _lock = self.zone_locks[zone].lock();

            let cell_value = self.table[idx].value(self.value_bits);
            if cell_value == 0 || self.table[idx].hashed_key(self.value_bits) == compacted_key {
                // Found our slot
                if *old_value == cell_value {
                    // Safety: we hold the zone lock, so no concurrent writes to this cell
                    unsafe {
                        let cell_ptr = self.table.as_ptr().add(idx) as *mut CompactHashCell;
                        (*cell_ptr).populate(
                            compacted_key,
                            new_value,
                            self.key_bits,
                            self.value_bits,
                        );
                    }
                    if *old_value == 0 {
                        self.size.fetch_add(1, Ordering::Relaxed);
                    }
                    set_successful = true;
                } else {
                    *old_value = cell_value;
                }
                drop(_lock);
                return set_successful;
            }

            drop(_lock);

            if step == 0 {
                step = self.second_hash(hc);
            }
            idx = (idx + step) % self.capacity;
            if idx == first_idx {
                // Matches C++ errx(EX_SOFTWARE, ...) — capacity exhausted is fatal
                panic!("compact hash table capacity exceeded (occupancy 100%)");
            }
        }
    }

    /// Direct compare-and-set at a specific index.
    /// Exact port of C++ `CompactHashTable::DirectCompareAndSet()`.
    pub fn direct_compare_and_set(
        &self,
        idx: usize,
        key: u64,
        new_value: HValue,
        old_value: &mut HValue,
    ) -> bool {
        let hc = murmurhash3(key);
        let compacted_key = (hc >> (32 + self.value_bits)) as HKey;

        let zone = idx % LOCK_ZONES;
        let _lock = self.zone_locks[zone].lock();

        let cell_value = self.table[idx].value(self.value_bits);
        if *old_value == cell_value {
            unsafe {
                let cell_ptr = self.table.as_ptr().add(idx) as *mut CompactHashCell;
                (*cell_ptr).populate(compacted_key, new_value, self.key_bits, self.value_bits);
            }
            if *old_value == 0 {
                self.size.fetch_add(1, Ordering::Relaxed);
            }
            true
        } else {
            *old_value = cell_value;
            false
        }
    }

    /// Get value counts across all cells (parallel).
    pub fn get_value_counts(&self) -> TaxonCounts {
        let mut counts = TaxonCounts::new();
        for i in 0..self.capacity {
            let val = self.table[i].value(self.value_bits) as u64;
            if val != 0 {
                *counts.entry(val).or_insert(0) += 1;
            }
        }
        counts
    }

    /// LINEAR_PROBING: always returns 1.
    #[inline]
    fn second_hash(&self, _first_hash: u64) -> usize {
        1
    }

    pub fn capacity(&self) -> usize {
        self.capacity
    }

    pub fn size(&self) -> usize {
        self.size.load(Ordering::Relaxed)
    }

    pub fn key_bits(&self) -> usize {
        self.key_bits
    }

    pub fn value_bits(&self) -> usize {
        self.value_bits
    }

    pub fn occupancy(&self) -> f64 {
        self.size() as f64 / self.capacity as f64
    }

    pub fn destroy(self) {}
}

// The CompactHashTable uses interior mutability via zone locks for CompareAndSet,
// so it can be shared across threads safely.
unsafe impl Sync for CompactHashTable {}
unsafe impl Send for CompactHashTable {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compact_hash_basic_operations() {
        let cht = CompactHashTable::new(1000, 22, 10);

        assert_eq!(cht.capacity(), 1000);
        assert_eq!(cht.size(), 0);
        assert_eq!(cht.key_bits(), 22);
        assert_eq!(cht.value_bits(), 10);

        // Insert via compare_and_set
        let mut old_val: u32 = 0;
        assert!(cht.compare_and_set(42, 5, &mut old_val));
        assert_eq!(cht.size(), 1);
        assert_eq!(cht.get(42), 5);

        // Update via compare_and_set
        old_val = 5;
        assert!(cht.compare_and_set(42, 7, &mut old_val));
        assert_eq!(cht.get(42), 7);

        // Failed compare_and_set (wrong old value)
        old_val = 99;
        assert!(!cht.compare_and_set(42, 10, &mut old_val));
        assert_eq!(old_val, 7); // old_value updated to current value
    }

    #[test]
    fn test_compact_hash_write_read_roundtrip() {
        let cht = CompactHashTable::new(100, 22, 10);

        // Insert some values
        let keys = [10u64, 20, 30, 42, 100, 200, 500];
        let vals = [1u32, 2, 3, 4, 5, 6, 7];
        for (&k, &v) in keys.iter().zip(vals.iter()) {
            let mut old = 0u32;
            cht.compare_and_set(k, v, &mut old);
        }

        let tmp_dir = tempfile::tempdir().unwrap();
        let filepath = tmp_dir.path().join("test_hash.k2d");

        // Write
        cht.write_table(filepath.to_str().unwrap()).unwrap();

        // Read back
        let cht2 = CompactHashTable::from_file(filepath.to_str().unwrap(), false).unwrap();

        assert_eq!(cht2.capacity(), cht.capacity());
        assert_eq!(cht2.size(), cht.size());
        assert_eq!(cht2.key_bits(), cht.key_bits());
        assert_eq!(cht2.value_bits(), cht.value_bits());

        // Verify all values
        for (&k, &v) in keys.iter().zip(vals.iter()) {
            assert_eq!(cht2.get(k), v, "Value mismatch for key {}", k);
        }
    }

    #[test]
    fn test_compact_hash_load_reference() {
        // Load the C++-generated reference hash table and verify it works
        let ref_path = format!("{}/tests/reference/hash.k2d", env!("CARGO_MANIFEST_DIR"));
        if !std::path::Path::new(&ref_path).exists() {
            eprintln!("Skipping: reference data not available");
            return;
        }
        let cht = CompactHashTable::from_file(&ref_path, false).unwrap();
        assert_eq!(cht.capacity(), 92526);
        assert!(cht.size() > 0);
        assert_eq!(cht.key_bits() + cht.value_bits(), 32);

        // Verify write+read roundtrip produces identical binary
        let tmp_dir = tempfile::tempdir().unwrap();
        let out_path = tmp_dir.path().join("roundtrip.k2d");
        cht.write_table(out_path.to_str().unwrap()).unwrap();
        let ref_bytes = std::fs::read(&ref_path).unwrap();
        let out_bytes = std::fs::read(&out_path).unwrap();
        assert_eq!(
            ref_bytes, out_bytes,
            "Write roundtrip of reference hash table differs"
        );
    }
}

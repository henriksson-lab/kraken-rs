use memmap2::{Mmap, MmapMut};
use std::fs::{File, OpenOptions};
use std::io;

/// Memory-mapped file wrapper using memmap2.
/// Replaces the C++ MMapFile class.
#[derive(Default)]
pub enum MMapFile {
    ReadOnly {
        mmap: Mmap,
        _file: File,
    },
    ReadWrite {
        mmap: MmapMut,
        _file: File,
    },
    #[default]
    Empty,
}

impl MMapFile {
    /// Create an empty (invalid) MMapFile.
    pub fn new() -> Self {
        MMapFile::Empty
    }

    /// Open a file for read-only memory mapping.
    pub fn open_read_only(filename: &str) -> io::Result<Self> {
        let file = File::open(filename)?;
        let mmap = unsafe { Mmap::map(&file)? };
        Ok(MMapFile::ReadOnly { mmap, _file: file })
    }

    /// Open or create a file for read-write memory mapping.
    pub fn open_read_write(filename: &str, size: Option<u64>) -> io::Result<Self> {
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(size.is_some())
            .open(filename)?;

        if let Some(sz) = size {
            file.set_len(sz)?;
        }

        let mmap = unsafe { MmapMut::map_mut(&file)? };
        Ok(MMapFile::ReadWrite { mmap, _file: file })
    }

    pub fn open_file(
        &mut self,
        filename: &str,
        read_only: bool,
        create_size: Option<u64>,
    ) -> io::Result<()> {
        *self = if read_only {
            Self::open_read_only(filename)?
        } else {
            Self::open_read_write(filename, create_size)?
        };
        Ok(())
    }

    /// Get a read-only pointer to the mapped data.
    pub fn as_slice(&self) -> &[u8] {
        match self {
            MMapFile::ReadOnly { mmap, .. } => &mmap[..],
            MMapFile::ReadWrite { mmap, .. } => &mmap[..],
            MMapFile::Empty => &[],
        }
    }

    /// Get a mutable pointer to the mapped data (read-write mode only).
    pub fn as_mut_slice(&mut self) -> Option<&mut [u8]> {
        match self {
            MMapFile::ReadWrite { mmap, .. } => Some(&mut mmap[..]),
            _ => None,
        }
    }

    pub fn fptr(&self) -> *const u8 {
        match self {
            MMapFile::ReadOnly { mmap, .. } => mmap.as_ptr(),
            MMapFile::ReadWrite { mmap, .. } => mmap.as_ptr(),
            MMapFile::Empty => std::ptr::null(),
        }
    }

    /// Get the file size.
    pub fn filesize(&self) -> usize {
        match self {
            MMapFile::ReadOnly { mmap, .. } => mmap.len(),
            MMapFile::ReadWrite { mmap, .. } => mmap.len(),
            MMapFile::Empty => 0,
        }
    }

    /// Load file into OS page cache (equivalent to C++ LoadFile).
    pub fn load_file(&self) {
        // Touch all pages to fault them into memory.
        // The madvise(WILLNEED) hint is the Rust/OS equivalent.
        if let MMapFile::ReadOnly { mmap, .. } = self {
            #[cfg(unix)]
            {
                let ptr = mmap.as_ptr();
                let len = mmap.len();
                unsafe {
                    libc::madvise(ptr as *mut libc::c_void, len, libc::MADV_WILLNEED);
                }
            }
            #[cfg(not(unix))]
            {
                // Fallback: touch pages sequentially
                let page_size = 4096;
                let mut _sum: u8 = 0;
                for i in (0..mmap.len()).step_by(page_size) {
                    _sum = _sum.wrapping_add(mmap[i]);
                }
            }
        }
    }

    /// Sync changes to disk (read-write mode only).
    pub fn sync(&self) -> io::Result<()> {
        match self {
            MMapFile::ReadWrite { mmap, .. } => mmap.flush(),
            _ => Ok(()),
        }
    }

    pub fn sync_file(&self) -> io::Result<()> {
        self.sync()
    }

    pub fn close_file(&mut self) -> io::Result<()> {
        if matches!(self, MMapFile::Empty) {
            return Ok(());
        }
        self.sync_file()?;
        *self = MMapFile::Empty;
        Ok(())
    }

    /// Check if this is a valid (opened) mapping.
    pub fn is_valid(&self) -> bool {
        !matches!(self, MMapFile::Empty)
    }
}

impl Drop for MMapFile {
    fn drop(&mut self) {
        let _ = self.sync_file();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_open_file_read_only_and_close_file() {
        let tmp = NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), b"abc").unwrap();

        let mut mmap = MMapFile::new();
        mmap.open_file(tmp.path().to_str().unwrap(), true, None)
            .unwrap();
        assert!(mmap.is_valid());
        assert_eq!(mmap.filesize(), 3);
        assert!(!mmap.fptr().is_null());

        mmap.close_file().unwrap();
        assert!(!mmap.is_valid());
        assert_eq!(mmap.filesize(), 0);
    }

    #[test]
    fn test_open_file_read_write_and_sync_file() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path().to_str().unwrap().to_string();

        let mut mmap = MMapFile::new();
        mmap.open_file(&path, false, Some(4)).unwrap();
        if let Some(buf) = mmap.as_mut_slice() {
            buf.copy_from_slice(b"rust");
        }
        mmap.sync_file().unwrap();
        mmap.close_file().unwrap();

        assert_eq!(std::fs::read(path).unwrap(), b"rust");
    }
}

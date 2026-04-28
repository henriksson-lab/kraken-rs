use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};

use crate::types::{Sequence, SequenceFormat};

/// Batched sequence reader for FASTA/FASTQ files.
/// Replaces the C++ BatchSequenceReader (which wraps kseq.h).
/// Uses a simple line-based parser to match kseq.h behavior exactly.
/// Supports gzip and bzip2 compressed input (auto-detected by file extension).
pub struct BatchSequenceReader {
    reader: BufReader<Box<dyn Read + Send>>,
    seqs: Vec<Sequence>,
    curr: usize,
    size: usize,
    file_format: SequenceFormat,
    /// If set, the next record's header line has already been read.
    peeked_line: Option<String>,
    /// Scratch buffer reused across `read_line` calls — avoids allocating
    /// a fresh String per FASTA/FASTQ line in the hot path.
    line_buf: String,
}

/// Open a file, auto-detecting compression from extension.
fn open_possibly_compressed(filename: &str) -> io::Result<Box<dyn Read + Send>> {
    let file = File::open(filename)?;
    if filename.ends_with(".gz") {
        Ok(Box::new(flate2::read::GzDecoder::new(file)))
    } else if filename.ends_with(".bz2") {
        Ok(Box::new(bzip2::read::BzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

impl BatchSequenceReader {
    /// Open a file for reading. If filename is None, reads from stdin.
    /// Automatically decompresses `.gz` and `.bz2` files.
    pub fn new(filename: Option<&str>) -> io::Result<Self> {
        let reader: Box<dyn Read + Send> = match filename {
            Some(f) => open_possibly_compressed(f)?,
            None => Box::new(io::stdin()),
        };

        Ok(BatchSequenceReader {
            reader: BufReader::new(reader),
            seqs: Vec::new(),
            curr: 0,
            size: 0,
            file_format: SequenceFormat::AutoDetect,
            peeked_line: None,
            line_buf: String::new(),
        })
    }

    /// Create a reader from any `Read` source (e.g., in-memory bytes).
    ///
    /// # Example
    /// ```
    /// use kraken2_rs::seq::BatchSequenceReader;
    /// let data = b">seq1\nACGT\n>seq2\nTTTT\n";
    /// let mut reader = BatchSequenceReader::from_reader(std::io::Cursor::new(data));
    /// reader.load_block(1000);
    /// let seq = reader.next_sequence().unwrap();
    /// assert_eq!(seq.header, "seq1");
    /// ```
    pub fn from_reader<R: Read + Send + 'static>(reader: R) -> Self {
        BatchSequenceReader {
            reader: BufReader::new(Box::new(reader)),
            seqs: Vec::new(),
            curr: 0,
            size: 0,
            file_format: SequenceFormat::AutoDetect,
            peeked_line: None,
            line_buf: String::new(),
        }
    }

    /// Load a block of sequences totaling at least `block_size` bases.
    /// Returns true if at least one sequence was loaded.
    pub fn load_block(&mut self, block_size: usize) -> bool {
        let mut total = 0;
        self.size = 0;
        self.curr = 0;

        while total < block_size {
            if self.size >= self.seqs.len() {
                self.seqs.push(Sequence::default());
            }
            // Take ownership of the slot temporarily so we can pass &mut Sequence
            // without simultaneously borrowing &mut self.
            let mut slot = std::mem::take(&mut self.seqs[self.size]);
            let ok = self.read_one_into(&mut slot);
            let seq_len = slot.seq.len();
            self.seqs[self.size] = slot;
            if !ok {
                break;
            }
            self.size += 1;
            total += seq_len;
        }

        self.size > 0
    }

    /// Load a batch of exactly `record_count` sequences.
    /// Returns true if at least one sequence was loaded.
    pub fn load_batch(&mut self, record_count: usize) -> bool {
        if self.seqs.len() < record_count {
            self.seqs.resize_with(record_count, Sequence::default);
        }
        self.curr = 0;
        self.size = 0;

        for i in 0..record_count {
            let mut slot = std::mem::take(&mut self.seqs[i]);
            let ok = self.read_one_into(&mut slot);
            self.seqs[i] = slot;
            if !ok {
                break;
            }
            self.size += 1;
        }

        self.size > 0
    }

    /// Get the next sequence from the current batch.
    pub fn next_sequence(&mut self) -> Option<&Sequence> {
        if self.size > 0 && self.curr < self.size {
            let idx = self.curr;
            self.curr += 1;
            Some(&self.seqs[idx])
        } else {
            None
        }
    }

    pub fn file_format(&self) -> SequenceFormat {
        self.file_format
    }

    /// Read one sequence (FASTA or FASTQ) from the stream.
    /// Mimics kseq.h behavior: header is text after @/> up to whitespace,
    /// comment is the rest of the header line.
    #[allow(dead_code)]
    fn read_one(&mut self) -> Option<Sequence> {
        let mut s = Sequence::default();
        if self.read_one_into(&mut s) {
            Some(s)
        } else {
            None
        }
    }

    /// Refill an existing `Sequence` from the stream, reusing its String
    /// buffers. Returns false on EOF or read error. Hot path during classify.
    fn read_one_into(&mut self, out: &mut Sequence) -> bool {
        out.header.clear();
        out.comment.clear();
        out.seq.clear();
        out.quals.clear();

        // Get the header line (starts with > or @). If we previously peeked
        // a header line, swap it into line_buf so we can parse it without
        // allocating.
        if let Some(peeked) = self.peeked_line.take() {
            self.line_buf = peeked;
        } else {
            loop {
                self.line_buf.clear();
                match self.reader.read_line(&mut self.line_buf) {
                    Ok(0) => return false,
                    Ok(_) => {
                        trim_eol_in_place(&mut self.line_buf);
                        if self.line_buf.starts_with('>') || self.line_buf.starts_with('@') {
                            break;
                        }
                        // Skip blank / non-header lines before the first header
                    }
                    Err(_) => return false,
                }
            }
        }

        let is_fastq = self.line_buf.starts_with('@');
        let format = if is_fastq {
            SequenceFormat::Fastq
        } else {
            SequenceFormat::Fasta
        };
        self.file_format = format;
        out.format = format;

        // Parse header: first word is name, rest is comment.
        // line_buf contains "@name comment..." or ">name comment..."
        let content = &self.line_buf[1..];
        match content.bytes().position(|b| b == b' ' || b == b'\t') {
            Some(pos) => {
                out.header.push_str(&content[..pos]);
                out.comment.push_str(content[pos..].trim_start());
            }
            None => {
                out.header.push_str(content);
            }
        }

        // Read sequence lines
        loop {
            self.line_buf.clear();
            match self.reader.read_line(&mut self.line_buf) {
                Ok(0) => break,
                Ok(_) => {
                    trim_eol_in_place(&mut self.line_buf);
                    let first = self.line_buf.as_bytes().first().copied();
                    if first == Some(b'>') || first == Some(b'@') || first == Some(b'+') {
                        if is_fastq && first == Some(b'+') {
                            // FASTQ quality separator line
                            break;
                        }
                        if first == Some(b'>') || first == Some(b'@') {
                            // Next record's header — stash for the next call.
                            // Move out of line_buf to avoid copy.
                            let stash = std::mem::take(&mut self.line_buf);
                            self.peeked_line = Some(stash);
                            break;
                        }
                    }
                    out.seq.push_str(&self.line_buf);
                }
                Err(_) => break,
            }
        }

        // Read quality lines (FASTQ only)
        if is_fastq {
            while out.quals.len() < out.seq.len() {
                self.line_buf.clear();
                match self.reader.read_line(&mut self.line_buf) {
                    Ok(0) => break,
                    Ok(_) => {
                        trim_eol_in_place(&mut self.line_buf);
                        out.quals.push_str(&self.line_buf);
                    }
                    Err(_) => break,
                }
            }
        }

        true
    }
}

/// Trim trailing CR/LF from a String in place — avoids the slice allocation
/// that `trim_end_matches` would produce when the caller already needs to
/// mutate the String.
fn trim_eol_in_place(s: &mut String) {
    let bytes = s.as_bytes();
    let mut end = bytes.len();
    while end > 0 {
        let b = bytes[end - 1];
        if b == b'\n' || b == b'\r' {
            end -= 1;
        } else {
            break;
        }
    }
    s.truncate(end);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_read_fasta() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        write!(
            tmp.as_file(),
            ">seq1 some comment\nACGTACGT\nGGGG\n>seq2\nTTTT\n"
        )
        .unwrap();

        let mut reader = BatchSequenceReader::new(Some(tmp.path().to_str().unwrap())).unwrap();
        assert!(reader.load_block(1000));

        let s1 = reader.next_sequence().unwrap();
        assert_eq!(s1.header, "seq1");
        assert_eq!(s1.comment, "some comment");
        assert_eq!(s1.seq, "ACGTACGTGGGG");
        assert_eq!(s1.format, SequenceFormat::Fasta);

        let s2 = reader.next_sequence().unwrap();
        assert_eq!(s2.header, "seq2");
        assert_eq!(s2.seq, "TTTT");

        assert!(reader.next_sequence().is_none());
    }

    #[test]
    fn test_read_fastq() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        write!(
            tmp.as_file(),
            "@read1 desc\nACGT\n+\nIIII\n@read2\nTTTT\n+\nJJJJ\n"
        )
        .unwrap();

        let mut reader = BatchSequenceReader::new(Some(tmp.path().to_str().unwrap())).unwrap();
        assert!(reader.load_block(1000));

        let s1 = reader.next_sequence().unwrap();
        assert_eq!(s1.header, "read1");
        assert_eq!(s1.comment, "desc");
        assert_eq!(s1.seq, "ACGT");
        assert_eq!(s1.quals, "IIII");
        assert_eq!(s1.format, SequenceFormat::Fastq);

        let s2 = reader.next_sequence().unwrap();
        assert_eq!(s2.header, "read2");
        assert_eq!(s2.seq, "TTTT");
        assert_eq!(s2.quals, "JJJJ");
    }

    #[test]
    fn test_read_test_fasta() {
        let data_dir = format!("{}/kraken2/data", env!("CARGO_MANIFEST_DIR"));
        let covid = format!("{}/COVID_19.fa", data_dir);
        if !std::path::Path::new(&covid).exists() {
            return; // skip if test data not available
        }
        let mut reader = BatchSequenceReader::new(Some(&covid)).unwrap();
        assert!(reader.load_block(100_000));
        let seq = reader.next_sequence().unwrap();
        assert!(!seq.header.is_empty());
        assert!(seq.seq.len() > 29000); // COVID genome is ~29.9kb
        assert_eq!(seq.format, SequenceFormat::Fasta);
    }
}

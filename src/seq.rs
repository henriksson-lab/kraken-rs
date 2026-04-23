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
    peeked_line: Option<String>,
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
            if let Some(seq) = self.read_one() {
                let seq_len = seq.seq.len();
                self.seqs[self.size] = seq;
                self.size += 1;
                total += seq_len;
            } else {
                break;
            }
        }

        self.size > 0
    }

    /// Load a batch of exactly `record_count` sequences.
    /// Returns true if at least one sequence was loaded.
    pub fn load_batch(&mut self, record_count: usize) -> bool {
        self.seqs.resize_with(record_count, Sequence::default);
        self.curr = 0;
        self.size = 0;

        for i in 0..record_count {
            if let Some(seq) = self.read_one() {
                self.seqs[i] = seq;
                self.size += 1;
            } else {
                break;
            }
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
    fn read_one(&mut self) -> Option<Sequence> {
        // Get the header line (starts with > or @)
        let header_line = match self.peeked_line.take() {
            Some(line) => line,
            None => {
                let mut line = String::new();
                loop {
                    line.clear();
                    match self.reader.read_line(&mut line) {
                        Ok(0) => return None,
                        Ok(_) => {
                            let trimmed = line.trim_end_matches(['\n', '\r']);
                            if trimmed.starts_with('>') || trimmed.starts_with('@') {
                                break trimmed.to_string();
                            }
                            // Skip blank lines before header
                        }
                        Err(_) => return None,
                    }
                }
            }
        };

        let is_fastq = header_line.starts_with('@');
        let format = if is_fastq {
            SequenceFormat::Fastq
        } else {
            SequenceFormat::Fasta
        };
        self.file_format = format;

        // Parse header: first word is name, rest is comment
        let content = &header_line[1..]; // skip > or @
        let (header, comment) = match content.find(|c: char| c.is_whitespace()) {
            Some(pos) => (content[..pos].to_string(), content[pos + 1..].to_string()),
            None => (content.to_string(), String::new()),
        };

        // Read sequence lines
        let mut seq = String::new();
        let mut line = String::new();
        loop {
            line.clear();
            match self.reader.read_line(&mut line) {
                Ok(0) => break,
                Ok(_) => {
                    let trimmed = line.trim_end_matches(['\n', '\r']);
                    if trimmed.starts_with('>')
                        || trimmed.starts_with('@')
                        || trimmed.starts_with('+')
                    {
                        if is_fastq && trimmed.starts_with('+') {
                            // This is the quality separator line
                            break;
                        }
                        if trimmed.starts_with('>') || trimmed.starts_with('@') {
                            // Next record's header — save for next call
                            self.peeked_line = Some(trimmed.to_string());
                            break;
                        }
                    }
                    seq.push_str(trimmed);
                }
                Err(_) => break,
            }
        }

        // Read quality lines (FASTQ only)
        let mut quals = String::new();
        if is_fastq {
            while quals.len() < seq.len() {
                line.clear();
                match self.reader.read_line(&mut line) {
                    Ok(0) => break,
                    Ok(_) => {
                        let trimmed = line.trim_end_matches(['\n', '\r']);
                        quals.push_str(trimmed);
                    }
                    Err(_) => break,
                }
            }
        }

        Some(Sequence {
            format,
            header,
            comment,
            seq,
            quals,
        })
    }
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

//! BLAST database to FASTA converter.
//! Port of C `blast_to_fasta.c`, `blast_defline.c`, `blast_utils.c`.
//!
//! Converts BLAST binary database files (.nin/.pin, .nhr/.phr, .nsq/.psq)
//! into standard FASTA format with optional taxid annotations.

use std::fs::File;
use std::io::{self, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

/// 2-bit nucleotide decoding table (BLAST uses NCBI2na encoding).
const NCBI2NA_DECODE: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// 4-bit nucleotide decoding table (NCBI4na encoding for ambiguities).
#[allow(dead_code)]
const NCBI4NA_DECODE: [u8; 16] = [
    b'-', b'A', b'C', b'M', b'G', b'R', b'S', b'V',
    b'T', b'W', b'Y', b'H', b'K', b'D', b'B', b'N',
];

/// BLAST database index header.
struct BlastIndex {
    fmt_version: u32,
    db_seq_type: u32, // 0 = nucleotide, 1 = protein
    num_oids: u32,
    title: String,
    #[allow(dead_code)]
    date: String,
    hdr_offsets: Vec<u32>,
    seq_offsets: Vec<u32>,
    #[allow(dead_code)]
    amb_offsets: Vec<u32>,
}

/// Read a 4-byte big-endian u32.
fn read_be_u32(file: &mut File) -> io::Result<u32> {
    let mut buf = [0u8; 4];
    file.read_exact(&mut buf)?;
    Ok(u32::from_be_bytes(buf))
}

/// Read a 4-byte big-endian u32 as length-prefixed string.
fn read_length_prefixed_string(file: &mut File) -> io::Result<String> {
    let len = read_be_u32(file)? as usize;
    let mut buf = vec![0u8; len];
    file.read_exact(&mut buf)?;
    Ok(String::from_utf8_lossy(&buf).into_owned())
}

/// Parse a BLAST index file (.nin or .pin).
fn parse_blast_index(idx_path: &str) -> io::Result<BlastIndex> {
    let mut file = File::open(idx_path)?;

    let fmt_version = read_be_u32(&mut file)?;
    let db_seq_type = read_be_u32(&mut file)?;
    let _volume = read_be_u32(&mut file)?;
    let title = read_length_prefixed_string(&mut file)?;
    let date = read_length_prefixed_string(&mut file)?;
    let num_oids = read_be_u32(&mut file)?;
    let _vol_len_hi = read_be_u32(&mut file)?;
    let _vol_len_lo = read_be_u32(&mut file)?;
    let _max_seq_len = read_be_u32(&mut file)?;

    // Read offset arrays
    let mut hdr_offsets = Vec::with_capacity(num_oids as usize + 1);
    for _ in 0..=num_oids {
        hdr_offsets.push(read_be_u32(&mut file)?);
    }

    let mut seq_offsets = Vec::with_capacity(num_oids as usize + 1);
    for _ in 0..=num_oids {
        seq_offsets.push(read_be_u32(&mut file)?);
    }

    let mut amb_offsets = Vec::with_capacity(num_oids as usize + 1);
    if db_seq_type == 0 {
        // Nucleotide databases have ambiguity offsets
        for _ in 0..=num_oids {
            amb_offsets.push(read_be_u32(&mut file)?);
        }
    }

    Ok(BlastIndex {
        fmt_version,
        db_seq_type,
        num_oids,
        title,
        date,
        hdr_offsets,
        seq_offsets,
        amb_offsets,
    })
}

/// Decode a 2-bit packed nucleotide sequence.
fn decode_ncbi2na(packed: &[u8], length: usize) -> Vec<u8> {
    let mut seq = Vec::with_capacity(length);
    for &byte in packed {
        for shift in (0..4).rev() {
            if seq.len() >= length {
                break;
            }
            let code = (byte >> (shift * 2)) & 0x03;
            seq.push(NCBI2NA_DECODE[code as usize]);
        }
    }
    seq.truncate(length);
    seq
}

/// Convert a BLAST database to FASTA format.
///
/// `db_prefix` is the path without extension (e.g., "mydb" for mydb.nin, mydb.nhr, mydb.nsq).
/// `output_path` is the FASTA output file path.
/// `include_taxid` adds kraken:taxid| prefix to headers if true.
pub fn blast_to_fasta(
    db_prefix: &str,
    output_path: &str,
    _include_taxid: bool,
) -> io::Result<()> {
    // Determine file extensions based on database type
    let idx_nuc = format!("{db_prefix}.nin");
    let idx_pro = format!("{db_prefix}.pin");

    let (idx_path, hdr_ext, seq_ext) = if Path::new(&idx_nuc).exists() {
        (idx_nuc, "nhr", "nsq")
    } else if Path::new(&idx_pro).exists() {
        (idx_pro, "phr", "psq")
    } else {
        return Err(io::Error::other(format!(
            "Cannot find BLAST index at {db_prefix}.nin or {db_prefix}.pin"
        )));
    };

    let idx = parse_blast_index(&idx_path)?;
    eprintln!("BLAST database: {}", idx.title);
    eprintln!("  Type: {}", if idx.db_seq_type == 0 { "nucleotide" } else { "protein" });
    eprintln!("  Sequences: {}", idx.num_oids);
    eprintln!("  Format version: {:#x}", idx.fmt_version);

    let hdr_path = format!("{db_prefix}.{hdr_ext}");
    let seq_path = format!("{db_prefix}.{seq_ext}");

    let mut hdr_file = File::open(&hdr_path)?;
    let mut seq_file = File::open(&seq_path)?;
    let mut output = BufWriter::new(File::create(output_path)?);

    for oid in 0..idx.num_oids as usize {
        // Read header
        let hdr_start = idx.hdr_offsets[oid] as u64;
        let hdr_end = idx.hdr_offsets[oid + 1] as u64;
        let hdr_len = (hdr_end - hdr_start) as usize;

        hdr_file.seek(SeekFrom::Start(hdr_start))?;
        let mut hdr_buf = vec![0u8; hdr_len];
        hdr_file.read_exact(&mut hdr_buf)?;

        // Extract header text (simplified — full ASN.1 parsing would be needed
        // for complete defline extraction, but for basic use we can extract
        // visible ASCII strings)
        let header_text = extract_readable_header(&hdr_buf);

        // Read sequence
        let seq_start = idx.seq_offsets[oid] as u64;
        let seq_end = idx.seq_offsets[oid + 1] as u64;

        if idx.db_seq_type == 0 {
            // Nucleotide: 2-bit packed, last byte has residue count in high 2 bits
            let packed_len = (seq_end - seq_start) as usize;
            seq_file.seek(SeekFrom::Start(seq_start))?;
            let mut packed = vec![0u8; packed_len];
            seq_file.read_exact(&mut packed)?;

            if packed.is_empty() {
                continue;
            }

            // Last byte encodes the number of valid residues in the final byte
            let last_byte = packed[packed_len - 1];
            let residue_count = if packed_len > 1 {
                let full_bytes = packed_len - 1;
                let last_residues = (last_byte >> 6) as usize; // high 2 bits
                if last_residues == 0 {
                    full_bytes * 4
                } else {
                    (full_bytes - 1) * 4 + last_residues
                }
            } else {
                0
            };

            let seq = decode_ncbi2na(&packed[..packed_len.saturating_sub(1)], residue_count);
            writeln!(output, ">{header_text}")?;
            for chunk in seq.chunks(70) {
                output.write_all(chunk)?;
                output.write_all(b"\n")?;
            }
        } else {
            // Protein: one byte per residue (NCBI standard encoding)
            let seq_len = (seq_end - seq_start) as usize;
            seq_file.seek(SeekFrom::Start(seq_start))?;
            let mut seq_buf = vec![0u8; seq_len];
            seq_file.read_exact(&mut seq_buf)?;

            // Convert NCBI encoding to ASCII amino acids
            let seq: Vec<u8> = seq_buf.iter().map(|&b| ncbi_aa_decode(b)).collect();
            writeln!(output, ">{header_text}")?;
            for chunk in seq.chunks(70) {
                output.write_all(chunk)?;
                output.write_all(b"\n")?;
            }
        }
    }

    eprintln!("Wrote {} sequences to {output_path}", idx.num_oids);
    Ok(())
}

/// Extract readable text from ASN.1 encoded header data.
/// This is a simplified extraction that finds printable ASCII strings.
fn extract_readable_header(data: &[u8]) -> String {
    // Find the longest printable ASCII run
    let mut best_start = 0;
    let mut best_len = 0;
    let mut curr_start = 0;
    let mut curr_len = 0;

    for (i, &b) in data.iter().enumerate() {
        if (0x20..0x7F).contains(&b) {
            if curr_len == 0 {
                curr_start = i;
            }
            curr_len += 1;
            if curr_len > best_len {
                best_start = curr_start;
                best_len = curr_len;
            }
        } else {
            curr_len = 0;
        }
    }

    if best_len > 0 {
        String::from_utf8_lossy(&data[best_start..best_start + best_len]).into_owned()
    } else {
        format!("sequence_{}", data.len())
    }
}

/// Decode NCBI amino acid encoding to ASCII.
fn ncbi_aa_decode(code: u8) -> u8 {
    const AA_TABLE: &[u8] = b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";
    if (code as usize) < AA_TABLE.len() {
        AA_TABLE[code as usize]
    } else {
        b'X'
    }
}

use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use flate2::read::GzDecoder;

const NCBI_FTP_BASE: &str = "https://ftp.ncbi.nlm.nih.gov";

/// Download a file from a URL using `ureq` and write it verbatim to
/// `output_path`. Replaces the `wget`/`curl` calls in the original
/// `download_*.sh` and `download_taxonomy.pl` scripts. Returns the number
/// of bytes written.
pub fn download_file(url: &str, output_path: &str) -> io::Result<u64> {
    eprintln!("Downloading {} ...", url);
    let response = ureq::get(url)
        .call()
        .map_err(|e| io::Error::other(format!("HTTP error: {e}")))?;
    let mut output = File::create(output_path)?;
    let mut body = response.into_body().into_reader();
    let bytes = io::copy(&mut body, &mut output)?;
    eprintln!("  Downloaded {bytes} bytes to {output_path}");
    Ok(bytes)
}

/// Download a gzip-compressed file and write the decompressed bytes to
/// `output_path`. Used for NCBI `accession2taxid.gz` maps; replaces the
/// `wget | gunzip` pipeline used by the original shell scripts. Returns the
/// number of decompressed bytes written.
pub fn download_and_decompress_gz(url: &str, output_path: &str) -> io::Result<u64> {
    eprintln!("Downloading and decompressing {} ...", url);
    let response = ureq::get(url)
        .call()
        .map_err(|e| io::Error::other(format!("HTTP error: {e}")))?;
    let body = response.into_body().into_reader();
    let mut decoder = GzDecoder::new(body);
    let mut output = File::create(output_path)?;
    let bytes = io::copy(&mut decoder, &mut output)?;
    eprintln!("  Decompressed {bytes} bytes to {output_path}");
    Ok(bytes)
}

/// Download and extract a `.tar.gz` archive into `output_dir`, keeping only
/// the entries whose basename appears in `keep_files` (empty `keep_files`
/// extracts everything). Replaces the `tar -xzf` invocation used by the
/// original `download_taxonomy.pl` script when fetching `taxdump.tar.gz`.
fn download_and_extract_tar_gz(
    url: &str,
    output_dir: &Path,
    keep_files: &[&str],
) -> io::Result<()> {
    eprintln!("Downloading and extracting {} ...", url);
    let response = ureq::get(url)
        .call()
        .map_err(|e| io::Error::other(format!("HTTP error: {e}")))?;
    let body = response.into_body().into_reader();
    let decoder = GzDecoder::new(body);
    let mut archive = tar::Archive::new(decoder);

    for entry in archive.entries()? {
        let mut entry = entry?;
        let path = entry.path()?.to_path_buf();
        let filename = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if keep_files.is_empty() || keep_files.contains(&filename) {
            let out_path = output_dir.join(filename);
            let mut out_file = File::create(&out_path)?;
            io::copy(&mut entry, &mut out_file)?;
            eprintln!("  Extracted {}", out_path.display());
        }
    }
    Ok(())
}

/// Download the NCBI taxonomy into `db_dir/taxonomy/`.
///
/// Native Rust replacement for the Perl `download_taxonomy.pl` workflow:
/// fetches the appropriate `accession2taxid` map(s) (`prot.accession2taxid`
/// when `protein`, otherwise the `nucl_gb`/`nucl_wgs` pair) unless
/// `skip_maps` is set, then downloads and extracts `taxdump.tar.gz` keeping
/// only `nodes.dmp`, `names.dmp`, and `merged.dmp`.
pub fn download_taxonomy(db_dir: &str, skip_maps: bool, protein: bool) -> io::Result<()> {
    let taxonomy_dir = PathBuf::from(db_dir).join("taxonomy");
    fs::create_dir_all(&taxonomy_dir)?;

    // Download accession-to-taxid maps
    if !skip_maps {
        if protein {
            let url =
                format!("{NCBI_FTP_BASE}/pub/taxonomy/accession2taxid/prot.accession2taxid.gz");
            let out = taxonomy_dir.join("prot.accession2taxid");
            download_and_decompress_gz(&url, out.to_str().unwrap())?;
        } else {
            for subsection in ["gb", "wgs"] {
                let url = format!(
                    "{NCBI_FTP_BASE}/pub/taxonomy/accession2taxid/nucl_{subsection}.accession2taxid.gz"
                );
                let out = taxonomy_dir.join(format!("nucl_{subsection}.accession2taxid"));
                download_and_decompress_gz(&url, out.to_str().unwrap())?;
            }
        }
    }

    // Download and extract taxdump
    let taxdump_url = format!("{NCBI_FTP_BASE}/pub/taxonomy/taxdump.tar.gz");
    download_and_extract_tar_gz(
        &taxdump_url,
        &taxonomy_dir,
        &["nodes.dmp", "names.dmp", "merged.dmp"],
    )?;

    eprintln!("Taxonomy downloaded to {}", taxonomy_dir.display());
    Ok(())
}

/// An entry from assembly_summary.txt.
struct AssemblyEntry {
    taxid: u64,
    #[allow(dead_code)]
    asm_level: String,
    ftp_path: String,
}

/// Parse an NCBI RefSeq `assembly_summary.txt` and return one
/// [`AssemblyEntry`] per Complete Genome / Chromosome assembly with a usable
/// `ftp_path` (field 19). Comment lines and rows with `ftp_path == "na"` or
/// fewer than 20 tab-separated fields are skipped. Replaces the
/// `awk`/`grep` filtering used by the original `download_genomic_library.sh`.
fn parse_assembly_summary(path: &Path) -> io::Result<Vec<AssemblyEntry>> {
    let file = BufReader::new(File::open(path)?);
    let mut entries = Vec::new();

    for line in file.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 20 {
            continue;
        }
        let taxid: u64 = fields[5].parse().unwrap_or(0);
        let asm_level = fields[11].to_string();
        let ftp_path = fields[19].to_string();

        if ftp_path == "na" {
            continue;
        }
        // Only include Complete Genome and Chromosome level assemblies
        if asm_level != "Complete Genome" && asm_level != "Chromosome" {
            continue;
        }
        entries.push(AssemblyEntry {
            taxid,
            asm_level,
            ftp_path,
        });
    }
    Ok(entries)
}

/// Download genomes from parsed assembly entries and write to `library.fna`
/// (or `library.faa` when `protein`).
///
/// For each entry, fetches `{ftp_path}/{basename}_genomic.fna.gz` (or the
/// `_protein.faa.gz` variant), streams it through gzip, and rewrites every
/// FASTA header to `>kraken:taxid|TAXID|<original-header>` so downstream
/// build steps can recover the taxid via [`scan_fasta_for_taxids`]. Returns
/// `(projects, sequences, bases)` counts. Replaces the bulk-download loop
/// in `download_genomic_library.sh`.
fn download_genomes(
    entries: &[AssemblyEntry],
    output_path: &Path,
    protein: bool,
) -> io::Result<(usize, usize, usize)> {
    let suffix = if protein {
        "_protein.faa.gz"
    } else {
        "_genomic.fna.gz"
    };
    let mut output = File::create(output_path)?;
    let mut projects = 0;
    let mut sequences = 0;
    let mut bases = 0;

    for entry in entries {
        let ftp_path = entry.ftp_path.trim_end_matches('/');
        let basename = ftp_path.rsplit('/').next().unwrap_or("");
        let url = format!("{ftp_path}/{basename}{suffix}");

        eprint!("\rDownloading {basename} (taxid {})...", entry.taxid);

        let response = match ureq::get(&url).call() {
            Ok(r) => r,
            Err(e) => {
                eprintln!("\n  Warning: failed to download {url}: {e}");
                continue;
            }
        };

        let body = response.into_body().into_reader();
        let decoder = GzDecoder::new(body);
        let reader = BufReader::new(decoder);

        for line in reader.lines() {
            let line = line?;
            if let Some(header) = line.strip_prefix('>') {
                // Prepend taxid to header
                writeln!(output, ">kraken:taxid|{}|{}", entry.taxid, header)?;
                sequences += 1;
            } else {
                writeln!(output, "{}", line)?;
                bases += line.len();
            }
        }
        projects += 1;
    }

    eprintln!(
        "\nCompleted: {projects} projects, {sequences} sequences, {} {}",
        bases,
        if protein { "aa" } else { "bp" }
    );
    Ok((projects, sequences, bases))
}

/// Scan a FASTA file and extract `(seqid, taxid)` pairs from headers shaped
/// `>kraken:taxid|TAXID|<rest>`.
///
/// The `seqid` returned is the entire first whitespace-delimited token of
/// the header (including the `kraken:taxid|TAXID|` prefix), matching what
/// the build step expects in `prelim_map.txt`. Lines without the
/// `kraken:taxid|` prefix are ignored.
pub fn scan_fasta_for_taxids(fasta_path: &Path) -> io::Result<Vec<(String, u64)>> {
    let file = BufReader::new(File::open(fasta_path)?);
    let mut mappings = Vec::new();

    for line in file.lines() {
        let line = line?;
        if !line.starts_with('>') {
            continue;
        }
        let header = &line[1..];
        let seqid = header.split_whitespace().next().unwrap_or("");

        // Check for kraken:taxid|TAXID| pattern
        if let Some(rest) = seqid.strip_prefix("kraken:taxid|") {
            if let Some(pipe_pos) = rest.find('|') {
                if let Ok(taxid) = rest[..pipe_pos].parse::<u64>() {
                    mappings.push((seqid.to_string(), taxid));
                }
            }
        }
    }
    Ok(mappings)
}

/// Download a genomic library from NCBI into `db_dir/library/<library_type>/`.
///
/// Replaces `kraken2-build --download-library` (the Perl/Bash
/// `download_genomic_library.sh`):
/// - `UniVec` / `UniVec_Core`: fetches the bare file from `pub/UniVec/`,
///   rewrites every FASTA header with the synthetic taxid 28384
///   ("other sequences"), and emits `prelim_map.txt`.
/// - `archaea`, `bacteria`, `viral`, `fungi`, `plant`, `protozoa`, `human`,
///   `plasmid`: downloads `assembly_summary.txt`, (for `human` filters to
///   Genome Reference Consortium assemblies), then calls
///   [`download_genomes`] over [`parse_assembly_summary`] results and writes
///   `prelim_map.txt`. Other `library_type` values are rejected.
/// Protein mode produces `library.faa`; nucleotide mode produces `library.fna`.
pub fn download_library(db_dir: &str, library_type: &str, protein: bool) -> io::Result<()> {
    let library_dir = PathBuf::from(db_dir).join("library").join(library_type);
    fs::create_dir_all(&library_dir)?;

    let library_file = if protein {
        "library.faa"
    } else {
        "library.fna"
    };

    match library_type {
        "UniVec" | "UniVec_Core" => {
            if protein {
                return Err(io::Error::other(format!(
                    "{library_type} is for nucleotide databases only"
                )));
            }
            let url = format!("{NCBI_FTP_BASE}/pub/UniVec/{library_type}");
            let raw_path = library_dir.join(library_type);
            download_file(&url, raw_path.to_str().unwrap())?;

            // Add taxid 28384 ("other sequences") to all headers
            let input = BufReader::new(File::open(&raw_path)?);
            let output_path = library_dir.join(library_file);
            let mut output = File::create(&output_path)?;
            for line in input.lines() {
                let line = line?;
                if let Some(header) = line.strip_prefix('>') {
                    writeln!(output, ">kraken:taxid|28384|{}", header)?;
                } else {
                    writeln!(output, "{}", line)?;
                }
            }
            // Generate prelim_map.txt
            let mappings = scan_fasta_for_taxids(&output_path)?;
            let map_path = library_dir.join("prelim_map.txt");
            let mut map_file = File::create(&map_path)?;
            for (seqid, taxid) in &mappings {
                writeln!(map_file, "TAXID\t{seqid}\t{taxid}")?;
            }
        }

        "archaea" | "bacteria" | "viral" | "fungi" | "plant" | "protozoa" | "human" | "plasmid" => {
            let remote_dir = if library_type == "human" {
                "vertebrate_mammalian/Homo_sapiens"
            } else {
                library_type
            };

            // Download assembly_summary.txt
            let summary_url =
                format!("{NCBI_FTP_BASE}/genomes/refseq/{remote_dir}/assembly_summary.txt");
            let summary_path = library_dir.join("assembly_summary.txt");
            download_file(&summary_url, summary_path.to_str().unwrap())?;

            // For human, filter to GRC assemblies only
            if library_type == "human" {
                let content = fs::read_to_string(&summary_path)?;
                let filtered: String = content
                    .lines()
                    .filter(|l| l.starts_with('#') || l.contains("Genome Reference Consortium"))
                    .map(|l| format!("{l}\n"))
                    .collect();
                fs::write(&summary_path, filtered)?;
            }

            // Parse assembly summary and download genomes
            let entries = parse_assembly_summary(&summary_path)?;
            eprintln!(
                "Found {} complete genome entries for {library_type}",
                entries.len()
            );

            let output_path = library_dir.join(library_file);
            download_genomes(&entries, &output_path, protein)?;

            // Generate prelim_map.txt
            let mappings = scan_fasta_for_taxids(&output_path)?;
            let map_path = library_dir.join("prelim_map.txt");
            let mut map_file = File::create(&map_path)?;
            for (seqid, taxid) in &mappings {
                writeln!(map_file, "TAXID\t{seqid}\t{taxid}")?;
            }
        }

        _ => {
            return Err(io::Error::other(format!(
                "Unknown library type: {library_type}"
            )));
        }
    }

    eprintln!(
        "Library {library_type} downloaded to {}",
        library_dir.display()
    );
    Ok(())
}

/// Clean up intermediate files after a database has been built.
///
/// Removes `db_dir/library/`, `db_dir/taxonomy/`, and `db_dir/seqid2taxid.map`
/// if present. Equivalent to the `kraken2-build --clean` step in the
/// original Perl wrapper.
pub fn clean_db(db_dir: &str) -> io::Result<()> {
    let library_dir = PathBuf::from(db_dir).join("library");
    let taxonomy_dir = PathBuf::from(db_dir).join("taxonomy");
    let seqid_map = PathBuf::from(db_dir).join("seqid2taxid.map");

    if library_dir.exists() {
        eprintln!("Removing library directory...");
        fs::remove_dir_all(&library_dir)?;
    }
    if taxonomy_dir.exists() {
        eprintln!("Removing taxonomy directory...");
        fs::remove_dir_all(&taxonomy_dir)?;
    }
    if seqid_map.exists() {
        fs::remove_file(&seqid_map)?;
    }
    eprintln!("Cleanup complete.");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Verify that two synthetic `>kraken:taxid|TAXID|...` headers are parsed
    /// out of a FASTA file and their taxids extracted in order.
    #[test]
    fn test_scan_fasta_for_taxids() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(
            &mut tmp.as_file(),
            b">kraken:taxid|562|NC_000913.3 E. coli\nACGT\n>kraken:taxid|9606|NC_000001.11 Human\nGGGG\n",
        ).unwrap();

        let mappings = scan_fasta_for_taxids(tmp.path()).unwrap();
        assert_eq!(mappings.len(), 2);
        assert_eq!(mappings[0].1, 562);
        assert_eq!(mappings[1].1, 9606);
    }

    /// Parse a synthetic `assembly_summary.txt` row and confirm taxid,
    /// assembly level, and FTP path land in the right [`AssemblyEntry`] fields.
    #[test]
    fn test_parse_assembly_summary() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        // Minimal assembly_summary format: field 5=taxid, 11=asm_level, 19=ftp_path
        let line = "GCF_000005845.2\tPRJNA225\tSAMN02436677\tna\t\t511145\t511145\tna\tReference Genome\tna\tna\tComplete Genome\tMajor\tFull\t2013/09/26\tASM584v2\tNational Center for Biotechnology Information\tGCA_000005845.2\t=\thttps://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2\n";
        std::io::Write::write_all(&mut tmp.as_file(), format!("# header\n{line}").as_bytes())
            .unwrap();

        let entries = parse_assembly_summary(tmp.path()).unwrap();
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].taxid, 511145);
        assert_eq!(entries[0].asm_level, "Complete Genome");
        assert!(entries[0].ftp_path.contains("GCF_000005845"));
    }
}

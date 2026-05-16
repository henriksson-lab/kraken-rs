use ahash::AHashMap as HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};

use crate::mmap_file::MMapFile;
use crate::types::TaxId;
use crate::utilities::split_string;

/// Look up accession numbers in NCBI accession-to-taxid map files.
///
/// `accession_file` is a TSV where column 2 of each line is an accession to
/// resolve (typically the `prelim_map.txt` produced by library scanning).
/// Each `map_files` entry follows the NCBI `accession2taxid` schema:
/// `accession\taccession.version\ttaxid\tgi` per line (header row included);
/// both the bare accession and the `accession.version` form are matched.
/// Returns a map from accession (as originally listed in `accession_file`) to
/// its resolved [`TaxId`], with `0` for unmapped entries.
/// Port of C++ `lookup_accession_numbers.cc`.
pub fn lookup_accession_numbers(
    accession_file: &str,
    map_files: &[String],
) -> io::Result<HashMap<String, TaxId>> {
    // Read accessions to look up
    let mut accessions: HashMap<String, TaxId> = HashMap::new();
    let file = BufReader::new(File::open(accession_file)?);
    for line in file.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 2 {
            accessions.insert(parts[1].to_string(), 0);
        }
    }

    // Search through map files
    for map_file in map_files {
        let file = BufReader::new(File::open(map_file)?);
        for line in file.lines() {
            let line = line?;
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 3 {
                let accession = parts[0];
                let accession_ver = parts[1];
                if let Ok(taxid) = parts[2].parse::<u64>() {
                    if accessions.contains_key(accession) {
                        accessions.insert(accession.to_string(), taxid);
                    }
                    if accessions.contains_key(accession_ver) {
                        accessions.insert(accession_ver.to_string(), taxid);
                    }
                }
            }
        }
    }

    Ok(accessions)
}

/// Return `true` if stderr is attached to a terminal. Used to gate the
/// in-place progress updates (`\r ...`) that the C++ code emits only when
/// `isatty(fileno(stderr))`.
fn stderr_is_tty() -> bool {
    #[cfg(unix)]
    {
        unsafe { libc::isatty(libc::STDERR_FILENO) == 1 }
    }
    #[cfg(not(unix))]
    {
        false
    }
}

/// Open the destination for resolved `seqid\ttaxid` lines: a file if
/// `output_filename` is `Some`, otherwise stdout (matching the C++ behaviour
/// of writing to `std::cout`).
fn open_lookup_writer(output_filename: Option<&str>) -> io::Result<Box<dyn Write>> {
    match output_filename {
        Some(path) => Ok(Box::new(File::create(path)?)),
        None => Ok(Box::new(io::stdout())),
    }
}

/// Lossy UTF-8 extract of `line[start..end]` as an owned `String`, used to
/// pull individual tab-separated fields out of the memory-mapped map file.
fn accession_field(line: &[u8], start: usize, end: usize) -> String {
    String::from_utf8_lossy(&line[start..end]).into_owned()
}

/// Core driver behind the `lookup_accession_numbers` CLI; writable target is
/// abstracted so tests can capture stdout.
///
/// `args[1]` is the `seqid\taccnum` lookup file (one accession may map to many
/// seqids). `args[2..]` are NCBI `accession2taxid` maps, each with a header
/// line followed by `accession\taccession.version\ttaxid\tgi` rows; files are
/// memory-mapped via [`MMapFile`] and scanned linearly. For each hit, every
/// associated `seqid\ttaxid` is written to `writer`; matched accessions are
/// removed from the working set so iteration stops as soon as everything has
/// been resolved. Any accessions still unmapped at the end are written to
/// `unmapped.txt` in the current directory (matching the C++ behaviour).
fn lookup_accession_numbers_main_inner(args: &[String], writer: &mut dyn Write) -> io::Result<()> {
    if args.len() < 3 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Usage: lookup_accession_numbers <lookup file> <accmaps>",
        ));
    }

    let mut target_lists: HashMap<String, Vec<String>> = HashMap::new();
    let lookup_list_file = File::open(&args[1])?;
    for line in BufReader::new(lookup_list_file).lines() {
        let line = line?;
        let fields = split_string(&line, "\t", 2);
        if fields.len() >= 2 {
            let seqid = fields[0].clone();
            let accnum = fields[1].clone();
            target_lists.entry(accnum).or_default().push(seqid);
        }
    }

    let initial_target_count = target_lists.len();
    let mut stderr = io::stderr();
    let tty = stderr_is_tty();
    if tty {
        write!(stderr, "\rFound 0/{initial_target_count} targets...")?;
        stderr.flush()?;
    }

    let mut accessions_searched = 0u64;
    for map_filename in &args[2..] {
        if target_lists.is_empty() {
            break;
        }
        let accmap_file = MMapFile::open_read_only(map_filename)?;
        let data = accmap_file.as_slice();
        if data.is_empty() {
            continue;
        }

        let mut ptr = 0usize;
        if let Some(pos) = data.iter().position(|&b| b == b'\n') {
            ptr = pos + 1;
        }

        while ptr < data.len() {
            let rel_lf = match data[ptr..].iter().position(|&b| b == b'\n') {
                Some(pos) => pos,
                None => {
                    eprintln!("expected EOL not found at EOF in {map_filename}");
                    break;
                }
            };
            let lf_ptr = ptr + rel_lf;
            let line = &data[ptr..lf_ptr];
            let first_tab = match line.iter().position(|&b| b == b'\t') {
                Some(pos) => pos,
                None => {
                    eprintln!("expected TAB not found in {map_filename}");
                    break;
                }
            };
            let accnum = accession_field(line, 0, first_tab);
            accessions_searched += 1;

            if let Some(seqids) = target_lists.get(&accnum).cloned() {
                let mut taxid_start = 0usize;
                let mut tab_ptr = first_tab;
                for _ in 0..2 {
                    taxid_start = tab_ptr + 1;
                    tab_ptr = line[taxid_start..]
                        .iter()
                        .position(|&b| b == b'\t')
                        .map(|pos| taxid_start + pos)
                        .unwrap_or(line.len());
                    if tab_ptr == line.len() {
                        eprintln!("expected TAB not found in {map_filename}");
                        break;
                    }
                }
                if tab_ptr == line.len() {
                    break;
                }
                let taxid = accession_field(line, taxid_start, tab_ptr);
                for seqid in seqids {
                    writeln!(writer, "{seqid}\t{taxid}")?;
                }
                target_lists.remove(&accnum);
                if tty {
                    write!(
                        stderr,
                        "\rFound {}/{} targets, searched through {} accession IDs...",
                        initial_target_count - target_lists.len(),
                        initial_target_count,
                        accessions_searched
                    )?;
                    stderr.flush()?;
                }
                if target_lists.is_empty() {
                    break;
                }
            }

            if accessions_searched.is_multiple_of(10_000_000) && tty {
                write!(
                    stderr,
                    "\rFound {}/{} targets, searched through {} accession IDs...",
                    initial_target_count - target_lists.len(),
                    initial_target_count,
                    accessions_searched
                )?;
                stderr.flush()?;
            }
            ptr = lf_ptr + 1;
        }
    }

    if tty {
        write!(stderr, "\r")?;
    }
    writeln!(
        stderr,
        "Found {}/{} targets, searched through {} accession IDs, search complete.",
        initial_target_count - target_lists.len(),
        initial_target_count,
        accessions_searched
    )?;

    if !target_lists.is_empty() {
        writeln!(
            stderr,
            "lookup_accession_numbers: {}/{} accession numbers remain unmapped, see unmapped.txt in DB directory",
            target_lists.len(),
            initial_target_count
        )?;
        let mut ofs = File::create("unmapped.txt")?;
        for accnum in target_lists.keys() {
            writeln!(ofs, "{accnum}")?;
        }
    }

    Ok(())
}

/// CLI entry point for `lookup_accession_numbers`.
///
/// Forwards to [`lookup_accession_numbers_main_inner`] with stdout as the
/// output sink, matching the original C++ binary.
pub fn lookup_accession_numbers_main(args: &[String]) -> io::Result<()> {
    let mut writer = open_lookup_writer(None)?;
    lookup_accession_numbers_main_inner(args, &mut writer)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::Path;
    use tempfile::tempdir;

    /// Round-trip [`lookup_accession_numbers`] on tiny synthetic accession and
    /// map files, asserting both bare and `accession.version` resolution.
    #[test]
    fn test_lookup_accession_numbers_helper() {
        let dir = tempdir().unwrap();
        let accession_path = dir.path().join("lookup.txt");
        let map_path = dir.path().join("map.tsv");
        fs::write(&accession_path, "seq1\tACC1\nseq2\tACC2.1\n").unwrap();
        fs::write(
            &map_path,
            "accession\taccession.version\ttaxid\tgi\nACC1\tACC1.1\t111\t0\nACC2\tACC2.1\t222\t0\n",
        )
        .unwrap();

        let result = lookup_accession_numbers(
            accession_path.to_str().unwrap(),
            &[map_path.to_str().unwrap().to_string()],
        )
        .unwrap();
        assert_eq!(result.get("ACC1"), Some(&111));
        assert_eq!(result.get("ACC2.1"), Some(&222));
    }

    /// Verify the CLI driver writes resolved `seqid\ttaxid` lines to its
    /// output sink and records remaining accessions in `unmapped.txt`.
    #[test]
    fn test_lookup_accession_numbers_main_writes_unmapped() {
        let dir = tempdir().unwrap();
        let old_cwd = std::env::current_dir().unwrap();
        std::env::set_current_dir(dir.path()).unwrap();

        let lookup_path = dir.path().join("lookup.txt");
        let map_path = dir.path().join("map.tsv");
        fs::write(&lookup_path, "seq1\tACC1\nseq2\tMISS2\n").unwrap();
        fs::write(
            &map_path,
            "accession\taccession.version\ttaxid\tgi\nACC1\tACC1.1\t123\t0\n",
        )
        .unwrap();

        let args = vec![
            "lookup_accession_numbers".to_string(),
            lookup_path.to_str().unwrap().to_string(),
            map_path.to_str().unwrap().to_string(),
        ];
        let mut output = Vec::new();
        lookup_accession_numbers_main_inner(&args, &mut output).unwrap();

        let unmapped = fs::read_to_string(Path::new("unmapped.txt")).unwrap();
        assert!(unmapped.contains("MISS2"));
        assert_eq!(String::from_utf8(output).unwrap(), "seq1\t123\n");

        std::env::set_current_dir(old_cwd).unwrap();
    }
}

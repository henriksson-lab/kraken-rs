use ahash::AHashMap as HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

use crate::types::TaxId;

/// Look up accession numbers in accession-to-taxid map files.
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

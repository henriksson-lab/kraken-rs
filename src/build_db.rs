use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, IsTerminal};

use crate::compact_hash::CompactHashTable;
use crate::hash::murmurhash3;
use crate::minimizer::MinimizerScanner;
use crate::seq::BatchSequenceReader;
use crate::taxonomy::{NCBITaxonomy, Taxonomy};
use crate::types::*;

/// Options for database building.
pub struct BuildOptions {
    pub id_to_taxon_map_filename: String,
    pub ncbi_taxonomy_directory: String,
    pub hashtable_filename: String,
    pub options_filename: String,
    pub taxonomy_filename: String,
    pub block_size: usize,
    pub subblock_size: usize,
    pub requested_bits_for_taxid: usize,
    pub num_threads: usize,
    pub input_is_protein: bool,
    pub k: usize,
    pub l: usize,
    pub capacity: usize,
    pub maximum_capacity: usize,
    pub spaced_seed_mask: u64,
    pub toggle_mask: u64,
    pub min_clear_hash_value: u64,
    pub deterministic_build: bool,
}

impl Default for BuildOptions {
    fn default() -> Self {
        BuildOptions {
            id_to_taxon_map_filename: String::new(),
            ncbi_taxonomy_directory: String::new(),
            hashtable_filename: String::new(),
            options_filename: String::new(),
            taxonomy_filename: String::new(),
            block_size: 10 * 1024 * 1024,
            subblock_size: 1024,
            requested_bits_for_taxid: 0,
            num_threads: 1,
            input_is_protein: false,
            k: 0,
            l: 0,
            capacity: 0,
            maximum_capacity: 0,
            spaced_seed_mask: DEFAULT_SPACED_SEED_MASK,
            toggle_mask: DEFAULT_TOGGLE_MASK,
            min_clear_hash_value: 0,
            deterministic_build: true,
        }
    }
}

/// Extract NCBI sequence IDs from a FASTA header.
/// Handles \x01 delimiters used in non-redundant databases.
pub fn extract_ncbi_sequence_ids(header: &str) -> Vec<String> {
    let mut list = Vec::new();
    let mut current = String::new();
    let mut in_id = true;

    for ch in header.bytes() {
        if ch == 0x01 {
            if !current.is_empty() {
                list.push(current.clone());
            }
            current.clear();
            in_id = true;
        } else if in_id && (ch as char).is_whitespace() {
            if !current.is_empty() {
                list.push(current.clone());
            }
            current.clear();
            in_id = false;
        } else if in_id {
            current.push(ch as char);
        }
    }
    if !current.is_empty() {
        list.push(current);
    }
    list
}

/// Read a seqid-to-taxon mapping file.
pub fn read_id_to_taxon_map(filename: &str) -> io::Result<BTreeMap<String, TaxId>> {
    let mut id_map = BTreeMap::new();
    let file = BufReader::new(File::open(filename)?);
    for line in file.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            if let Ok(taxid) = parts[1].parse::<u64>() {
                if taxid != 0 {
                    id_map.insert(parts[0].to_string(), taxid);
                }
            }
        }
    }
    Ok(id_map)
}

/// Generate taxonomy from NCBI dumps and a seqid-to-taxon map.
pub fn generate_taxonomy(
    ncbi_dir: &str,
    id_map: &BTreeMap<String, TaxId>,
    output_filename: &str,
) -> io::Result<()> {
    let nodes_file = format!("{}/nodes.dmp", ncbi_dir);
    let names_file = format!("{}/names.dmp", ncbi_dir);

    let mut ncbi_tax = NCBITaxonomy::new(&nodes_file, &names_file)?;
    for &taxid in id_map.values() {
        if taxid != 0 {
            ncbi_tax.mark_node(taxid);
        }
    }
    ncbi_tax.convert_to_kraken_taxonomy(output_filename)
}

/// Set the LCA for a minimizer in the hash table.
/// Retries CompareAndSet in a loop until success.
/// Used in the deterministic build path.
pub fn set_minimizer_lca(
    hash: &CompactHashTable,
    minimizer: u64,
    taxid: HValue,
    taxonomy: &Taxonomy,
) {
    let mut old_value: HValue = 0;
    let mut new_value = taxid;
    while !hash.compare_and_set(minimizer, new_value, &mut old_value) {
        new_value = taxonomy.lowest_common_ancestor(new_value as u64, old_value as u64) as HValue;
    }
}

/// Deterministic sequence processing for a single sequence.
/// Exact port of C++ `ProcessSequence()` from `build_db.cc:324-433`.
/// Uses block/subblock decomposition with safe-prefix parallel insertion.
fn process_sequence_deterministic(
    seq: &str,
    taxid: HValue,
    hash: &CompactHashTable,
    taxonomy: &Taxonomy,
    opts: &BuildOptions,
) {
    use parking_lot::Mutex;
    use std::collections::BTreeSet;

    const SET_CT: usize = 256;

    for j in (0..seq.len()).step_by(opts.block_size) {
        let block_start = j;
        let block_finish = (j + opts.block_size + opts.k - 1).min(seq.len());

        // Gather minimizers from subblocks into zone-partitioned sets
        let minimizer_sets: Vec<Mutex<BTreeSet<u64>>> =
            (0..SET_CT).map(|_| Mutex::new(BTreeSet::new())).collect();

        let subblocks: Vec<(usize, usize)> = (block_start..block_finish)
            .step_by(opts.subblock_size)
            .map(|i| {
                let sub_finish = (i + opts.subblock_size + opts.k - 1).min(block_finish);
                (i, sub_finish)
            })
            .collect();

        // Process subblocks (could use rayon here for parallelism)
        for &(sub_start, sub_finish) in &subblocks {
            let mut scanner = MinimizerScanner::new(
                opts.k as isize,
                opts.l as isize,
                opts.spaced_seed_mask,
                !opts.input_is_protein,
                opts.toggle_mask,
                CURRENT_REVCOM_VERSION,
            );
            scanner.load_sequence(seq, sub_start, sub_finish);
            while let Some(minimizer) = scanner.next_minimizer() {
                if scanner.is_ambiguous() {
                    continue;
                }
                let hc = murmurhash3(minimizer);
                if opts.min_clear_hash_value != 0 && hc < opts.min_clear_hash_value {
                    continue;
                }
                let zone = (hc % SET_CT as u64) as usize;
                minimizer_sets[zone].lock().insert(minimizer);
            }
        }

        // Combine sets into a single sorted list
        let mut minimizer_list: Vec<u64> = Vec::new();
        for set_mutex in &minimizer_sets {
            let set = set_mutex.lock();
            let mut sorted: Vec<u64> = set.iter().copied().collect();
            sorted.sort();
            minimizer_list.extend(sorted);
        }

        let mm_ct = minimizer_list.len();
        let mut index_list = vec![0usize; mm_ct];
        let mut insertion_list = vec![true; mm_ct];

        // Loop to enforce deterministic order
        let mut offset = 0;
        while offset < mm_ct {
            // Gather insertion point information
            for i in offset..mm_ct {
                if insertion_list[i] {
                    match hash.find_index(minimizer_list[i]) {
                        Some(idx) => {
                            insertion_list[i] = false;
                            index_list[i] = idx;
                        }
                        None => {
                            insertion_list[i] = true;
                            // index_list[i] left as-is (not used for non-insertions)
                        }
                    }
                }
            }

            // Determine safe prefix: insertions with unique insertion points
            let mut novel_insertion_points = BTreeSet::new();
            let mut safe_ct = offset;
            while safe_ct < mm_ct {
                if insertion_list[safe_ct] {
                    if novel_insertion_points.contains(&index_list[safe_ct]) {
                        break;
                    }
                    novel_insertion_points.insert(index_list[safe_ct]);
                }
                safe_ct += 1;
            }

            // Insert safe prefix
            for i in offset..safe_ct {
                set_minimizer_lca(hash, minimizer_list[i], taxid, taxonomy);
            }

            offset = safe_ct;
        }
    }
}

/// Fast (nondeterministic) sequence processing.
fn process_sequence_fast(
    seq: &str,
    taxid: HValue,
    hash: &CompactHashTable,
    taxonomy: &Taxonomy,
    scanner: &mut MinimizerScanner,
    min_clear_hash_value: u64,
) {
    scanner.load_sequence(seq, 0, usize::MAX);
    while let Some(minimizer) = scanner.next_minimizer() {
        if scanner.is_ambiguous() {
            continue;
        }
        if min_clear_hash_value != 0 && murmurhash3(minimizer) < min_clear_hash_value {
            continue;
        }
        let mut existing_taxid: HValue = 0;
        let mut new_taxid = taxid;
        while !hash.compare_and_set(minimizer, new_taxid, &mut existing_taxid) {
            new_taxid =
                taxonomy.lowest_common_ancestor(new_taxid as u64, existing_taxid as u64) as HValue;
        }
    }
}

fn process_sequences_fast(
    opts: &BuildOptions,
    id_to_taxon_map: &BTreeMap<String, TaxId>,
    kraken_index: &CompactHashTable,
    taxonomy: &Taxonomy,
) -> io::Result<()> {
    let mut processed_seq_ct = 0usize;
    let mut processed_ch_ct = 0usize;
    let mut scanner = MinimizerScanner::new(
        opts.k as isize,
        opts.l as isize,
        opts.spaced_seed_mask,
        !opts.input_is_protein,
        opts.toggle_mask,
        CURRENT_REVCOM_VERSION,
    );
    let mut reader = BatchSequenceReader::new(None)?;

    while reader.load_block(opts.block_size) {
        while let Some(seq) = reader.next_sequence() {
            let mut sequence = seq.clone();
            let all_sequence_ids = extract_ncbi_sequence_ids(&sequence.header);
            let mut taxid: TaxId = 0;
            for seqid in &all_sequence_ids {
                if let Some(&ext_taxid) = id_to_taxon_map.get(seqid) {
                    if ext_taxid == 0 {
                        continue;
                    }
                    taxid =
                        taxonomy.lowest_common_ancestor(taxid, taxonomy.get_internal_id(ext_taxid));
                }
            }
            if taxid != 0 {
                if opts.input_is_protein && !sequence.seq.ends_with('*') {
                    sequence.seq.push('*');
                }
                process_sequence_fast(
                    &sequence.seq,
                    taxid as HValue,
                    kraken_index,
                    taxonomy,
                    &mut scanner,
                    opts.min_clear_hash_value,
                );
                processed_seq_ct += 1;
                processed_ch_ct += sequence.seq.len();
            }
        }
        if io::stderr().is_terminal() {
            eprint!(
                "\rProcessed {} sequences ({} {})...",
                processed_seq_ct,
                processed_ch_ct,
                if opts.input_is_protein { "aa" } else { "bp" }
            );
        }
    }
    if io::stderr().is_terminal() {
        eprint!("\r");
    }
    eprintln!(
        "Completed processing of {} sequences, {} {}",
        processed_seq_ct,
        processed_ch_ct,
        if opts.input_is_protein { "aa" } else { "bp" }
    );
    Ok(())
}

fn process_sequences(
    opts: &BuildOptions,
    id_to_taxon_map: &BTreeMap<String, TaxId>,
    kraken_index: &CompactHashTable,
    taxonomy: &Taxonomy,
) -> io::Result<()> {
    let mut processed_seq_ct = 0usize;
    let mut processed_ch_ct = 0usize;
    let mut reader = BatchSequenceReader::new(None)?;

    while reader.load_block(opts.block_size) {
        while let Some(sequence) = reader.next_sequence() {
            let all_sequence_ids = extract_ncbi_sequence_ids(&sequence.header);
            let mut taxid: TaxId = 0;
            for seqid in &all_sequence_ids {
                if let Some(&ext_taxid) = id_to_taxon_map.get(seqid) {
                    if ext_taxid == 0 {
                        continue;
                    }
                    taxid =
                        taxonomy.lowest_common_ancestor(taxid, taxonomy.get_internal_id(ext_taxid));
                }
            }
            if taxid != 0 {
                let mut seq = sequence.seq.clone();
                if opts.input_is_protein && !seq.ends_with('*') {
                    seq.push('*');
                }
                process_sequence_deterministic(&seq, taxid as HValue, kraken_index, taxonomy, opts);
                processed_seq_ct += 1;
                processed_ch_ct += seq.len();
            }
        }
        if io::stderr().is_terminal() {
            eprint!(
                "\rProcessed {} sequences ({} {})...",
                processed_seq_ct,
                processed_ch_ct,
                if opts.input_is_protein { "aa" } else { "bp" }
            );
        }
    }
    if io::stderr().is_terminal() {
        eprint!("\r");
    }
    eprintln!(
        "Completed processing of {} sequences, {} {}",
        processed_seq_ct,
        processed_ch_ct,
        if opts.input_is_protein { "aa" } else { "bp" }
    );
    Ok(())
}

fn parse_command_line(args: &[String], opts: &mut BuildOptions) -> io::Result<()> {
    let mut i = 1usize;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "-?" => return Err(io::Error::other(usage(0))),
            "-B" => {
                i += 1;
                let sig = args
                    .get(i)
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing -B value"))?
                    .parse::<i64>()
                    .map_err(|_| {
                        io::Error::new(io::ErrorKind::InvalidInput, "must have positive block size")
                    })?;
                if sig < 1 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "must have positive block size",
                    ));
                }
                opts.block_size = sig as usize;
            }
            "-b" => {
                i += 1;
                let sig = args
                    .get(i)
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing -b value"))?
                    .parse::<i64>()
                    .map_err(|_| {
                        io::Error::new(
                            io::ErrorKind::InvalidInput,
                            "must have positive subblock size",
                        )
                    })?;
                if sig < 1 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "must have positive subblock size",
                    ));
                }
                opts.subblock_size = sig as usize;
            }
            "-r" => {
                i += 1;
                let sig = args
                    .get(i)
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing -r value"))?
                    .parse::<i64>()
                    .map_err(|_| {
                        io::Error::new(
                            io::ErrorKind::InvalidInput,
                            "can't have negative bit storage",
                        )
                    })?;
                if sig < 0 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "can't have negative bit storage",
                    ));
                }
                if sig > 31 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "can't have more than 31 bits of storage for taxid",
                    ));
                }
                opts.requested_bits_for_taxid = sig as usize;
            }
            "-p" => {
                i += 1;
                let sig = args
                    .get(i)
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing -p value"))?
                    .parse::<i64>()
                    .map_err(|_| {
                        io::Error::new(
                            io::ErrorKind::InvalidInput,
                            "can't have negative number of threads",
                        )
                    })?;
                if sig < 1 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "can't have negative number of threads",
                    ));
                }
                opts.num_threads = sig as usize;
            }
            "-H" => {
                i += 1;
                opts.hashtable_filename = args.get(i).cloned().unwrap_or_default();
            }
            "-m" => {
                i += 1;
                opts.id_to_taxon_map_filename = args.get(i).cloned().unwrap_or_default();
            }
            "-n" => {
                i += 1;
                opts.ncbi_taxonomy_directory = args.get(i).cloned().unwrap_or_default();
            }
            "-o" => {
                i += 1;
                opts.options_filename = args.get(i).cloned().unwrap_or_default();
            }
            "-t" => {
                i += 1;
                opts.taxonomy_filename = args.get(i).cloned().unwrap_or_default();
            }
            "-S" => {
                i += 1;
                opts.spaced_seed_mask =
                    u64::from_str_radix(args.get(i).map(String::as_str).unwrap_or(""), 2).map_err(
                        |_| io::Error::new(io::ErrorKind::InvalidInput, "invalid spaced seed mask"),
                    )?;
            }
            "-T" => {
                i += 1;
                opts.toggle_mask =
                    u64::from_str_radix(args.get(i).map(String::as_str).unwrap_or(""), 2).map_err(
                        |_| io::Error::new(io::ErrorKind::InvalidInput, "invalid toggle mask"),
                    )?;
            }
            "-k" => {
                i += 1;
                let sig = args
                    .get(i)
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing -k value"))?
                    .parse::<i64>()
                    .map_err(|_| {
                        io::Error::new(io::ErrorKind::InvalidInput, "k must be positive integer")
                    })?;
                if sig < 1 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "k must be positive integer",
                    ));
                }
                opts.k = sig as usize;
            }
            "-l" => {
                i += 1;
                let sig = args
                    .get(i)
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing -l value"))?
                    .parse::<i64>()
                    .map_err(|_| {
                        io::Error::new(io::ErrorKind::InvalidInput, "l must be positive integer")
                    })?;
                if sig < 1 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "l must be positive integer",
                    ));
                }
                if sig > 31 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "l must be no more than 31",
                    ));
                }
                opts.l = sig as usize;
            }
            "-c" => {
                i += 1;
                let sig = args
                    .get(i)
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing -c value"))?
                    .parse::<i64>()
                    .map_err(|_| {
                        io::Error::new(
                            io::ErrorKind::InvalidInput,
                            "capacity must be positive integer",
                        )
                    })?;
                if sig < 1 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "capacity must be positive integer",
                    ));
                }
                opts.capacity = sig as usize;
            }
            "-M" => {
                i += 1;
                let sig = args
                    .get(i)
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing -M value"))?
                    .parse::<i64>()
                    .map_err(|_| {
                        io::Error::new(
                            io::ErrorKind::InvalidInput,
                            "max capacity must be positive integer",
                        )
                    })?;
                if sig < 1 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "max capacity must be positive integer",
                    ));
                }
                opts.maximum_capacity = sig as usize;
            }
            "-F" => opts.deterministic_build = false,
            "-X" => opts.input_is_protein = true,
            _ => {}
        }
        i += 1;
    }

    if opts.spaced_seed_mask != DEFAULT_SPACED_SEED_MASK {
        crate::utilities::expand_spaced_seed_mask(
            &mut opts.spaced_seed_mask,
            if opts.input_is_protein {
                BITS_PER_CHAR_PRO as i32
            } else {
                BITS_PER_CHAR_DNA as i32
            },
        );
    }
    if opts.hashtable_filename.is_empty()
        || opts.id_to_taxon_map_filename.is_empty()
        || opts.ncbi_taxonomy_directory.is_empty()
        || opts.options_filename.is_empty()
        || opts.taxonomy_filename.is_empty()
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "missing mandatory filename parameter",
        ));
    }
    if opts.k == 0 || opts.l == 0 || opts.capacity == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "missing mandatory integer parameter",
        ));
    }
    if opts.k < opts.l {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "k cannot be less than l",
        ));
    }
    if opts.block_size < opts.subblock_size {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "block size cannot be less than subblock size",
        ));
    }
    if opts.maximum_capacity > opts.capacity {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "maximum capacity option shouldn't specify larger capacity than normal",
        ));
    }

    Ok(())
}

fn usage(_exit_code: i32) -> String {
    [
        "Usage: build_db <options>",
        "",
        "Options (*mandatory):",
        "* -H FILENAME   Kraken 2 hash table filename",
        "* -m FILENAME   Sequence ID to taxon map filename",
        "* -t FILENAME   Kraken 2 taxonomy filename",
        "* -n DIR        NCBI taxonomy directory name",
        "* -o FILENAME   Kraken 2 options filename",
        "* -k INT        Set length of k-mers",
        "* -l INT        Set length of minimizers",
        "* -c INT        Set capacity of hash table",
        "  -M INT        Set maximum capacity of hash table (MiniKraken)",
        "  -S BITSTRING  Spaced seed mask",
        "  -T BITSTRING  Minimizer toggle mask",
        "  -X            Input seqs. are proteins",
        "  -p INT        Number of threads",
        "  -F            Use fast, nondeterministic building method",
        "  -B INT        Read block size",
        "  -b INT        Read subblock size",
        "  -r INT        Bit storage requested for taxid",
    ]
    .join("\n")
}

pub fn build_db_main(args: &[String]) -> io::Result<()> {
    let mut opts = BuildOptions {
        spaced_seed_mask: DEFAULT_SPACED_SEED_MASK,
        toggle_mask: DEFAULT_TOGGLE_MASK,
        input_is_protein: false,
        num_threads: 1,
        block_size: 10 * 1024 * 1024,
        subblock_size: 1024,
        requested_bits_for_taxid: 0,
        min_clear_hash_value: 0,
        maximum_capacity: 0,
        deterministic_build: true,
        ..BuildOptions::default()
    };
    parse_command_line(args, &mut opts)?;

    let id_map = read_id_to_taxon_map(&opts.id_to_taxon_map_filename)?;
    generate_taxonomy(
        &opts.ncbi_taxonomy_directory,
        &id_map,
        &opts.taxonomy_filename,
    )?;
    eprintln!("Taxonomy parsed and converted.");

    let mut taxonomy = Taxonomy::from_file(&opts.taxonomy_filename, false)?;
    taxonomy.generate_external_to_internal_id_map();
    let mut bits_needed_for_value = 1usize;
    while (1usize << bits_needed_for_value) < taxonomy.node_count() {
        bits_needed_for_value += 1;
    }
    if opts.requested_bits_for_taxid > 0 && bits_needed_for_value > opts.requested_bits_for_taxid {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "more bits required for storing taxid",
        ));
    }

    let bits_for_taxid = bits_needed_for_value.max(opts.requested_bits_for_taxid);
    let mut actual_capacity = opts.capacity;
    if opts.maximum_capacity > 0 {
        let frac = opts.maximum_capacity as f64 / opts.capacity as f64;
        if frac > 1.0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "maximum capacity larger than requested capacity",
            ));
        }
        opts.min_clear_hash_value = ((1.0 - frac) * u64::MAX as f64) as u64;
        actual_capacity = opts.maximum_capacity;
    }

    let kraken_index = CompactHashTable::new(actual_capacity, 32 - bits_for_taxid, bits_for_taxid);
    eprintln!(
        "CHT created with {} bits reserved for taxid.",
        bits_for_taxid
    );

    if opts.deterministic_build {
        process_sequences(&opts, &id_map, &kraken_index, &taxonomy)?;
    } else {
        process_sequences_fast(&opts, &id_map, &kraken_index, &taxonomy)?;
    }

    eprint!("Writing data to disk... ");
    kraken_index.write_table(&opts.hashtable_filename)?;

    let mut index_opts = IndexOptions::new();
    index_opts.k = opts.k;
    index_opts.l = opts.l;
    index_opts.spaced_seed_mask = opts.spaced_seed_mask;
    index_opts.toggle_mask = opts.toggle_mask;
    index_opts.dna_db = !opts.input_is_protein;
    index_opts.minimum_acceptable_hash_value = opts.min_clear_hash_value;
    index_opts.revcom_version = CURRENT_REVCOM_VERSION;
    index_opts.db_version = 0;
    index_opts.db_type = 0;
    index_opts.write_to_file(&opts.options_filename)?;
    eprintln!(" complete.");
    Ok(())
}

/// Build the database from input sequences read from stdin.
pub fn build_database(opts: &mut BuildOptions) -> io::Result<()> {
    // Read ID-to-taxon map
    let id_map = read_id_to_taxon_map(&opts.id_to_taxon_map_filename)?;

    // Generate taxonomy
    generate_taxonomy(
        &opts.ncbi_taxonomy_directory,
        &id_map,
        &opts.taxonomy_filename,
    )?;
    eprintln!("Taxonomy parsed and converted.");

    // Load taxonomy
    let mut taxonomy = Taxonomy::from_file(&opts.taxonomy_filename, false)?;
    taxonomy.generate_external_to_internal_id_map();

    // Compute bits needed for taxid
    let mut bits_needed: usize = 1;
    while (1usize << bits_needed) < taxonomy.node_count() {
        bits_needed += 1;
    }
    let bits_for_taxid = bits_needed.max(opts.requested_bits_for_taxid);

    // Handle maximum capacity (subsampling)
    let actual_capacity = if opts.maximum_capacity > 0 {
        let frac = opts.maximum_capacity as f64 / opts.capacity as f64;
        opts.min_clear_hash_value = ((1.0 - frac) * u64::MAX as f64) as u64;
        opts.maximum_capacity
    } else {
        opts.capacity
    };

    let hash = CompactHashTable::new(actual_capacity, 32 - bits_for_taxid, bits_for_taxid);
    eprintln!(
        "CHT created with {} bits reserved for taxid.",
        bits_for_taxid
    );

    if opts.deterministic_build {
        process_sequences(opts, &id_map, &hash, &taxonomy)?;
    } else {
        process_sequences_fast(opts, &id_map, &hash, &taxonomy)?;
    }

    // Write outputs
    eprint!("Writing data to disk... ");
    hash.write_table(&opts.hashtable_filename)?;

    let mut index_opts = IndexOptions::new();
    index_opts.k = opts.k;
    index_opts.l = opts.l;
    index_opts.spaced_seed_mask = opts.spaced_seed_mask;
    index_opts.toggle_mask = opts.toggle_mask;
    index_opts.dna_db = !opts.input_is_protein;
    index_opts.minimum_acceptable_hash_value = opts.min_clear_hash_value;
    index_opts.revcom_version = CURRENT_REVCOM_VERSION;

    index_opts
        .write_to_file(&opts.options_filename)
        .expect("Failed to write index options");
    eprintln!(" complete.");

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_db_parse_command_line() {
        let args = vec![
            "build_db".to_string(),
            "-H".to_string(),
            "hash.k2d".to_string(),
            "-m".to_string(),
            "map.txt".to_string(),
            "-t".to_string(),
            "taxo.k2d".to_string(),
            "-n".to_string(),
            "taxonomy".to_string(),
            "-o".to_string(),
            "opts.k2d".to_string(),
            "-k".to_string(),
            "35".to_string(),
            "-l".to_string(),
            "31".to_string(),
            "-c".to_string(),
            "1000".to_string(),
            "-F".to_string(),
            "-X".to_string(),
        ];
        let mut opts = BuildOptions::default();
        parse_command_line(&args, &mut opts).unwrap();
        assert_eq!(opts.hashtable_filename, "hash.k2d");
        assert_eq!(opts.id_to_taxon_map_filename, "map.txt");
        assert_eq!(opts.taxonomy_filename, "taxo.k2d");
        assert_eq!(opts.ncbi_taxonomy_directory, "taxonomy");
        assert_eq!(opts.options_filename, "opts.k2d");
        assert_eq!(opts.k, 35);
        assert_eq!(opts.l, 31);
        assert_eq!(opts.capacity, 1000);
        assert!(!opts.deterministic_build);
        assert!(opts.input_is_protein);
    }

    #[test]
    fn test_build_db_usage() {
        let text = usage(0);
        assert!(text.contains("Usage: build_db <options>"));
        assert!(text.contains("* -H FILENAME"));
        assert!(text.contains("-F            Use fast, nondeterministic building method"));
    }
}

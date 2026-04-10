use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

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
    use std::collections::BTreeSet;
    use parking_lot::Mutex;

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
            new_taxid = taxonomy.lowest_common_ancestor(new_taxid as u64, existing_taxid as u64) as HValue;
        }
    }
}

/// Build the database from input sequences read from stdin.
pub fn build_database(opts: &mut BuildOptions) -> io::Result<()> {
    // Read ID-to-taxon map
    let id_map = read_id_to_taxon_map(&opts.id_to_taxon_map_filename)?;

    // Generate taxonomy
    generate_taxonomy(&opts.ncbi_taxonomy_directory, &id_map, &opts.taxonomy_filename)?;
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
    eprintln!("CHT created with {} bits reserved for taxid.", bits_for_taxid);

    // Process sequences from stdin
    let mut reader = BatchSequenceReader::new(None)?;
    let mut scanner = MinimizerScanner::new(
        opts.k as isize,
        opts.l as isize,
        opts.spaced_seed_mask,
        !opts.input_is_protein,
        opts.toggle_mask,
        CURRENT_REVCOM_VERSION,
    );

    let mut processed_seq_ct = 0usize;
    let mut processed_ch_ct = 0usize;

    while reader.load_block(opts.block_size) {
        while let Some(seq) = reader.next_sequence() {
            let all_ids = extract_ncbi_sequence_ids(&seq.header);
            let mut taxid: TaxId = 0;
            for seqid in &all_ids {
                if let Some(&ext_taxid) = id_map.get(seqid) {
                    if ext_taxid != 0 {
                        taxid = taxonomy.lowest_common_ancestor(
                            taxid,
                            taxonomy.get_internal_id(ext_taxid),
                        );
                    }
                }
            }
            if taxid != 0 {
                let mut seq_str = seq.seq.clone();
                if opts.input_is_protein && !seq_str.ends_with('*') {
                    seq_str.push('*');
                }
                if opts.deterministic_build {
                    process_sequence_deterministic(
                        &seq_str, taxid as HValue, &hash, &taxonomy, opts,
                    );
                } else {
                    process_sequence_fast(
                        &seq_str, taxid as HValue, &hash, &taxonomy,
                        &mut scanner, opts.min_clear_hash_value,
                    );
                }
                processed_seq_ct += 1;
                processed_ch_ct += seq_str.len();
            }
        }
        eprint!("\rProcessed {} sequences ({} {})...",
            processed_seq_ct, processed_ch_ct,
            if opts.input_is_protein { "aa" } else { "bp" });
    }
    eprintln!("\nCompleted processing of {} sequences, {} {}",
        processed_seq_ct, processed_ch_ct,
        if opts.input_is_protein { "aa" } else { "bp" });

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

    index_opts.write_to_file(&opts.options_filename)
        .expect("Failed to write index options");
    eprintln!(" complete.");

    Ok(())
}

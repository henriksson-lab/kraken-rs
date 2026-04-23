use std::io;

use crate::aa_translate::translate_to_all_frames;
use crate::compact_hash::CompactHashTable;
use crate::hash::murmurhash3;
use crate::minimizer::MinimizerScanner;
use crate::taxonomy::{generate_taxonomy_libtax, Taxonomy};
use crate::types::IndexOptions;
use crate::utilities::expand_spaced_seed_mask;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TaxonomyNodeFields {
    pub parent_id: u64,
    pub first_child: u64,
    pub child_count: u64,
    pub name_offset: u64,
    pub rank_offset: u64,
    pub external_id: u64,
    pub godparent_id: u64,
}

pub struct ScannerContext {
    pub scanner: MinimizerScanner,
    pub seq_storage: String,
}

pub fn ffi_scanner_new(
    k: isize,
    l: isize,
    spaced_seed_mask: u64,
    dna_sequence: bool,
    toggle_mask: u64,
    revcom_version: i32,
) -> MinimizerScanner {
    MinimizerScanner::new(
        k,
        l,
        spaced_seed_mask,
        dna_sequence,
        toggle_mask,
        revcom_version,
    )
}

pub fn ffi_scanner_destroy(_scanner: MinimizerScanner) {}

pub fn ffi_scanner_context_new(
    k: isize,
    l: isize,
    spaced_seed_mask: u64,
    dna_sequence: bool,
    toggle_mask: u64,
    revcom_version: i32,
) -> ScannerContext {
    ScannerContext {
        scanner: MinimizerScanner::new(
            k,
            l,
            spaced_seed_mask,
            dna_sequence,
            toggle_mask,
            revcom_version,
        ),
        seq_storage: String::new(),
    }
}

pub fn ffi_scanner_context_destroy(_context: ScannerContext) {}

pub fn ffi_scanner_context_load_sequence(
    context: &mut ScannerContext,
    seq: &str,
    start: usize,
    finish: usize,
) {
    context.seq_storage.clear();
    context.seq_storage.push_str(seq);
    context
        .scanner
        .load_sequence(&context.seq_storage, start, finish);
}

pub fn ffi_scanner_context_next_minimizer(context: &mut ScannerContext) -> Option<u64> {
    context.scanner.next_minimizer()
}

pub fn ffi_scanner_context_is_ambiguous(context: &ScannerContext) -> bool {
    context.scanner.is_ambiguous()
}

pub fn ffi_cht_new(capacity: usize, key_bits: usize, value_bits: usize) -> CompactHashTable {
    CompactHashTable::new(capacity, key_bits, value_bits)
}

pub fn ffi_cht_load(filename: &str, memory_mapping: bool) -> io::Result<CompactHashTable> {
    CompactHashTable::from_file(filename, memory_mapping)
}

pub fn ffi_cht_destroy(_cht: CompactHashTable) {}

pub fn ffi_cht_get(cht: &CompactHashTable, key: u64) -> u32 {
    cht.get(key)
}

pub fn ffi_cht_find_index(cht: &CompactHashTable, key: u64) -> Option<usize> {
    cht.find_index(key)
}

pub fn ffi_cht_compare_and_set(
    cht: &CompactHashTable,
    key: u64,
    new_value: u32,
    old_value: &mut u32,
) -> bool {
    cht.compare_and_set(key, new_value, old_value)
}

pub fn ffi_cht_direct_compare_and_set(
    cht: &CompactHashTable,
    idx: usize,
    key: u64,
    new_value: u32,
    old_value: &mut u32,
) -> bool {
    cht.direct_compare_and_set(idx, key, new_value, old_value)
}

pub fn ffi_cht_write(cht: &CompactHashTable, filename: &str) -> io::Result<()> {
    cht.write_table(filename)
}

pub fn ffi_cht_capacity(cht: &CompactHashTable) -> usize {
    cht.capacity()
}

pub fn ffi_cht_size(cht: &CompactHashTable) -> usize {
    cht.size()
}

pub fn ffi_cht_key_bits(cht: &CompactHashTable) -> usize {
    cht.key_bits()
}

pub fn ffi_cht_value_bits(cht: &CompactHashTable) -> usize {
    cht.value_bits()
}

pub fn ffi_murmurhash3(key: u64) -> u64 {
    murmurhash3(key)
}

pub fn ffi_expand_spaced_seed_mask(mask: &mut u64, bit_expansion_factor: i32) {
    expand_spaced_seed_mask(mask, bit_expansion_factor);
}

pub fn ffi_translate_to_all_frames(dna: &str) -> Vec<String> {
    translate_to_all_frames(dna)
}

pub fn ffi_free_string(_s: String) {}

pub fn ffi_read_index_options(filename: &str, opts_out: &mut IndexOptions) -> bool {
    match IndexOptions::read_from_file(filename) {
        Ok(opts) => {
            *opts_out = opts;
            true
        }
        Err(_) => false,
    }
}

pub fn ffi_write_index_options(filename: &str, opts: &IndexOptions) -> bool {
    opts.write_to_file(filename).is_ok()
}

pub fn ffi_sizeof_index_options() -> usize {
    std::mem::size_of::<IndexOptions>()
}

pub fn ffi_taxonomy_load(filename: &str, memory_mapping: bool) -> io::Result<Taxonomy> {
    let mut tax = Taxonomy::from_file(filename, memory_mapping)?;
    tax.generate_external_to_internal_id_map();
    Ok(tax)
}

pub fn ffi_taxonomy_destroy(_taxonomy: Taxonomy) {}

pub fn ffi_taxonomy_node_count(taxonomy: &Taxonomy) -> u64 {
    taxonomy.node_count() as u64
}

pub fn ffi_taxonomy_get_internal_id(taxonomy: &Taxonomy, external_id: u64) -> u64 {
    taxonomy.get_internal_id(external_id)
}

pub fn ffi_taxonomy_lca(taxonomy: &Taxonomy, a: u64, b: u64) -> u64 {
    taxonomy.lowest_common_ancestor(a, b)
}

pub fn ffi_taxonomy_is_a_ancestor_of_b(taxonomy: &Taxonomy, a: u64, b: u64) -> bool {
    taxonomy.is_a_ancestor_of_b(a, b)
}

pub fn ffi_taxonomy_get_node(taxonomy: &Taxonomy, internal_id: u64) -> TaxonomyNodeFields {
    let node = taxonomy.node(internal_id);
    TaxonomyNodeFields {
        parent_id: node.parent_id,
        first_child: node.first_child,
        child_count: node.child_count,
        name_offset: node.name_offset,
        rank_offset: node.rank_offset,
        external_id: node.external_id,
        godparent_id: node.godparent_id,
    }
}

pub fn ffi_taxonomy_name_data(taxonomy: &Taxonomy) -> &[u8] {
    taxonomy.name_data()
}

pub fn ffi_taxonomy_rank_data(taxonomy: &Taxonomy) -> &[u8] {
    taxonomy.rank_data()
}

pub fn ffi_taxonomy_write_to_disk(taxonomy: &Taxonomy, filename: &str) -> io::Result<()> {
    taxonomy.write_to_disk(filename)
}

pub fn ffi_generate_taxonomy(
    nodes_filename: &str,
    names_filename: &str,
    seqid2taxid_filename: &str,
    output_filename: &str,
) -> io::Result<()> {
    generate_taxonomy_libtax(
        names_filename,
        nodes_filename,
        seqid2taxid_filename,
        output_filename,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{CURRENT_REVCOM_VERSION, DEFAULT_SPACED_SEED_MASK, DEFAULT_TOGGLE_MASK};

    #[test]
    fn test_ffi_hash_and_mask_helpers() {
        assert_eq!(ffi_murmurhash3(42), murmurhash3(42));

        let mut mask = 0b010110u64;
        ffi_expand_spaced_seed_mask(&mut mask, 2);
        assert_eq!(mask, 0b001100111100);
    }

    #[test]
    fn test_ffi_translate_to_all_frames_helper() {
        let frames = ffi_translate_to_all_frames("AAAGGGCCCTTT");
        assert_eq!(frames.len(), 6);
        assert_eq!(frames[0], "KGP");
        assert!(frames.iter().all(|frame| !frame.is_empty()));

        ffi_free_string(frames[0].clone());
    }

    #[test]
    fn test_ffi_index_options_helpers() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().to_str().unwrap();

        let mut opts = IndexOptions::new();
        opts.k = 35;
        opts.l = 31;
        opts.spaced_seed_mask = DEFAULT_SPACED_SEED_MASK;
        opts.toggle_mask = DEFAULT_TOGGLE_MASK;
        opts.dna_db = true;
        opts.minimum_acceptable_hash_value = 17;
        opts.revcom_version = CURRENT_REVCOM_VERSION;
        opts.db_version = 2;
        opts.db_type = 1;

        assert!(ffi_write_index_options(path, &opts));
        assert_eq!(
            ffi_sizeof_index_options(),
            std::mem::size_of::<IndexOptions>()
        );

        let mut loaded = IndexOptions::new();
        assert!(ffi_read_index_options(path, &mut loaded));
        assert_eq!(loaded.k, opts.k);
        assert_eq!(loaded.l, opts.l);
        assert_eq!(loaded.spaced_seed_mask, opts.spaced_seed_mask);
        assert_eq!(loaded.toggle_mask, opts.toggle_mask);
        assert_eq!(loaded.dna_db, opts.dna_db);
        assert_eq!(
            loaded.minimum_acceptable_hash_value,
            opts.minimum_acceptable_hash_value
        );
        assert_eq!(loaded.revcom_version, opts.revcom_version);
        assert_eq!(loaded.db_version, opts.db_version);
        assert_eq!(loaded.db_type, opts.db_type);
    }

    #[test]
    fn test_ffi_taxonomy_helpers() {
        let ref_path = format!("{}/tests/reference/taxo.k2d", env!("CARGO_MANIFEST_DIR"));
        if !std::path::Path::new(&ref_path).exists() {
            eprintln!("Skipping: reference data not available");
            return;
        }

        let tax = ffi_taxonomy_load(&ref_path, false).unwrap();
        assert!(ffi_taxonomy_node_count(&tax) > 1);
        assert_eq!(ffi_taxonomy_get_internal_id(&tax, 1), 1);
        assert_eq!(ffi_taxonomy_lca(&tax, 1, 1), 1);
        assert!(ffi_taxonomy_is_a_ancestor_of_b(&tax, 1, 1));

        let root = ffi_taxonomy_get_node(&tax, 1);
        assert_eq!(root.external_id, 1);
        assert!(!ffi_taxonomy_name_data(&tax).is_empty());
        assert!(!ffi_taxonomy_rank_data(&tax).is_empty());
    }

    #[test]
    fn test_ffi_generate_taxonomy() {
        let data_dir = format!("{}/kraken2/data", env!("CARGO_MANIFEST_DIR"));
        let nodes_file = format!("{}/nodes.dmp", data_dir);
        let names_file = format!("{}/names.dmp", data_dir);
        if !std::path::Path::new(&nodes_file).exists()
            || !std::path::Path::new(&names_file).exists()
        {
            eprintln!("Skipping taxonomy generation test: source data not available");
            return;
        }

        let tmp_dir = tempfile::tempdir().unwrap();
        let seqid2taxid = tmp_dir.path().join("seqid2taxid.map");
        let out = tmp_dir.path().join("taxo.k2d");
        std::fs::write(&seqid2taxid, "seq1\t2697049\nseq2\t694009\n").unwrap();

        ffi_generate_taxonomy(
            &nodes_file,
            &names_file,
            seqid2taxid.to_str().unwrap(),
            out.to_str().unwrap(),
        )
        .unwrap();

        let tax = ffi_taxonomy_load(out.to_str().unwrap(), false).unwrap();
        assert!(ffi_taxonomy_node_count(&tax) > 1);
    }

    #[test]
    fn test_ffi_scanner_context_round_trip() {
        let mut context = ffi_scanner_context_new(
            10,
            5,
            DEFAULT_SPACED_SEED_MASK,
            true,
            DEFAULT_TOGGLE_MASK,
            CURRENT_REVCOM_VERSION,
        );
        ffi_scanner_context_load_sequence(&mut context, "ACGTACGTACGTACGT", 0, usize::MAX);

        let mut minimizers = Vec::new();
        while let Some(minimizer) = ffi_scanner_context_next_minimizer(&mut context) {
            minimizers.push((minimizer, ffi_scanner_context_is_ambiguous(&context)));
        }

        assert!(!minimizers.is_empty());
        assert!(minimizers.iter().any(|(_, ambiguous)| !ambiguous));
    }

    #[test]
    fn test_ffi_compact_hash_helpers() {
        let cht = ffi_cht_new(128, 22, 10);
        assert_eq!(ffi_cht_capacity(&cht), 128);
        assert_eq!(ffi_cht_size(&cht), 0);
        assert_eq!(ffi_cht_key_bits(&cht), 22);
        assert_eq!(ffi_cht_value_bits(&cht), 10);

        let mut old = 0u32;
        assert!(ffi_cht_compare_and_set(&cht, 42, 7, &mut old));
        assert_eq!(ffi_cht_get(&cht, 42), 7);
        let idx = ffi_cht_find_index(&cht, 42).unwrap();

        old = 7;
        assert!(ffi_cht_direct_compare_and_set(&cht, idx, 42, 9, &mut old));
        assert_eq!(ffi_cht_get(&cht, 42), 9);

        let tmp_dir = tempfile::tempdir().unwrap();
        let out = tmp_dir.path().join("hash.k2d");
        ffi_cht_write(&cht, out.to_str().unwrap()).unwrap();
        let loaded = ffi_cht_load(out.to_str().unwrap(), false).unwrap();
        assert_eq!(ffi_cht_get(&loaded, 42), 9);
    }
}

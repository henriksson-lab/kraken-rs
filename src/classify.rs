use ahash::AHashMap as HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};

use rayon::prelude::*;

use crate::aa_translate::translate_to_all_frames;
use crate::compact_hash::CompactHashTable;
use crate::hash::murmurhash3;
use crate::minimizer::MinimizerScanner;
use crate::readcounts::TaxonCounters;
use crate::reports;
use crate::seq::BatchSequenceReader;
use crate::taxonomy::Taxonomy;
use crate::types::*;

/// Special taxon IDs used as markers in the taxa vector.
const MATE_PAIR_BORDER_TAXON: TaxId = TaxId::MAX;
const READING_FRAME_BORDER_TAXON: TaxId = TaxId::MAX - 1;
const AMBIGUOUS_SPAN_TAXON: TaxId = TaxId::MAX - 2;

/// Classification options.
#[derive(Clone)]
pub struct ClassifyOptions {
    pub confidence_threshold: f64,
    pub minimum_quality_score: u8,
    pub minimum_hit_groups: i64,
    pub paired_end_processing: bool,
    pub single_file_pairs: bool,
    pub use_translated_search: bool,
    pub print_scientific_name: bool,
    pub quick_mode: bool,
    pub report_filename: String,
    pub report_kmer_data: bool,
    pub report_zero_counts: bool,
    pub mpa_style_report: bool,
    pub use_memory_mapping: bool,
    pub num_threads: usize,
}

impl Default for ClassifyOptions {
    fn default() -> Self {
        ClassifyOptions {
            confidence_threshold: 0.0,
            minimum_quality_score: 0,
            minimum_hit_groups: 0,
            paired_end_processing: false,
            single_file_pairs: false,
            use_translated_search: false,
            print_scientific_name: false,
            quick_mode: false,
            report_filename: String::new(),
            report_kmer_data: false,
            report_zero_counts: false,
            mpa_style_report: false,
            use_memory_mapping: false,
            num_threads: 1,
        }
    }
}

/// Classification statistics.
#[derive(Default)]
pub struct ClassificationStats {
    pub total_sequences: u64,
    pub total_classified: u64,
    pub total_unclassified: u64,
}

/// Resolve the classification tree using LTR scoring and confidence walk-up.
/// Exact port of C++ `ResolveTree()`.
pub fn resolve_tree(
    hit_counts: &HashMap<TaxId, u32>,
    taxonomy: &Taxonomy,
    total_minimizers: usize,
    confidence_threshold: f64,
) -> TaxId {
    let mut max_taxon: TaxId = 0;
    let mut max_score: u32 = 0;
    let required_score = (confidence_threshold * total_minimizers as f64).ceil() as u32;

    // Phase 1: Find taxon with highest LTR score
    for &taxon in hit_counts.keys() {
        let mut score: u32 = 0;
        for (&taxon2, &count2) in hit_counts {
            if taxonomy.is_a_ancestor_of_b(taxon2, taxon) {
                score += count2;
            }
        }
        if score > max_score {
            max_score = score;
            max_taxon = taxon;
        } else if score == max_score {
            max_taxon = taxonomy.lowest_common_ancestor(max_taxon, taxon);
        }
    }

    // Phase 2: Confidence walk-up
    max_score = hit_counts.get(&max_taxon).copied().unwrap_or(0);
    while max_taxon != 0 && max_score < required_score {
        max_score = 0;
        for (&taxon, &count) in hit_counts {
            if taxonomy.is_a_ancestor_of_b(max_taxon, taxon) {
                max_score += count;
            }
        }
        if max_score >= required_score {
            return max_taxon;
        }
        max_taxon = taxonomy.node(max_taxon).parent_id;
    }

    max_taxon
}

/// Format the hit list string from the taxa vector.
/// Exact port of C++ `AddHitlistString()`.
pub fn add_hitlist_string(taxa: &[TaxId], taxonomy: &Taxonomy) -> String {
    if taxa.is_empty() {
        return "0:0".to_string();
    }

    let mut result = String::new();
    let mut last_code = taxa[0];
    let mut code_count = 1u32;

    for i in 1..taxa.len() {
        let code = taxa[i];
        if code == last_code {
            code_count += 1;
        } else {
            append_hitlist_entry(&mut result, last_code, code_count, taxonomy, true);
            code_count = 1;
            last_code = code;
        }
    }
    // Last entry (no trailing space for non-border entries)
    append_hitlist_entry(&mut result, last_code, code_count, taxonomy, false);

    result
}

fn append_hitlist_entry(
    result: &mut String,
    code: TaxId,
    count: u32,
    taxonomy: &Taxonomy,
    with_trailing_space: bool,
) {
    use std::fmt::Write as FmtWrite;
    if code != MATE_PAIR_BORDER_TAXON && code != READING_FRAME_BORDER_TAXON {
        if code == AMBIGUOUS_SPAN_TAXON {
            let _ = write!(result, "A:{}", count);
        } else {
            let ext_code = taxonomy.node(code).external_id;
            let _ = write!(result, "{}:{}", ext_code, count);
        }
        if with_trailing_space {
            result.push(' ');
        }
    } else {
        if code == MATE_PAIR_BORDER_TAXON {
            result.push_str("|:|");
        } else {
            result.push_str("-:-");
        }
        if with_trailing_space {
            result.push(' ');
        }
    }
}

/// Trim /1 or /2 from paired-end read names.
fn trim_pair_info(id: &str) -> &str {
    let bytes = id.as_bytes();
    if bytes.len() > 2 && bytes[bytes.len() - 2] == b'/' {
        let last = bytes[bytes.len() - 1];
        if last == b'1' || last == b'2' {
            return &id[..id.len() - 2];
        }
    }
    id
}

/// Mask low-quality bases in a FASTQ sequence.
fn mask_low_quality_bases(seq: &mut Sequence, min_quality: u8) {
    if seq.format != SequenceFormat::Fastq || seq.quals.is_empty() {
        return;
    }
    let mut seq_bytes = std::mem::take(&mut seq.seq).into_bytes();
    let qual_bytes = seq.quals.as_bytes();
    for i in 0..seq_bytes.len().min(qual_bytes.len()) {
        if qual_bytes[i].wrapping_sub(b'!') < min_quality {
            seq_bytes[i] = b'x';
        }
    }
    // Safety: we only replace ASCII bases with ASCII 'x', so UTF-8 validity is preserved
    seq.seq = String::from_utf8(seq_bytes)
        .unwrap_or_else(|e| String::from_utf8_lossy(e.as_bytes()).into_owned());
}

/// Classify a single sequence (or pair).
/// Returns the taxonomy call (0 = unclassified).
pub fn classify_sequence(
    dna: &Sequence,
    dna2: Option<&Sequence>,
    hash: &CompactHashTable,
    taxonomy: &Taxonomy,
    idx_opts: &IndexOptions,
    opts: &ClassifyOptions,
    scanner: &mut MinimizerScanner,
    taxa: &mut Vec<TaxId>,
    hit_counts: &mut HashMap<TaxId, u32>,
    tx_frames: &mut Vec<String>,
    curr_taxon_counts: &mut TaxonCounters,
    output_buf: &mut String,
) -> TaxId {
    taxa.clear();
    hit_counts.clear();
    output_buf.clear();
    let frame_ct = if opts.use_translated_search { 6 } else { 1 };
    let mut minimizer_hit_groups: i64 = 0;
    let mut call: TaxId = 0;
    let mut quick_exit = false;

    'mate_loop: for mate_num in 0..2 {
        if mate_num == 1 && !opts.paired_end_processing {
            break;
        }

        let seq_str = if mate_num == 0 {
            &dna.seq
        } else {
            &dna2.unwrap().seq
        };

        if opts.use_translated_search {
            *tx_frames = translate_to_all_frames(seq_str);
        }

        for frame_idx in 0..frame_ct {
            let scan_seq = if opts.use_translated_search {
                &tx_frames[frame_idx]
            } else {
                seq_str
            };

            scanner.load_sequence(scan_seq, 0, usize::MAX);
            let mut last_minimizer: u64 = u64::MAX;
            let mut last_taxon: TaxId = TAXID_MAX;

            while let Some(minimizer) = scanner.next_minimizer() {
                let taxon;
                if scanner.is_ambiguous() {
                    taxon = AMBIGUOUS_SPAN_TAXON;
                } else {
                    if minimizer != last_minimizer {
                        let mut skip_lookup = false;
                        if idx_opts.minimum_acceptable_hash_value != 0
                            && murmurhash3(minimizer) < idx_opts.minimum_acceptable_hash_value {
                                skip_lookup = true;
                            }
                        let t = if !skip_lookup {
                            hash.get(minimizer) as TaxId
                        } else {
                            0
                        };
                        last_taxon = t;
                        last_minimizer = minimizer;

                        if t != 0 {
                            minimizer_hit_groups += 1;
                            if !opts.report_filename.is_empty() {
                                curr_taxon_counts.entry(t).or_default().add_kmer(minimizer);
                            }
                        }
                    } else {
                        // Same minimizer as before
                    }
                    taxon = last_taxon;

                    if taxon != 0 {
                        if opts.quick_mode && minimizer_hit_groups >= opts.minimum_hit_groups {
                            call = taxon;
                            quick_exit = true;
                            break 'mate_loop;
                        }
                        *hit_counts.entry(taxon).or_insert(0) += 1;
                    }
                }
                taxa.push(taxon);
            }

            if opts.use_translated_search && frame_idx != 5 {
                taxa.push(READING_FRAME_BORDER_TAXON);
            }
        }

        if opts.paired_end_processing && mate_num == 0 {
            taxa.push(MATE_PAIR_BORDER_TAXON);
        }
    }

    // Calculate total k-mers (minus markers)
    let mut total_kmers = taxa.len();
    if opts.paired_end_processing {
        total_kmers = total_kmers.saturating_sub(1);
    }
    if opts.use_translated_search {
        let frame_markers = if opts.paired_end_processing { 4 } else { 2 };
        total_kmers = total_kmers.saturating_sub(frame_markers);
    }

    if !quick_exit {
        call = resolve_tree(hit_counts, taxonomy, total_kmers, opts.confidence_threshold);
    }

    // Void call if too few hit groups
    if call != 0 && minimizer_hit_groups < opts.minimum_hit_groups {
        call = 0;
    }

    if call != 0 && !opts.report_filename.is_empty() {
        curr_taxon_counts.entry(call).or_default().increment_read_count();
    }

    // Build Kraken output line into reusable buffer
    use std::fmt::Write as FmtWrite;
    let koss = output_buf;
    koss.push_str(if call != 0 { "C\t" } else { "U\t" });

    if !opts.paired_end_processing {
        koss.push_str(&dna.header);
    } else {
        koss.push_str(trim_pair_info(&dna.header));
    }
    koss.push('\t');

    let ext_call = taxonomy.node(call).external_id;
    if opts.print_scientific_name {
        if call != 0 {
            let name = taxonomy.name_at_offset(taxonomy.node(call).name_offset);
            let _ = write!(koss, "{} (taxid {})", name, ext_call);
        } else {
            let _ = write!(koss, "unclassified (taxid {})", ext_call);
        }
    } else {
        let _ = write!(koss, "{}", ext_call);
    }
    koss.push('\t');

    if !opts.paired_end_processing {
        let _ = write!(koss, "{}", dna.seq.len());
    } else {
        let len2 = dna2.map(|d| d.seq.len()).unwrap_or(0);
        let _ = write!(koss, "{}|{}", dna.seq.len(), len2);
    }
    koss.push('\t');

    if opts.quick_mode && quick_exit {
        let _ = write!(koss, "{}:Q", ext_call);
    } else if taxa.is_empty() {
        koss.push_str("0:0");
    } else {
        koss.push_str(&add_hitlist_string(taxa, taxonomy));
    }
    koss.push('\n');

    call
}

/// Helper to write classification output for a single sequence.
fn process_output(
    call: TaxId,
    kraken_str: &str,
    seq: &Sequence,
    taxonomy: &Taxonomy,
    kraken_out: &mut Option<BufWriter<File>>,
    classified_out: &mut Option<BufWriter<File>>,
    unclassified_out: &mut Option<BufWriter<File>>,
    stats: &mut ClassificationStats,
) -> io::Result<()> {
    if call != 0 {
        stats.total_classified += 1;
    } else {
        stats.total_unclassified += 1;
    }

    if let Some(ref mut out) = kraken_out {
        write!(out, "{}", kraken_str)?;
    }

    if call != 0 {
        if let Some(ref mut out) = classified_out {
            let ext = taxonomy.node(call).external_id;
            if seq.format == SequenceFormat::Fastq {
                writeln!(out, "@{} kraken:taxid|{}", seq.header, ext)?;
                writeln!(out, "{}", seq.seq)?;
                writeln!(out, "+")?;
                writeln!(out, "{}", seq.quals)?;
            } else {
                writeln!(out, ">{} kraken:taxid|{}", seq.header, ext)?;
                writeln!(out, "{}", seq.seq)?;
            }
        }
    } else if let Some(ref mut out) = unclassified_out {
        if seq.format == SequenceFormat::Fastq {
            writeln!(out, "@{}", seq.header)?;
            writeln!(out, "{}", seq.seq)?;
            writeln!(out, "+")?;
            writeln!(out, "{}", seq.quals)?;
        } else {
            writeln!(out, ">{}", seq.header)?;
            writeln!(out, "{}", seq.seq)?;
        }
    }
    Ok(())
}

/// Run the full classification pipeline on input files.
pub fn run_classify(
    input_files: &[String],
    hash_filename: &str,
    taxonomy_filename: &str,
    opts_filename: &str,
    opts: &ClassifyOptions,
    kraken_output: Option<&str>,
    classified_output: Option<&str>,
    unclassified_output: Option<&str>,
    report_filename: Option<&str>,
) -> io::Result<ClassificationStats> {
    // Load index options
    let idx_opts = IndexOptions::read_from_file(opts_filename)
        .expect("Failed to read index options");

    // Load taxonomy
    let mut taxonomy = Taxonomy::from_file(taxonomy_filename, opts.use_memory_mapping)?;
    taxonomy.generate_external_to_internal_id_map();

    // Load hash table
    let hash = CompactHashTable::from_file(hash_filename, opts.use_memory_mapping)?;

    let mut classify_opts = opts.clone();
    classify_opts.use_translated_search = !idx_opts.dna_db;
    if let Some(rf) = report_filename {
        classify_opts.report_filename = rf.to_string();
    }

    // Open output streams
    let mut kraken_out: Option<BufWriter<File>> = kraken_output
        .map(|f| BufWriter::new(File::create(f).expect("Failed to open kraken output")));
    let mut classified_out: Option<BufWriter<File>> = classified_output
        .map(|f| BufWriter::new(File::create(f).expect("Failed to open classified output")));
    let mut unclassified_out: Option<BufWriter<File>> = unclassified_output
        .map(|f| BufWriter::new(File::create(f).expect("Failed to open unclassified output")));

    let mut stats = ClassificationStats::default();
    let mut taxon_counters = TaxonCounters::new();

    let mut scanner = MinimizerScanner::new(
        idx_opts.k as isize,
        idx_opts.l as isize,
        idx_opts.spaced_seed_mask,
        idx_opts.dna_db,
        idx_opts.toggle_mask,
        idx_opts.revcom_version,
    );
    let mut taxa = Vec::new();
    let mut hit_counts = HashMap::new();
    let mut tx_frames = Vec::new();
    let mut output_buf = String::with_capacity(512);

    if classify_opts.paired_end_processing && !classify_opts.single_file_pairs && input_files.len() >= 2 {
        // Paired-end from two files
        let mut reader1 = BatchSequenceReader::new(Some(&input_files[0]))?;
        let mut reader2 = BatchSequenceReader::new(Some(&input_files[1]))?;

        while reader1.load_batch(10000) {
            reader2.load_batch(10000);
            loop {
                let seq1 = reader1.next_sequence().cloned();
                let seq2 = reader2.next_sequence().cloned();
                match (seq1, seq2) {
                    (Some(mut s1), Some(mut s2)) => {
                        stats.total_sequences += 1;
                        if classify_opts.minimum_quality_score > 0 {
                            mask_low_quality_bases(&mut s1, classify_opts.minimum_quality_score);
                            mask_low_quality_bases(&mut s2, classify_opts.minimum_quality_score);
                        }
                        let call = classify_sequence(
                            &s1, Some(&s2), &hash, &taxonomy, &idx_opts, &classify_opts,
                            &mut scanner, &mut taxa, &mut hit_counts, &mut tx_frames,
                            &mut taxon_counters, &mut output_buf,
                        );
                        process_output(call, &output_buf, &s1, &taxonomy,
                            &mut kraken_out, &mut classified_out, &mut unclassified_out,
                            &mut stats)?;
                    }
                    _ => break,
                }
            }
        }
    } else if classify_opts.paired_end_processing && classify_opts.single_file_pairs {
        // Paired-end from single interleaved file
        let filename = if input_files[0] == "-" { None } else { Some(input_files[0].as_str()) };
        let mut reader = BatchSequenceReader::new(filename)?;

        while reader.load_batch(20000) {
            loop {
                let seq1 = reader.next_sequence().cloned();
                let seq2 = reader.next_sequence().cloned();
                match (seq1, seq2) {
                    (Some(mut s1), Some(mut s2)) => {
                        stats.total_sequences += 1;
                        if classify_opts.minimum_quality_score > 0 {
                            mask_low_quality_bases(&mut s1, classify_opts.minimum_quality_score);
                            mask_low_quality_bases(&mut s2, classify_opts.minimum_quality_score);
                        }
                        let call = classify_sequence(
                            &s1, Some(&s2), &hash, &taxonomy, &idx_opts, &classify_opts,
                            &mut scanner, &mut taxa, &mut hit_counts, &mut tx_frames,
                            &mut taxon_counters, &mut output_buf,
                        );
                        process_output(call, &output_buf, &s1, &taxonomy,
                            &mut kraken_out, &mut classified_out, &mut unclassified_out,
                            &mut stats)?;
                    }
                    _ => break,
                }
            }
        }
    } else {
        // Single-end processing
        if classify_opts.num_threads > 1 {
            // Multi-threaded: C++-style parallel pipeline.
            // Each thread: lock reader -> read batch -> unlock -> classify -> queue output.
            // Output is flushed in block_id order to preserve input ordering.
            use std::sync::atomic::{AtomicU64, Ordering as AtomOrd};
            use std::collections::BinaryHeap;
            use std::cmp::Reverse;

            struct OutputBlock {
                block_id: u64,
                results: Vec<(TaxId, String, Sequence)>,
                local_counters: TaxonCounters,
                local_stats: ClassificationStats,
            }
            impl PartialEq for OutputBlock { fn eq(&self, o: &Self) -> bool { self.block_id == o.block_id } }
            impl Eq for OutputBlock {}
            impl PartialOrd for OutputBlock { fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> { Some(self.cmp(o)) } }
            impl Ord for OutputBlock { fn cmp(&self, o: &Self) -> std::cmp::Ordering { o.block_id.cmp(&self.block_id) } } // min-heap

            for input_file in input_files {
                let filename = if input_file == "-" { None } else { Some(input_file.as_str()) };
                let reader = parking_lot::Mutex::new(BatchSequenceReader::new(filename)?);
                let next_block_id = AtomicU64::new(0);
                let next_output_id = parking_lot::Mutex::new(0u64);
                let output_queue: parking_lot::Mutex<BinaryHeap<OutputBlock>> =
                    parking_lot::Mutex::new(BinaryHeap::new());

                // Wrap mutable output state in a mutex
                let output_state: parking_lot::Mutex<(
                    &mut Option<BufWriter<File>>,
                    &mut Option<BufWriter<File>>,
                    &mut Option<BufWriter<File>>,
                    &mut ClassificationStats,
                    &mut TaxonCounters,
                )> = parking_lot::Mutex::new((
                    &mut kraken_out, &mut classified_out, &mut unclassified_out,
                    &mut stats, &mut taxon_counters,
                ));

                crossbeam::scope(|s| {
                    for _ in 0..classify_opts.num_threads {
                        s.spawn(|_| {
                            let mut scanner = MinimizerScanner::new(
                                idx_opts.k as isize, idx_opts.l as isize,
                                idx_opts.spaced_seed_mask, idx_opts.dna_db,
                                idx_opts.toggle_mask, idx_opts.revcom_version,
                            );
                            let mut taxa = Vec::new();
                            let mut hit_counts = HashMap::new();
                            let mut tx_frames = Vec::new();
                            let mut local_buf = String::with_capacity(512);

                            loop {
                                // Read a batch (critical section)
                                let (batch, block_id) = {
                                    let mut rdr = reader.lock();
                                    if !rdr.load_block(3 * 1024 * 1024) {
                                        break;
                                    }
                                    let bid = next_block_id.fetch_add(1, AtomOrd::SeqCst);
                                    let mut batch = Vec::new();
                                    while let Some(seq) = rdr.next_sequence() {
                                        batch.push(seq.clone());
                                    }
                                    (batch, bid)
                                };

                                // Classify batch (no locks held)
                                let mut results = Vec::with_capacity(batch.len());
                                let mut local_counters = TaxonCounters::new();
                                let mut local_stats = ClassificationStats::default();

                                for mut seq in batch {
                                    if classify_opts.minimum_quality_score > 0 {
                                        mask_low_quality_bases(&mut seq, classify_opts.minimum_quality_score);
                                    }
                                    let call = classify_sequence(
                                        &seq, None, &hash, &taxonomy, &idx_opts, &classify_opts,
                                        &mut scanner, &mut taxa, &mut hit_counts, &mut tx_frames,
                                        &mut local_counters, &mut local_buf,
                                    );
                                    local_stats.total_sequences += 1;
                                    if call != 0 { local_stats.total_classified += 1; }
                                    else { local_stats.total_unclassified += 1; }
                                    results.push((call, local_buf.clone(), seq));
                                }

                                // Queue output block
                                let block = OutputBlock {
                                    block_id, results, local_counters, local_stats,
                                };
                                {
                                    let mut q = output_queue.lock();
                                    q.push(block);
                                }

                                // Try to flush output in order
                                loop {
                                    let mut q = output_queue.lock();
                                    let mut next_out = next_output_id.lock();
                                    if q.peek().map(|b| b.block_id) == Some(*next_out) {
                                        let blk = q.pop().unwrap();
                                        *next_out += 1;
                                        drop(q);
                                        drop(next_out);

                                        // Write output (holds output_state lock)
                                        let mut state = output_state.lock();
                                        let (ko, co, uo, st, tc) = &mut *state;
                                        st.total_sequences += blk.local_stats.total_sequences;
                                        st.total_classified += blk.local_stats.total_classified;
                                        st.total_unclassified += blk.local_stats.total_unclassified;
                                        for (taxid, counter) in blk.local_counters {
                                            tc.entry(taxid).or_default().merge(&counter);
                                        }
                                        for (call, kraken_str, seq) in &blk.results {
                                            let _ = process_output(*call, kraken_str, seq, &taxonomy,
                                                ko, co, uo, &mut ClassificationStats::default());
                                        }
                                    } else {
                                        break;
                                    }
                                }
                            }
                        });
                    }
                }).expect("Thread panicked");
            }
        } else {
            // Single-threaded processing
            for input_file in input_files {
                let filename = if input_file == "-" { None } else { Some(input_file.as_str()) };
                let mut reader = BatchSequenceReader::new(filename)?;

                while reader.load_block(3 * 1024 * 1024) {
                    while let Some(seq) = reader.next_sequence() {
                        let mut seq = seq.clone();
                        stats.total_sequences += 1;

                        if classify_opts.minimum_quality_score > 0 {
                            mask_low_quality_bases(&mut seq, classify_opts.minimum_quality_score);
                        }

                        let call = classify_sequence(
                            &seq, None, &hash, &taxonomy, &idx_opts, &classify_opts,
                            &mut scanner, &mut taxa, &mut hit_counts, &mut tx_frames,
                            &mut taxon_counters, &mut output_buf,
                        );
                        process_output(call, &output_buf, &seq, &taxonomy,
                            &mut kraken_out, &mut classified_out, &mut unclassified_out,
                            &mut stats)?;
                    }
                }
            }
        }
    }

    // Write report if requested
    if let Some(rf) = report_filename {
        if classify_opts.mpa_style_report {
            reports::report_mpa_style(rf, classify_opts.report_zero_counts,
                &taxonomy, &taxon_counters)?;
        } else {
            reports::report_kraken_style(rf, classify_opts.report_zero_counts,
                classify_opts.report_kmer_data, &taxonomy, &taxon_counters,
                stats.total_sequences, stats.total_unclassified)?;
        }
    }

    Ok(stats)
}

/// Result of classifying a single sequence.
#[derive(Debug, Clone)]
pub struct ClassificationResult {
    /// Taxonomy call (internal ID, 0 = unclassified).
    pub call: TaxId,
    /// External (NCBI) taxonomy ID.
    pub external_id: u64,
    /// Whether the sequence was classified.
    pub classified: bool,
    /// Kraken-format output line.
    pub kraken_line: String,
}

/// Classify a batch of sequences in memory, without file I/O.
///
/// This is the primary library API for classification. It takes loaded database
/// components and a slice of sequences, and returns classification results.
///
/// # Example
/// ```no_run
/// use kraken2::classify::{ClassifyOptions, classify_sequences_in_memory};
/// use kraken2::compact_hash::CompactHashTable;
/// use kraken2::taxonomy::Taxonomy;
/// use kraken2::types::{IndexOptions, Sequence};
///
/// let idx_opts = IndexOptions::read_from_file("db/opts.k2d").unwrap();
/// let mut taxonomy = Taxonomy::from_file("db/taxo.k2d", false).unwrap();
/// taxonomy.generate_external_to_internal_id_map();
/// let hash = CompactHashTable::from_file("db/hash.k2d", false).unwrap();
///
/// let seq = Sequence {
///     header: "read1".to_string(),
///     seq: "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_string(),
///     ..Default::default()
/// };
///
/// let opts = ClassifyOptions::default();
/// let results = classify_sequences_in_memory(
///     &[seq], &hash, &taxonomy, &idx_opts, &opts,
/// );
/// for r in &results {
///     println!("{}: taxid={}", if r.classified { "C" } else { "U" }, r.external_id);
/// }
/// ```
pub fn classify_sequences_in_memory(
    sequences: &[Sequence],
    hash: &CompactHashTable,
    taxonomy: &Taxonomy,
    idx_opts: &IndexOptions,
    opts: &ClassifyOptions,
) -> Vec<ClassificationResult> {
    let mut scanner = MinimizerScanner::new(
        idx_opts.k as isize,
        idx_opts.l as isize,
        idx_opts.spaced_seed_mask,
        idx_opts.dna_db,
        idx_opts.toggle_mask,
        idx_opts.revcom_version,
    );
    let mut taxa = Vec::new();
    let mut hit_counts = HashMap::new();
    let mut tx_frames = Vec::new();
    let mut taxon_counters = TaxonCounters::new();
    let mut output_buf = String::with_capacity(512);

    let mut results = Vec::with_capacity(sequences.len());

    for seq in sequences {
        let call = classify_sequence(
            seq, None, hash, taxonomy, idx_opts, opts,
            &mut scanner, &mut taxa, &mut hit_counts, &mut tx_frames,
            &mut taxon_counters, &mut output_buf,
        );

        let external_id = taxonomy.node(call).external_id;
        results.push(ClassificationResult {
            call,
            external_id,
            classified: call != 0,
            kraken_line: output_buf.clone(),
        });
    }

    results
}

/// A loaded Kraken 2 database ready for classification.
///
/// Load once, classify many times. This is the recommended high-level API.
///
/// # Example
/// ```no_run
/// use kraken2::classify::{ClassifyDb, ClassifyOptions};
/// use kraken2::types::Sequence;
///
/// let db = ClassifyDb::from_directory("path/to/db").unwrap();
/// let opts = ClassifyOptions::default();
///
/// // Classify one sequence
/// let seq = Sequence {
///     header: "read1".to_string(),
///     seq: "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_string(),
///     ..Default::default()
/// };
/// let results = db.classify(&[seq], &opts);
///
/// // Classify another batch with the same loaded database
/// let results2 = db.classify(&[/* more sequences */], &opts);
/// ```
pub struct ClassifyDb {
    pub idx_opts: IndexOptions,
    pub taxonomy: Taxonomy,
    pub hash: CompactHashTable,
}

impl ClassifyDb {
    /// Load a database from a directory containing hash.k2d, taxo.k2d, opts.k2d.
    pub fn from_directory(db_dir: &str) -> io::Result<Self> {
        Self::from_files(
            &format!("{}/hash.k2d", db_dir),
            &format!("{}/taxo.k2d", db_dir),
            &format!("{}/opts.k2d", db_dir),
            false,
        )
    }

    /// Load a database from individual file paths.
    pub fn from_files(
        hash_path: &str,
        taxonomy_path: &str,
        opts_path: &str,
        memory_mapping: bool,
    ) -> io::Result<Self> {
        let idx_opts = IndexOptions::read_from_file(opts_path)?;
        let mut taxonomy = Taxonomy::from_file(taxonomy_path, memory_mapping)?;
        taxonomy.generate_external_to_internal_id_map();
        let hash = CompactHashTable::from_file(hash_path, memory_mapping)?;
        Ok(ClassifyDb { idx_opts, taxonomy, hash })
    }

    /// Classify a batch of sequences. Can be called repeatedly.
    pub fn classify(
        &self,
        sequences: &[Sequence],
        opts: &ClassifyOptions,
    ) -> Vec<ClassificationResult> {
        classify_sequences_in_memory(sequences, &self.hash, &self.taxonomy, &self.idx_opts, opts)
    }

    /// Classify a single sequence. Convenience wrapper.
    pub fn classify_one(
        &self,
        seq: &Sequence,
        opts: &ClassifyOptions,
    ) -> ClassificationResult {
        self.classify(std::slice::from_ref(seq), opts).into_iter().next().unwrap()
    }
}

/// Run classification in daemon mode.
/// The daemon forks into the background, creates FIFOs for IPC,
/// and loops processing classification requests until STOP is received.
/// Port of C++ `ClassifyDaemon()`.
#[cfg(unix)]
pub fn run_daemon(
    hash_filename: &str,
    taxonomy_filename: &str,
    opts_filename: &str,
    opts: &ClassifyOptions,
) -> io::Result<()> {
    use std::io::BufRead;

    // Load database once
    let idx_opts = IndexOptions::read_from_file(opts_filename)?;
    let mut taxonomy = Taxonomy::from_file(taxonomy_filename, opts.use_memory_mapping)?;
    taxonomy.generate_external_to_internal_id_map();
    let hash = std::sync::Arc::new(CompactHashTable::from_file(hash_filename, opts.use_memory_mapping)?);

    eprintln!("Database loaded. Daemon ready.");
    eprintln!("PID: {}", std::process::id());

    // Write PID file
    let mut pid_file = File::create("/tmp/classify.pid")?;
    writeln!(pid_file, "{}", std::process::id())?;

    let stdin = io::stdin();
    let mut line = String::new();

    loop {
        line.clear();
        if stdin.lock().read_line(&mut line)? == 0 {
            break;
        }
        let trimmed = line.trim();

        if trimmed == "PING" {
            eprintln!("OK");
            continue;
        }
        if trimmed == "STOP" {
            eprintln!("OK");
            break;
        }

        // Treat the line as a list of input files to classify
        let files: Vec<String> = trimmed.split_whitespace()
            .map(|s| s.to_string())
            .collect();

        if files.is_empty() {
            continue;
        }

        let mut classify_opts = opts.clone();
        classify_opts.use_translated_search = !idx_opts.dna_db;

        let mut stats = ClassificationStats::default();
        let mut taxon_counters = TaxonCounters::new();

        let mut scanner = MinimizerScanner::new(
            idx_opts.k as isize, idx_opts.l as isize,
            idx_opts.spaced_seed_mask, idx_opts.dna_db,
            idx_opts.toggle_mask, idx_opts.revcom_version,
        );
        let mut taxa = Vec::new();
        let mut hit_counts = HashMap::new();
        let mut tx_frames = Vec::new();
        let mut output_buf = String::with_capacity(512);

        for input_file in &files {
            let filename = if input_file == "-" { None } else { Some(input_file.as_str()) };
            let mut reader = crate::seq::BatchSequenceReader::new(filename)?;

            while reader.load_block(3 * 1024 * 1024) {
                while let Some(seq) = reader.next_sequence() {
                    let seq = seq.clone();
                    stats.total_sequences += 1;

                    let call = classify_sequence(
                        &seq, None, &hash, &taxonomy, &idx_opts, &classify_opts,
                        &mut scanner, &mut taxa, &mut hit_counts, &mut tx_frames,
                        &mut taxon_counters, &mut output_buf,
                    );

                    if call != 0 {
                        stats.total_classified += 1;
                    }
                    print!("{}", output_buf);
                }
            }
        }

        eprintln!(
            "{} sequences ({:.2}% classified)",
            stats.total_sequences,
            if stats.total_sequences > 0 {
                100.0 * stats.total_classified as f64 / stats.total_sequences as f64
            } else {
                0.0
            }
        );
        println!("DONE");
    }

    // Cleanup
    let _ = std::fs::remove_file("/tmp/classify.pid");
    Ok(())
}

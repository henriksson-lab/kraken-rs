use ahash::AHashMap as HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::time::{Duration, Instant};

#[cfg(unix)]
use std::ffi::CString;
#[cfg(unix)]
use std::os::fd::RawFd;

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

type IndexData = (IndexOptions, Taxonomy, CompactHashTable);

#[derive(Clone, Default)]
struct Options {
    index_filename: String,
    taxonomy_filename: String,
    options_filename: String,
    report_filename: String,
    classified_output_filename: String,
    unclassified_output_filename: String,
    kraken_output_filename: String,
    mpa_style_report: bool,
    report_kmer_data: bool,
    quick_mode: bool,
    report_zero_counts: bool,
    use_translated_search: bool,
    print_scientific_name: bool,
    confidence_threshold: f64,
    num_threads: usize,
    paired_end_processing: bool,
    single_file_pairs: bool,
    minimum_quality_score: u8,
    minimum_hit_groups: i64,
    use_memory_mapping: bool,
    match_input_order: bool,
    filenames: Vec<String>,
    daemon_mode: bool,
}

impl Options {
    fn reset(&mut self) {
        self.quick_mode = false;
        self.confidence_threshold = 0.0;
        self.paired_end_processing = false;
        self.single_file_pairs = false;
        self.num_threads = 1;
        self.mpa_style_report = false;
        self.report_kmer_data = false;
        self.report_zero_counts = false;
        self.use_translated_search = false;
        self.print_scientific_name = false;
        self.minimum_quality_score = 0;
        self.minimum_hit_groups = 0;
        self.use_memory_mapping = false;
        self.daemon_mode = false;
        self.match_input_order = false;

        self.index_filename.clear();
        self.taxonomy_filename.clear();
        self.options_filename.clear();
        self.report_filename.clear();
        self.classified_output_filename.clear();
        self.unclassified_output_filename.clear();
        self.kraken_output_filename.clear();
        self.filenames.clear();
    }
}

#[derive(Default)]
struct OutputStreamData {
    initialized: bool,
    printing_sequences: bool,
    classified_output1: Option<BufWriter<File>>,
    classified_output2: Option<BufWriter<File>>,
    unclassified_output1: Option<BufWriter<File>>,
    unclassified_output2: Option<BufWriter<File>>,
    kraken_output: Option<BufWriter<File>>,
}

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
    pub total_bases: u64,
    pub total_classified: u64,
    pub total_unclassified: u64,
}

fn load_index(opts: &mut Options) -> io::Result<IndexData> {
    let idx_opts = IndexOptions::read_from_file(&opts.options_filename)?;
    opts.use_translated_search = !idx_opts.dna_db;
    let mut taxonomy = Taxonomy::from_file(&opts.taxonomy_filename, opts.use_memory_mapping)?;
    taxonomy.generate_external_to_internal_id_map();
    let hash = CompactHashTable::from_file(&opts.index_filename, opts.use_memory_mapping)?;
    Ok((idx_opts, taxonomy, hash))
}

fn tokenize_string(argv: &mut Vec<String>, s: &str) {
    argv.clear();
    argv.extend(
        s.split(' ')
            .filter(|t| !t.is_empty())
            .map(|t| t.to_string()),
    );
}

fn report_stats(start: Duration, end: Duration, stats: &ClassificationStats) -> String {
    let seconds = end.saturating_sub(start).as_secs_f64();
    let total_unclassified = stats.total_sequences.saturating_sub(stats.total_classified);
    let mbp = stats.total_bases as f64 / 1.0e6;
    let kseq_per_m = if seconds > 0.0 {
        stats.total_sequences as f64 / 1.0e3 / (seconds / 60.0)
    } else {
        0.0
    };
    let mbp_per_m = if seconds > 0.0 {
        mbp / (seconds / 60.0)
    } else {
        0.0
    };
    let classified_pct = if stats.total_sequences > 0 {
        stats.total_classified as f64 * 100.0 / stats.total_sequences as f64
    } else {
        0.0
    };
    let unclassified_pct = if stats.total_sequences > 0 {
        total_unclassified as f64 * 100.0 / stats.total_sequences as f64
    } else {
        0.0
    };

    format!(
        "{seq} sequences ({mbp:.2} Mbp) processed in {seconds:.3}s ({kseq:.1} Kseq/m, {mbp_m:.2} Mbp/m).\n  {classed} sequences classified ({classed_pct:.2}%)\n  {unclassed} sequences unclassified ({unclassed_pct:.2}%)\n",
        seq = stats.total_sequences,
        kseq = kseq_per_m,
        mbp_m = mbp_per_m,
        classed = stats.total_classified,
        classed_pct = classified_pct,
        unclassed = total_unclassified,
        unclassed_pct = unclassified_pct,
    )
}

fn open_output_stream(filename: &str) -> io::Result<BufWriter<File>> {
    File::create(filename).map(BufWriter::new)
}

fn initialize_outputs(
    opts: &Options,
    outputs: &mut OutputStreamData,
    _format: SequenceFormat,
) -> io::Result<()> {
    if outputs.initialized {
        return Ok(());
    }

    if !opts.classified_output_filename.is_empty() {
        if opts.paired_end_processing {
            let fields = crate::utilities::split_string(&opts.classified_output_filename, "#", 3);
            if fields.len() < 2 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Paired filename format missing # character: {}",
                        opts.classified_output_filename
                    ),
                ));
            } else if fields.len() > 2 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Paired filename format has >1 # character: {}",
                        opts.classified_output_filename
                    ),
                ));
            }
            outputs.classified_output1 = Some(open_output_stream(
                &(fields[0].clone() + "_1" + &fields[1]),
            )?);
            outputs.classified_output2 = Some(open_output_stream(
                &(fields[0].clone() + "_2" + &fields[1]),
            )?);
        } else {
            outputs.classified_output1 =
                Some(open_output_stream(&opts.classified_output_filename)?);
        }
        outputs.printing_sequences = true;
    }

    if !opts.unclassified_output_filename.is_empty() {
        if opts.paired_end_processing {
            let fields = crate::utilities::split_string(&opts.unclassified_output_filename, "#", 3);
            if fields.len() < 2 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Paired filename format missing # character: {}",
                        opts.unclassified_output_filename
                    ),
                ));
            } else if fields.len() > 2 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Paired filename format has >1 # character: {}",
                        opts.unclassified_output_filename
                    ),
                ));
            }
            outputs.unclassified_output1 = Some(open_output_stream(
                &(fields[0].clone() + "_1" + &fields[1]),
            )?);
            outputs.unclassified_output2 = Some(open_output_stream(
                &(fields[0].clone() + "_2" + &fields[1]),
            )?);
        } else {
            outputs.unclassified_output1 =
                Some(open_output_stream(&opts.unclassified_output_filename)?);
        }
        outputs.printing_sequences = true;
    }

    if !opts.kraken_output_filename.is_empty() && opts.kraken_output_filename != "-" {
        outputs.kraken_output = Some(open_output_stream(&opts.kraken_output_filename)?);
    }

    outputs.initialized = true;
    Ok(())
}

fn parse_command_line(args: &[String], opts: &mut Options) -> io::Result<()> {
    opts.reset();
    let mut i = 1usize;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "-?" => return Err(io::Error::new(io::ErrorKind::InvalidInput, usage(0))),
            "-H" => {
                i += 1;
                opts.index_filename = args.get(i).cloned().unwrap_or_default();
            }
            "-t" => {
                i += 1;
                opts.taxonomy_filename = args.get(i).cloned().unwrap_or_default();
            }
            "-o" => {
                i += 1;
                opts.options_filename = args.get(i).cloned().unwrap_or_default();
            }
            "-T" => {
                i += 1;
                opts.confidence_threshold =
                    args.get(i).and_then(|s| s.parse().ok()).unwrap_or(-1.0);
                if !(0.0..=1.0).contains(&opts.confidence_threshold) {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "confidence threshold must be in [0, 1]",
                    ));
                }
            }
            "-q" => opts.quick_mode = true,
            "-p" => {
                i += 1;
                opts.num_threads = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(0);
                if opts.num_threads < 1 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "number of threads can't be less than 1",
                    ));
                }
            }
            "-g" => {
                i += 1;
                opts.minimum_hit_groups = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(0);
            }
            "-P" => opts.paired_end_processing = true,
            "-S" => {
                opts.paired_end_processing = true;
                opts.single_file_pairs = true;
            }
            "-m" => opts.mpa_style_report = true,
            "-K" => opts.report_kmer_data = true,
            "-R" => {
                i += 1;
                opts.report_filename = args.get(i).cloned().unwrap_or_default();
            }
            "-z" => opts.report_zero_counts = true,
            "-C" => {
                i += 1;
                opts.classified_output_filename = args.get(i).cloned().unwrap_or_default();
            }
            "-U" => {
                i += 1;
                opts.unclassified_output_filename = args.get(i).cloned().unwrap_or_default();
            }
            "-O" => {
                i += 1;
                opts.kraken_output_filename = args.get(i).cloned().unwrap_or_default();
            }
            "-n" => opts.print_scientific_name = true,
            "-Q" => {
                i += 1;
                opts.minimum_quality_score = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(0);
            }
            "-M" => opts.use_memory_mapping = true,
            "-D" => opts.daemon_mode = true,
            other if other.starts_with('-') => {}
            other => opts.filenames.push(other.to_string()),
        }
        i += 1;
    }

    if opts.index_filename.is_empty()
        || opts.taxonomy_filename.is_empty()
        || opts.options_filename.is_empty()
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "mandatory filename missing",
        ));
    }
    if opts.mpa_style_report && opts.report_filename.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "-m requires -R be used",
        ));
    }
    Ok(())
}

fn usage(_exit_code: i32) -> String {
    [
        "Usage: classify [options] <fasta/fastq file(s)>",
        "",
        "Options: (*mandatory)",
        "* -H filename      Kraken 2 index filename",
        "* -t filename      Kraken 2 taxonomy filename",
        "* -o filename      Kraken 2 options filename",
        "  -q               Quick mode",
        "  -M               Use memory mapping to access hash & taxonomy",
        "  -T NUM           Confidence score threshold (def. 0)",
        "  -p NUM           Number of threads (def. 1)",
        "  -Q NUM           Minimum quality score (FASTQ only, def. 0)",
        "  -P               Process pairs of reads",
        "  -S               Process pairs with mates in same file",
        "  -R filename      Print report to filename",
        "  -m               In comb. w/ -R, use mpa-style report",
        "  -z               In comb. w/ -R, report taxa w/ 0 count",
        "  -n               Print scientific name instead of taxid in Kraken output",
        "  -g NUM           Minimum number of hit groups needed for call",
        "  -C filename      Filename/format to have classified sequences",
        "  -U filename      Filename/format to have unclassified sequences",
        "  -O filename      Output file for normal Kraken output",
        "  -K               In comb. w/ -R, provide minimizer information in report",
        "  -D               Start a daemon, this options is intended to be used with wrappers",
    ]
    .join("\n")
}

#[cfg(unix)]
fn remove_blocking(fd: RawFd) -> io::Result<()> {
    let flags = unsafe { libc::fcntl(fd, libc::F_GETFL, 0) };
    if flags < 0 {
        return Err(io::Error::last_os_error());
    }
    if unsafe { libc::fcntl(fd, libc::F_SETFL, flags & !libc::O_NONBLOCK) } < 0 {
        return Err(io::Error::last_os_error());
    }
    Ok(())
}

#[cfg(unix)]
fn mkfifo_if_needed(path: &str) -> io::Result<()> {
    let c_path = CString::new(path).unwrap();
    let rc = unsafe { libc::mkfifo(c_path.as_ptr(), libc::S_IRWXU) };
    if rc == 0 || io::Error::last_os_error().kind() == io::ErrorKind::AlreadyExists {
        Ok(())
    } else {
        Err(io::Error::last_os_error())
    }
}

#[cfg(unix)]
fn open_fd(path: &CString, flags: i32) -> io::Result<RawFd> {
    let fd = unsafe { libc::open(path.as_ptr(), flags) };
    if fd < 0 {
        Err(io::Error::last_os_error())
    } else {
        Ok(fd)
    }
}

#[cfg(unix)]
fn dup2_checked(src: RawFd, dst: RawFd) -> io::Result<()> {
    if unsafe { libc::dup2(src, dst) } < 0 {
        Err(io::Error::last_os_error())
    } else {
        Ok(())
    }
}

#[cfg(unix)]
fn daemonize() -> io::Result<i32> {
    let mut pid = unsafe { libc::fork() };
    if pid < 0 {
        return Err(io::Error::last_os_error());
    }

    if pid == 0 {
        if unsafe { libc::setsid() } == -1 {
            return Err(io::Error::last_os_error());
        }
    } else {
        std::process::exit(0);
    }

    pid = unsafe { libc::fork() };
    if pid < 0 {
        return Err(io::Error::last_os_error());
    }
    if pid != 0 {
        std::process::exit(0);
    }

    mkfifo_if_needed("/tmp/classify_stdin")?;
    mkfifo_if_needed("/tmp/classify_stdout")?;

    let stdin_name = CString::new("/tmp/classify_stdin").unwrap();
    let stdout_name = CString::new("/tmp/classify_stdout").unwrap();
    let read_fd = unsafe { libc::open(stdin_name.as_ptr(), libc::O_RDONLY | libc::O_NONBLOCK) };
    if read_fd < 0 {
        return Err(io::Error::last_os_error());
    }
    let dummy_fd_1 = unsafe { libc::open(stdin_name.as_ptr(), libc::O_WRONLY) };
    if dummy_fd_1 < 0 {
        return Err(io::Error::last_os_error());
    }
    let dummy_fd_2 = unsafe { libc::open(stdout_name.as_ptr(), libc::O_RDONLY | libc::O_NONBLOCK) };
    if dummy_fd_2 < 0 {
        return Err(io::Error::last_os_error());
    }
    let write_fd = unsafe { libc::open(stdout_name.as_ptr(), libc::O_WRONLY) };
    if write_fd < 0 {
        return Err(io::Error::last_os_error());
    }

    remove_blocking(read_fd)?;
    remove_blocking(dummy_fd_2)?;

    for fd in 0..2 {
        unsafe {
            libc::close(fd);
        }
    }
    if unsafe { libc::dup2(read_fd, 0) } < 0
        || unsafe { libc::dup2(write_fd, 1) } < 0
        || unsafe { libc::dup2(write_fd, 2) } < 0
    {
        return Err(io::Error::last_os_error());
    }

    Ok(0)
}

#[cfg(unix)]
fn open_fifos(opts: &Options, pid: libc::pid_t) -> io::Result<()> {
    let stdin_filename = format!("/tmp/classify_{}_stdin", pid);
    let stdout_filename = format!("/tmp/classify_{}_stdout", pid);

    mkfifo_if_needed(&stdin_filename)?;
    mkfifo_if_needed(&stdout_filename)?;

    let stdin_c = CString::new(stdin_filename.clone()).unwrap();
    let stdout_c = CString::new(stdout_filename.clone()).unwrap();

    let read_fd;
    // If we are expecting input from stdin, open the FIFO in blocking mode and
    // wait for the input.
    if opts.filenames.is_empty() {
        read_fd = open_fd(&stdin_c, libc::O_RDONLY)?;
    } else {
        read_fd = open_fd(&stdin_c, libc::O_RDONLY | libc::O_NONBLOCK)?;
        // Keep this open so the daemon does not receive EOF when the external
        // process closes its end of the fifo.
        let _dummy_fd_1 = open_fd(&stdin_c, libc::O_WRONLY)?;
    }

    let mut dummy_fd_2 = -1;
    // If we are outputting to a file, open the write FIFO in non-blocking mode
    // and keep a read end open so that we do not block the process.
    if !opts.kraken_output_filename.is_empty() {
        dummy_fd_2 = open_fd(&stdout_c, libc::O_RDONLY | libc::O_NONBLOCK)?;
    }
    let write_fd = open_fd(&stdout_c, libc::O_WRONLY)?;

    if opts.kraken_output_filename.is_empty() {
        remove_blocking(read_fd)?;
    } else {
        remove_blocking(dummy_fd_2)?;
    }

    for fd in 0..3 {
        unsafe {
            libc::close(fd);
        }
    }
    dup2_checked(read_fd, 0)?;
    dup2_checked(write_fd, 1)?;
    dup2_checked(write_fd, 2)?;
    Ok(())
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
    let mut result = String::new();
    append_hitlist_string(taxa, taxonomy, &mut result);
    result
}

/// Same as [`add_hitlist_string`] but appends directly to an existing buffer,
/// avoiding the per-read String allocation.
pub fn append_hitlist_string(taxa: &[TaxId], taxonomy: &Taxonomy, out: &mut String) {
    if taxa.is_empty() {
        out.push_str("0:0");
        return;
    }

    let mut last_code = taxa[0];
    let mut code_count = 1u32;

    for i in 1..taxa.len() {
        let code = taxa[i];
        if code == last_code {
            code_count += 1;
        } else {
            append_hitlist_entry(out, last_code, code_count, taxonomy, true);
            code_count = 1;
            last_code = code;
        }
    }
    append_hitlist_entry(out, last_code, code_count, taxonomy, false);
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
                            && murmurhash3(minimizer) < idx_opts.minimum_acceptable_hash_value
                        {
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

    call = resolve_tree(hit_counts, taxonomy, total_kmers, opts.confidence_threshold);

    // Void call if too few hit groups
    if call != 0 && minimizer_hit_groups < opts.minimum_hit_groups {
        call = 0;
    }

    if call != 0 && !opts.report_filename.is_empty() {
        curr_taxon_counts
            .entry(call)
            .or_default()
            .increment_read_count();
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

    if opts.quick_mode {
        let _ = write!(koss, "{}:Q", ext_call);
    } else {
        append_hitlist_string(taxa, taxonomy, koss);
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

struct OutputHandles<'a> {
    kraken_output: &'a mut Option<BufWriter<File>>,
    classified_output: &'a mut Option<BufWriter<File>>,
    unclassified_output: &'a mut Option<BufWriter<File>>,
}

enum ClassificationInput {
    Single(Sequence),
    Pair(Sequence, Sequence),
}

fn process_files(
    filename1: Option<&str>,
    filename2: Option<&str>,
    hash: &CompactHashTable,
    taxonomy: &Taxonomy,
    idx_opts: &IndexOptions,
    classify_opts: &ClassifyOptions,
    stats: &mut ClassificationStats,
    outputs: &mut OutputHandles<'_>,
    taxon_counters: &mut TaxonCounters,
) -> io::Result<()> {
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

    if classify_opts.num_threads > 1 {
        use std::collections::BinaryHeap;
        use std::sync::atomic::{AtomicU64, Ordering as AtomOrd};

        struct OutputBlock {
            block_id: u64,
            results: Vec<(TaxId, String, Sequence)>,
            local_counters: TaxonCounters,
            local_stats: ClassificationStats,
        }
        impl PartialEq for OutputBlock {
            fn eq(&self, o: &Self) -> bool {
                self.block_id == o.block_id
            }
        }
        impl Eq for OutputBlock {}
        impl PartialOrd for OutputBlock {
            fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
                Some(self.cmp(o))
            }
        }
        impl Ord for OutputBlock {
            fn cmp(&self, o: &Self) -> std::cmp::Ordering {
                o.block_id.cmp(&self.block_id)
            }
        }

        let next_block_id = AtomicU64::new(0);
        let next_output_id = parking_lot::Mutex::new(0u64);
        let output_queue: parking_lot::Mutex<BinaryHeap<OutputBlock>> =
            parking_lot::Mutex::new(BinaryHeap::new());

        #[allow(clippy::type_complexity)]
        let output_state: parking_lot::Mutex<(
            &mut Option<BufWriter<File>>,
            &mut Option<BufWriter<File>>,
            &mut Option<BufWriter<File>>,
            &mut ClassificationStats,
            &mut TaxonCounters,
        )> = parking_lot::Mutex::new((
            outputs.kraken_output,
            outputs.classified_output,
            outputs.unclassified_output,
            stats,
            taxon_counters,
        ));

        if classify_opts.paired_end_processing && !classify_opts.single_file_pairs {
            let readers = parking_lot::Mutex::new((
                BatchSequenceReader::new(filename1)?,
                BatchSequenceReader::new(filename2)?,
            ));

            crossbeam::scope(|s| {
                for _ in 0..classify_opts.num_threads {
                    s.spawn(|_| {
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
                        let mut local_buf = String::with_capacity(512);

                        loop {
                            let (batch, block_id) = {
                                let mut readers = readers.lock();
                                if !readers.0.load_batch(10000) {
                                    break;
                                }
                                readers.1.load_batch(10000);
                                let bid = next_block_id.fetch_add(1, AtomOrd::SeqCst);
                                let mut batch = Vec::new();
                                loop {
                                    let seq1 = readers.0.next_sequence().cloned();
                                    let seq2 = readers.1.next_sequence().cloned();
                                    match (seq1, seq2) {
                                        (Some(s1), Some(s2)) => {
                                            batch.push(ClassificationInput::Pair(s1, s2));
                                        }
                                        _ => break,
                                    }
                                }
                                (batch, bid)
                            };

                            let mut results = Vec::with_capacity(batch.len());
                            let mut local_counters = TaxonCounters::new();
                            let mut local_stats = ClassificationStats::default();

                            for input in batch {
                                if let ClassificationInput::Pair(mut s1, mut s2) = input {
                                    if classify_opts.minimum_quality_score > 0 {
                                        mask_low_quality_bases(
                                            &mut s1,
                                            classify_opts.minimum_quality_score,
                                        );
                                        mask_low_quality_bases(
                                            &mut s2,
                                            classify_opts.minimum_quality_score,
                                        );
                                    }
                                    let call = classify_sequence(
                                        &s1,
                                        Some(&s2),
                                        hash,
                                        taxonomy,
                                        idx_opts,
                                        classify_opts,
                                        &mut scanner,
                                        &mut taxa,
                                        &mut hit_counts,
                                        &mut tx_frames,
                                        &mut local_counters,
                                        &mut local_buf,
                                    );
                                    local_stats.total_sequences += 1;
                                    local_stats.total_bases +=
                                        (s1.seq.len() + s2.seq.len()) as u64;
                                    if call != 0 {
                                        local_stats.total_classified += 1;
                                    } else {
                                        local_stats.total_unclassified += 1;
                                    }
                                    results.push((call, local_buf.clone(), s1));
                                }
                            }

                            let block = OutputBlock {
                                block_id,
                                results,
                                local_counters,
                                local_stats,
                            };
                            {
                                let mut q = output_queue.lock();
                                q.push(block);
                            }

                            loop {
                                let mut q = output_queue.lock();
                                let mut next_out = next_output_id.lock();
                                if q.peek().map(|b| b.block_id) == Some(*next_out) {
                                    let blk = q.pop().unwrap();
                                    *next_out += 1;
                                    drop(q);
                                    drop(next_out);

                                    let mut state = output_state.lock();
                                    let (ko, co, uo, st, tc) = &mut *state;
                                    st.total_sequences += blk.local_stats.total_sequences;
                                    st.total_bases += blk.local_stats.total_bases;
                                    st.total_classified += blk.local_stats.total_classified;
                                    st.total_unclassified += blk.local_stats.total_unclassified;
                                    for (taxid, counter) in blk.local_counters {
                                        tc.entry(taxid).or_default().merge(&counter);
                                    }
                                    for (call, kraken_str, seq) in &blk.results {
                                        let _ = process_output(
                                            *call,
                                            kraken_str,
                                            seq,
                                            taxonomy,
                                            ko,
                                            co,
                                            uo,
                                            &mut ClassificationStats::default(),
                                        );
                                    }
                                } else {
                                    break;
                                }
                            }
                        }
                    });
                }
            })
            .expect("Thread panicked");
        } else if classify_opts.paired_end_processing && classify_opts.single_file_pairs {
            let reader = parking_lot::Mutex::new(BatchSequenceReader::new(filename1)?);

            crossbeam::scope(|s| {
                for _ in 0..classify_opts.num_threads {
                    s.spawn(|_| {
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
                        let mut local_buf = String::with_capacity(512);

                        loop {
                            let (batch, block_id) = {
                                let mut rdr = reader.lock();
                                if !rdr.load_batch(20000) {
                                    break;
                                }
                                let bid = next_block_id.fetch_add(1, AtomOrd::SeqCst);
                                let mut batch = Vec::new();
                                loop {
                                    let seq1 = rdr.next_sequence().cloned();
                                    let seq2 = rdr.next_sequence().cloned();
                                    match (seq1, seq2) {
                                        (Some(s1), Some(s2)) => {
                                            batch.push(ClassificationInput::Pair(s1, s2));
                                        }
                                        _ => break,
                                    }
                                }
                                (batch, bid)
                            };

                            let mut results = Vec::with_capacity(batch.len());
                            let mut local_counters = TaxonCounters::new();
                            let mut local_stats = ClassificationStats::default();

                            for input in batch {
                                if let ClassificationInput::Pair(mut s1, mut s2) = input {
                                    if classify_opts.minimum_quality_score > 0 {
                                        mask_low_quality_bases(
                                            &mut s1,
                                            classify_opts.minimum_quality_score,
                                        );
                                        mask_low_quality_bases(
                                            &mut s2,
                                            classify_opts.minimum_quality_score,
                                        );
                                    }
                                    let call = classify_sequence(
                                        &s1,
                                        Some(&s2),
                                        hash,
                                        taxonomy,
                                        idx_opts,
                                        classify_opts,
                                        &mut scanner,
                                        &mut taxa,
                                        &mut hit_counts,
                                        &mut tx_frames,
                                        &mut local_counters,
                                        &mut local_buf,
                                    );
                                    local_stats.total_sequences += 1;
                                    local_stats.total_bases +=
                                        (s1.seq.len() + s2.seq.len()) as u64;
                                    if call != 0 {
                                        local_stats.total_classified += 1;
                                    } else {
                                        local_stats.total_unclassified += 1;
                                    }
                                    results.push((call, local_buf.clone(), s1));
                                }
                            }

                            let block = OutputBlock {
                                block_id,
                                results,
                                local_counters,
                                local_stats,
                            };
                            {
                                let mut q = output_queue.lock();
                                q.push(block);
                            }

                            loop {
                                let mut q = output_queue.lock();
                                let mut next_out = next_output_id.lock();
                                if q.peek().map(|b| b.block_id) == Some(*next_out) {
                                    let blk = q.pop().unwrap();
                                    *next_out += 1;
                                    drop(q);
                                    drop(next_out);

                                    let mut state = output_state.lock();
                                    let (ko, co, uo, st, tc) = &mut *state;
                                    st.total_sequences += blk.local_stats.total_sequences;
                                    st.total_bases += blk.local_stats.total_bases;
                                    st.total_classified += blk.local_stats.total_classified;
                                    st.total_unclassified += blk.local_stats.total_unclassified;
                                    for (taxid, counter) in blk.local_counters {
                                        tc.entry(taxid).or_default().merge(&counter);
                                    }
                                    for (call, kraken_str, seq) in &blk.results {
                                        let _ = process_output(
                                            *call,
                                            kraken_str,
                                            seq,
                                            taxonomy,
                                            ko,
                                            co,
                                            uo,
                                            &mut ClassificationStats::default(),
                                        );
                                    }
                                } else {
                                    break;
                                }
                            }
                        }
                    });
                }
            })
            .expect("Thread panicked");
        } else {
            let reader = parking_lot::Mutex::new(BatchSequenceReader::new(filename1)?);

            crossbeam::scope(|s| {
                for _ in 0..classify_opts.num_threads {
                    s.spawn(|_| {
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
                        let mut local_buf = String::with_capacity(512);

                        loop {
                            let (batch, block_id) = {
                                let mut rdr = reader.lock();
                                if !rdr.load_block(3 * 1024 * 1024) {
                                    break;
                                }
                                let bid = next_block_id.fetch_add(1, AtomOrd::SeqCst);
                                let mut batch = Vec::new();
                                while let Some(seq) = rdr.next_sequence() {
                                    batch.push(ClassificationInput::Single(seq.clone()));
                                }
                                (batch, bid)
                            };

                            let mut results = Vec::with_capacity(batch.len());
                            let mut local_counters = TaxonCounters::new();
                            let mut local_stats = ClassificationStats::default();

                            for input in batch {
                                if let ClassificationInput::Single(mut seq) = input {
                                    if classify_opts.minimum_quality_score > 0 {
                                        mask_low_quality_bases(
                                            &mut seq,
                                            classify_opts.minimum_quality_score,
                                        );
                                    }
                                    let call = classify_sequence(
                                        &seq,
                                        None,
                                        hash,
                                        taxonomy,
                                        idx_opts,
                                        classify_opts,
                                        &mut scanner,
                                        &mut taxa,
                                        &mut hit_counts,
                                        &mut tx_frames,
                                        &mut local_counters,
                                        &mut local_buf,
                                    );
                                    local_stats.total_sequences += 1;
                                    local_stats.total_bases += seq.seq.len() as u64;
                                    if call != 0 {
                                        local_stats.total_classified += 1;
                                    } else {
                                        local_stats.total_unclassified += 1;
                                    }
                                    results.push((call, local_buf.clone(), seq));
                                }
                            }

                            let block = OutputBlock {
                                block_id,
                                results,
                                local_counters,
                                local_stats,
                            };
                            {
                                let mut q = output_queue.lock();
                                q.push(block);
                            }

                            loop {
                                let mut q = output_queue.lock();
                                let mut next_out = next_output_id.lock();
                                if q.peek().map(|b| b.block_id) == Some(*next_out) {
                                    let blk = q.pop().unwrap();
                                    *next_out += 1;
                                    drop(q);
                                    drop(next_out);

                                    let mut state = output_state.lock();
                                    let (ko, co, uo, st, tc) = &mut *state;
                                    st.total_sequences += blk.local_stats.total_sequences;
                                    st.total_bases += blk.local_stats.total_bases;
                                    st.total_classified += blk.local_stats.total_classified;
                                    st.total_unclassified += blk.local_stats.total_unclassified;
                                    for (taxid, counter) in blk.local_counters {
                                        tc.entry(taxid).or_default().merge(&counter);
                                    }
                                    for (call, kraken_str, seq) in &blk.results {
                                        let _ = process_output(
                                            *call,
                                            kraken_str,
                                            seq,
                                            taxonomy,
                                            ko,
                                            co,
                                            uo,
                                            &mut ClassificationStats::default(),
                                        );
                                    }
                                } else {
                                    break;
                                }
                            }
                        }
                    });
                }
            })
            .expect("Thread panicked");
        }
    } else if classify_opts.paired_end_processing && !classify_opts.single_file_pairs {
        let mut reader1 = BatchSequenceReader::new(filename1)?;
        let mut reader2 = BatchSequenceReader::new(filename2)?;

        while reader1.load_batch(10000) {
            reader2.load_batch(10000);
            loop {
                let seq1 = reader1.next_sequence().cloned();
                let seq2 = reader2.next_sequence().cloned();
                match (seq1, seq2) {
                    (Some(mut s1), Some(mut s2)) => {
                        stats.total_sequences += 1;
                        stats.total_bases += (s1.seq.len() + s2.seq.len()) as u64;
                        if classify_opts.minimum_quality_score > 0 {
                            mask_low_quality_bases(&mut s1, classify_opts.minimum_quality_score);
                            mask_low_quality_bases(&mut s2, classify_opts.minimum_quality_score);
                        }
                        let call = classify_sequence(
                            &s1,
                            Some(&s2),
                            hash,
                            taxonomy,
                            idx_opts,
                            classify_opts,
                            &mut scanner,
                            &mut taxa,
                            &mut hit_counts,
                            &mut tx_frames,
                            taxon_counters,
                            &mut output_buf,
                        );
                        process_output(
                            call,
                            &output_buf,
                            &s1,
                            taxonomy,
                            outputs.kraken_output,
                            outputs.classified_output,
                            outputs.unclassified_output,
                            stats,
                        )?;
                    }
                    _ => break,
                }
            }
        }
    } else if classify_opts.paired_end_processing && classify_opts.single_file_pairs {
        let mut reader = BatchSequenceReader::new(filename1)?;

        while reader.load_batch(20000) {
            loop {
                let seq1 = reader.next_sequence().cloned();
                let seq2 = reader.next_sequence().cloned();
                match (seq1, seq2) {
                    (Some(mut s1), Some(mut s2)) => {
                        stats.total_sequences += 1;
                        stats.total_bases += (s1.seq.len() + s2.seq.len()) as u64;
                        if classify_opts.minimum_quality_score > 0 {
                            mask_low_quality_bases(&mut s1, classify_opts.minimum_quality_score);
                            mask_low_quality_bases(&mut s2, classify_opts.minimum_quality_score);
                        }
                        let call = classify_sequence(
                            &s1,
                            Some(&s2),
                            hash,
                            taxonomy,
                            idx_opts,
                            classify_opts,
                            &mut scanner,
                            &mut taxa,
                            &mut hit_counts,
                            &mut tx_frames,
                            taxon_counters,
                            &mut output_buf,
                        );
                        process_output(
                            call,
                            &output_buf,
                            &s1,
                            taxonomy,
                            outputs.kraken_output,
                            outputs.classified_output,
                            outputs.unclassified_output,
                            stats,
                        )?;
                    }
                    _ => break,
                }
            }
        }
    } else {
        let mut reader = BatchSequenceReader::new(filename1)?;

        while reader.load_block(3 * 1024 * 1024) {
            while let Some(seq) = reader.next_sequence() {
                let mut seq = seq.clone();
                stats.total_sequences += 1;
                stats.total_bases += seq.seq.len() as u64;

                if classify_opts.minimum_quality_score > 0 {
                    mask_low_quality_bases(&mut seq, classify_opts.minimum_quality_score);
                }

                let call = classify_sequence(
                    &seq,
                    None,
                    hash,
                    taxonomy,
                    idx_opts,
                    classify_opts,
                    &mut scanner,
                    &mut taxa,
                    &mut hit_counts,
                    &mut tx_frames,
                    taxon_counters,
                    &mut output_buf,
                );
                process_output(
                    call,
                    &output_buf,
                    &seq,
                    taxonomy,
                    outputs.kraken_output,
                    outputs.classified_output,
                    outputs.unclassified_output,
                    stats,
                )?;
            }
        }
    }

    Ok(())
}

/// Run the full classification pipeline on input files.
fn classify_with_index_data(
    input_files: &[String],
    hash: &CompactHashTable,
    taxonomy: &Taxonomy,
    idx_opts: &IndexOptions,
    opts: &ClassifyOptions,
    kraken_output: Option<&str>,
    classified_output: Option<&str>,
    unclassified_output: Option<&str>,
    report_filename: Option<&str>,
) -> io::Result<ClassificationStats> {
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
    let mut outputs = OutputHandles {
        kraken_output: &mut kraken_out,
        classified_output: &mut classified_out,
        unclassified_output: &mut unclassified_out,
    };

    if input_files.is_empty() {
        if classify_opts.paired_end_processing && !classify_opts.single_file_pairs {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "paired end processing used with no files specified",
            ));
        }
        process_files(
            None,
            None,
            hash,
            taxonomy,
            idx_opts,
            &classify_opts,
            &mut stats,
            &mut outputs,
            &mut taxon_counters,
        )?;
    } else {
        let mut i = 0usize;
        while i < input_files.len() {
            if classify_opts.paired_end_processing && !classify_opts.single_file_pairs {
                if i + 1 == input_files.len() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "paired end processing used with unpaired file",
                    ));
                }
                process_files(
                    Some(input_files[i].as_str()),
                    Some(input_files[i + 1].as_str()),
                    hash,
                    taxonomy,
                    idx_opts,
                    &classify_opts,
                    &mut stats,
                    &mut outputs,
                    &mut taxon_counters,
                )?;
                i += 1;
            } else {
                let filename1 = if input_files[i] == "-" {
                    None
                } else {
                    Some(input_files[i].as_str())
                };
                process_files(
                    filename1,
                    None,
                    hash,
                    taxonomy,
                    idx_opts,
                    &classify_opts,
                    &mut stats,
                    &mut outputs,
                    &mut taxon_counters,
                )?;
            }
            i += 1;
        }
    }

    if let Some(rf) = report_filename {
        if classify_opts.mpa_style_report {
            reports::report_mpa_style(
                rf,
                classify_opts.report_zero_counts,
                taxonomy,
                &taxon_counters,
            )?;
        } else {
            reports::report_kraken_style(
                rf,
                classify_opts.report_zero_counts,
                classify_opts.report_kmer_data,
                taxonomy,
                &taxon_counters,
                stats.total_sequences,
                stats.total_unclassified,
            )?;
        }
    }

    Ok(stats)
}

fn classify(opts: &Options, index_data: &mut IndexData) -> io::Result<()> {
    let started = Instant::now();
    let start = Duration::from_secs(0);
    let idx_opts = &index_data.0;
    let taxonomy = &index_data.1;
    let hash = &index_data.2;

    let classify_opts = ClassifyOptions {
        confidence_threshold: opts.confidence_threshold,
        minimum_quality_score: opts.minimum_quality_score,
        minimum_hit_groups: opts.minimum_hit_groups,
        paired_end_processing: opts.paired_end_processing,
        single_file_pairs: opts.single_file_pairs,
        use_translated_search: opts.use_translated_search,
        print_scientific_name: opts.print_scientific_name,
        quick_mode: opts.quick_mode,
        report_filename: opts.report_filename.clone(),
        report_kmer_data: opts.report_kmer_data,
        report_zero_counts: opts.report_zero_counts,
        mpa_style_report: opts.mpa_style_report,
        use_memory_mapping: opts.use_memory_mapping,
        num_threads: opts.num_threads,
    };

    let stats = classify_with_index_data(
        &opts.filenames,
        hash,
        taxonomy,
        idx_opts,
        &classify_opts,
        if opts.kraken_output_filename.is_empty() || opts.kraken_output_filename == "-" {
            None
        } else {
            Some(opts.kraken_output_filename.as_str())
        },
        if opts.classified_output_filename.is_empty() {
            None
        } else {
            Some(opts.classified_output_filename.as_str())
        },
        if opts.unclassified_output_filename.is_empty() {
            None
        } else {
            Some(opts.unclassified_output_filename.as_str())
        },
        if opts.report_filename.is_empty() {
            None
        } else {
            Some(opts.report_filename.as_str())
        },
    )?;
    let end = started.elapsed();
    eprintln!("{}", report_stats(start, end, &stats));
    Ok(())
}

pub fn classify_main(args: &[String]) -> io::Result<()> {
    let mut opts = Options::default();
    parse_command_line(args, &mut opts)?;

    if opts.daemon_mode {
        #[cfg(unix)]
        {
            let classify_opts = ClassifyOptions {
                confidence_threshold: opts.confidence_threshold,
                minimum_quality_score: opts.minimum_quality_score,
                minimum_hit_groups: opts.minimum_hit_groups,
                paired_end_processing: opts.paired_end_processing,
                single_file_pairs: opts.single_file_pairs,
                use_translated_search: false,
                print_scientific_name: opts.print_scientific_name,
                quick_mode: opts.quick_mode,
                report_filename: opts.report_filename.clone(),
                report_kmer_data: opts.report_kmer_data,
                report_zero_counts: opts.report_zero_counts,
                mpa_style_report: opts.mpa_style_report,
                use_memory_mapping: opts.use_memory_mapping,
                num_threads: opts.num_threads,
            };
            return run_daemon(
                &opts.index_filename,
                &opts.taxonomy_filename,
                &opts.options_filename,
                &classify_opts,
            );
        }
        #[cfg(not(unix))]
        {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "daemon mode is only supported on Unix",
            ));
        }
    }

    let mut index_data = load_index(&mut opts)?;
    classify(&opts, &mut index_data)
}

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
    let idx_opts =
        IndexOptions::read_from_file(opts_filename).expect("Failed to read index options");

    // Load taxonomy
    let mut taxonomy = Taxonomy::from_file(taxonomy_filename, opts.use_memory_mapping)?;
    taxonomy.generate_external_to_internal_id_map();

    // Load hash table
    let hash = CompactHashTable::from_file(hash_filename, opts.use_memory_mapping)?;

    classify_with_index_data(
        input_files,
        &hash,
        &taxonomy,
        &idx_opts,
        opts,
        kraken_output,
        classified_output,
        unclassified_output,
        report_filename,
    )
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
/// use kraken2_rs::classify::{ClassifyOptions, classify_sequences_in_memory};
/// use kraken2_rs::compact_hash::CompactHashTable;
/// use kraken2_rs::taxonomy::Taxonomy;
/// use kraken2_rs::types::{IndexOptions, Sequence};
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
            seq,
            None,
            hash,
            taxonomy,
            idx_opts,
            opts,
            &mut scanner,
            &mut taxa,
            &mut hit_counts,
            &mut tx_frames,
            &mut taxon_counters,
            &mut output_buf,
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
/// use kraken2_rs::classify::{ClassifyDb, ClassifyOptions};
/// use kraken2_rs::types::Sequence;
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
        Ok(ClassifyDb {
            idx_opts,
            taxonomy,
            hash,
        })
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
    pub fn classify_one(&self, seq: &Sequence, opts: &ClassifyOptions) -> ClassificationResult {
        self.classify(std::slice::from_ref(seq), opts)
            .into_iter()
            .next()
            .unwrap()
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
    let hash = std::sync::Arc::new(CompactHashTable::from_file(
        hash_filename,
        opts.use_memory_mapping,
    )?);

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
        let files: Vec<String> = trimmed.split_whitespace().map(|s| s.to_string()).collect();

        if files.is_empty() {
            continue;
        }

        let mut classify_opts = opts.clone();
        classify_opts.use_translated_search = !idx_opts.dna_db;

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

        for input_file in &files {
            let filename = if input_file == "-" {
                None
            } else {
                Some(input_file.as_str())
            };
            let mut reader = crate::seq::BatchSequenceReader::new(filename)?;

            while reader.load_block(3 * 1024 * 1024) {
                while let Some(seq) = reader.next_sequence() {
                    let seq = seq.clone();
                    stats.total_sequences += 1;

                    let call = classify_sequence(
                        &seq,
                        None,
                        &hash,
                        &taxonomy,
                        &idx_opts,
                        &classify_opts,
                        &mut scanner,
                        &mut taxa,
                        &mut hit_counts,
                        &mut tx_frames,
                        &mut taxon_counters,
                        &mut output_buf,
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    fn test_data_dir() -> String {
        std::env::var("KRAKEN2_TEST_DATA_DIR")
            .unwrap_or_else(|_| "/tmp/kraken2_test_data".to_string())
    }

    #[test]
    fn test_tokenize_string() {
        let mut argv = Vec::new();
        tokenize_string(&mut argv, "classify -H db/hash.k2d -q reads.fq");
        assert_eq!(
            argv,
            vec!["classify", "-H", "db/hash.k2d", "-q", "reads.fq"]
        );
    }

    #[test]
    fn test_parse_command_line() {
        let args = vec![
            "classify".to_string(),
            "-H".to_string(),
            "hash".to_string(),
            "-t".to_string(),
            "taxo".to_string(),
            "-o".to_string(),
            "opts".to_string(),
            "-P".to_string(),
            "-C".to_string(),
            "classified#.fq".to_string(),
            "-R".to_string(),
            "report.txt".to_string(),
            "reads_1.fq".to_string(),
            "reads_2.fq".to_string(),
        ];
        let mut opts = Options::default();
        parse_command_line(&args, &mut opts).unwrap();
        assert_eq!(opts.index_filename, "hash");
        assert!(opts.paired_end_processing);
        assert_eq!(opts.classified_output_filename, "classified#.fq");
        assert_eq!(opts.report_filename, "report.txt");
        assert_eq!(opts.filenames, vec!["reads_1.fq", "reads_2.fq"]);
    }

    #[test]
    fn test_initialize_outputs() {
        let dir =
            std::env::temp_dir().join(format!("kraken2_classify_outputs_{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();

        let mut opts = Options::default();
        opts.paired_end_processing = true;
        opts.classified_output_filename = dir.join("classified#.fq").display().to_string();
        opts.unclassified_output_filename = dir.join("unclassified#.fq").display().to_string();
        opts.kraken_output_filename = dir.join("kraken.out").display().to_string();

        let mut outputs = OutputStreamData::default();
        initialize_outputs(&opts, &mut outputs, SequenceFormat::Fastq).unwrap();
        assert!(outputs.initialized);
        assert!(outputs.printing_sequences);
        assert!(outputs.classified_output1.is_some());
        assert!(outputs.classified_output2.is_some());
        assert!(outputs.unclassified_output1.is_some());
        assert!(outputs.unclassified_output2.is_some());
        assert!(outputs.kraken_output.is_some());

        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_report_stats() {
        let stats = ClassificationStats {
            total_sequences: 10,
            total_bases: 2_000_000,
            total_classified: 7,
            total_unclassified: 3,
        };
        let text = report_stats(Duration::from_secs(1), Duration::from_secs(3), &stats);
        assert!(text.contains("10 sequences"));
        assert!(text.contains("2.00 Mbp"));
        assert!(text.contains("7 sequences classified"));
        assert!(text.contains("3 sequences unclassified"));
    }

    #[test]
    fn test_load_index() {
        let data_dir = std::path::PathBuf::from(test_data_dir());
        if !data_dir.exists() {
            return;
        }
        let mut opts = Options {
            index_filename: data_dir.join("hash.k2d").display().to_string(),
            taxonomy_filename: data_dir.join("taxo.k2d").display().to_string(),
            options_filename: data_dir.join("opts.k2d").display().to_string(),
            ..Default::default()
        };
        let (idx_opts, taxonomy, hash) = load_index(&mut opts).unwrap();
        assert_eq!(idx_opts.k, 35);
        assert!(taxonomy.node_count() > 0);
        assert!(hash.capacity() > 0);
    }

    #[test]
    fn test_paired_multithread_matches_single_thread_stats_and_output() {
        let root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let ref_dir = root.join("tests/reference");
        if !ref_dir.exists() {
            return;
        }

        let mut taxonomy =
            Taxonomy::from_file(ref_dir.join("taxo.k2d").to_str().unwrap(), false).unwrap();
        taxonomy.generate_external_to_internal_id_map();
        let hash = CompactHashTable::from_file(ref_dir.join("hash.k2d").to_str().unwrap(), false)
            .unwrap();
        let idx_opts = IndexOptions::read_from_file(ref_dir.join("opts.k2d").to_str().unwrap())
            .unwrap();

        let dir = tempfile::tempdir().unwrap();
        let r1 = dir.path().join("reads_1.fq");
        let r2 = dir.path().join("reads_2.fq");
        fs::write(
            &r1,
            "@r1/1\nAAAA\n+\nIIII\n@r2/1\nTTTT\n+\nIIII\n@r3/1\nCCCC\n+\nIIII\n@r4/1\nGGGG\n+\nIIII\n",
        )
        .unwrap();
        fs::write(
            &r2,
            "@r1/2\nAAAA\n+\nIIII\n@r2/2\nTTTT\n+\nIIII\n@r3/2\nCCCC\n+\nIIII\n@r4/2\nGGGG\n+\nIIII\n",
        )
        .unwrap();

        let common = ClassifyOptions {
            paired_end_processing: true,
            ..Default::default()
        };

        let out1 = dir.path().join("single.out");
        let stats_single = classify_with_index_data(
            &[
                r1.to_str().unwrap().to_string(),
                r2.to_str().unwrap().to_string(),
            ],
            &hash,
            &taxonomy,
            &idx_opts,
            &ClassifyOptions {
                num_threads: 1,
                ..common.clone()
            },
            Some(out1.to_str().unwrap()),
            None,
            None,
            None,
        )
        .unwrap();

        let out4 = dir.path().join("multi.out");
        let stats_multi = classify_with_index_data(
            &[
                r1.to_str().unwrap().to_string(),
                r2.to_str().unwrap().to_string(),
            ],
            &hash,
            &taxonomy,
            &idx_opts,
            &ClassifyOptions {
                num_threads: 4,
                ..common
            },
            Some(out4.to_str().unwrap()),
            None,
            None,
            None,
        )
        .unwrap();

        assert_eq!(stats_single.total_sequences, 4);
        assert_eq!(stats_single.total_sequences, stats_multi.total_sequences);
        assert_eq!(stats_single.total_bases, stats_multi.total_bases);
        assert_eq!(stats_single.total_classified, stats_multi.total_classified);
        assert_eq!(stats_single.total_unclassified, stats_multi.total_unclassified);
        assert_eq!(fs::read_to_string(out1).unwrap(), fs::read_to_string(out4).unwrap());
    }
}

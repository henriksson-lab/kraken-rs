use clap::{Parser, Subcommand};
use std::process;

use kraken2::classify::{self, ClassifyOptions};
use kraken2::build_db::{self, BuildOptions};
use kraken2::estimate::{self, EstimateOptions};
use kraken2::types::*;
use kraken2::utilities::expand_spaced_seed_mask;

#[derive(Parser)]
#[command(name = "kraken2", about = "Taxonomic sequence classification system")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Classify sequences against a Kraken 2 database
    Classify {
        /// Database directory
        #[arg(long, short = 'd')]
        db: String,

        /// Number of threads
        #[arg(long, short = 'p', default_value = "1")]
        threads: usize,

        /// Quick mode (use first hit group)
        #[arg(long)]
        quick: bool,

        /// Confidence threshold (0.0-1.0)
        #[arg(long, short = 'T', default_value = "0.0")]
        confidence: f64,

        /// Minimum hit groups for classification
        #[arg(long, short = 'g', default_value = "0")]
        minimum_hit_groups: i64,

        /// Paired-end processing (two input files)
        #[arg(long, short = 'P')]
        paired: bool,

        /// Single-file paired-end (interleaved mates in one file)
        #[arg(long, short = 'S')]
        single_file_pairs: bool,

        /// Output classified sequences to file
        #[arg(long, short = 'C')]
        classified_out: Option<String>,

        /// Output unclassified sequences to file
        #[arg(long, short = 'U')]
        unclassified_out: Option<String>,

        /// Kraken output file (- to suppress)
        #[arg(long, short = 'O')]
        output: Option<String>,

        /// Report file
        #[arg(long, short = 'R')]
        report: Option<String>,

        /// MPA-style report
        #[arg(long)]
        use_mpa_style: bool,

        /// Report zero-count taxa
        #[arg(long)]
        report_zero_counts: bool,

        /// Include minimizer data in report
        #[arg(long)]
        report_minimizer_data: bool,

        /// Use memory mapping for database
        #[arg(long)]
        memory_mapping: bool,

        /// Print scientific names instead of taxids
        #[arg(long)]
        use_names: bool,

        /// Minimum base quality score (FASTQ)
        #[arg(long, short = 'Q', default_value = "0")]
        minimum_base_quality: u8,

        /// Run as daemon (persistent database, reads commands from stdin)
        #[arg(long, short = 'D')]
        daemon: bool,

        /// Input files
        input: Vec<String>,
    },

    /// Build a Kraken 2 database
    Build {
        /// Hash table output file
        #[arg(long, short = 'H')]
        hashtable: String,

        /// Taxonomy output file
        #[arg(long, short = 't')]
        taxonomy: String,

        /// Options output file
        #[arg(long, short = 'o')]
        options: String,

        /// Sequence ID to taxon map file
        #[arg(long, short = 'm')]
        id_map: String,

        /// NCBI taxonomy directory
        #[arg(long, short = 'n')]
        ncbi_taxonomy_dir: String,

        /// K-mer length
        #[arg(long, short = 'k')]
        kmer_len: usize,

        /// Minimizer length
        #[arg(long, short = 'l')]
        minimizer_len: usize,

        /// Hash table capacity
        #[arg(long, short = 'c')]
        capacity: usize,

        /// Maximum capacity (for subsampled DBs)
        #[arg(long, short = 'M', default_value = "0")]
        max_capacity: usize,

        /// Number of threads
        #[arg(long, short = 'p', default_value = "1")]
        threads: usize,

        /// Input is protein sequences
        #[arg(long, short = 'X')]
        protein: bool,

        /// Fast (nondeterministic) building
        #[arg(long, short = 'F')]
        fast_build: bool,

        /// Block size
        #[arg(long, short = 'B', default_value = "10485760")]
        block_size: usize,

        /// Subblock size
        #[arg(long, short = 'b', default_value = "1024")]
        subblock_size: usize,

        /// Bits for taxid storage
        #[arg(long, short = 'r', default_value = "0")]
        taxid_bits: usize,

        /// Spaced seed mask (binary string)
        #[arg(long, short = 'S')]
        spaced_seed: Option<String>,

        /// Toggle mask (binary string)
        #[arg(long)]
        toggle_mask: Option<String>,
    },

    /// Estimate hash table capacity needed
    Estimate {
        /// K-mer length
        #[arg(long, short = 'k')]
        kmer_len: usize,

        /// Minimizer length
        #[arg(long, short = 'l')]
        minimizer_len: usize,

        /// Max qualifying hash sections (1-1024)
        #[arg(long, short = 'n', default_value = "4")]
        n: usize,

        /// Input is protein
        #[arg(long, short = 'X')]
        protein: bool,

        /// Number of threads
        #[arg(long, short = 'p', default_value = "1")]
        threads: usize,

        /// Read block size
        #[arg(long, short = 'B', default_value = "31457280")]
        block_size: usize,

        /// Spaced seed mask (binary string)
        #[arg(long, short = 'S')]
        spaced_seed: Option<String>,

        /// Toggle mask (binary string)
        #[arg(long)]
        toggle_mask: Option<String>,
    },

    /// Download taxonomy or library data from NCBI
    Download {
        /// Database directory
        #[arg(long, short = 'd')]
        db: String,

        /// Download taxonomy (nodes.dmp, names.dmp)
        #[arg(long)]
        taxonomy: bool,

        /// Download a library (archaea, bacteria, viral, fungi, plant, protozoa, human, plasmid, UniVec, UniVec_Core)
        #[arg(long)]
        library: Option<String>,

        /// Protein database
        #[arg(long)]
        protein: bool,

        /// Skip accession-to-taxid map downloads
        #[arg(long)]
        skip_maps: bool,
    },

    /// Clean up intermediate files after database build
    Clean {
        /// Database directory
        #[arg(long, short = 'd')]
        db: String,
    },

    /// Inspect a Kraken 2 database
    Inspect {
        /// Database directory
        #[arg(long, short = 'd')]
        db: String,

        /// Use memory mapping
        #[arg(long)]
        memory_mapping: bool,

        /// Skip counts (just show DB stats)
        #[arg(long)]
        skip_counts: bool,

        /// MPA-style output
        #[arg(long)]
        use_mpa_style: bool,

        /// Report zero-count taxa
        #[arg(long)]
        report_zero_counts: bool,

        /// Number of threads
        #[arg(long, short = 'p', default_value = "1")]
        threads: usize,
    },

    /// Convert a BLAST database to FASTA format
    Blast2fasta {
        /// BLAST database prefix (without .nin/.pin extension)
        #[arg(long, short = 'i')]
        input: String,

        /// Output FASTA file
        #[arg(long, short = 'o')]
        output: String,

        /// Include taxid in headers
        #[arg(long)]
        include_taxid: bool,
    },
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Classify {
            db, threads, quick, confidence, minimum_hit_groups,
            paired, single_file_pairs, classified_out, unclassified_out, output, report,
            use_mpa_style, report_zero_counts, report_minimizer_data,
            memory_mapping, use_names, minimum_base_quality, daemon, input,
        } => {
            let hash_file = format!("{}/hash.k2d", db);
            let taxo_file = format!("{}/taxo.k2d", db);
            let opts_file = format!("{}/opts.k2d", db);

            let opts = ClassifyOptions {
                confidence_threshold: confidence,
                minimum_quality_score: minimum_base_quality,
                minimum_hit_groups,
                paired_end_processing: paired || single_file_pairs,
                single_file_pairs,
                print_scientific_name: use_names,
                quick_mode: quick,
                report_kmer_data: report_minimizer_data,
                report_zero_counts,
                mpa_style_report: use_mpa_style,
                use_memory_mapping: memory_mapping,
                num_threads: threads,
                ..Default::default()
            };

            if daemon {
                #[cfg(unix)]
                {
                    if let Err(e) = classify::run_daemon(&hash_file, &taxo_file, &opts_file, &opts) {
                        eprintln!("Daemon error: {e}");
                        process::exit(1);
                    }
                    return;
                }
                #[cfg(not(unix))]
                {
                    eprintln!("Daemon mode is only supported on Unix");
                    process::exit(1);
                }
            }

            match classify::run_classify(
                &input, &hash_file, &taxo_file, &opts_file, &opts,
                output.as_deref(), classified_out.as_deref(),
                unclassified_out.as_deref(), report.as_deref(),
            ) {
                Ok(stats) => {
                    eprintln!(
                        "{} sequences ({:.2}% classified)",
                        stats.total_sequences,
                        if stats.total_sequences > 0 {
                            100.0 * stats.total_classified as f64 / stats.total_sequences as f64
                        } else {
                            0.0
                        }
                    );
                }
                Err(e) => {
                    eprintln!("Error: {}", e);
                    process::exit(1);
                }
            }
        }

        Commands::Build {
            hashtable, taxonomy, options, id_map, ncbi_taxonomy_dir,
            kmer_len, minimizer_len, capacity, max_capacity, threads,
            protein, fast_build, block_size, subblock_size, taxid_bits,
            spaced_seed, toggle_mask,
        } => {
            let mut spaced_seed_mask = DEFAULT_SPACED_SEED_MASK;
            if let Some(ref s) = spaced_seed {
                spaced_seed_mask = u64::from_str_radix(s, 2).expect("Invalid spaced seed mask");
                let bits = if protein { BITS_PER_CHAR_PRO } else { BITS_PER_CHAR_DNA };
                expand_spaced_seed_mask(&mut spaced_seed_mask, bits as i32);
            }
            let mut toggle = DEFAULT_TOGGLE_MASK;
            if let Some(ref t) = toggle_mask {
                toggle = u64::from_str_radix(t, 2).expect("Invalid toggle mask");
            }

            let mut opts = BuildOptions {
                id_to_taxon_map_filename: id_map,
                ncbi_taxonomy_directory: ncbi_taxonomy_dir,
                hashtable_filename: hashtable,
                options_filename: options,
                taxonomy_filename: taxonomy,
                block_size,
                subblock_size,
                requested_bits_for_taxid: taxid_bits,
                num_threads: threads,
                input_is_protein: protein,
                k: kmer_len,
                l: minimizer_len,
                capacity,
                maximum_capacity: max_capacity,
                spaced_seed_mask,
                toggle_mask: toggle,
                deterministic_build: !fast_build,
                ..Default::default()
            };

            if let Err(e) = build_db::build_database(&mut opts) {
                eprintln!("Error building database: {}", e);
                process::exit(1);
            }
        }

        Commands::Estimate {
            kmer_len, minimizer_len, n, protein, threads, block_size,
            spaced_seed, toggle_mask,
        } => {
            let mut spaced_seed_mask = DEFAULT_SPACED_SEED_MASK;
            if let Some(ref s) = spaced_seed {
                spaced_seed_mask = u64::from_str_radix(s, 2).expect("Invalid spaced seed mask");
                let bits = if protein { BITS_PER_CHAR_PRO } else { BITS_PER_CHAR_DNA };
                expand_spaced_seed_mask(&mut spaced_seed_mask, bits as i32);
            }
            let mut toggle = DEFAULT_TOGGLE_MASK;
            if let Some(ref t) = toggle_mask {
                toggle = u64::from_str_radix(t, 2).expect("Invalid toggle mask");
            }

            let opts = EstimateOptions {
                k: kmer_len,
                l: minimizer_len,
                n,
                input_is_protein: protein,
                threads,
                block_size,
                spaced_seed_mask,
                toggle_mask: toggle,
            };

            match estimate::estimate_capacity(&opts) {
                Ok(capacity) => println!("{}", capacity),
                Err(e) => {
                    eprintln!("Error estimating capacity: {}", e);
                    process::exit(1);
                }
            }
        }

        Commands::Inspect {
            db, memory_mapping, skip_counts, use_mpa_style, report_zero_counts, threads: _,
        } => {
            let hash_file = format!("{}/hash.k2d", db);
            let taxo_file = format!("{}/taxo.k2d", db);
            let opts_file = format!("{}/opts.k2d", db);

            // Load taxonomy
            let mut taxonomy = match kraken2::taxonomy::Taxonomy::from_file(&taxo_file, memory_mapping) {
                Ok(t) => t,
                Err(e) => {
                    eprintln!("Error loading taxonomy: {}", e);
                    process::exit(1);
                }
            };
            taxonomy.generate_external_to_internal_id_map();

            // Load hash table
            let hash = match kraken2::compact_hash::CompactHashTable::from_file(&hash_file, memory_mapping) {
                Ok(h) => h,
                Err(e) => {
                    eprintln!("Error loading hash table: {}", e);
                    process::exit(1);
                }
            };

            if skip_counts {
                println!("Database options:");
                if let Ok(idx_opts) = kraken2::types::IndexOptions::read_from_file(&opts_file) {
                    println!("  k = {}", idx_opts.k);
                    println!("  l = {}", idx_opts.l);
                    println!("  DNA database: {}", idx_opts.dna_db);
                }
                println!("  Capacity: {}", hash.capacity());
                println!("  Size: {}", hash.size());
                println!("  Occupancy: {:.4}", hash.occupancy());
            } else {
                // Get value counts and produce report
                let value_counts = hash.get_value_counts();
                let mut taxon_counters = kraken2::readcounts::TaxonCounters::new();
                for (&taxid, &count) in &value_counts {
                    let counter = taxon_counters.entry(taxid).or_default();
                    counter.n_reads = count;
                }

                let total = value_counts.values().sum::<u64>();
                if use_mpa_style {
                    let _ = kraken2::reports::report_mpa_style(
                        "/dev/stdout", report_zero_counts, &taxonomy, &taxon_counters);
                } else {
                    let _ = kraken2::reports::report_kraken_style(
                        "/dev/stdout", report_zero_counts, false,
                        &taxonomy, &taxon_counters, total, 0);
                }
            }
        }

        Commands::Download { db, taxonomy, library, protein, skip_maps } => {
            if taxonomy {
                if let Err(e) = kraken2::download::download_taxonomy(&db, skip_maps, protein) {
                    eprintln!("Error downloading taxonomy: {e}");
                    process::exit(1);
                }
            }
            if let Some(ref lib_type) = library {
                if let Err(e) = kraken2::download::download_library(&db, lib_type, protein) {
                    eprintln!("Error downloading library: {e}");
                    process::exit(1);
                }
            }
            if !taxonomy && library.is_none() {
                eprintln!("Specify --taxonomy and/or --library TYPE");
                process::exit(1);
            }
        }

        Commands::Clean { db } => {
            if let Err(e) = kraken2::download::clean_db(&db) {
                eprintln!("Error cleaning database: {e}");
                process::exit(1);
            }
        }

        Commands::Blast2fasta { input, output, include_taxid } => {
            if let Err(e) = kraken2::blast::blast_to_fasta(&input, &output, include_taxid) {
                eprintln!("Error converting BLAST database: {e}");
                process::exit(1);
            }
        }
    }
}

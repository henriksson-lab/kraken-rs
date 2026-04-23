use clap::{Parser, Subcommand};
use std::process;

use kraken2::build_db;
use kraken2::blast;
use kraken2::classify;
use kraken2::dump_table;
use kraken2::dust;
use kraken2::estimate;
use kraken2::lookup;
use kraken2::mmtest;

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

    /// Dump a Kraken 2 table using the translated dump_table entrypoint
    DumpTable {
        #[arg(trailing_var_arg = true, allow_hyphen_values = true)]
        args: Vec<String>,
    },

    /// Look up accession numbers using the translated standalone tool
    LookupAccessionNumbers {
        #[arg(trailing_var_arg = true, allow_hyphen_values = true)]
        args: Vec<String>,
    },

    /// Mask low-complexity regions using the translated k2mask entrypoint
    K2mask {
        #[arg(trailing_var_arg = true, allow_hyphen_values = true)]
        args: Vec<String>,
    },

    /// Run the translated minimizer scanner smoke-test tool
    Mmtest,
}

fn push_flag(args: &mut Vec<String>, flag: &str) {
    args.push(flag.to_string());
}

fn push_option<T: ToString>(args: &mut Vec<String>, flag: &str, value: T) {
    args.push(flag.to_string());
    args.push(value.to_string());
}

fn classify_args(
    db: String,
    threads: usize,
    quick: bool,
    confidence: f64,
    minimum_hit_groups: i64,
    paired: bool,
    single_file_pairs: bool,
    classified_out: Option<String>,
    unclassified_out: Option<String>,
    output: Option<String>,
    report: Option<String>,
    use_mpa_style: bool,
    report_zero_counts: bool,
    report_minimizer_data: bool,
    memory_mapping: bool,
    use_names: bool,
    minimum_base_quality: u8,
    daemon: bool,
    input: Vec<String>,
) -> Vec<String> {
    let mut args = vec!["classify".to_string()];
    push_option(&mut args, "-H", format!("{db}/hash.k2d"));
    push_option(&mut args, "-t", format!("{db}/taxo.k2d"));
    push_option(&mut args, "-o", format!("{db}/opts.k2d"));
    push_option(&mut args, "-p", threads);
    push_option(&mut args, "-T", confidence);
    push_option(&mut args, "-g", minimum_hit_groups);
    push_option(&mut args, "-Q", minimum_base_quality);
    if quick {
        push_flag(&mut args, "-q");
    }
    if paired {
        push_flag(&mut args, "-P");
    }
    if single_file_pairs {
        push_flag(&mut args, "-S");
    }
    if use_mpa_style {
        push_flag(&mut args, "-m");
    }
    if report_zero_counts {
        push_flag(&mut args, "-z");
    }
    if report_minimizer_data {
        push_flag(&mut args, "-K");
    }
    if memory_mapping {
        push_flag(&mut args, "-M");
    }
    if use_names {
        push_flag(&mut args, "-n");
    }
    if daemon {
        push_flag(&mut args, "-D");
    }
    if let Some(path) = report {
        push_option(&mut args, "-R", path);
    }
    if let Some(path) = classified_out {
        push_option(&mut args, "-C", path);
    }
    if let Some(path) = unclassified_out {
        push_option(&mut args, "-U", path);
    }
    if let Some(path) = output {
        push_option(&mut args, "-O", path);
    }
    args.extend(input);
    args
}

#[allow(clippy::too_many_arguments)]
fn build_db_args(
    hashtable: String,
    taxonomy: String,
    options: String,
    id_map: String,
    ncbi_taxonomy_dir: String,
    kmer_len: usize,
    minimizer_len: usize,
    capacity: usize,
    max_capacity: usize,
    threads: usize,
    protein: bool,
    fast_build: bool,
    block_size: usize,
    subblock_size: usize,
    taxid_bits: usize,
    spaced_seed: Option<String>,
    toggle_mask: Option<String>,
) -> Vec<String> {
    let mut args = vec!["build_db".to_string()];
    push_option(&mut args, "-H", hashtable);
    push_option(&mut args, "-t", taxonomy);
    push_option(&mut args, "-o", options);
    push_option(&mut args, "-m", id_map);
    push_option(&mut args, "-n", ncbi_taxonomy_dir);
    push_option(&mut args, "-k", kmer_len);
    push_option(&mut args, "-l", minimizer_len);
    push_option(&mut args, "-c", capacity);
    push_option(&mut args, "-p", threads);
    push_option(&mut args, "-B", block_size);
    push_option(&mut args, "-b", subblock_size);
    push_option(&mut args, "-r", taxid_bits);
    if max_capacity > 0 {
        push_option(&mut args, "-M", max_capacity);
    }
    if protein {
        push_flag(&mut args, "-X");
    }
    if fast_build {
        push_flag(&mut args, "-F");
    }
    if let Some(mask) = spaced_seed {
        push_option(&mut args, "-S", mask);
    }
    if let Some(mask) = toggle_mask {
        push_option(&mut args, "-T", mask);
    }
    args
}

fn estimate_args(
    kmer_len: usize,
    minimizer_len: usize,
    n: usize,
    protein: bool,
    threads: usize,
    block_size: usize,
    spaced_seed: Option<String>,
    toggle_mask: Option<String>,
) -> Vec<String> {
    let mut args = vec!["estimate_capacity".to_string()];
    push_option(&mut args, "-k", kmer_len);
    push_option(&mut args, "-l", minimizer_len);
    push_option(&mut args, "-n", n);
    push_option(&mut args, "-p", threads);
    push_option(&mut args, "-B", block_size);
    if protein {
        push_flag(&mut args, "-X");
    }
    if let Some(mask) = spaced_seed {
        push_option(&mut args, "-S", mask);
    }
    if let Some(mask) = toggle_mask {
        push_option(&mut args, "-T", mask);
    }
    args
}

fn blast_to_fasta_args(input: String, output: String, include_taxid: bool) -> Vec<String> {
    let mut args = vec!["blast_to_fasta".to_string()];
    push_option(&mut args, "-o", output);
    if include_taxid {
        push_flag(&mut args, "-t");
    }
    args.push(input);
    args
}

fn passthrough_args(prog: &str, args: Vec<String>) -> Vec<String> {
    let mut full_args = vec![prog.to_string()];
    full_args.extend(args);
    full_args
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Classify {
            db,
            threads,
            quick,
            confidence,
            minimum_hit_groups,
            paired,
            single_file_pairs,
            classified_out,
            unclassified_out,
            output,
            report,
            use_mpa_style,
            report_zero_counts,
            report_minimizer_data,
            memory_mapping,
            use_names,
            minimum_base_quality,
            daemon,
            input,
        } => {
            let args = classify_args(
                db,
                threads,
                quick,
                confidence,
                minimum_hit_groups,
                paired,
                single_file_pairs,
                classified_out,
                unclassified_out,
                output,
                report,
                use_mpa_style,
                report_zero_counts,
                report_minimizer_data,
                memory_mapping,
                use_names,
                minimum_base_quality,
                daemon,
                input,
            );
            if let Err(e) = classify::classify_main(&args) {
                eprintln!("Error: {e}");
                process::exit(1);
            }
        }

        Commands::Build {
            hashtable,
            taxonomy,
            options,
            id_map,
            ncbi_taxonomy_dir,
            kmer_len,
            minimizer_len,
            capacity,
            max_capacity,
            threads,
            protein,
            fast_build,
            block_size,
            subblock_size,
            taxid_bits,
            spaced_seed,
            toggle_mask,
        } => {
            let args = build_db_args(
                hashtable,
                taxonomy,
                options,
                id_map,
                ncbi_taxonomy_dir,
                kmer_len,
                minimizer_len,
                capacity,
                max_capacity,
                threads,
                protein,
                fast_build,
                block_size,
                subblock_size,
                taxid_bits,
                spaced_seed,
                toggle_mask,
            );
            if let Err(e) = build_db::build_db_main(&args) {
                eprintln!("Error building database: {}", e);
                process::exit(1);
            }
        }

        Commands::Estimate {
            kmer_len,
            minimizer_len,
            n,
            protein,
            threads,
            block_size,
            spaced_seed,
            toggle_mask,
        } => {
            let args = estimate_args(
                kmer_len,
                minimizer_len,
                n,
                protein,
                threads,
                block_size,
                spaced_seed,
                toggle_mask,
            );
            match estimate::estimate_capacity_main(&args) {
                Ok(capacity) => println!("{}", capacity),
                Err(e) => {
                    eprintln!("Error: {}", e);
                    process::exit(1);
                }
            }
        }

        Commands::Inspect {
            db,
            memory_mapping,
            skip_counts,
            use_mpa_style,
            report_zero_counts,
            threads: _,
        } => {
            let hash_file = format!("{}/hash.k2d", db);
            let taxo_file = format!("{}/taxo.k2d", db);
            let opts_file = format!("{}/opts.k2d", db);

            // Load taxonomy
            let mut taxonomy =
                match kraken2::taxonomy::Taxonomy::from_file(&taxo_file, memory_mapping) {
                    Ok(t) => t,
                    Err(e) => {
                        eprintln!("Error loading taxonomy: {}", e);
                        process::exit(1);
                    }
                };
            taxonomy.generate_external_to_internal_id_map();

            // Load hash table
            let hash = match kraken2::compact_hash::CompactHashTable::from_file(
                &hash_file,
                memory_mapping,
            ) {
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
                        "/dev/stdout",
                        report_zero_counts,
                        &taxonomy,
                        &taxon_counters,
                    );
                } else {
                    let _ = kraken2::reports::report_kraken_style(
                        "/dev/stdout",
                        report_zero_counts,
                        false,
                        &taxonomy,
                        &taxon_counters,
                        total,
                        0,
                    );
                }
            }
        }

        Commands::Download {
            db,
            taxonomy,
            library,
            protein,
            skip_maps,
        } => {
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

        Commands::Blast2fasta {
            input,
            output,
            include_taxid,
        } => {
            let args = blast_to_fasta_args(input, output, include_taxid);
            if let Err(e) = blast::blast_to_fasta_main(&args) {
                eprintln!("Error converting BLAST database: {e}");
                process::exit(1);
            }
        }

        Commands::DumpTable { args } => {
            let args = passthrough_args("dump_table", args);
            if let Err(e) = dump_table::dump_table_main(&args) {
                eprintln!("Error dumping table: {e}");
                process::exit(1);
            }
        }

        Commands::LookupAccessionNumbers { args } => {
            let args = passthrough_args("lookup_accession_numbers", args);
            if let Err(e) = lookup::lookup_accession_numbers_main(&args) {
                eprintln!("Error looking up accession numbers: {e}");
                process::exit(1);
            }
        }

        Commands::K2mask { args } => {
            let args = passthrough_args("k2mask", args);
            if let Err(e) = dust::k2mask_main(&args) {
                eprintln!("Error masking low-complexity sequence: {e}");
                process::exit(1);
            }
        }

        Commands::Mmtest => {
            print!("{}", mmtest::mmtest_main());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_classify_args_routes_to_translated_flags() {
        let args = classify_args(
            "db".to_string(),
            4,
            true,
            0.2,
            3,
            true,
            false,
            Some("classified#.fq".to_string()),
            Some("unclassified#.fq".to_string()),
            Some("-".to_string()),
            Some("report.txt".to_string()),
            true,
            true,
            true,
            true,
            true,
            7,
            false,
            vec!["reads.fq".to_string()],
        );
        assert_eq!(args[0], "classify");
        assert!(args.contains(&"-H".to_string()));
        assert!(args.contains(&"-m".to_string()));
        assert!(args.contains(&"-K".to_string()));
        assert!(args.contains(&"reads.fq".to_string()));
    }

    #[test]
    fn test_build_db_args_routes_to_translated_flags() {
        let args = build_db_args(
            "hash.k2d".to_string(),
            "taxo.k2d".to_string(),
            "opts.k2d".to_string(),
            "seqid.map".to_string(),
            "taxonomy".to_string(),
            35,
            31,
            1000,
            500,
            2,
            true,
            true,
            4096,
            512,
            12,
            Some("101".to_string()),
            Some("11".to_string()),
        );
        assert_eq!(args[0], "build_db");
        assert!(args.contains(&"-F".to_string()));
        assert!(args.contains(&"-X".to_string()));
        assert!(args.contains(&"-M".to_string()));
    }

    #[test]
    fn test_estimate_args_routes_to_translated_flags() {
        let args = estimate_args(
            35,
            31,
            4,
            true,
            2,
            4096,
            Some("101".to_string()),
            Some("11".to_string()),
        );
        assert_eq!(args[0], "estimate_capacity");
        assert!(args.contains(&"-X".to_string()));
        assert!(args.contains(&"-S".to_string()));
        assert!(args.contains(&"-T".to_string()));
    }

    #[test]
    fn test_blast_to_fasta_args_routes_to_translated_flags() {
        let args =
            blast_to_fasta_args("db/core_nt.00".to_string(), "out.fna".to_string(), true);
        assert_eq!(args[0], "blast_to_fasta");
        assert!(args.contains(&"-o".to_string()));
        assert!(args.contains(&"-t".to_string()));
        assert!(args.contains(&"db/core_nt.00".to_string()));
    }

    #[test]
    fn test_passthrough_args_prefixes_program_name() {
        let args = passthrough_args("dump_table", vec!["-H".to_string(), "hash".to_string()]);
        assert_eq!(args, vec!["dump_table", "-H", "hash"]);
    }
}

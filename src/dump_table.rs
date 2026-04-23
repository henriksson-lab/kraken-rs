use std::io;

use crate::compact_hash::CompactHashTable;
use crate::readcounts::TaxonCounters;
use crate::reports;
use crate::taxonomy::Taxonomy;
use crate::types::IndexOptions;

#[derive(Clone, Default)]
pub struct Options {
    pub hashtable_filename: String,
    pub taxonomy_filename: String,
    pub options_filename: String,
    pub output_filename: String,
    pub use_mpa_style: bool,
    pub report_zeros: bool,
    pub skip_counts: bool,
    pub memory_mapping: bool,
    pub num_threads: usize,
}

pub fn mask2str(mask: u64, digits: i32) -> String {
    let mut s = String::new();
    for i in (0..digits).rev() {
        s.push(if ((mask >> i) & 1) != 0 { '1' } else { '0' });
    }
    s
}

pub fn parse_command_line(args: &[String], opts: &mut Options) -> io::Result<()> {
    let mut i = 1usize;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "-?" => return Err(io::Error::new(io::ErrorKind::InvalidInput, usage(0))),
            "-H" => {
                i += 1;
                opts.hashtable_filename = args.get(i).cloned().unwrap_or_default();
            }
            "-t" => {
                i += 1;
                opts.taxonomy_filename = args.get(i).cloned().unwrap_or_default();
            }
            "-o" => {
                i += 1;
                opts.options_filename = args.get(i).cloned().unwrap_or_default();
            }
            "-z" => opts.report_zeros = true,
            "-O" => {
                i += 1;
                opts.output_filename = args.get(i).cloned().unwrap_or_default();
            }
            "-m" => opts.use_mpa_style = true,
            "-s" => opts.skip_counts = true,
            "-p" => {
                i += 1;
                opts.num_threads = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(1);
            }
            "-M" => opts.memory_mapping = true,
            _ => {}
        }
        i += 1;
    }

    if opts.output_filename.is_empty() {
        opts.output_filename = "/dev/fd/1".to_string();
    }
    if opts.num_threads == 0 {
        opts.num_threads = 1;
    }

    if opts.hashtable_filename.is_empty()
        || opts.taxonomy_filename.is_empty()
        || opts.options_filename.is_empty()
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "missing mandatory filename parameter",
        ));
    }

    Ok(())
}

pub fn usage(_exit_code: i32) -> String {
    [
        "Usage: dump_table <options>",
        "",
        "Options (*mandatory):",
        "* -H FILENAME   Kraken 2 hash table filename",
        "* -t FILENAME   Kraken 2 taxonomy filename",
        "* -o FILENAME   Kraken 2 database options filename",
        "  -O FILENAME   Output filename (def: /dev/fd/1)",
        "  -m            Use MPA style output instead of Kraken 2 output",
        "  -M            Use memory mapping to access hash & taxonomy",
        "  -s            Skip reporting minimizer counts, just show DB stats",
        "  -z            Report taxa with zero counts",
    ]
    .join("\n")
}

pub fn dump_table_main(args: &[String]) -> io::Result<()> {
    let mut opts = Options {
        output_filename: "/dev/fd/1".to_string(),
        num_threads: 1,
        ..Default::default()
    };
    parse_command_line(args, &mut opts)?;

    let kraken_index = CompactHashTable::from_file(&opts.hashtable_filename, opts.memory_mapping)?;
    let mut taxonomy = Taxonomy::from_file(&opts.taxonomy_filename, opts.memory_mapping)?;
    taxonomy.generate_external_to_internal_id_map();
    let idx_opts = IndexOptions::read_from_file(&opts.options_filename)?;

    println!(
        "# Database options: {} db, k = {}, l = {}",
        if idx_opts.dna_db {
            "nucleotide"
        } else {
            "protein"
        },
        idx_opts.k,
        idx_opts.l
    );
    println!(
        "# Spaced mask = {}",
        mask2str(
            idx_opts.spaced_seed_mask,
            (idx_opts.l
                * if idx_opts.dna_db {
                    crate::types::BITS_PER_CHAR_DNA as usize
                } else {
                    crate::types::BITS_PER_CHAR_PRO as usize
                }) as i32
        )
    );
    println!("# Toggle mask = {}", mask2str(idx_opts.toggle_mask, 64));
    println!("# Total taxonomy nodes: {}", taxonomy.node_count());
    println!("# Table size: {}", kraken_index.size());
    println!("# Table capacity: {}", kraken_index.capacity());
    println!(
        "# Min clear hash value = {}",
        idx_opts.minimum_acceptable_hash_value
    );
    if idx_opts.revcom_version != crate::types::CURRENT_REVCOM_VERSION {
        println!("# Built with outdated revcom version");
    }

    if !opts.skip_counts {
        let taxid_counts = kraken_index.get_value_counts();
        let mut total_seqs = 0u64;
        let mut taxid_counters = TaxonCounters::new();
        for (taxid, count) in taxid_counts {
            total_seqs += u64::from(count);
            let counter = taxid_counters.entry(taxid).or_default();
            counter.n_reads = u64::from(count);
        }

        if opts.use_mpa_style {
            reports::report_mpa_style(
                &opts.output_filename,
                opts.report_zeros,
                &taxonomy,
                &taxid_counters,
            )?;
        } else {
            reports::report_kraken_style(
                &opts.output_filename,
                opts.report_zeros,
                false,
                &taxonomy,
                &taxid_counters,
                total_seqs,
                0,
            )?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mask2str() {
        assert_eq!(mask2str(0b1011, 4), "1011");
        assert_eq!(mask2str(0, 3), "000");
    }

    #[test]
    fn test_parse_command_line() {
        let args = vec![
            "dump_table".to_string(),
            "-H".to_string(),
            "hash".to_string(),
            "-t".to_string(),
            "taxo".to_string(),
            "-o".to_string(),
            "opts".to_string(),
            "-m".to_string(),
            "-z".to_string(),
        ];
        let mut opts = Options::default();
        parse_command_line(&args, &mut opts).unwrap();
        assert_eq!(opts.hashtable_filename, "hash");
        assert_eq!(opts.output_filename, "/dev/fd/1");
        assert!(opts.use_mpa_style);
        assert!(opts.report_zeros);
    }

    #[test]
    fn test_usage() {
        let text = usage(0);
        assert!(text.contains("Usage: dump_table <options>"));
        assert!(text.contains("-H FILENAME"));
    }
}

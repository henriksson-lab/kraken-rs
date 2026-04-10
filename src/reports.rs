use std::fs::File;
use std::io::{self, Write, BufWriter};

use crate::readcounts::{ReadCounter, TaxonCounters};
use crate::taxonomy::Taxonomy;
use crate::types::TaxonCounts;

/// Compute clade counts by propagating call counts up the taxonomy tree.
pub fn get_clade_counts(tax: &Taxonomy, call_counts: &TaxonCounts) -> TaxonCounts {
    let mut clade_counts = TaxonCounts::new();
    for (&taxid, &count) in call_counts {
        let mut tid = taxid;
        while tid != 0 {
            *clade_counts.entry(tid).or_insert(0) += count;
            tid = tax.node(tid).parent_id;
        }
    }
    clade_counts
}

/// Compute clade counters by propagating call counters up the taxonomy tree.
pub fn get_clade_counters(tax: &Taxonomy, call_counters: &TaxonCounters) -> TaxonCounters {
    let mut clade_counters = TaxonCounters::new();
    for (&taxid, counter) in call_counters {
        let mut tid = taxid;
        while tid != 0 {
            clade_counters.entry(tid).or_default().merge(counter);
            tid = tax.node(tid).parent_id;
        }
    }
    clade_counters
}

fn rank_to_mpa_code(rank: &str) -> Option<char> {
    match rank {
        "superkingdom" | "domain" => Some('d'),
        "kingdom" => Some('k'),
        "phylum" => Some('p'),
        "class" => Some('c'),
        "order" => Some('o'),
        "family" => Some('f'),
        "genus" => Some('g'),
        "species" => Some('s'),
        _ => None,
    }
}

fn rank_to_kraken_code(rank: &str) -> Option<char> {
    match rank {
        "superkingdom" | "domain" => Some('D'),
        "kingdom" => Some('K'),
        "phylum" => Some('P'),
        "class" => Some('C'),
        "order" => Some('O'),
        "family" => Some('F'),
        "genus" => Some('G'),
        "species" => Some('S'),
        _ => None,
    }
}

/// MPA-style report DFS.
fn mpa_report_dfs(
    taxid: u64,
    ofs: &mut dyn Write,
    report_zeros: bool,
    taxonomy: &Taxonomy,
    clade_counts: &TaxonCounts,
    taxonomy_names: &mut Vec<String>,
) {
    let clade_count = clade_counts.get(&taxid).copied().unwrap_or(0);
    if !report_zeros && clade_count == 0 {
        return;
    }

    let node = taxonomy.node(taxid);
    let rank = taxonomy.rank_at_offset(node.rank_offset);
    let rank_code = rank_to_mpa_code(rank);

    if let Some(code) = rank_code {
        let name = taxonomy.name_at_offset(node.name_offset);
        let entry = format!("{}__{}",  code, name);
        taxonomy_names.push(entry);

        let taxonomy_line = taxonomy_names.join("|");
        writeln!(ofs, "{}\t{}", taxonomy_line, clade_count).ok();
    }

    if node.child_count != 0 {
        let mut children: Vec<u64> = (0..node.child_count)
            .map(|i| node.first_child + i)
            .collect();
        children.sort_by(|a, b| {
            let ca = clade_counts.get(a).copied().unwrap_or(0);
            let cb = clade_counts.get(b).copied().unwrap_or(0);
            cb.cmp(&ca)
        });
        for child in children {
            mpa_report_dfs(child, ofs, report_zeros, taxonomy, clade_counts, taxonomy_names);
        }
    }

    if rank_code.is_some() {
        taxonomy_names.pop();
    }
}

/// Generate an MPA-style report.
pub fn report_mpa_style(
    filename: &str,
    report_zeros: bool,
    taxonomy: &Taxonomy,
    call_counters: &TaxonCounters,
) -> io::Result<()> {
    let mut call_counts = TaxonCounts::new();
    for (&taxid, counter) in call_counters {
        call_counts.insert(taxid, counter.read_count());
    }
    let clade_counts = get_clade_counts(taxonomy, &call_counts);

    let file = File::create(filename)?;
    let mut ofs = BufWriter::new(file);
    let mut taxonomy_names = Vec::new();
    mpa_report_dfs(1, &mut ofs, report_zeros, taxonomy, &clade_counts, &mut taxonomy_names);
    Ok(())
}

/// Print a single line of Kraken-style report.
fn print_kraken_style_report_line(
    ofs: &mut dyn Write,
    report_kmer_data: bool,
    total_seqs: u64,
    clade_counter: &ReadCounter,
    taxon_counter: &ReadCounter,
    rank_str: &str,
    taxid: u32,
    sci_name: &str,
    depth: usize,
) {
    let pct = if total_seqs > 0 {
        100.0 * clade_counter.read_count() as f64 / total_seqs as f64
    } else {
        0.0
    };

    // Match C++ sprintf("%6.2f", ...) exactly
    write!(ofs, "{:6.2}\t{}\t{}\t", pct, clade_counter.read_count(), taxon_counter.read_count()).ok();
    if report_kmer_data {
        write!(ofs, "{}\t{}\t", clade_counter.kmer_count(), clade_counter.distinct_kmer_count()).ok();
    }
    write!(ofs, "{}\t{}\t", rank_str, taxid).ok();
    for _ in 0..depth {
        write!(ofs, "  ").ok();
    }
    writeln!(ofs, "{}", sci_name).ok();
}

/// Kraken-style report DFS.
fn kraken_report_dfs(
    taxid: u64,
    ofs: &mut dyn Write,
    report_zeros: bool,
    report_kmer_data: bool,
    taxonomy: &Taxonomy,
    clade_counters: &TaxonCounters,
    call_counters: &TaxonCounters,
    total_seqs: u64,
    mut rank_code: char,
    mut rank_depth: i32,
    depth: usize,
) {
    let default_counter = ReadCounter::new();
    let clade_counter = clade_counters.get(&taxid).unwrap_or(&default_counter);
    let call_counter = call_counters.get(&taxid).unwrap_or(&default_counter);

    if !report_zeros && clade_counter.read_count() == 0 {
        return;
    }

    let node = taxonomy.node(taxid);
    let rank = taxonomy.rank_at_offset(node.rank_offset);

    if let Some(code) = rank_to_kraken_code(rank) {
        rank_code = code;
        rank_depth = 0;
    } else {
        rank_depth += 1;
    }

    let rank_str = if rank_depth != 0 {
        format!("{}{}", rank_code, rank_depth)
    } else {
        rank_code.to_string()
    };

    let name = taxonomy.name_at_offset(node.name_offset);

    print_kraken_style_report_line(
        ofs, report_kmer_data, total_seqs,
        clade_counter, call_counter,
        &rank_str, node.external_id as u32, name, depth,
    );

    if node.child_count != 0 {
        let mut children: Vec<u64> = (0..node.child_count)
            .map(|i| node.first_child + i)
            .collect();
        children.sort_by(|a, b| {
            let has_a = clade_counters.contains_key(a);
            let has_b = clade_counters.contains_key(b);
            if !has_a && !has_b {
                return std::cmp::Ordering::Equal;
            }
            if !has_a {
                return std::cmp::Ordering::Greater;
            }
            if !has_b {
                return std::cmp::Ordering::Less;
            }
            let ca = clade_counters[a].read_count();
            let cb = clade_counters[b].read_count();
            cb.cmp(&ca)
        });
        for child in children {
            kraken_report_dfs(
                child, ofs, report_zeros, report_kmer_data, taxonomy,
                clade_counters, call_counters, total_seqs,
                rank_code, rank_depth, depth + 1,
            );
        }
    }
}

/// Generate a Kraken-style report.
pub fn report_kraken_style(
    filename: &str,
    report_zeros: bool,
    report_kmer_data: bool,
    taxonomy: &Taxonomy,
    call_counters: &TaxonCounters,
    total_seqs: u64,
    total_unclassified: u64,
) -> io::Result<()> {
    let clade_counters = get_clade_counters(taxonomy, call_counters);

    let file = File::create(filename)?;
    let mut ofs = BufWriter::new(file);
    let rank_code = 'R';

    // Special handling of unclassified
    if total_unclassified != 0 || report_zeros {
        let rc = ReadCounter { n_reads: total_unclassified, ..Default::default() };
        print_kraken_style_report_line(
            &mut ofs, report_kmer_data, total_seqs,
            &rc, &rc, "U", 0, "unclassified", 0,
        );
    }

    kraken_report_dfs(
        1, &mut ofs, report_zeros, report_kmer_data, taxonomy,
        &clade_counters, call_counters, total_seqs,
        rank_code, -1, 0,
    );
    Ok(())
}

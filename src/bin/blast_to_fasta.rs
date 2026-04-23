use std::process;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if let Err(e) = kraken2::blast::blast_to_fasta_main(&args) {
        eprintln!("{e}");
        process::exit(1);
    }
}

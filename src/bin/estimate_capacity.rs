use std::process;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match kraken2::estimate::estimate_capacity_main(&args) {
        Ok(capacity) => println!("{capacity}"),
        Err(e) => {
            eprintln!("{e}");
            process::exit(1);
        }
    }
}

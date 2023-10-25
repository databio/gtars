use std::{fs::File, io::Write};
use std::path::Path;

use clap::{arg, Command};
use bedtools_rs::create_count_map;

const DEFAULT_OUTPUT: &str = "output.bed";

fn build_parser() -> Command {
    Command::new("bedtools-rs")
        .version("0.1.0")
        .author("Nathan LeRoy")
        .about("A universe pruning tool for use in genomic machine learning pipelines.")
        .arg(arg!(--universe <VALUE> "Path to the universe file we want to prune.").required(true))
        .arg(arg!(--data <VALUE> "Path to the training data. This should be a folder of bed files.").required(true))
        .arg(arg!(--output <VALUE> "Path to the output file (the pruned universe).").required(false))

}

fn main() {
    let app = build_parser();
    let matches = app.get_matches();
    let universe = matches.get_one::<String>("universe").expect("Universe path is required");
    let data = matches.get_one::<String>("data").expect("Data path is required");

    let default_output = format!("{}/{}", data, DEFAULT_OUTPUT);

    let output = matches.get_one::<String>("output").unwrap_or(&default_output);

    let cnt_map = create_count_map(data, universe).unwrap();

    // write the output
    let file_path = Path::new(output);
    let mut file = File::create(file_path).unwrap();
    for (region, _) in cnt_map {
        let line = format!("{}\t{}\t{}\n", region.chr, region.start, region.end);
        file.write_all(line.as_bytes()).unwrap();
    }

}
use std::io;
use std::io::Write;

use clap::{ArgMatches, Command, Arg};

use super::*;
use crate::common::models::region_set::RegionSet;
use crate::tokenizers::TreeTokenizer;

pub fn make_tokenization_cli() -> Command {
    Command::new(consts::TOKENIZE_CMD)
        .author("Nathan LeRoy")
        .about("Tokenize a bed file into a specific vocabulary.")
        .arg(
            Arg::new("bed")
                .long("bed")
                .short('b')
                .help("Path to the bed file we want to tokenize.")
                .required(true)
        )
        .arg(
            Arg::new("universe")
                .long("universe")
                .short('u')
                .help("Path to the universe file we want to use.")
                .required(true)
        )
        
}

pub mod handlers {
    
    use std::path::Path;

    use super::*;

    pub fn tokenize_bed_file(matches: &ArgMatches) {
        let bed = matches
            .get_one::<String>("bed")
            .expect("Bed file path is required");

        let universe = matches
            .get_one::<String>("universe")
            .expect("Universe path is required");

        // core logic/algorithm here
        let universe = Path::new(&universe);
        let tokenizer = TreeTokenizer::from(universe);

        let bed = Path::new(&bed);
        let regions = RegionSet::try_from(bed).expect("Failed to read bed file");
        
        let mut stdout = io::stdout().lock();

        let tokenized_regions = tokenizer.tokenize_region_set(&regions).expect("Could not tokenize region set.");

        for tokenized_region in tokenized_regions.into_iter() {
            let chr = tokenized_region.chr;
            let start = tokenized_region.start;
            let end = tokenized_region.end;
            
            
            let line = format!("{}\t{}\t{}\n", chr, start, end);

            // push to stdout
            stdout.write_all(line.as_bytes()).unwrap();
        }
    }
}
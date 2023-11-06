use std::io;
use std::io::Write;

use clap::{arg, ArgMatches, Command};
use vocab::create_count_map;

pub mod interfaces {
    
    use super::*;

    pub fn make_prune_cli() -> Command {
        Command::new(crate::consts::PRUNE_CMD)
        .author("Nathan LeRoy")
        .about("A universe pruning tool to remove regions below a minimum count threshold.")
        .arg(arg!(--universe <VALUE> "Path to the universe file we want to prune.").required(true))
        .arg(arg!(--data <VALUE> "Path to the training data. This should be a folder of bed files.").required(true))
        .arg(arg!(--"min-count" <VALUE> "Minimum number of overlaps required to keep a region.").required(false))
    }

    pub fn make_uniwig_cli() -> Command {
        Command::new(crate::consts::UNIWIG_CMD)
            .author("Nathan LeRoy")
            .about("Given a set of bed files, we want to produce 2 wiggle files: one is the track of start coordinates, the other is of end coordinates.")
            // uniwig bedfile stepsize smoothSize variableformat
            .arg(arg!(<BEDFILE> "Path to the bed file.").required(true))
            .arg(arg!(<STEPSIZE> "Step size.").required(false))
            .arg(arg!(<SMOOTH_SIZE> "Smooth size.").required(false))
            .arg(arg!(<VARIABLE_FORMAT> "Variable format.").required(false))
    }
}

pub mod functions {
    
    use super::*;

    pub fn prune_universe(matches: &ArgMatches) {
        let universe = matches
            .get_one::<String>("universe")
            .expect("Universe path is required");
    
        let data = matches
            .get_one::<String>("data")
            .expect("Data path is required");
    
        let min_count = matches
            .get_one::<u32>("min-count")
            .unwrap_or(&vocab::consts::DEFAULT_MIN_COUNT)
            .to_owned();
    
        let cnt_map = create_count_map(data, universe).unwrap();
    
        let mut stdout = io::stdout().lock();
    
        for (region, cnt) in cnt_map {
            if cnt < min_count {
                // skip this region
                continue;
            }
            let line = format!("{}\t{}\t{}\n", region.chr, region.start, region.end);
    
            // push to stdout
            stdout.write_all(line.as_bytes()).unwrap();
        }
    }

    pub fn uniwig(_matches: &ArgMatches) {
        println!("uniwig");
    }    
}

use std::io;
use std::io::Write;

use super::*;
use clap::{arg, ArgMatches, Command};

pub fn make_prune_cli() -> Command {
    Command::new(consts::PRUNE_CMD)
        .author("Nathan LeRoy")
        .about("A universe pruning tool to remove regions below a minimum count threshold.")
        .arg(arg!(--universe <VALUE> "Path to the universe file we want to prune.").required(true))
        .arg(
            arg!(--data <VALUE> "Path to the training data. This should be a folder of bed files.")
                .required(true),
        )
        .arg(
            arg!(--"min-count" <VALUE> "Minimum number of overlaps required to keep a region.")
                .required(false),
        )
}

pub mod handlers {

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
            .unwrap_or(&consts::DEFAULT_MIN_COUNT)
            .to_owned();

        // core logic/algorithm here
        let cnt_map = create_count_map(data, universe).unwrap();

        let mut stdout = io::stdout().lock();

        for (region, cnt) in cnt_map {
            if cnt < min_count {
                // skip this region
                continue;
            }

            let chr = region.chr;
            let start = region.start;
            let end = region.end;

            let line = format!("{}\t{}\t{}\n", chr, start, end);

            // push to stdout
            stdout.write_all(line.as_bytes()).unwrap();
        }
    }
}

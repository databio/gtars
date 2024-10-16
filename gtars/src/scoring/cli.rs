use std::collections::HashSet;
use std::io::BufRead;
use std::path::PathBuf;

use anyhow::Result;
use clap::{arg, Arg, ArgMatches, Command};

use super::*;
use crate::scoring::{region_scoring_from_fragments, ConsensusSet, FragmentFileGlob};

pub fn make_fscoring_cli() -> Command {
    Command::new(consts::FSCORING_CMD)
        .author("Nathan LeRoy")
        .about("Create a scoring matrix for a set of fragment files over a consensus peak set.")
        .arg(Arg::new("fragments"))
        .arg(Arg::new("consensus"))
        .arg(arg!(--mode <mode>))
        .arg(arg!(--output <output>))
        .arg(arg!(--whitelist <whitelist>))
}

pub mod handlers {

    use std::str::FromStr;

    use consts::DEFAULT_SCORING_MODE;

    use crate::common::utils::get_dynamic_reader;

    use super::*;

    pub fn region_fragment_scoring(matches: &ArgMatches) -> Result<()> {
        // get arguments from CLI
        let fragments = matches
            .get_one::<String>("fragments")
            .expect("A path to fragment files is required.");

        let consensus = matches
            .get_one::<String>("consensus")
            .expect("A path to a mapping file is required.");

        let default_out = consts::DEFAULT_OUT.to_string();
        let output = matches.get_one::<String>("output").unwrap_or(&default_out);
        let mode = match matches.get_one::<String>("mode") {
            Some(mode) => ScoringMode::from_str(mode),
            None => Ok(DEFAULT_SCORING_MODE),
        };
        let mode = mode.unwrap_or(DEFAULT_SCORING_MODE);

        let whitelist = matches.get_one::<String>("whitelist");

        // coerce arguments to types
        let mut fragments = FragmentFileGlob::new(fragments)?;
        let consensus = PathBuf::from(consensus);
        let consensus = ConsensusSet::new(consensus)?;

        let whitelist = match whitelist {
            Some(whitelist) => {
                // open whitelist and read to HashSet<String>
                let whitelist = PathBuf::from(whitelist);
                let reader = get_dynamic_reader(&whitelist)?;
                let mut whitelist: HashSet<String> = HashSet::new();
                for line in reader.lines() {
                    let line = line?;
                    if !whitelist.contains(&line) {
                        whitelist.insert(line);
                    }
                }
                Some(whitelist)
            }
            None => None,
        };

        region_scoring_from_fragments(
            &mut fragments,
            &consensus,
            output,
            whitelist.as_ref(),
            mode,
        )?;

        Ok(())
    }
}

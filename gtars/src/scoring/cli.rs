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
}

pub mod handlers {

    use std::str::FromStr;

    use consts::DEFAULT_SCORING_MODE;

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
            Some(mode) => {
                let supplied_mode = ScoringMode::from_str(mode);
                match supplied_mode {
                    Ok(mode) => mode,
                    Err(_err) => anyhow::bail!("Unknown scoring mode supplied: {}", mode),
                }
            }
            None => DEFAULT_SCORING_MODE,
        };

        // coerce arguments to types
        let mut fragments = FragmentFileGlob::new(fragments)?;
        let consensus = PathBuf::from(consensus);
        let consensus = ConsensusSet::new(consensus)?;

        let count_mat = region_scoring_from_fragments(&mut fragments, &consensus, mode)?;

        count_mat.write_to_file(output)?;

        Ok(())
    }
}

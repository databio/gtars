use std::path::PathBuf;
use std::str::FromStr;

use anyhow::Result;
use clap::ArgMatches;

use gtars_scoring::consts;
use gtars_scoring::{ConsensusSet, FragmentFileGlob, ScoringMode, region_scoring_from_fragments};

pub fn run_scoring(matches: &ArgMatches) -> Result<()> {
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
        None => consts::DEFAULT_SCORING_MODE,
    };

    // coerce arguments to types
    let mut fragments = FragmentFileGlob::new(fragments)?;
    let consensus = PathBuf::from(consensus);
    let consensus = ConsensusSet::new(consensus)?;

    let count_mat = region_scoring_from_fragments(&mut fragments, &consensus, mode)?;

    count_mat.write_to_file(output)?;

    Ok(())
}

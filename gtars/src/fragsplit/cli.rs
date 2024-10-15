use anyhow::Result;
use clap::{arg, Arg, ArgMatches, Command};

use super::*;
use crate::fragsplit::{pseudobulk_fragment_files, BarcodeToClusterMap};

pub fn make_fragsplit_cli() -> Command {
    Command::new(consts::FRAGSPLIT_CMD)
        .author("Nathan LeRoy")
        .about("Split fragment files into pseudobulks based on cluster labels.")
        .arg(Arg::new("fragments"))
        .arg(Arg::new("mapping"))
        .arg(arg!(--output <output>))
}

pub mod handlers {

    use std::path::Path;

    use super::*;

    pub fn split_fragment_files(matches: &ArgMatches) -> Result<()> {
        let fragments = matches
            .get_one::<String>("fragments")
            .expect("A path to fragment files is required.");

        let mapping = matches
            .get_one::<String>("mapping")
            .expect("A path to a mapping file is required.");

        let default_out = consts::DEFAULT_OUT.to_string();
        let output = matches.get_one::<String>("output").unwrap_or(&default_out);

        let fragments = Path::new(fragments);
        let mapping = &BarcodeToClusterMap::from_file(Path::new(mapping))?;
        let output = Path::new(output);

        pseudobulk_fragment_files(fragments, mapping, output)?;

        Ok(())
    }
}

use std::path::Path;

use anyhow::Result;
use clap::ArgMatches;

use gtars_fragsplit::consts::*;
use gtars_fragsplit::{BarcodeToClusterMap, pseudobulk_fragment_files};

pub fn run_fragsplit(matches: &ArgMatches) -> Result<()> {
    let fragments = matches
        .get_one::<String>("fragments")
        .expect("A path to fragment files is required.");

    let mapping = matches
        .get_one::<String>("mapping")
        .expect("A path to a mapping file is required.");

    let default_out = DEFAULT_OUT.to_string();
    let output = matches.get_one::<String>("output").unwrap_or(&default_out);

    let fragments = Path::new(fragments);
    let mapping = &BarcodeToClusterMap::from_file(Path::new(mapping))?;
    let output = Path::new(output);

    pseudobulk_fragment_files(fragments, mapping, output)?;

    Ok(())
}

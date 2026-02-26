use std::path::PathBuf;
use std::str::FromStr;

use anyhow::Result;
use clap::ArgMatches;

use gtars_scoring::{
    ConsensusSet, FragmentFileGlob, ScoringMode, barcode_scoring_from_fragments,
    region_scoring_from_fragments, write_sparse_counts_to_mtx,
};

use super::cli::{DEFAULT_OUT, DEFAULT_SCORING_MODE};

pub fn run_scoring(matches: &ArgMatches) -> Result<()> {
    // Check if barcode mode is enabled
    let barcode_mode = matches.get_flag("barcode");

    if barcode_mode {
        // Barcode-based scoring for single-cell data
        let fragment_file = matches
            .get_one::<String>("fragments")
            .expect("A path to fragment file is required.");

        let consensus_file = matches
            .get_one::<String>("consensus")
            .expect("A path to consensus peaks is required.");

        let default_output = "output".to_string();
        let output_prefix = matches
            .get_one::<String>("output")
            .unwrap_or(&default_output);

        // Load consensus peaks
        let consensus = ConsensusSet::new(PathBuf::from(consensus_file))?;

        // Count fragments by barcode (single-pass, sparse)
        let barcode_counts = barcode_scoring_from_fragments(fragment_file, &consensus)?;

        // Write directly to Matrix Market format
        write_sparse_counts_to_mtx(&barcode_counts, consensus.len(), output_prefix)?;

        println!(
            "Created {} cells × {} peaks sparse matrix",
            barcode_counts.len(),
            consensus.len()
        );
        println!(
            "Output files: {}_matrix.mtx.gz, {}_barcodes.tsv.gz, {}_features.tsv.gz",
            output_prefix, output_prefix, output_prefix
        );
    } else {
        // Original file-based scoring
        let fragments = matches
            .get_one::<String>("fragments")
            .expect("A path to fragment files is required.");

        let consensus = matches
            .get_one::<String>("consensus")
            .expect("A path to a mapping file is required.");

        let default_out = DEFAULT_OUT.to_string();
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
    }

    Ok(())
}

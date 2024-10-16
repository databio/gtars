use std::collections::HashSet;
use std::io::BufRead;
use std::str::FromStr;

use crate::common::models::Fragment;
use crate::common::utils::get_dynamic_reader;
use crate::fragsplit::utils::remove_all_extensions;
use crate::scoring::counts::CountMatrix;
use crate::scoring::files::FragmentFileGlob;
use crate::scoring::files::{ConsensusSet, FindOverlaps};

use anyhow::Result;
use indicatif::{ProgressBar, ProgressStyle};

type BarcodeWhiteList = HashSet<String>;

pub fn region_scoring_from_fragments(
    fragments: &mut FragmentFileGlob,
    consensus: &ConsensusSet,
    outfile: &str,
    barcode_whitelist: Option<&BarcodeWhiteList>,
) -> Result<()> {
    let binding = HashSet::new();
    let barcode_whitelist = barcode_whitelist.unwrap_or(&binding);

    let rows = fragments.len();
    let cols = consensus.len();

    let mut count_mat: CountMatrix<u32> = CountMatrix::new(rows, cols);

    let spinner = ProgressBar::new_spinner();
    spinner.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed}] {msg} ({per_sec})")
            .unwrap()
            .tick_strings(&["-", "\\", "|", "/"]),
    );

    spinner.set_message("Processing fragment files...");

    let mut processed_reads: u64 = 0;

    for (file_num, file) in fragments.into_iter().enumerate() {

        let reader = get_dynamic_reader(&file)?;
        let file_path = file.as_path();
        let file_stem = remove_all_extensions(file_path);

        for line in reader.lines() {
            let line = line?;
            let fragment = Fragment::from_str(&line)?;

            let whitelist_check_value = format!("{file_stem}+{}", fragment.barcode);
            // skip anything not in the whitelist
            if !barcode_whitelist.contains(&whitelist_check_value) {
                continue;
            }
            let olaps = consensus.find_overlaps(&fragment.into());
            if olaps.is_some() {
                let olaps = olaps.unwrap();
                for olap in olaps {
                    count_mat.increment(file_num, olap.1 as usize);
                }
            }

            // update the spinner
            processed_reads += 1;
            if processed_reads % 10_000 == 0 {
                spinner.set_message(format!("Processed {} reads", processed_reads));
            }
            spinner.inc(1);
        }
    }

    // write to a file
    count_mat.write_to_file(outfile)?;

    Ok(())
}
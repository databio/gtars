use std::io::BufRead;
use std::str::FromStr;

use crate::common::models::{Fragment, Region};
use crate::common::utils::get_dynamic_reader;
use crate::fragsplit::utils::remove_all_extensions;
use crate::scoring::consts::{END_SHIFT, START_SHIFT};
use crate::scoring::counts::CountMatrix;
use crate::scoring::files::FragmentFileGlob;
use crate::scoring::files::{ConsensusSet, FindOverlaps};
use crate::scoring::scoring_modes::ScoringMode;

use anyhow::Result;
use indicatif::{ProgressBar, ProgressStyle};

pub fn region_scoring_from_fragments(
    fragments: &mut FragmentFileGlob,
    consensus: &ConsensusSet,
    scoring_mode: ScoringMode,
) -> Result<CountMatrix<u32>> {
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

    spinner.set_message("Processing file...");

    let mut processed_reads: u64 = 0;
    let mut total_overlaps: u64 = 0;
    let total_fragments = fragments.len();

    for (file_num, file) in fragments.into_iter().enumerate() {
        let reader = get_dynamic_reader(&file)?;
        let file_path = file.as_path();
        let file_stem = remove_all_extensions(file_path);

        for line in reader.lines() {
            let line = line?;

            // convert to fragment and then get new positions of start and end
            let fragment = Fragment::from_str(&line)?;

            match scoring_mode {
                ScoringMode::Atac => {
                    let new_start = fragment.start + START_SHIFT as u32;
                    let new_end = fragment.end - END_SHIFT as u32;

                    let start_region = Region {
                        chr: fragment.chr,
                        start: new_start,
                        end: new_start + 1,
                    };

                    let olaps = consensus.find_overlaps(&start_region);
                    if let Some(olaps) = olaps {
                        total_overlaps += olaps.len() as u64;
                        for olap in olaps {
                            count_mat.increment(file_num, olap.1 as usize);
                        }
                    }

                    let end_region = Region {
                        // take from start_region to avoid a clone
                        chr: start_region.chr,
                        start: new_end,
                        end: new_end - 1,
                    };

                    let olaps = consensus.find_overlaps(&end_region);
                    if let Some(olaps) = olaps {
                        total_overlaps += olaps.len() as u64;
                        for olap in olaps {
                            count_mat.increment(file_num, olap.1 as usize);
                        }
                    }
                }
                ScoringMode::Chip => {
                    let olaps = consensus.find_overlaps(&fragment.into());
                    if let Some(olaps) = olaps {
                        total_overlaps += olaps.len() as u64;
                        for olap in olaps {
                            count_mat.increment(file_num, olap.1 as usize);
                        }
                    }
                }
            }

            // update the spinner
            processed_reads += 1;
            if processed_reads % 10_000 == 0 {
                spinner.set_message(format!(
                    "{file_stem} ({}/{total_fragments}) | {} overlaps | Processed {} reads",
                    file_num + 1,
                    total_overlaps,
                    processed_reads
                ));
            }
            spinner.inc(1);
        }
    }

    Ok(count_mat)
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[fixture]
    fn path_to_fragment_files() -> &'static str {
        "tests/data/fragments/region_scoring/*.bed.gz"
    }

    #[fixture]
    fn consensus_set() -> &'static str {
        "tests/data/consensus/consensus1.bed"
    }

    #[fixture]
    fn output_file() -> &'static str {
        "tests/data/out/region_scoring_count.csv.gz"
    }

    #[rstest]
    fn test_region_scoring_from_fragments_atac(
        path_to_fragment_files: &str,
        consensus_set: &str,
        output_file: &str,
    ) {
        let mut fragments = FragmentFileGlob::new(path_to_fragment_files).unwrap();
        let consensus = ConsensusSet::new(consensus_set.into()).unwrap();

        let res = region_scoring_from_fragments(&mut fragments, &consensus, ScoringMode::Atac);
        assert_eq!(res.is_ok(), true);

        let count_mat = res.unwrap();
        assert_eq!(count_mat.cols == 4, true);
        assert_eq!(count_mat.rows == 2, true);

        // Matrix should look like:
        // 2   2   1   3
        // 4   1   3   1
        // assert this is true
        assert_eq!(*count_mat.get(0, 0).unwrap(), 2);
        assert_eq!(*count_mat.get(0, 1).unwrap(), 2);
        assert_eq!(*count_mat.get(0, 2).unwrap(), 1);
        assert_eq!(*count_mat.get(0, 3).unwrap(), 3);

        assert_eq!(*count_mat.get(1, 0).unwrap(), 4);
        assert_eq!(*count_mat.get(1, 1).unwrap(), 1);
        assert_eq!(*count_mat.get(1, 2).unwrap(), 3);
        assert_eq!(*count_mat.get(1, 3).unwrap(), 1);

        let res = count_mat.write_to_file(output_file);
        assert_eq!(res.is_ok(), true);
    }
}

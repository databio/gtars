use std::io::BufRead;
use std::str::FromStr;

use crate::common::models::Fragment;
use crate::common::utils::get_dynamic_reader;
use crate::scoring::files::FragmentFileGlob;
use crate::scoring::files::ConsensusSet;
use crate::scoring::counts::CountMatrix;

use anyhow::Result;

pub fn region_scoring_from_fragments(fragments: &mut FragmentFileGlob, consensus: &ConsensusSet) -> Result<()> {

    let rows = fragments.len();
    let cols = consensus.len();

    let mut count_mat: CountMatrix<u32> = CountMatrix::new(rows, cols);

    for file in fragments.into_iter() {
        let reader = get_dynamic_reader(&file)?;
        for line in reader.lines() {
            let line = line?;
            let fragment = Fragment::from_str(&line)?;
        }
    }
    Ok(())
}

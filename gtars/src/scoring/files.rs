use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::Result;
use glob::glob;
use rust_lapper::{Interval, Lapper};

use crate::common::utils::{extract_regions_from_bed_file, generate_region_to_id_map};

pub struct FragmentFileGlob {
    curr: usize,
    files: Vec<PathBuf>,
}

#[allow(dead_code)]
pub struct ConsensusSet {
    overlap_trees: HashMap<String, Lapper<u32, u32>>,
}

impl FragmentFileGlob {
    pub fn new(pattern: &str) -> Result<Self> {
        let files = glob(pattern)?;
        let files = files
            .map(|f| match f {
                Ok(path) => Ok(path),
                Err(_) => anyhow::bail!(format!("Error reading file entry: {:?}", f)),
            })
            .collect::<Result<Vec<_>>>()?;
        let curr = 0_usize;
        Ok(FragmentFileGlob { files, curr })
    }
}

impl Iterator for FragmentFileGlob {
    type Item = PathBuf;
    fn next(&mut self) -> Option<Self::Item> {
        let result = self.files.get(self.curr).cloned();
        self.curr += 1;
        result
    }
}

impl ConsensusSet {
    pub fn new(path: PathBuf) -> Result<Self> {
        let regions = extract_regions_from_bed_file(&path)?;

        let mut trees: HashMap<String, Lapper<u32, u32>> = HashMap::new();
        let mut intervals: HashMap<String, Vec<Interval<u32, u32>>> = HashMap::new();

        let region_to_id_map = generate_region_to_id_map(&regions);

        for region in regions.iter() {
            // create interval
            let interval = Interval {
                start: region.start,
                stop: region.end,
                val: *region_to_id_map.get(region).unwrap(),
            };

            // use chr to get the vector of intervals
            let chr_intervals = intervals.entry(region.chr.clone()).or_default();

            // push interval to vector
            chr_intervals.push(interval);
        }

        // build the tree
        for (chr, chr_intervals) in intervals.into_iter() {
            let lapper: Lapper<u32, u32> = Lapper::new(chr_intervals);
            trees.insert(chr.to_string(), lapper);
        }

        Ok(ConsensusSet {
            overlap_trees: trees,
        })
    }
}
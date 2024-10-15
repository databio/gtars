use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::Result;
use glob::glob;
use rust_lapper::{Interval, Lapper};

use crate::common::models::Region;
use crate::common::utils::{extract_regions_from_bed_file, generate_region_to_id_map};

#[allow(unused)]
pub struct OverlapResult(Region, pub(crate) u32);

pub trait FindOverlaps {
    fn find_overlaps(&self, region: &Region) -> Option<Vec<OverlapResult>>;
}

pub struct FragmentFileGlob {
    curr: usize,
    files: Vec<PathBuf>,
}

pub struct ConsensusSet {
    len: usize,
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

    pub fn len(&self) -> usize {
        self.files.len()
    }

    pub fn is_empty(&self) -> bool {
        self.files.is_empty()
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
        let len = regions.len();

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

        // build the trees
        for (chr, chr_intervals) in intervals.into_iter() {
            let lapper: Lapper<u32, u32> = Lapper::new(chr_intervals);
            trees.insert(chr.to_string(), lapper);
        }

        Ok(ConsensusSet {
            overlap_trees: trees,
            len,
        })
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
}

impl FindOverlaps for ConsensusSet {
    fn find_overlaps(&self, region: &Region) -> Option<Vec<OverlapResult>> {
        let tree = self.overlap_trees.get(&region.chr);
        if tree.is_none() {
            None
        } else {
            let olaps = tree.unwrap().find(region.start, region.end);
            let olaps = olaps
                .into_iter()
                .map(|olap| {
                    OverlapResult(
                        Region {
                            chr: region.chr.clone(),
                            start: region.start,
                            end: region.end,
                        },
                        olap.val,
                    )
                })
                .collect();

            Some(olaps)
        }
    }
}

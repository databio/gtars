use std::collections::HashMap;
use std::path::Path;

use anyhow::Result;
use polars::prelude::*;
use rust_lapper::{Interval, Lapper};

use crate::common::consts::{UNKNOWN_CHR, UNKNOWN_END, UNKNOWN_START};
use crate::common::models::{Region, RegionSet, TokenizedRegionSet, Universe};
use crate::common::utils::extract_regions_from_bed_file;
use crate::tokenizers::traits::Tokenizer;

pub struct TreeTokenizer {
    pub universe: Universe,
    pub tree: HashMap<String, Lapper<u32, ()>>,
}

impl TryFrom<&Path> for TreeTokenizer {
    type Error = anyhow::Error;
    ///
    /// # Arguments
    /// - `value` - the path to the bed file
    ///
    /// # Returns
    /// A new TreeTokenizer
    fn try_from(value: &Path) -> Result<Self> {
        let universe = Universe::try_from(value)?;

        let mut tree: HashMap<String, Lapper<u32, ()>> = HashMap::new();
        let mut intervals: HashMap<String, Vec<Interval<u32, ()>>> = HashMap::new();

        for region in universe.regions.iter() {
            // create interval
            let interval = Interval {
                start: region.start,
                stop: region.end,
                val: (),
            };

            // use chr to get the vector of intervals
            let chr_intervals = intervals.entry(region.chr.to_owned()).or_default();

            // push interval to vector
            chr_intervals.push(interval);
        }

        for (chr, chr_intervals) in intervals.iter() {
            let lapper: Lapper<u32, ()> = Lapper::new(chr_intervals.to_owned());
            tree.insert(chr.to_string(), lapper);
        }

        Ok(TreeTokenizer { universe, tree })
    }
}

impl Tokenizer for TreeTokenizer {
    fn tokenize_region(&self, region: &Region) -> TokenizedRegionSet {
        let lapper = self.tree.get(&region.chr);
        match lapper {
            Some(lapper) => {
                let intervals = lapper.find(region.start, region.end);
                let regions: Vec<Region> = intervals
                    .map(|interval| Region {
                        chr: region.chr.to_owned(),
                        start: interval.start,
                        end: interval.stop,
                    })
                    .collect();

                if regions.is_empty() {
                    let regions = vec![self.unknown_token()];
                    return TokenizedRegionSet {
                        regions,
                        universe: &self.universe,
                    };
                }

                TokenizedRegionSet {
                    regions,
                    universe: &self.universe,
                }
            }
            None => TokenizedRegionSet {
                regions: vec![Region {
                    chr: UNKNOWN_CHR.to_string(),
                    start: UNKNOWN_START as u32,
                    end: UNKNOWN_END as u32,
                }],
                universe: &self.universe,
            },
        }
    }

    fn tokenize_region_set(&self, region_set: &RegionSet) -> Option<TokenizedRegionSet> {
        let chrs = region_set.chrs();
        let starts = region_set.starts();
        let ends = region_set.ends();

        let mut tokenized_regions: Vec<Region> = Vec::new();

        for i in 0..region_set.len() {
            let chr: String;
            let start: u32;
            let end: u32;

            if let AnyValue::Utf8(v) = chrs.get(i).unwrap() {
                chr = v.to_string();
            } else {
                println!(
                    "chr column must be of type Utf8, instead found {:?}",
                    chrs.get(i).unwrap()
                );
                return None;
            }

            if let AnyValue::UInt32(v) = starts.get(i).unwrap() {
                start = v;
            } else {
                println!(
                    "start column must be of type UInt32, instead found {:?}",
                    starts.get(i).unwrap()
                );
                return None;
            }

            if let AnyValue::UInt32(v) = ends.get(i).unwrap() {
                end = v;
            } else {
                println!(
                    "end column must be of type UInt32, instead found {:?}",
                    ends.get(i).unwrap()
                );
                return None;
            }

            let lapper = self.tree.get(&chr);
            match lapper {
                Some(tree) => {
                    let intervals = tree.find(start, end);

                    let regions: Vec<Region> = intervals
                        .map(|interval| Region {
                            chr: chr.to_owned(),
                            start: interval.start,
                            end: interval.stop,
                        })
                        .collect();

                    if regions.is_empty() {
                        tokenized_regions.push(self.unknown_token());
                        continue;
                    }

                    tokenized_regions.extend(regions);
                }
                None => {
                    tokenized_regions.push(self.unknown_token());
                }
            }
        }

        Some(TokenizedRegionSet {
            regions: tokenized_regions,
            universe: &self.universe,
        })
    }
}

impl TreeTokenizer {
    pub fn tokenize_bed_file(&self, bed_file: &Path) -> Option<TokenizedRegionSet> {
        let regions = extract_regions_from_bed_file(bed_file);
        match regions {
            Ok(regions) => {
                let rs = RegionSet::from(regions);
                self.tokenize_region_set(&rs)
            }
            Err(e) => panic!("Error reading bedfile: {}", e),
        }
    }
}

use std::collections::HashMap;
use std::path::Path;

use anyhow::Result;
use polars::prelude::*;
use rust_lapper::{Interval, Lapper};

use crate::common::consts::special_tokens::*;
use crate::common::models::{Region, RegionSet, TokenizedRegionSet, Universe};
use crate::common::utils::extract_regions_from_bed_file;
use crate::tokenizers::traits::{Pad, SpecialTokens, Tokenizer};

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
        let mut universe = Universe::try_from(value)?;

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

        // add special tokens to the universe
        // unk
        universe.insert_token(&Region {
            chr: UNKNOWN_CHR.to_string(),
            start: UNKNOWN_START as u32,
            end: UNKNOWN_END as u32,
        });

        // pad
        universe.insert_token(&Region {
            chr: PAD_CHR.to_string(),
            start: PAD_START as u32,
            end: PAD_END as u32,
        });

        // mask
        universe.insert_token(&Region {
            chr: MASK_CHR.to_string(),
            start: MASK_START as u32,
            end: MASK_END as u32,
        });

        // eos
        universe.insert_token(&Region {
            chr: EOS_CHR.to_string(),
            start: EOS_START as u32,
            end: EOS_END as u32,
        });

        // bos
        universe.insert_token(&Region {
            chr: BOS_CHR.to_string(),
            start: BOS_START as u32,
            end: BOS_END as u32,
        });

        // cls
        universe.insert_token(&Region {
            chr: CLS_CHR.to_string(),
            start: CLS_START as u32,
            end: CLS_END as u32,
        });

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

    fn tokenize_region_set(&self, region_set: &RegionSet) -> TokenizedRegionSet {
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
                panic!(
                    "chr column must be of type Utf8, instead found {:?}",
                    chrs.get(i).unwrap()
                );
            }

            if let AnyValue::UInt32(v) = starts.get(i).unwrap() {
                start = v;
            } else {
                panic!(
                    "start column must be of type UInt32, instead found {:?}",
                    starts.get(i).unwrap()
                );
            }

            if let AnyValue::UInt32(v) = ends.get(i).unwrap() {
                end = v;
            } else {
                panic!(
                    "end column must be of type UInt32, instead found {:?}",
                    ends.get(i).unwrap()
                );
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

        TokenizedRegionSet {
            regions: tokenized_regions,
            universe: &self.universe,
        }
    }

    fn vocab_size(&self) -> usize {
        self.universe.len()
    }

    fn convert_token_to_id(&self, token: &crate::common::models::TokenizedRegion) -> u32 {
        self.universe
            .convert_region_to_id(&Region::from(token))
            .unwrap()
    }

    fn convert_id_to_token(&self, id: u32) -> crate::common::models::TokenizedRegion {
        // TODO: this is a slow, naive approach, but I'll put it here anyways...
        let regions: Vec<Region> = self
            .universe
            .region_to_id
            .iter()
            .filter_map(|(key, &val)| if val == id { Some(key.clone()) } else { None })
            .collect();

        regions.0
    }
}

impl SpecialTokens for TreeTokenizer {
    fn unknown_token(&self) -> Region {
        Region {
            chr: UNKNOWN_CHR.to_string(),
            start: UNKNOWN_START as u32,
            end: UNKNOWN_END as u32,
        }
    }

    fn padding_token(&self) -> Region {
        Region {
            chr: PAD_CHR.to_string(),
            start: PAD_START as u32,
            end: PAD_END as u32,
        }
    }

    fn mask_token(&self) -> Region {
        Region {
            chr: MASK_CHR.to_string(),
            start: MASK_START as u32,
            end: MASK_END as u32,
        }
    }

    fn cls_token(&self) -> Region {
        Region {
            chr: CLS_CHR.to_string(),
            start: CLS_START as u32,
            end: CLS_END as u32,
        }
    }

    fn bos_token(&self) -> Region {
        Region {
            chr: BOS_CHR.to_string(),
            start: BOS_START as u32,
            end: BOS_END as u32,
        }
    }

    fn eos_token(&self) -> Region {
        Region {
            chr: EOS_CHR.to_string(),
            start: EOS_START as u32,
            end: EOS_END as u32,
        }
    }

    fn sep_token(&self) -> Region {
        Region {
            chr: SEP_CHR.to_string(),
            start: SEP_START as u32,
            end: SEP_END as u32,
        }
    }
}

impl TreeTokenizer {
    pub fn tokenize_bed_file(&self, bed_file: &Path) -> Result<TokenizedRegionSet> {
        let regions = extract_regions_from_bed_file(bed_file)?;
        let rs = RegionSet::from(regions);

        Ok(self.tokenize_region_set(&rs))
    }
}

// use default implementation
impl Pad for TreeTokenizer {}

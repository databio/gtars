use std::collections::HashMap;
use std::path::Path;

use anyhow::Result;
use rust_lapper::{Interval, Lapper};

use crate::common::consts::special_tokens::*;
use crate::common::models::{Region, RegionSet, TokenizedRegionSet, Universe};
use crate::common::utils::extract_regions_from_bed_file;
use crate::tokenizers::traits::{Pad, SpecialTokens, Tokenizer};

pub struct TreeTokenizer {
    pub universe: Universe,
    tree: HashMap<String, Lapper<u32, u32>>,
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

        // sep
        universe.insert_token(&Region {
            chr: SEP_CHR.to_string(),
            start: SEP_START as u32,
            end: SEP_END as u32,
        });

        let mut tree: HashMap<String, Lapper<u32, u32>> = HashMap::new();
        let mut intervals: HashMap<String, Vec<Interval<u32, u32>>> = HashMap::new();

        for region in universe.regions.iter() {
            // create interval
            let interval = Interval {
                start: region.start,
                stop: region.end,
                val: universe.convert_region_to_id(region).unwrap(),
            };

            // use chr to get the vector of intervals
            let chr_intervals = intervals.entry(region.chr.to_owned()).or_default();

            // push interval to vector
            chr_intervals.push(interval);
        }

        for (chr, chr_intervals) in intervals.iter() {
            let lapper: Lapper<u32, u32> = Lapper::new(chr_intervals.to_owned());
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
                let ids: Vec<u32> = intervals.map(|interval| interval.val).collect();

                if ids.is_empty() {
                    let ids = vec![self.unknown_token_id()];
                    return TokenizedRegionSet {
                        ids,
                        universe: &self.universe,
                    };
                }

                TokenizedRegionSet {
                    ids,
                    universe: &self.universe,
                }
            }
            None => TokenizedRegionSet {
                ids: vec![self.unknown_token_id()],
                universe: &self.universe,
            },
        }
    }

    fn tokenize_region_set(&self, region_set: &RegionSet) -> TokenizedRegionSet {
        let mut tokenized_regions: Vec<u32> = Vec::new();

        for region in region_set {
            let lapper = self.tree.get(&region.chr);

            match lapper {
                Some(tree) => {
                    let intervals = tree.find(region.start, region.end);

                    let regions: Vec<u32> = intervals.map(|interval| interval.val).collect();

                    if regions.is_empty() {
                        tokenized_regions.push(self.unknown_token_id());
                        continue;
                    }

                    tokenized_regions.extend(regions);
                }
                None => {
                    tokenized_regions.push(self.unknown_token_id());
                }
            }
        }

        TokenizedRegionSet {
            ids: tokenized_regions,
            universe: &self.universe,
        }
    }

    fn vocab_size(&self) -> usize {
        self.universe.len()
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

    fn unknown_token_id(&self) -> u32 {
        self.universe
            .convert_region_to_id(&self.unknown_token())
            .unwrap()
    }

    fn padding_token_id(&self) -> u32 {
        self.universe
            .convert_region_to_id(&self.padding_token())
            .unwrap()
    }

    fn mask_token_id(&self) -> u32 {
        self.universe
            .convert_region_to_id(&self.mask_token())
            .unwrap()
    }

    fn cls_token_id(&self) -> u32 {
        self.universe
            .convert_region_to_id(&self.cls_token())
            .unwrap()
    }

    fn bos_token_id(&self) -> u32 {
        self.universe
            .convert_region_to_id(&self.bos_token())
            .unwrap()
    }

    fn eos_token_id(&self) -> u32 {
        self.universe
            .convert_region_to_id(&self.eos_token())
            .unwrap()
    }

    fn sep_token_id(&self) -> u32 {
        self.universe
            .convert_region_to_id(&self.sep_token())
            .unwrap()
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

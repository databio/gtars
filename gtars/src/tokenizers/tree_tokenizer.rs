use std::collections::HashMap;
use std::fs::read_to_string;
use std::path::Path;

use anyhow::Result;
use rust_lapper::Lapper;

use crate::common::consts::special_tokens::*;
use crate::common::models::{Region, RegionSet, TokenizedRegionSet, Universe};
use crate::common::utils::{create_interval_tree_from_universe, extract_regions_from_bed_file};
use crate::tokenizers::config::TokenizerConfig;
use crate::tokenizers::traits::{Pad, SpecialTokens, Tokenizer};

pub struct TreeTokenizer {
    pub universe: Universe,
    tree: HashMap<String, Lapper<u32, u32>>,
    secondary_trees: Option<Vec<HashMap<String, Lapper<u32, u32>>>>,
    exclude_ranges: Option<HashMap<String, Lapper<u32, u32>>>,
}

impl TryFrom<&Path> for TreeTokenizer {
    type Error = anyhow::Error;
    ///
    /// # Arguments
    /// - `value` - the path to the tokenizer config file (a TOML) or bed file
    ///
    /// # Returns
    /// A new TreeTokenizer
    fn try_from(value: &Path) -> Result<Self> {
        // detect file type... if ends in toml assume toml otherwise assume its a bed file
        // and just build the universe + tree from that and move on.
        //
        // This maintains backwards compatibility with the old way of creating tokenizers from bed files
        // and allows for the new way of creating tokenizers from toml files
        let file_extension = value.extension().unwrap().to_str().unwrap();

        let (mut universe, tree, secondary_trees, exclude_ranges) = match file_extension {
            // parse config file
            "toml" => {
                let toml_str = read_to_string(value)?;
                let config: TokenizerConfig = toml::from_str(&toml_str)?;

                // universe path is relative to the config file
                let universe_path = value.parent().unwrap().join(&config.universe);

                // create initial universe from the *required* universe field
                let mut universe = Universe::try_from(Path::new(&universe_path))?;

                let tree = create_interval_tree_from_universe(&universe);

                // create secondary trees if they exist
                let secondary_trees = match config.hierarchical_universes {
                    Some(hierarchical_universes) => {
                        let mut secondary_trees = Vec::new();
                        for hierarchical_universe in hierarchical_universes {
                            let hierarchical_universe_path =
                                value.parent().unwrap().join(&hierarchical_universe);

                            let hierarchical_universe_regions =
                                extract_regions_from_bed_file(&hierarchical_universe_path)?;

                            for region in hierarchical_universe_regions {
                                universe.insert_token(&region);
                            }

                            let hierarchical_tree = create_interval_tree_from_universe(&universe);

                            secondary_trees.push(hierarchical_tree);
                        }

                        Some(secondary_trees)
                    }
                    None => None,
                };

                // create exclude ranges if they exist
                let exclude_ranges = match config.exclude_ranges {
                    Some(exclude_ranges) => {
                        let exclude_ranges_path = value.parent().unwrap().join(exclude_ranges);

                        // universe gets discarded since its not conasidered a part of the tokenizers universe
                        let exclude_ranges_universe =
                            Universe::try_from(exclude_ranges_path.as_path())?;

                        let exclude_ranges_map =
                            create_interval_tree_from_universe(&exclude_ranges_universe);

                        Some(exclude_ranges_map)
                    }

                    None => None,
                };

                (universe, tree, secondary_trees, exclude_ranges)
            }
            // else assume its a bed file
            _ => {
                let regions = extract_regions_from_bed_file(value)?;
                let universe = Universe::from(regions);
                let tree = create_interval_tree_from_universe(&universe);
                (universe, tree, None, None)
            }
        };

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

        Ok(TreeTokenizer {
            universe,
            tree,
            secondary_trees,
            exclude_ranges,
        })
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

    fn get_universe(&self) -> &Universe {
        &self.universe
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

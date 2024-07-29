use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};
use rust_lapper::{Interval, Lapper};

use crate::common::consts::special_tokens::*;
use crate::common::models::{Region, RegionSet, TokenizedRegionSet, Universe};
use crate::common::utils::{create_interval_tree_from_universe, extract_regions_from_bed_file};
use crate::tokenizers::config::TokenizerConfig;
use crate::tokenizers::traits::{Pad, SpecialTokens, Tokenizer};

///
/// The TreeTokenizer is a basic tokenizer that can "tokenize" genomic regions
/// into a known universe (or vocabulary). This is especially useful as a
/// pre-processor for machine learning pipelines
pub struct TreeTokenizer {
    pub universe: Universe,
    config: TokenizerConfig,
    tree: HashMap<String, Lapper<u32, u32>>,
    secondary_trees: Option<Vec<HashMap<String, Lapper<u32, u32>>>>,
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

        let (config, mut universe, tree, secondary_trees, _exclude_ranges) = match file_extension {
            // parse config file
            "toml" => {
                let config = TokenizerConfig::try_from(value).with_context(|| {
                    format!(
                        "Invalid tokenizer configuration found for file: {}",
                        value.to_str().unwrap()
                    )
                })?;

                if config.universes.is_empty() {
                    anyhow::bail!(
                        "You must have at least one universe in your universe list. Found zero."
                    )
                }

                let primary_universe = &config.universes[0];
                let other_universes = match config.universes.len() {
                    1 => None,
                    _ => Some(&config.universes[1..]),
                };

                // universe path is relative to the config file
                let universe_path = value.parent().unwrap().join(primary_universe);

                // create initial universe from the *required* universe field
                let mut universe = Universe::try_from(Path::new(&universe_path))?;

                let tree = create_interval_tree_from_universe(&universe);

                // create secondary trees if they exist
                let secondary_trees = match other_universes {
                    Some(hierarchical_universes) => {
                        let mut secondary_trees = Vec::new();
                        for hierarchical_universe in hierarchical_universes {
                            let mut hierarchical_tree: HashMap<String, Lapper<u32, u32>> =
                                HashMap::new();

                            let hierarchical_universe_path =
                                value.parent().unwrap().join(hierarchical_universe);

                            let hierarchical_universe_regions =
                                extract_regions_from_bed_file(&hierarchical_universe_path)?;

                            let mut intervals: HashMap<String, Vec<Interval<u32, u32>>> =
                                HashMap::new();
                            for region in hierarchical_universe_regions {
                                universe.insert_token(&region);
                                let interval = Interval {
                                    start: region.start,
                                    stop: region.end,
                                    val: universe.convert_region_to_id(&region).unwrap(),
                                };

                                intervals
                                    .entry(region.chr.clone())
                                    .or_default()
                                    .push(interval);
                            }

                            for (chr, chr_intervals) in intervals.iter() {
                                let lapper: Lapper<u32, u32> =
                                    Lapper::new(chr_intervals.to_owned());
                                hierarchical_tree.insert(chr.to_string(), lapper);
                            }

                            secondary_trees.push(hierarchical_tree);
                        }

                        Some(secondary_trees)
                    }
                    None => None,
                };

                // create exclude ranges if they exist
                let exclude_ranges = match &config.exclude_ranges {
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

                (config, universe, tree, secondary_trees, exclude_ranges)
            }
            // else assume its a bed file
            _ => {
                let regions = extract_regions_from_bed_file(value)?;
                let universe = Universe::from(regions);
                let tree = create_interval_tree_from_universe(&universe);

                let universe_as_path = Path::new(value).file_name().unwrap();
                let universe_as_path = universe_as_path.to_string_lossy().to_string();

                let config =
                    TokenizerConfig::new(Some("tree".to_string()), vec![universe_as_path], None);
                (config, universe, tree, None, None)
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
            config,
            universe,
            tree,
            secondary_trees,
        })
    }
}

impl Tokenizer for TreeTokenizer {
    fn tokenize_region(&self, region: &Region) -> TokenizedRegionSet {
        let lapper = self.tree.get(&region.chr);

        match lapper {
            Some(lapper) => {
                let intervals = lapper.find(region.start, region.end);
                let mut ids: Vec<u32> = intervals.map(|interval| interval.val).collect();

                // tokenized to nothing... check secondary trees
                if ids.is_empty() {
                    // oh, we have no secondary trees, just return the unknown token
                    if self.secondary_trees.is_none() {
                        ids = vec![self.unknown_token_id()];
                    // iterate over secondary trees and check if the region is in any of them
                    } else {
                        for s_tree in self.secondary_trees.as_ref().unwrap() {
                            // default to unknown token
                            ids = vec![self.unknown_token_id()];

                            let s_lapper = s_tree.get(&region.chr);
                            if s_lapper.is_none() {
                                continue;
                            }
                            // get overlapped intervals -- map to regions
                            let intervals = s_lapper.unwrap().find(region.start, region.end);
                            let regions: Vec<u32> =
                                intervals.map(|interval| interval.val).collect();

                            // a hit
                            if !regions.is_empty() {
                                ids = regions;
                                break;
                            }
                        }
                    }
                }

                TokenizedRegionSet {
                    ids,
                    universe: &self.universe,
                }
            }
            // primary universe didnt have that chromosome/contig/seqname
            // so, check secondary trees
            None => {
                let mut ids = Vec::new();
                // oh, we have no secondary trees, just return the unknown token
                if self.secondary_trees.is_none() {
                    ids = vec![self.unknown_token_id()];
                // iterate over secondary trees and check if the region is in any of them
                } else {
                    for s_tree in self.secondary_trees.as_ref().unwrap() {
                        // default to unknown token
                        ids = vec![self.unknown_token_id()];

                        let s_lapper = s_tree.get(&region.chr);
                        if s_lapper.is_none() {
                            continue;
                        }

                        // get overlapped intervals -- map to regions
                        let intervals = s_lapper.unwrap().find(region.start, region.end);
                        let regions: Vec<u32> = intervals.map(|interval| interval.val).collect();

                        // a hit
                        if !regions.is_empty() {
                            ids = regions;
                            break;
                        } else {
                            ids = vec![self.unknown_token_id()];
                        }
                    }
                }

                TokenizedRegionSet {
                    ids,
                    universe: &self.universe,
                }
            }
        }
    }

    fn tokenize_region_set(&self, region_set: &RegionSet) -> TokenizedRegionSet {
        let mut tokenized_regions: Vec<u32> = Vec::new();

        for region in region_set {
            let tokenized_region = self.tokenize_region(region);
            tokenized_regions.extend(tokenized_region.ids);
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

    fn export(&self, path: &Path) -> Result<()> {
        let toml_str = toml::to_string(&self.config)?;
        std::fs::write(path, toml_str)?;
        Ok(())
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

#[cfg(test)]
mod tests {

    use crate::common::models::{Region, RegionSet};
    use crate::tokenizers::traits::SpecialTokens;
    use std::path::Path;

    use super::*;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[fixture]
    fn path_to_bed_file() -> &'static str {
        "tests/data/peaks.bed"
    }

    #[fixture]
    fn path_to_config_file() -> &'static str {
        "tests/data/tokenizer.toml"
    }

    #[fixture]
    fn path_to_bad_config_file() -> &'static str {
        "tests/data/tokenizer_bad.toml"
    }

    #[fixture]
    fn path_to_tokenize_bed_file() -> &'static str {
        "tests/data/to_tokenize.bed"
    }

    #[rstest]
    fn test_create_tokenizer_from_bed(path_to_bed_file: &str) {
        let tokenizer = TreeTokenizer::try_from(Path::new(path_to_bed_file)).unwrap();
        assert_eq!(tokenizer.vocab_size(), 32); // 25 regions + 7 special tokens
    }

    #[rstest]
    fn test_create_tokenizer_from_config(path_to_config_file: &str) {
        let tokenizer = TreeTokenizer::try_from(Path::new(path_to_config_file)).unwrap();
        assert_eq!(tokenizer.vocab_size(), 56); // 25 regions in main universe + 24 in hierarchical + 7 special tokens
    }

    #[rstest]
    #[should_panic]
    fn test_bad_config_file(path_to_bad_config_file: &str) {
        let tokenizer = TreeTokenizer::try_from(Path::new(path_to_bad_config_file));
        let _tokenizer = tokenizer.unwrap();
    }

    #[rstest]
    fn test_get_special_token_ids(path_to_bed_file: &str) {
        let tokenizer = TreeTokenizer::try_from(Path::new(path_to_bed_file)).unwrap();
        let unk_id = tokenizer.unknown_token_id();
        let pad_id = tokenizer.padding_token_id();
        let mask_id = tokenizer.mask_token_id();
        let eos_id = tokenizer.eos_token_id();
        let bos_id = tokenizer.bos_token_id();
        let cls_id = tokenizer.cls_token_id();
        let sep_id = tokenizer.sep_token_id();

        assert_eq!(unk_id, 25);
        assert_eq!(pad_id, 26);
        assert_eq!(mask_id, 27);
        assert_eq!(eos_id, 28);
        assert_eq!(bos_id, 29);
        assert_eq!(cls_id, 30);
        assert_eq!(sep_id, 31);
    }

    #[rstest]
    fn test_tokenize_bed_file(path_to_bed_file: &str, path_to_tokenize_bed_file: &str) {
        let tokenizer = TreeTokenizer::try_from(Path::new(path_to_bed_file)).unwrap();
        let rs = RegionSet::try_from(Path::new(path_to_tokenize_bed_file)).unwrap();
        let tokenized_regions = tokenizer.tokenize_region_set(&rs);

        println!("{}", tokenized_regions.len());
        assert_eq!(tokenized_regions.len(), 4);

        // last should be the unknown token
        let unknown_token = tokenizer
            .universe
            .convert_id_to_region(tokenized_regions[3])
            .unwrap();
        assert!(unknown_token.chr == "chrUNK");
    }

    #[rstest]
    fn test_hierarchical_universe_hit(path_to_config_file: &str) {
        let tokenizer = TreeTokenizer::try_from(Path::new(path_to_config_file)).unwrap();
        let res = tokenizer.tokenize_region(&Region {
            chr: "chr1".to_string(),
            start: 100,
            end: 200,
        });
        assert_eq!(res.len(), 1);

        // check the id, it should be len(primary_universe) + 1 (since its chr1)
        assert_eq!(res.ids, vec![25]);

        let res = res.into_region_vec();
        let region = &res[0];

        assert_eq!(region.chr, "chr1");
        assert_eq!(region.start, 0);
        assert_eq!(region.end, 248_956_422);
    }

    #[rstest]
    fn test_hierarchical_universe_no_hit(path_to_config_file: &str) {
        let tokenizer = TreeTokenizer::try_from(Path::new(path_to_config_file)).unwrap();
        let res = tokenizer.tokenize_region(&Region {
            chr: "chrFOO".to_string(),
            start: 100,
            end: 200,
        });
        assert_eq!(res.len(), 1);

        // check the id, it should be the id of the UNK token
        assert_eq!(res.ids, vec![49]);

        let res = res.into_region_vec();
        let region = &res[0];

        assert_eq!(region.chr, "chrUNK");
        assert_eq!(region.start, 0);
        assert_eq!(region.end, 0);
    }
}

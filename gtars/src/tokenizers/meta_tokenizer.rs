use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use anyhow::{Context, Result};
use rust_lapper::{Interval, Lapper};

use crate::common::consts::special_tokens::*;
use crate::common::models::{Region, Universe};
use crate::common::utils::get_dynamic_reader;
use crate::tokenizers::TokenizerConfig;

///
/// The MetaTokenizer is a TreeTokenizer that implements the concept of meta-tokens. Meta
/// tokens are a way to reduce the size of the vocabulary for genomic interval-based
/// machine learning models.
///
/// In brief, meta-tokens are tokens that represent *clusters* of genomic intervals.
pub struct MetaTokenizer {
    pub universe: Universe,
    region_to_metatoken: HashMap<Region, Region>,
    tree: HashMap<String, Lapper<u32, u32>>,
    secondary_trees: Option<Vec<HashMap<String, Lapper<u32, u32>>>>,
}

impl TryFrom<&Path> for MetaTokenizer {
    type Error = anyhow::Error;

    ///
    /// # Arguments
    /// - `value` - the path to the tokenizer config file (a TOML) or bed file
    ///
    /// # Returns
    /// A new TreeTokenizer
    fn try_from(value: &Path) -> Result<Self, Self::Error> {
        let config = TokenizerConfig::new(value).with_context(|| {
            format!(
                "Invalid tokenizer configuration found for file: {}",
                value.to_str().unwrap()
            )
        })?;

        // verify the config is good
        if config.universes.is_empty() {
            anyhow::bail!("You must have at least one universe in your universe list. Found zero.")
        }

        // get priority list of universes
        let primary_universe = &config.universes[0];
        let other_universes = match config.universes.len() {
            1 => None,
            _ => Some(&config.universes[1..]),
        };

        let primary_universe = value.parent().unwrap().join(primary_universe);

        // parse first universe
        let reader = get_dynamic_reader(Path::new(&primary_universe))?;
        let mut universe = Universe::default();
        let mut intervals: HashMap<String, Vec<Interval<u32, u32>>> = HashMap::new();
        let mut region_to_metatoken: HashMap<Region, Region> = HashMap::new();

        let mut seen_metatokens: HashMap<String, u32> = HashMap::new();

        for line in reader.lines() {
            let line = line?;

            let fields: Vec<&str> = line.split('\t').collect();

            // check length of fields
            if fields.len() < 4 {
                anyhow::bail!("BED file line does not have at least 4 fields: {}", line);
            }

            // parse the fields
            let chr = fields[0];
            let start = fields[1].parse::<u32>().with_context(|| {
                format!("Failed to parse start position in BED file line: {}", line)
            })?;

            let end = fields[2].parse::<u32>().with_context(|| {
                format!("Failed to parse end position in BED file line: {}", line)
            })?;
            
            // why is primary_ being prepended to the metatoken id?
            // - this is a way to ensure that the metatoken id is unique,
            // imagine a secondary universe that has the same metatoken id
            let meta_id = format!("primary_{}", fields[3]);

            // get the id for the metatoken if we've seen it before
            // else create a new id and insert it into the hashmap
            let meta_id = match seen_metatokens.get(&meta_id) {
                Some(id) => *id,
                None => {
                    let id = seen_metatokens.len() as u32;
                    seen_metatokens.insert(meta_id, id);
                    id
                }
            };            

            // construct the actual region
            let region = Region {
                chr: chr.to_string(),
                start,
                end,
            };

            // construct the mapped meta token
            let meta_region = Region {
                chr: format!("chrM{}", meta_id),
                start: 0,
                end: 0,
            };

            // update the universe with the metatoken
            universe.insert_token(&meta_region);

            // insert a region into the appropriate list
            let ilist = intervals.entry(region.chr.clone()).or_default();
            ilist.push(Interval {
                start: region.start,
                stop: region.end,
                val: universe.convert_region_to_id(&meta_region).unwrap(),
            });

            // insert the region into the meta token map
            region_to_metatoken.insert(region, meta_region);
        }

        let mut tree: HashMap<String, Lapper<u32, u32>> = HashMap::new();

        for (chr, chr_intervals) in intervals.into_iter() {
            let lapper: Lapper<u32, u32> = Lapper::new(chr_intervals);
            tree.insert(chr, lapper);
        }

        let secondary_trees = match other_universes {
            None => None,
            Some(other_universes) => {
                let mut secondary_trees = Vec::new();

                for (u_num, other_universe) in other_universes.iter().enumerate() {

                    let other_universe = value.parent().unwrap().join(other_universe);

                    let reader = get_dynamic_reader(Path::new(&other_universe))?;
                    let mut intervals: HashMap<String, Vec<Interval<u32, u32>>> = HashMap::new();

                    for line in reader.lines() {
                        let line = line?;

                        let fields: Vec<&str> = line.split('\t').collect();

                        // check length of fields
                        if fields.len() < 4 {
                            anyhow::bail!(
                                "BED file line does not have at least 4 fields: {}",
                                line
                            );
                        }

                        // parse the fields
                        let chr = fields[0];
                        let start = fields[1].parse::<u32>().with_context(|| {
                            format!("Failed to parse start position in BED file line: {}", line)
                        })?;

                        let end = fields[2].parse::<u32>().with_context(|| {
                            format!("Failed to parse end position in BED file line: {}", line)
                        })?;

                        let meta_id = format!("secondary_{}_{}", u_num, fields[3]);

                        let meta_id = match seen_metatokens.get(&meta_id) {
                            Some(id) => *id,
                            None => {
                                let id = seen_metatokens.len() as u32;
                                seen_metatokens.insert(meta_id, id);
                                id
                            }
                        };

                        // construct the actual region
                        let region = Region {
                            chr: chr.to_string(),
                            start,
                            end,
                        };

                        // extract meta region id
                        let meta_region = Region {
                            chr: format!("chrM{}", meta_id),
                            start: 0,
                            end: 0,
                        };

                        // update the universe with the metatoken
                        universe.insert_token(&meta_region);

                        // insert a region into the appropriate list
                        let ilist = intervals.entry(region.chr.clone()).or_default();
                        ilist.push(Interval {
                            start: region.start,
                            stop: region.end,
                            val: universe.convert_region_to_id(&meta_region).unwrap(),
                        });

                        // insert the region into the meta token map
                        region_to_metatoken.insert(region, meta_region);
                    }

                    let mut tree: HashMap<String, Lapper<u32, u32>> = HashMap::new();

                    for (chr, chr_intervals) in intervals.into_iter() {
                        let lapper: Lapper<u32, u32> = Lapper::new(chr_intervals);
                        tree.insert(chr, lapper);
                    }

                    secondary_trees.push(tree);
                }

                Some(secondary_trees)
            }
        };

        // now we can insert the special tokens
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

        Ok(MetaTokenizer {
            universe,
            region_to_metatoken,
            tree,
            secondary_trees,
        })
    }
}


// tests
#[cfg(test)]
mod tests {

    use super::*;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[fixture]
    fn path_to_config_file() -> &'static str {
        "tests/data/tokenizer.meta.toml"
    }

    #[fixture]
    fn path_to_tokenize_bed_file() -> &'static str {
        "tests/data/to_tokenize.bed"
    }

    #[rstest]
    fn test_create_tokenizer(path_to_config_file: &str) {
        let tokenizer = MetaTokenizer::try_from(Path::new(path_to_config_file)).unwrap();
        assert_eq!(tokenizer.universe.len(), 27);
    }
}
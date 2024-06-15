use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use anyhow::{Context, Result};
use rust_lapper::{Interval, Lapper};

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

        // parse first universe
        let reader = get_dynamic_reader(Path::new(primary_universe))?;
        let mut universe = Universe::default();
        let mut intervals: HashMap<String, Vec<Interval<u32, u32>>> = HashMap::new();
        let mut region_to_metatoken: HashMap<Region, Region> = HashMap::new();

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

            let meta_id = fields[3]
                .parse::<u32>()
                .with_context(|| format!("Failed to parse meta ID in BED file line: {}", line))?;

            // construct the actual region
            let region = Region {
                chr: chr.to_string(),
                start,
                end,
            };

            // construct the mapped meta token
            let meta_region = Region {
                chr: meta_id.to_string(),
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

                for other_universe in other_universes {
                    let reader = get_dynamic_reader(Path::new(other_universe))?;
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

                        let meta_id = fields[3].parse::<u32>().with_context(|| {
                            format!("Failed to parse meta ID in BED file line: {}", line)
                        })?;

                        // construct the actual region
                        let region = Region {
                            chr: chr.to_string(),
                            start,
                            end,
                        };

                        // TODO: this is actually not right... secondary universes
                        // shouldnt have to be aware of others. So, this might
                        // be the wrong meta token id.
                        // we need to keep track of meta tokens that
                        // already exist and increment from there.
                        // construct the mapped meta token
                        let meta_region = Region {
                            chr: meta_id.to_string(),
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

        Ok(MetaTokenizer {
            universe,
            region_to_metatoken,
            tree,
            secondary_trees,
        })
    }
}

use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use anyhow::Result;

use crate::common::consts::special_tokens::{PAD_CHR, PAD_END, PAD_START};
use crate::common::models::region::Region;
use crate::common::models::tokenized_region::TokenizedRegion;
use crate::common::models::universe::Universe;
use crate::io::write_tokens_to_gtok;

use super::RegionSet;

pub struct TokenizedRegionSet<'a> {
    pub ids: Vec<u32>,
    pub universe: &'a Universe,
}

impl Into<RegionSet> for TokenizedRegionSet<'_> {
    fn into(self) -> RegionSet {
        let regions: Vec<Region> = self
            .ids
            .iter()
            .map(|id| self.universe.regions[*id as usize])
            .collect();

        RegionSet::from(regions)
    }
}

impl<'a> IntoIterator for &'a TokenizedRegionSet<'_> {
    type Item = TokenizedRegion;
    type IntoIter = std::vec::IntoIter<TokenizedRegion>;

    fn into_iter(self) -> Self::IntoIter {
        let mut tokenized_regions = Vec::with_capacity(self.regions.len());
        for region in self.regions.iter() {
            // TODO: is unwrapping here the smartest thing in the world?
            let id = self.universe.convert_region_to_id(region).unwrap();

            let tokenized_region: TokenizedRegion = TokenizedRegion {
                chr: region.chr.to_owned(),
                start: region.start,
                end: region.end,
                id,
            };

            tokenized_regions.push(tokenized_region);
        }
        tokenized_regions.into_iter()
    }
}

impl<'a> TokenizedRegionSet<'a> {
    ///
    /// Create a new TokenizedRegionSet. the TokenizedRegionSet takes
    /// a reference to a Universe.
    ///
    /// # Arguments
    /// * `regions` - A vector of regions
    /// * `universe` - A reference to a Universe
    ///
    pub fn new(regions: Vec<Region>, universe: &'a Universe) -> Self {
        TokenizedRegionSet { regions, universe }
    }

    ///
    /// Write a TokenizedRegionSet to a BED file
    ///
    /// # Arguments
    /// * `path` - A PathBuf to write the BED file to
    ///
    pub fn to_bed_file(&self, path: &PathBuf) -> Result<()> {
        let mut file = File::create(path)?;
        for region in self.regions.iter() {
            let line = format!(
                "{}\t{}\t{}\n",
                region.chr.to_owned(),
                region.start,
                region.end
            );
            file.write_all(line.as_bytes())?;
        }
        Ok(())
    }

    ///
    /// Write a TokenizedRegionSet to a .gtok file
    /// * `path` - A PathBuf to write the .gtok file to
    ///
    pub fn to_gtok_file(&self, path: &str) -> Result<()> {
        let tokens = self.to_region_ids();
        write_tokens_to_gtok(path, &tokens)?;
        Ok(())
    }

    pub fn ids(&self) -> Vec<u32> {
        self.ids
    }

    pub fn into_region_set(self) -> RegionSet {
        self.into()
    }
}

impl<'a> TokenizedRegionSet<'a> {
    pub fn from(regions: Vec<Region>, universe: &'a Universe) -> Self {
        TokenizedRegionSet { regions, universe }
    }

    pub fn len(&self) -> usize {
        self.ids.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }
}

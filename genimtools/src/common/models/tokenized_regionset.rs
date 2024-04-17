use std::path::PathBuf;

use anyhow::Result;

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

impl Into<Vec<u8>> for TokenizedRegionSet<'_> {
    fn into(self) -> Vec<u8> {
        let mut bit_vector: Vec<u8> = Vec::with_capacity(self.universe.len());

        for id in self.ids {
            bit_vector[id as usize] = 1;
        }

        bit_vector
    }
}

impl<'a> IntoIterator for &'a TokenizedRegionSet<'_> {
    type Item = TokenizedRegion<'a>;
    type IntoIter = std::vec::IntoIter<TokenizedRegion<'a>>;

    fn into_iter(self) -> Self::IntoIter {
        let mut tokenized_regions = Vec::with_capacity(self.ids.len());
        for id in self.ids.iter() {
            let tokenized_region: TokenizedRegion = TokenizedRegion {
                universe: &self.universe,
                id: *id,
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
    pub fn new(ids: Vec<u32>, universe: &'a Universe) -> Self {
        TokenizedRegionSet { ids, universe }
    }

    ///
    /// Write a TokenizedRegionSet to a BED file
    ///
    /// # Arguments
    /// * `path` - A PathBuf to write the BED file to
    ///
    pub fn to_bed_file(&self, path: &PathBuf) -> Result<()> {
        self.into_region_set().to_bed(path)
    }

    ///
    /// Write a TokenizedRegionSet to a .gtok file
    /// * `path` - A PathBuf to write the .gtok file to
    ///
    pub fn to_gtok_file(&self, path: &str) -> Result<()> {
        let tokens = self.ids;
        write_tokens_to_gtok(path, &tokens)?;
        Ok(())
    }

    pub fn ids(&self) -> Vec<u32> {
        self.ids
    }

    pub fn into_region_set(self) -> RegionSet {
        self.into()
    }

    pub fn into_bit_vector(self) -> Vec<u8> {
        self.into()
    }
}

impl<'a> TokenizedRegionSet<'a> {
    pub fn len(&self) -> usize {
        self.ids.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }
}

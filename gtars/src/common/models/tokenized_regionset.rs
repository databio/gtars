use std::fs::File;
use std::io::Write;
use std::ops::Index;
use std::path::PathBuf;

use anyhow::Result;

use crate::common::models::region::Region;
use crate::common::models::tokenized_region::TokenizedRegion;
use crate::common::models::universe::Universe;
use crate::io::write_tokens_to_gtok;

use super::RegionSet;

///
/// A tokenized region set is a specific representation of a region set. It is
/// two things: 1) a list of ids, and 2) a pointer to a Universe. The ids correspond
/// to the regions in that universe this [TokenizedRegionSet] represents.
pub struct TokenizedRegionSet<'a> {
    pub ids: Vec<u32>,
    pub universe: &'a Universe,
}

impl From<TokenizedRegionSet<'_>> for RegionSet {
    fn from(val: TokenizedRegionSet<'_>) -> Self {
        let regions: Vec<Region> = val
            .ids
            .iter()
            .map(|id| val.universe.regions[*id as usize].clone())
            .collect();

        RegionSet::from(regions)
    }
}

impl From<TokenizedRegionSet<'_>> for Vec<u8> {
    fn from(val: TokenizedRegionSet<'_>) -> Self {
        let mut bit_vector: Vec<u8> = Vec::with_capacity(val.universe.len());

        for id in val.ids {
            bit_vector[id as usize] = 1;
        }

        bit_vector
    }
}

impl From<TokenizedRegionSet<'_>> for Vec<Region> {
    fn from(value: TokenizedRegionSet<'_>) -> Self {
        value
            .ids
            .iter()
            .map(|id| value.universe.id_to_region[id].to_owned())
            .collect()
    }
}

impl<'a> Index<usize> for TokenizedRegionSet<'a> {
    type Output = u32;

    fn index(&self, index: usize) -> &Self::Output {
        &self.ids[index]
    }
}

impl<'a> IntoIterator for &'a TokenizedRegionSet<'_> {
    type Item = TokenizedRegion<'a>;
    type IntoIter = std::vec::IntoIter<TokenizedRegion<'a>>;

    fn into_iter(self) -> Self::IntoIter {
        let mut tokenized_regions = Vec::with_capacity(self.ids.len());
        for id in self.ids.iter() {
            let tokenized_region: TokenizedRegion = TokenizedRegion {
                universe: self.universe,
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
        let mut file = File::create(path)?;
        for id in self.ids.iter() {
            let r = self.universe.id_to_region.get(id).unwrap();
            let line = format!("{}\t{}\t{}\n", r.chr, r.start, r.end);
            file.write_all(line.as_bytes())?;
        }

        Ok(())
    }

    ///
    /// Write a TokenizedRegionSet to a .gtok file
    /// * `path` - A PathBuf to write the .gtok file to
    ///
    pub fn to_gtok_file(&self, path: &str) -> Result<()> {
        let tokens = &self.ids;
        write_tokens_to_gtok(path, tokens)?;
        Ok(())
    }

    ///
    /// Get the tokenized regions as a vector of ids
    /// * Returns a vector of u32
    pub fn ids(&self) -> Vec<u32> {
        self.ids.clone()
    }

    ///
    /// Get the tokenized regions into a dedicated RegionSet
    /// * Returns a RegionSet
    pub fn into_region_set(self) -> RegionSet {
        self.into()
    }

    ///
    /// Get the tokenized regions as a vector of u8
    /// * Returns a vector of u8
    pub fn into_bit_vector(self) -> Vec<u8> {
        self.into()
    }

    ///
    /// Get the tokenized regions as a vector of regions
    /// * Returns a vector of regions
    pub fn into_region_vec(self) -> Vec<Region> {
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

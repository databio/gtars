use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use crate::common::consts::{PAD_CHR, PAD_END, PAD_START};
use crate::common::models::region::Region;
use crate::common::models::tokenized_region::TokenizedRegion;
use crate::common::models::universe::Universe;

pub struct TokenizedRegionSet<'a> {
    pub regions: Vec<Region>,
    pub universe: &'a Universe,
}

impl<'a> IntoIterator for &'a TokenizedRegionSet<'_> {
    type Item = TokenizedRegion;
    type IntoIter = std::vec::IntoIter<TokenizedRegion>;

    fn into_iter(self) -> Self::IntoIter {
        let mut tokenized_regions = Vec::with_capacity(self.regions.len());
        for region in self.regions.iter() {
            let id = self.universe.convert_region_to_id(region);

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
    pub fn to_bed_file(&self, path: &PathBuf) -> Result<(), Box<dyn Error>> {
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
    /// Convert a TokenizedRegionSet to a vector of region IDs
    ///
    pub fn to_region_ids(&self) -> Vec<u32> {
        let mut region_ids = Vec::new();
        for region in &self.regions {
            region_ids.push(self.universe.convert_chr_start_end_to_id(
                &region.chr,
                region.start,
                region.end,
            ));
        }
        region_ids
    }

    ///
    /// Pad a tokenized region set
    ///
    pub fn pad(&mut self, len: usize) {
        // this is wrong: the padding token might not be in the universe
        let pad_region = Region {
            chr: PAD_CHR.to_string(),
            start: PAD_START as u32,
            end: PAD_END as u32,
        };

        while self.regions.len() < len {
            self.regions.push(pad_region.clone());
        }
    }
}

impl<'a> TokenizedRegionSet<'a> {
    pub fn from(regions: Vec<Region>, universe: &'a Universe) -> Self {
        TokenizedRegionSet { regions, universe }
    }

    pub fn len(&self) -> usize {
        self.regions.len()
    }

    pub fn is_empty(&self) -> bool {
        self.regions.is_empty()
    }
}

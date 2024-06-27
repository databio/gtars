use crate::common::utils::extract_regions_from_bed_file;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use anyhow::Result;

use crate::common::models::Region;

pub struct RegionSet {
    pub regions: Vec<Region>,
}

pub struct RegionSetIterator<'a> {
    region_set: &'a RegionSet,
    index: usize,
}

impl TryFrom<&Path> for RegionSet {
    type Error = anyhow::Error;

    ///
    /// Create a new [RegionSet] from a bed file.
    ///
    /// # Arguments:
    /// - value: path to bed file on disk.
    fn try_from(value: &Path) -> Result<Self> {
        let regions = extract_regions_from_bed_file(value)?;
        Ok(RegionSet { regions })
    }
}

impl From<Vec<Region>> for RegionSet {
    fn from(regions: Vec<Region>) -> Self {
        RegionSet { regions }
    }
}

impl From<&[u8]> for RegionSet {
    fn from(value: &[u8]) -> Self {
        let region_str = String::from_utf8_lossy(value);
        let regions: Vec<Region> = region_str
            .split('\n')
            .map(|line| {
                let parts = line.split('\t').collect::<Vec<&str>>();

                let chr = parts[0].to_string();
                let start = parts[1].parse::<u32>().unwrap();
                let end = parts[2].parse::<u32>().unwrap();

                Region { chr, start, end }
            })
            .collect();

        RegionSet { regions }
    }
}

impl<'a> Iterator for RegionSetIterator<'a> {
    type Item = &'a Region;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.region_set.regions.len() {
            let region = &self.region_set.regions[self.index];
            self.index += 1;
            Some(region)
        } else {
            None
        }
    }
}

impl<'a> IntoIterator for &'a RegionSet {
    type Item = &'a Region;
    type IntoIter = RegionSetIterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        RegionSetIterator {
            region_set: self,
            index: 0,
        }
    }
}

impl RegionSet {
    ///
    /// Dump a regionset to disk
    ///
    /// # Arguments
    /// - path: the path to the file to dump to
    pub fn to_bed(&self, path: &Path) -> Result<()> {
        let mut file = File::create(path)?;
        // is there a better way to do this?
        for region in self.regions.iter() {
            let chr = region.chr.clone();
            let start = region.start;
            let end = region.end;

            let line = format!("{}\t{}\t{}\n", chr, start, end);
            file.write_all(line.as_bytes())?;
        }

        Ok(())
    }

    pub fn len(&self) -> usize {
        self.regions.len()
    }

    pub fn is_empty(&self) -> bool {
        self.regions.is_empty()
    }
}

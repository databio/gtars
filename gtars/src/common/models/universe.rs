use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};

use crate::common::models::region::Region;
use crate::common::utils::{
    extract_regions_from_bed_file, generate_id_to_region_map, generate_region_to_id_map,
};

#[derive(Clone, Eq, PartialEq, Default)]
pub struct Universe {
    pub regions: Vec<Region>,
    pub region_to_id: HashMap<Region, u32>,
    pub id_to_region: HashMap<u32, Region>,
    pub source_chrom_to_id: HashMap<String, u16>,
}

impl Universe {
    pub fn insert_token(&mut self, region: &Region) {
        let new_id = self.region_to_id.len();
        self.region_to_id.insert(region.to_owned(), new_id as u32);
        self.id_to_region.insert(new_id as u32, region.to_owned());
        self.regions.push(region.to_owned());
    }

    pub fn insert_chrom(&mut self, region: &Region) {
        if !self.source_chrom_to_id.contains_key(&region.chr) {
            let new_id = self.source_chrom_to_id.len();
            self.source_chrom_to_id.insert(region.chr.clone(), new_id as u16);
        }
    }


    pub fn convert_region_to_id(&self, region: &Region) -> Option<u32> {
        let id = self.region_to_id.get(region);
        id.map(|id| id.to_owned())
    }

    pub fn convert_chr_start_end_to_id(&self, chr: &str, start: u32, end: u32) -> Option<u32> {
        let region = Region {
            chr: chr.to_string(),
            start,
            end,
        };
        self.convert_region_to_id(&region)
    }

    pub fn convert_id_to_region(&self, id: u32) -> Option<Region> {
        self.id_to_region.get(&id).cloned()
    }

    pub fn convert_source_chrom_to_id(&self, chr: &str) -> Option<u16> {
        self.source_chrom_to_id.get(chr).copied()
    }

    pub fn len(&self) -> usize {
        self.region_to_id.len()
    }

    pub fn num_chroms(&self) -> usize {
        self.source_chrom_to_id.len()
    }

    pub fn is_empty(&self) -> bool {
        self.region_to_id.len() == 0
    }

    pub fn contains_region(&self, region: &Region) -> bool {
        self.region_to_id.contains_key(region)
    }
}

impl From<Vec<Region>> for Universe {
    fn from(value: Vec<Region>) -> Self {
        // make a copy of the regions
        let regions = value;

        // create the region to id map and add the Unk token if it doesn't exist
        let region_to_id = generate_region_to_id_map(&regions);
        let id_to_region = generate_id_to_region_map(&regions);

        let mut current_chrom_id = 0_u16;
        let mut source_chrom_to_id: HashMap<String, u16> = HashMap::new();
        for region in regions.iter() {
            source_chrom_to_id.entry(region.chr.clone()).or_insert_with(|| {
                let old_id = current_chrom_id;
                current_chrom_id += 1;
                old_id
            });
        }

        Universe {
            regions,
            region_to_id,
            id_to_region,
            source_chrom_to_id,
        }
    }
}

impl TryFrom<&Path> for Universe {
    type Error = anyhow::Error;

    fn try_from(value: &Path) -> Result<Self> {
        let regions = extract_regions_from_bed_file(value)
            .with_context(|| "There was an error reading the universe file!")?;

        let region_to_id = generate_region_to_id_map(&regions);
        let id_to_region = generate_id_to_region_map(&regions);

        let mut current_chrom_id = 0_u16;
        let mut source_chrom_to_id: HashMap<String, u16> = HashMap::new();
        for region in regions.iter() {
            source_chrom_to_id.entry(region.chr.clone()).or_insert_with(|| {
                let old_id = current_chrom_id;
                current_chrom_id += 1;
                old_id
            });
        }

        Ok(Universe {
            regions,
            region_to_id,
            id_to_region,
            source_chrom_to_id,
        })
    }
}

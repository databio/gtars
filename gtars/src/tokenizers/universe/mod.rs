use std::collections::HashMap;
use std::path::Path;

use crate::common::models::region::Region;
use crate::common::models::region_set::RegionSet;
use crate::common::utils::{generate_id_to_region_map, generate_region_to_id_map};

use super::utils::special_tokens::SpecialTokens;

#[derive(Clone, Eq, PartialEq, Default)]
pub struct Universe {
    pub regions: Vec<Region>,
    pub region_to_id: HashMap<Region, u32>,
    pub id_to_region: HashMap<u32, Region>,
}

impl Universe {
    pub fn add_token_to_universe(&mut self, region: &Region) {
        let new_id = self.region_to_id.len();
        self.region_to_id.insert(region.to_owned(), new_id as u32);
        self.id_to_region.insert(new_id as u32, region.to_owned());
        self.regions.push(region.to_owned());
    }

    pub fn convert_region_to_id(&self, region: &Region) -> Option<u32> {
        let id = self.region_to_id.get(region);
        id.map(|id| id.to_owned())
    }

    pub fn convert_id_to_region(&self, id: u32) -> Option<Region> {
        self.id_to_region.get(&id).cloned()
    }

    pub fn len(&self) -> usize {
        self.region_to_id.len()
    }

    pub fn is_empty(&self) -> bool {
        self.region_to_id.len() == 0
    }

    pub fn contains_region(&self, region: &Region) -> bool {
        self.region_to_id.contains_key(region)
    }

    pub fn add_special_tokens(&mut self, special_tokens: &SpecialTokens) {
        let special_tokens_vec: Vec<Region> = special_tokens.into();
        for token in special_tokens_vec.iter() {
            self.add_token_to_universe(token);
        }
    }
}

impl From<Vec<Region>> for Universe {
    fn from(value: Vec<Region>) -> Self {
        // make a copy of the regions
        let regions = value;

        // create the region to id map and add the Unk token if it doesn't exist
        let region_to_id = generate_region_to_id_map(&regions);
        let id_to_region = generate_id_to_region_map(&regions);

        Universe {
            regions,
            region_to_id,
            id_to_region,
        }
    }
}

impl TryFrom<&Path> for Universe {
    type Error = std::io::Error;

    fn try_from(value: &Path) -> Result<Self, Self::Error> {
        let regions = RegionSet::try_from(value)
            .map_err(|err| std::io::Error::new(std::io::ErrorKind::Other, err))?
            .regions;

        let region_to_id = generate_region_to_id_map(&regions);
        let id_to_region = generate_id_to_region_map(&regions);

        Ok(Universe {
            regions,
            region_to_id,
            id_to_region,
        })
    }
}

pub mod utils;

use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use thiserror::Error;

use utils::UniverseFileType;

use crate::common::models::region::Region;
use crate::common::utils::{generate_id_to_region_map, generate_region_to_id_map, get_dynamic_reader};

use super::utils::special_tokens::SpecialTokens;

#[derive(Debug, Error)]
pub enum UniverseError {
    #[error("Could not determine the universe type from the file")]
    UnknownUniverseType,
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error("Error parsing line: {0}")]
    ParsingError(String),
}

#[derive(Clone, Eq, PartialEq, Default)]
pub struct Universe {
    pub regions: Vec<Region>,
    pub region_to_id: HashMap<Region, u32>,
    pub id_to_region: HashMap<u32, Region>,
    pub ordered: bool
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

    pub fn is_ordered(&self) -> bool {
        self.ordered
    }

    pub fn with_ordered(mut self, ordered: bool) -> Self {
        self.ordered = ordered;
        self
    }

}

impl TryFrom<&Path> for Universe {
    type Error = UniverseError;

    fn try_from(value: &Path) -> Result<Self, Self::Error> {
        
        let reader = get_dynamic_reader(value).map_err(|err| 
            UniverseError::Io(err.downcast().unwrap())
        )?;
        
        let mut lines = reader.lines().peekable();
        let first_line = lines.peek()
            .ok_or(UniverseError::UnknownUniverseType)?
            .as_ref()
            .map_err(|_| UniverseError::UnknownUniverseType)?;

        let (ordered, regions) = match UniverseFileType::from(first_line) {
            UniverseFileType::BedThree => {
                let mut regions = vec![];
                for line in lines {
                    let line = line?;
                    let parts: Vec<&str> = line.split('\t').collect();

                    if parts.len() == 3 {
                        let region = Region {
                            chr: parts[0].to_string(),
                            start: parts[1].parse().map_err(|_| UniverseError::ParsingError(line.clone()))?,
                            end: parts[2].parse().map_err(|_| UniverseError::ParsingError(line.clone()))?,
                            rest: None
                        };
                        regions.push(region);
                    } else {
                        return Err(UniverseError::ParsingError(line));
                    }
                }

                (false, regions)
            }
            // contains a score, so order the regions
            UniverseFileType::BedFivePlus => {

                let mut regions = vec![];
                for line in lines {
                    let line = line?;
                    let parts: Vec<&str> = line.split('\t').collect();

                    if parts.len() >= 5 {
                        let region = Region {
                            chr: parts[0].to_string(),
                            start: parts[1].parse().map_err(|_| UniverseError::ParsingError(line.clone()))?,
                            end: parts[2].parse().map_err(|_| UniverseError::ParsingError(line.clone()))?,
                            rest: Some(parts[4..].join("\t")),
                        };
                        regions.push(region);
                    } else {
                        return Err(UniverseError::ParsingError(line));
                    }
                }
                
                // sort by the score in the 5th column
                // sort high to low, this enables us to then order the regions
                // post processing by score
                regions.sort_by(|a, b| {
                    let score_a: f32 = a.rest.as_ref().unwrap().split_whitespace().next().unwrap().parse().unwrap();
                    let score_b: f32 = b.rest.as_ref().unwrap().split_whitespace().next().unwrap().parse().unwrap();

                    score_b.partial_cmp(&score_a).unwrap()
                });

                (true, regions)
            }

            UniverseFileType::Unknown => {
                return Err(UniverseError::UnknownUniverseType);
            }
        };

        let region_to_id = generate_region_to_id_map(&regions);
        let id_to_region = generate_id_to_region_map(&regions);

        Ok(Universe {
            regions,
            region_to_id,
            id_to_region,
            ordered
        })
    }
}

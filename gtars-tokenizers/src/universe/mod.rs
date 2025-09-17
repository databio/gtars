//!
//! This module defines the `Universe` struct, which represents a collection of regions.
//! It also provides methods to convert between region strings and their corresponding IDs.
//!
//! The universe is at the core of all tokenizers. It can be thought of as the oracle of the tokenizer, providing the necessary
//! information to convert between region strings and their corresponding IDs.
//!
pub mod utils;

use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use thiserror::Error;

use utils::UniverseFileType;

use gtars_core::utils::{
    generate_id_to_region_string_map, generate_region_string_to_id_map, get_dynamic_reader,
};

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

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Universe {
    pub regions: Vec<String>,
    pub region_to_id: HashMap<String, u32>,
    pub id_to_region: HashMap<u32, String>,
    pub names: Option<HashMap<String, String>>,
    pub scores: Option<HashMap<String, f64>>,
    pub special_tokens: Option<Vec<String>>,
}

impl Universe {
    ///
    /// Add a new region to the universe. Regions must be in the format `chr:start-end`.
    ///
    /// # Arguments:
    /// - `region`: the region string to add
    ///
    pub fn add_token_to_universe(&mut self, region: &str) {
        let new_id = self.region_to_id.len();
        self.region_to_id.insert(region.to_owned(), new_id as u32);
        self.id_to_region.insert(new_id as u32, region.to_owned());
        self.regions.push(region.to_owned());
    }

    ///
    /// Convert a region string to its corresponding ID.
    /// # Arguments:
    /// - `region`: the region string to convert
    /// # Returns:
    /// - `Option<u32>`: the ID corresponding to the region string, or None if it doesn't exist
    ///
    pub fn convert_token_to_id(&self, region: &str) -> Option<u32> {
        let id = self.region_to_id.get(region);
        id.map(|id| id.to_owned())
    }

    ///
    /// Convert an ID to its corresponding region string.
    ///
    /// # Arguments:
    /// - `id`: the ID to convert
    /// # Returns:
    /// - `Option<String>`: the region string corresponding to the ID, or None if it doesn't exist
    ///
    pub fn convert_id_to_token(&self, id: u32) -> Option<String> {
        self.id_to_region.get(&id).cloned()
    }

    ///
    /// Get the number of regions in the universe.
    ///
    /// # Returns:
    /// - `usize`: the number of regions in the universe
    ///
    pub fn len(&self) -> usize {
        self.region_to_id.len()
    }

    ///
    /// Check if the universe is empty.
    ///
    pub fn is_empty(&self) -> bool {
        self.region_to_id.len() == 0
    }

    ///
    /// Check if a region exists in the universe.
    /// # Arguments:
    /// - `region`: the region string to check
    /// # Returns:
    /// - `bool`: true if the region exists, false otherwise
    ///     
    pub fn contains_region(&self, region: &String) -> bool {
        self.region_to_id.contains_key(region)
    }

    ///
    /// Add special tokens to the universe. This will also add the special tokens
    /// to the universe as new regions.
    ///
    pub fn add_special_tokens(&mut self, special_tokens: &SpecialTokens) {
        let special_tokens_vec: Vec<String> = special_tokens.into();
        self.special_tokens = Some(special_tokens_vec.clone());
        for token in special_tokens_vec.iter() {
            self.add_token_to_universe(token);
        }
    }
}

impl TryFrom<&Path> for Universe {
    type Error = UniverseError;

    fn try_from(value: &Path) -> Result<Self, Self::Error> {
        let reader =
            get_dynamic_reader(value).map_err(|err| UniverseError::Io(err.downcast().unwrap()))?;

        let mut lines = reader.lines().peekable();
        let first_line = lines
            .peek()
            .ok_or(UniverseError::UnknownUniverseType)?
            .as_ref()
            .map_err(|_| UniverseError::UnknownUniverseType)?;

        let (regions, names, scores) = match UniverseFileType::from(first_line) {
            UniverseFileType::BedThree => {
                let mut regions = vec![];
                for line in lines {
                    let line = line?;
                    let parts: Vec<&str> = line.split_whitespace().collect();

                    if parts.len() == 3 {
                        let region = format!("{}:{}-{}", parts[0], parts[1], parts[2]);
                        regions.push(region);
                    } else {
                        return Err(UniverseError::ParsingError(line));
                    }
                }

                (regions, None, None)
            }
            // contains a score, so order the regions
            UniverseFileType::BedFivePlus => {
                let mut regions = vec![];
                let mut names = HashMap::new();
                let mut scores = HashMap::new();
                for line in lines {
                    let line = line?;
                    let parts: Vec<&str> = line.split('\t').collect();

                    if parts.len() >= 5 {
                        let region = format!("{}:{}-{}", parts[0], parts[1], parts[2]);
                        regions.push(region.clone());

                        let name = parts[3].to_string();
                        let score: f64 = parts[4].trim().parse().unwrap();

                        names.insert(region.clone(), name);
                        scores.insert(region, score);
                    } else {
                        return Err(UniverseError::ParsingError(line));
                    }
                }

                (regions, Some(names), Some(scores))
            }

            UniverseFileType::Unknown => {
                return Err(UniverseError::UnknownUniverseType);
            }
        };

        let region_to_id = generate_region_string_to_id_map(&regions);
        let id_to_region = generate_id_to_region_string_map(&regions);

        Ok(Universe {
            regions,
            region_to_id,
            id_to_region,
            names,
            scores,
            special_tokens: None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[fixture]
    fn vanilla_peaks_path() -> String {
        "../tests/data/tokenizers/peaks.bed".to_string()
    }

    #[fixture]
    fn ordered_peaks_path() -> String {
        "../tests/data/tokenizers/peaks.scored.bed".to_string()
    }

    #[rstest]
    fn test_create_vanilla_universe(vanilla_peaks_path: String) {
        let universe = Universe::try_from(Path::new(&vanilla_peaks_path)).unwrap();
        assert_eq!(universe.len(), 25); // there are 25 regions in the file
        assert_eq!(universe.scores.is_none(), true);
    }

    #[rstest]
    fn test_create_ordered_universe(ordered_peaks_path: String) {
        let universe = Universe::try_from(Path::new(&ordered_peaks_path)).unwrap();
        assert_eq!(universe.len(), 25); // there are 25 regions in the file
        assert_eq!(universe.names.is_some(), true);
        assert_eq!(universe.scores.is_some(), true);
    }
}

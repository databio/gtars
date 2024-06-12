use serde::{Deserialize, Serialize};

#[derive(Deserialize, Serialize, Debug, PartialEq)]
pub struct TokenizerConfig {
    pub universe: String,
    pub hierarchical_universes: Option<Vec<String>>,
    pub exclude_ranges: Option<String>,
}

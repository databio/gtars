use std::fs::read_to_string;
use std::path::Path;

use anyhow::Result;
use serde::{Deserialize, Serialize};

#[derive(Deserialize, Serialize, Debug, PartialEq)]
pub struct TokenizerConfig {
    pub tokenizer_type: Option<String>,
    pub universes: Vec<String>,
    pub exclude_ranges: Option<String>,
}

impl TokenizerConfig {
    ///
    /// Create a new tokenizer config.
    ///
    /// # Arguments
    /// - path: Path to the config file (a .toml) file.
    pub fn try_from(path: &Path) -> Result<TokenizerConfig> {
        let toml_str = read_to_string(path)?;
        let config: TokenizerConfig = toml::from_str(&toml_str)?;

        Ok(config)
    }

    pub fn new(
        tokenizer_type: Option<String>,
        universes: Vec<String>,
        exclude_ranges: Option<String>,
    ) -> TokenizerConfig {
        TokenizerConfig {
            tokenizer_type,
            universes,
            exclude_ranges,
        }
    }
}

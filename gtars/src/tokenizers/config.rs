use serde::{Deserialize, Serialize};

#[derive(Deserialize, Serialize, Debug, PartialEq )]
pub struct TokenizerConfig {
    pub universe: String,
    pub excluderanges: Option<String>,
}
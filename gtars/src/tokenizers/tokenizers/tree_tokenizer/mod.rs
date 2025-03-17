use std::collections::HashMap;

use rust_lapper::{Interval, Lapper};

use crate::tokenizers::config::{self, TokenizerConfig};

use super::SpecialTokens;

pub struct TreeTokenizer {
    /// The core interval tree. Actually, its **many** interval trees. The hash-map will map chrom names
    /// to an interval tree for querying. The hash-map lookup should be constant time (O(1)), while
    /// the interval tree is [reported to be NlogN](https://academic.oup.com/bioinformatics/article/29/1/1/273289?login=false)
    tree: HashMap<String, Lapper<u32, u32>>,
    special_tokens: SpecialTokens
}

impl From<TokenizerConfig> for TreeTokenizer {
    fn from(value: TokenizerConfig) -> Self {
        // we start with the default, then will replace as they exist in the config
        let special_tokens = SpecialTokens::default();
        if let Some(config_special_tokens) = value.special_tokens {
            // if let Some(config_eos) = config_special_tokens.eos
        }
    }
}
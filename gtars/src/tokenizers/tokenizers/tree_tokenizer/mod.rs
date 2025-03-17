use std::collections::HashMap;
use std::path::Path;

use rust_lapper::Lapper;

use super::SpecialTokens;
use super::Tokenizer;
use super::TokenizerError;
use crate::common::models::Region;
use crate::tokenizers::config::{SpecialToken, TokenizerConfig};
use crate::tokenizers::tokens::TokenizedRegionSet;
use crate::tokenizers::universe::Universe;
use crate::tokenizers::utils::create_interval_tree_from_universe;

pub struct TreeTokenizer {
    /// The core interval tree. Actually, its **many** interval trees. The hash-map will map chrom names
    /// to an interval tree for querying. The hash-map lookup should be constant time (O(1)), while
    /// the interval tree is [reported to be NlogN](https://academic.oup.com/bioinformatics/article/29/1/1/273289?login=false)
    tree: HashMap<String, Lapper<u32, u32>>,
    universe: Universe,
    special_tokens: SpecialTokens,
}

impl TryFrom<TokenizerConfig> for TreeTokenizer {
    type Error = TokenizerError;

    fn try_from(value: TokenizerConfig) -> Result<Self, Self::Error> {
        let universe = Path::new(&value.universe);
        let mut universe = Universe::try_from(universe)?;
        let tree = create_interval_tree_from_universe(&universe);

        // we start with the default, then will replace as they exist in the config
        let mut special_tokens = SpecialTokens::default();
        if let Some(config_special_tokens) = value.special_tokens {
            for token in config_special_tokens {
                let token_type = token.name;
                let token_value = token.token;

                let parts = token_value.split(':').collect::<Vec<&str>>();
                if parts.len() != 2 {
                    return Err(TokenizerError::InvalidSpecialTokenConfig);
                }

                let chr = parts[0].to_string();
                let start_end = parts[1].split('-').collect::<Vec<&str>>();
                if start_end.len() != 2 {
                    return Err(TokenizerError::InvalidSpecialTokenConfig);
                }
                let start = start_end[0]
                    .parse::<u32>()
                    .map_err(|_| TokenizerError::InvalidSpecialTokenConfig)?;
                let end = start_end[1]
                    .parse::<u32>()
                    .map_err(|_| TokenizerError::InvalidSpecialTokenConfig)?;

                let token_value = Region {
                    chr,
                    start,
                    end,
                    rest: None,
                };

                match token_type {
                    SpecialToken::Unk => special_tokens.unk = token_value,
                    SpecialToken::Pad => special_tokens.pad = token_value,
                    SpecialToken::Mask => special_tokens.mask = token_value,
                    SpecialToken::Cls => special_tokens.cls = token_value,
                    SpecialToken::Sep => special_tokens.sep = token_value,
                    SpecialToken::Eos => special_tokens.eos = token_value,
                    SpecialToken::Bos => special_tokens.bos = token_value,
                }
            }
        }

        // insert all new special tokens into the universe
        let s_tokens: Vec<Region> = special_tokens.clone().into();
        for s_token in s_tokens.iter() {
            universe.insert_token(s_token);
        }

        Ok(TreeTokenizer {
            tree,
            universe,
            special_tokens,
        })
    }
}

impl Tokenizer for TreeTokenizer {
    fn tokenize<T: Into<Vec<Region>>>(
        &self,
        regions: T,
    ) -> Result<TokenizedRegionSet, TokenizerError> {
    }

    fn token_to_id(&self, token: &Region) -> Option<u32> {
        self.universe.convert_region_to_id(token)
    }

    fn id_to_token(&self, id: u32) -> Option<Region> {
        self.universe.convert_id_to_region(id)
    }

    fn get_vocab_size(&self) -> usize {
        self.universe.len()
    }

    fn save(
        &self,
        folder: &Path,
        prefix: Option<&str>,
    ) -> Result<Vec<std::path::PathBuf>, TokenizerError> {
        todo!()
    }
}

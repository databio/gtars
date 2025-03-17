pub mod bits_tree;

use std::path::Path;

use bits_tree::BitsTree;
use thiserror::Error;

use crate::common::models::Region;

use super::config::{TokenizerConfig, TokenizerConfigError, TokenizerType};
use super::tokens::TokenizedRegionSet;
use super::universe::Universe;
use super::utils::prepare_universe_and_special_tokens;
use super::utils::special_tokens::SpecialTokens;

#[derive(Error, Debug)]
pub enum TokenizerError {
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error("Invalid special token configuration")]
    InvalidSpecialTokenConfig,
    #[error(transparent)]
    Config(#[from] TokenizerConfigError),
}

pub trait GTokenize {
    /// Tokenize the given sequence into multiple underlying `Token`. The `offsets` on the `Token`
    /// are expected to be relative to the given sequence.
    fn tokenize(
        &self,
        regions: &[Region],
    ) -> Result<TokenizedRegionSet, TokenizerError>;
    /// Find the ID associated to a string token
    fn token_to_id(&self, token: &Region) -> Option<u32>;
    /// Find the string token associated to an ID
    fn id_to_token(&self, id: u32) -> Option<Region>;
    /// Retrieve the size of the vocabulary
    fn get_vocab_size(&self) -> usize;
    /// Retrieve the universe -- this is here to
    /// enforce that the tokenizer has a universe
    fn get_universe(&self) -> &Universe;
}

pub struct Tokenizer {
    core: Box<dyn GTokenize>,
    special_tokens: SpecialTokens
}

impl Tokenizer {
    ///
    /// Create a new tokenizer with the given core GTokenizer
    /// 
    pub fn new(core: Box<dyn GTokenize>) -> Self {
        let special_tokens = SpecialTokens::default();
        Tokenizer { core, special_tokens }
    }

    ///
    /// Add special tokens to the tokenizer
    /// 
    pub fn with_special_tokens(self, special_tokens: SpecialTokens) -> Self {
        Tokenizer { core: self.core, special_tokens }
    }

    ///
    /// Create a new tokenizer from a config file
    /// 
    pub fn from_config<P: AsRef<Path>>(cfg_path: P) -> Result<Self, TokenizerError> {
        let config = TokenizerConfig::try_from(cfg_path.as_ref())?;

        let config_path = cfg_path.as_ref().parent().unwrap();
        let universe_path = config_path.join(&config.universe);
        let special_tokens = match config.special_tokens {
            Some(tokens) => SpecialTokens::from(tokens),
            None => SpecialTokens::default(),
        };

        let (universe, special_tokens) = prepare_universe_and_special_tokens(universe_path, special_tokens)?;

        let core = match config.tokenizer_type {
            Some(TokenizerType::BitsTree) => {
                Box::new(BitsTree::from(universe))
            },
            // default to BitsTree if no tokenizer type is specified (yes, this means we only support BitsTree for now)
            None => {
                Box::new(BitsTree::from(universe))
            },
        };

        Ok(Tokenizer {
            core,
            special_tokens
        })
    }

    pub fn tokenize(&self, regions: &[Region]) -> Result<TokenizedRegionSet, TokenizerError> {
        let mut res = self.core.tokenize(regions)?;
        if res.is_empty() {
            res.ids.push(self.core.token_to_id(&self.special_tokens.unk).unwrap());
        }
        Ok(res)
    }

    pub fn tokenize_batch(&self, regions: &[&[Region]]) -> Result<Vec<TokenizedRegionSet>, TokenizerError> {
        let mut results = Vec::new();
        for region_set in regions {
            let mut res = self.core.tokenize(region_set)?;
            if res.is_empty() {
                res.ids.push(self.core.token_to_id(&self.special_tokens.unk).unwrap());
            }
            results.push(res);
        }
        Ok(results)
    }

}
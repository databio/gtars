pub mod tree_tokenizer;

use thiserror::Error;

use crate::common::models::Region;

use super::tokens::TokenizedRegionSet;
use super::universe::Universe;

#[derive(Error, Debug)]
pub enum TokenizerError {
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error("Invalid special token configuration")]
    InvalidSpecialTokenConfig,
}

pub trait GTokenize {
    /// Tokenize the given sequence into multiple underlying `Token`. The `offsets` on the `Token`
    /// are expected to be relative to the given sequence.
    fn tokenize<T: Into<Vec<Region>>>(
        &self,
        regions: T,
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

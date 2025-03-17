pub mod tree_tokenizer;

use thiserror::Error;

use crate::common::models::Region;
use crate::tokenizers::utils::padding::PaddingParams;
use crate::tokenizers::utils::truncation::TruncationParams;

use super::tokens::TokenizedRegionSet;

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
}

pub struct Tokenizer<T: GTokenize> {
    core: T,
    padding: Option<PaddingParams>,
    truncation: Option<TruncationParams>,
}

impl<T> Tokenizer<T>
where
    T: GTokenize,
{
    pub fn new(tokenizer: T) -> Self {
        Tokenizer {
            core: tokenizer,
            padding: None,
            truncation: None,
        }
    }

    pub fn with_padding(mut self, padding: PaddingParams) -> Self {
        self.padding = Some(padding);
        self
    }

    pub fn get_padding_params(&self) -> Option<&PaddingParams> {
        self.padding.as_ref()
    }

    pub fn with_truncation(mut self, truncation: TruncationParams) -> Self {
        self.truncation = Some(truncation);
        self
    }

    pub fn get_truncation_params(&self) -> Option<&TruncationParams> {
        self.truncation.as_ref()
    }
}
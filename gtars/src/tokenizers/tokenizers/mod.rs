pub mod tree_tokenizer;

use std::{
    default,
    path::{Path, PathBuf},
};

use thiserror::Error;

use crate::common::models::Region;

use super::tokens::TokenizedRegionSet;

#[derive(Error, Debug)]
pub enum TokenizerError {
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error("Invalid special token configuration")]
    InvalidSpecialTokenConfig,
}

pub trait Tokenizer {
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

#[derive(Clone, Debug)]
pub struct SpecialTokens {
    pub unk: Region,
    pub pad: Region,
    pub mask: Region,
    pub cls: Region,
    pub eos: Region,
    pub bos: Region,
    pub sep: Region,
}

impl Default for SpecialTokens {
    fn default() -> Self {
        SpecialTokens {
            unk: Region {
                chr: "chrUNK".to_string(),
                start: 0,
                end: 0,
                rest: None,
            },
            pad: Region {
                chr: "chrPAD".to_string(),
                start: 0,
                end: 0,
                rest: None,
            },
            mask: Region {
                chr: "chrMASK".to_string(),
                start: 0,
                end: 0,
                rest: None,
            },
            cls: Region {
                chr: "chrCLS".to_string(),
                start: 0,
                end: 0,
                rest: None,
            },
            eos: Region {
                chr: "chrEOS".to_string(),
                start: 0,
                end: 0,
                rest: None,
            },
            bos: Region {
                chr: "chrBOS".to_string(),
                start: 0,
                end: 0,
                rest: None,
            },
            sep: Region {
                chr: "chrSEP".to_string(),
                start: 0,
                end: 0,
                rest: None,
            },
        }
    }
}

impl From<SpecialTokens> for Vec<Region> {
    fn from(val: SpecialTokens) -> Self {
        vec![
            val.unk, val.pad, val.mask, val.cls, val.eos, val.bos, val.sep,
        ]
    }
}

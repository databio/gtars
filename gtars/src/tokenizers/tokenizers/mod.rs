pub mod tree_tokenizer;

use std::{default, path::{Path, PathBuf}};

use thiserror::Error;

use crate::common::models::Region;

use super::tokens::TokenizedRegionSet;

#[derive(Error, Debug)]
pub enum TokenizerError {}

pub trait Tokenizer {
    /// Tokenize the given sequence into multiple underlying `Token`. The `offsets` on the `Token`
    /// are expected to be relative to the given sequence.
    fn tokenize(&self, sequence: &str) -> Result<TokenizedRegionSet, TokenizerError>;
    /// Find the ID associated to a string token
    fn token_to_id(&self, token: &str) -> Option<u32>;
    /// Find the string token associated to an ID
    fn id_to_token(&self, id: u32) -> Option<String>;
    /// Retrieve the size of the vocabulary
    fn get_vocab_size(&self) -> usize;
    /// Save the current `Model` in the given folder, using the given `prefix` for the various
    /// files that need to be saved.
    fn save(&self, folder: &Path, prefix: Option<&str>) -> Result<Vec<PathBuf>, TokenizerError>;
}

pub struct SpecialTokens {
    pub unk: Region,
    pub pad: Region,
    pub mask: Region,
    pub cls: Region,
    pub eos: Region,
    pub bos: Region,
    pub sep: Region
}

impl Default for SpecialTokens {
    fn default() -> Self {
        SpecialTokens {
            unk: Region {
                chr: "chrUNK".to_string(),
                start: 0,
                end: 0,
                rest: None
            },
            pad: Region {
                chr: "chrPAD".to_string(),
                start: 0,
                end: 0,
                rest: None
            },
            mask: Region {
                chr: "chrMASK".to_string(),
                start: 0,
                end: 0,
                rest: None
            },
            cls: Region {
                chr: "chrCLS".to_string(),
                start: 0,
                end: 0,
                rest: None
            },
            eos: Region {
                chr: "chrEOS".to_string(),
                start: 0,
                end: 0,
                rest: None
            },
            bos: Region {
                chr: "chrBOS".to_string(),
                start: 0,
                end: 0,
                rest: None
            },
            sep: Region {
                chr: "chrSEP".to_string(),
                start: 0,
                end: 0,
                rest: None
            }
        }
    }
}
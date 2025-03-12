use std::fs::read_to_string;
use std::path::Path;

use thiserror::Error;

use serde::{Deserialize, Serialize};
use std::ffi::OsStr;

#[derive(Deserialize, Serialize, Debug, PartialEq)]
pub struct SpecialTokensConfig {
    pub unk: Option<String>,
    pub pad: Option<String>,
    pub bos: Option<String>,
    pub eos: Option<String>,
    pub cls: Option<String>,
    pub sep: Option<String>,
    pub mask: Option<String>,
}

#[derive(Deserialize, Serialize, Debug, PartialEq)]
pub struct TokenizerConfig {
    pub universe: String,
    pub special_tokens: Option<SpecialTokensConfig>,
}

#[derive(Debug)]
pub enum TokenizerInputFileType {
    Toml,
    Bed,
    BedGz,
}

#[derive(Error, Debug)]
pub enum TokenizerConfigError {
    #[error("Missing or invalid file extension in tokenizer config file. It must be `toml`, `bed` or `bed.gz`")]
    InvalidFileType,
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error(transparent)]
    Toml(#[from] toml::de::Error),
}

pub type TokenizerConfigResult<T> = std::result::Result<T, TokenizerConfigError>;

impl TokenizerInputFileType {
    pub fn from_path(path: &Path) -> TokenizerConfigResult<Self> {
        // first extension, might be ".gz"
        let ext1 = path
            .extension()
            .and_then(OsStr::to_str)
            .ok_or(TokenizerConfigError::InvalidFileType)?;

        if ext1 == "gz" {
            // if it ends with .gz, we look at the file stem’s extension:
            // e.g. "universe.bed.gz" → first strip .gz → "universe.bed"
            // so now we check if *that* ends with ".bed"
            let file_stem = path
                .file_stem()
                .ok_or(TokenizerConfigError::InvalidFileType)?;
            let ext2 = Path::new(file_stem)
                .extension()
                .and_then(OsStr::to_str)
                .ok_or(TokenizerConfigError::InvalidFileType)?;

            if ext2 == "bed" {
                Ok(TokenizerInputFileType::BedGz)
            } else {
                Err(TokenizerConfigError::InvalidFileType)
            }
        } else if ext1 == "toml" {
            Ok(TokenizerInputFileType::Toml)
        } else if ext1 == "bed" {
            Ok(TokenizerInputFileType::Bed)
        } else {
            Err(TokenizerConfigError::InvalidFileType)
        }
    }
}


impl TokenizerConfig {
    ///
    /// Create a new tokenizer config.
    ///
    /// # Arguments
    /// - path: Path to the config file (a .toml) file.
    pub fn try_from(path: &Path) -> TokenizerConfigResult<TokenizerConfig> {
        let toml_str = read_to_string(path)?;
        let config: TokenizerConfig = toml::from_str(&toml_str)?;

        Ok(config)
    }

    pub fn new(universe: String, special_tokens: Option<SpecialTokensConfig>) -> TokenizerConfig {
        TokenizerConfig {
            universe,
            special_tokens,
        }
    }
}

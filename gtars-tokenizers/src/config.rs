use std::fs::read_to_string;
use std::path::Path;

use thiserror::Error;

use serde::{Deserialize, Serialize};
use std::ffi::OsStr;

#[derive(Deserialize, Serialize, Debug, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum SpecialToken {
    Unk,
    Pad,
    Mask,
    Cls,
    Bos,
    Eos,
    Sep,
}

#[derive(Deserialize, Serialize, Debug, PartialEq)]
pub struct SpecialTokenAssignment {
    pub name: SpecialToken,
    pub token: String, // must be valid chr:start-end
}

#[derive(Deserialize, Serialize, Debug, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum TokenizerType {
    #[serde(rename = "bits")]
    Bits,
    #[serde(rename = "ailist")]
    AIList,
}

#[derive(Deserialize, Serialize, Debug, PartialEq)]
pub struct TokenizerConfig {
    pub universe: String,
    pub special_tokens: Option<Vec<SpecialTokenAssignment>>,
    pub tokenizer_type: Option<TokenizerType>,
}

#[derive(Debug)]
pub enum TokenizerInputFileType {
    Toml,
    Bed,
    BedGz,
}

#[derive(Error, Debug)]
pub enum TokenizerConfigError {
    #[error(
        "Missing or invalid file extension in tokenizer config file. It must be `toml`, `bed` or `bed.gz`"
    )]
    InvalidFileType,
    #[error("Invalid tokenizer type in config file")]
    InvalidTokenizerType,
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error(transparent)]
    Toml(#[from] toml::de::Error),
}

pub type TokenizerConfigResult<T> = std::result::Result<T, TokenizerConfigError>;

impl TokenizerInputFileType {
    ///
    /// Determine the type of the tokenizer input file based on its extension.
    /// # Arguments
    /// * `path` - A reference to a `Path` object representing the file path.
    /// # Returns
    /// * `TokenizerInputFileType` - An enum representing the type of the tokenizer input file.
    ///
    pub fn from_path(path: &Path) -> TokenizerConfigResult<Self> {
        match path.extension().and_then(OsStr::to_str) {
            Some("gz") => {
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
            }
            Some("toml") => Ok(TokenizerInputFileType::Toml),
            Some("bed") => Ok(TokenizerInputFileType::Bed),
            _ => Err(TokenizerConfigError::InvalidFileType),
        }
    }
}

impl TryFrom<&Path> for TokenizerConfig {
    type Error = TokenizerConfigError;

    fn try_from(path: &Path) -> Result<Self, Self::Error> {
        let toml_str = read_to_string(path)?;
        let config = toml::from_str(&toml_str)?;
        Ok(config)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use pretty_assertions::assert_eq;
    use rstest::rstest;

    use std::path::PathBuf;

    #[rstest]
    fn test_try_from_toml() {
        let path = PathBuf::from("../tests/data/tokenizers/tokenizer.toml");
        let result = TokenizerConfig::try_from(path.as_path());
        assert_eq!(result.is_ok(), true);
    }

    #[rstest]
    fn test_from_path_for_toml_extension() {
        let path = PathBuf::from("dummy.toml");
        let file_type = TokenizerInputFileType::from_path(path.as_path());
        assert_eq!(matches!(file_type, Ok(TokenizerInputFileType::Toml)), true);
    }

    #[rstest]
    fn test_from_path_for_invalid_extension() {
        let path = PathBuf::from("invalid.xyz");
        let file_type = TokenizerInputFileType::from_path(&path);
        assert_eq!(file_type.is_err(), true);
    }

    #[rstest]
    fn test_get_universe_name() {
        let path = PathBuf::from("../tests/data/tokenizers/tokenizer.toml");
        let result = TokenizerConfig::try_from(path.as_path()).unwrap();

        assert_eq!(result.universe, "peaks.bed.gz");
    }

    #[rstest]
    fn test_get_special_tokens() {
        let path = PathBuf::from("../tests/data/tokenizers/tokenizer_custom_specials.toml");
        let result = TokenizerConfig::try_from(path.as_path()).unwrap();
        let special_tokens = result.special_tokens;

        assert_eq!(special_tokens.is_some(), true);
    }
}

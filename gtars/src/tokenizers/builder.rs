use std::path::Path;

use anyhow::{Error, Result};

use super::{
    MetaTokenizer,
    Tokenizer,
    // FragmentTokenizer,
    TokenizerConfig,
    TreeTokenizer,
};

pub struct TokenizerBuilder;

impl TokenizerBuilder {
    pub fn from_toml(path: &Path) -> Result<Box<dyn Tokenizer>> {
        let config = TokenizerConfig::try_from(path)?;
        if let Some(tokenizer_type) = config.tokenizer_type {
            match tokenizer_type.as_str() {
                "tree" => Ok(Box::new(TreeTokenizer::try_from(path)?)),
                "meta" => Ok(Box::new(MetaTokenizer::try_from(path)?)),
                _ => Err(Error::msg("Tokenizer type not supported")),
            }
        } else {
            println!("No tokenizer type found in config file. Instantiating a default TreeTokenizer. Note that this may lead to unexpected behavior.");
            Ok(Box::new(TreeTokenizer::try_from(path)?))
        }
    }
}

#[cfg(test)]
mod tests {

    use std::path::Path;

    use super::*;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[fixture]
    fn path_to_bed_file() -> &'static str {
        "tests/data/peaks.bed"
    }

    #[fixture]
    fn path_to_config_file() -> &'static str {
        "tests/data/tokenizer.toml"
    }

    #[rstest]
    fn test_from_toml(path_to_config_file: &str) {
        let path = Path::new(path_to_config_file);
        let tokenizer = TokenizerBuilder::from_toml(path).unwrap();
        assert_eq!(tokenizer.vocab_size(), 56);
    }
}

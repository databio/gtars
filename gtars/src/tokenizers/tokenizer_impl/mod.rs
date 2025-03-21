pub mod bits_tree;

use std::collections::HashMap;
use std::path::Path;

use bits_tree::BitsTree;
use thiserror::Error;

use crate::common::models::Region;

use super::config::{TokenizerConfig, TokenizerConfigError, TokenizerInputFileType, TokenizerType};
use super::encoding::{Encoding, BatchEncoding};
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
    #[error("Universe error: {0}")]
    UniverseError(#[from] crate::tokenizers::universe::UniverseError),
    #[error(transparent)]
    Anyhow(#[from] anyhow::Error),
}


#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Token {
    pub id: u32,
    pub value: String,
}

impl Token {
    pub fn new(id: u32, value: String) -> Self {
        Self { id, value }
    }
}

pub trait GTokenize: Send + Sync {
    /// Tokenize the given sequence into multiple underlying `Token`. The `offsets` on the `Token`
    /// are expected to be relative to the given sequence.
    fn tokenize(&self, regions: &[Region]) -> Result<Vec<Token>, TokenizerError>;
    /// Find the ID associated to a string token
    fn token_to_id(&self, token: &str) -> Option<u32>;
    /// Find the string token associated to an ID
    fn id_to_token(&self, id: u32) -> Option<String>;
    /// Retrieve the size of the vocabulary
    fn get_vocab_size(&self) -> usize;
    /// Retrieve the universe -- this is here to
    /// enforce that the tokenizer has a universe
    fn get_universe(&self) -> &Universe;
    /// Retrieve the entire vocabulary mapping (token -> id)
    fn get_vocab(&self) -> HashMap<String, u32>;
}

pub struct Tokenizer {
    core: Box<dyn GTokenize + Send + Sync>,
    special_tokens: SpecialTokens,
}

impl Tokenizer {
    ///
    /// Create a new tokenizer with the given core GTokenizer
    ///
    pub fn new(core: Box<dyn GTokenize>) -> Self {
        let special_tokens = SpecialTokens::default();
        Tokenizer {
            core,
            special_tokens,
        }
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

        let (universe, special_tokens) =
            prepare_universe_and_special_tokens(universe_path, special_tokens)?;

        let core = match config.tokenizer_type {
            Some(TokenizerType::BitsTree) => Box::new(BitsTree::from(universe)),
            // default to BitsTree if no tokenizer type is specified (yes, this means we only support BitsTree for now)
            None => Box::new(BitsTree::from(universe)),
        };

        Ok(Tokenizer {
            core,
            special_tokens,
        })
    }

    ///
    /// Create a new tokenizer from a bed file
    ///
    pub fn from_bed<P: AsRef<Path>>(bed_path: P) -> Result<Self, TokenizerError> {
        let special_tokens = SpecialTokens::default();
        let (universe, special_tokens) =
            prepare_universe_and_special_tokens(bed_path.as_ref(), special_tokens)?;
        let core = Box::new(BitsTree::from(universe));

        Ok(Tokenizer {
            core,
            special_tokens,
        })
    }

    ///
    /// Create a new tokenizer from a file, automatically detecting the type
    ///
    pub fn from_auto<P: AsRef<Path>>(path: P) -> Result<Self, TokenizerError> {
        let file_type = TokenizerInputFileType::from_path(path.as_ref())?;
        match file_type {
            TokenizerInputFileType::Toml => Tokenizer::from_config(path),
            TokenizerInputFileType::Bed => Tokenizer::from_bed(path),
            TokenizerInputFileType::BedGz => Tokenizer::from_bed(path),
        }
    }

    pub fn tokenize(&self, regions: &[Region]) -> Result<Vec<String>, TokenizerError> {
        let tokenized = self.core.tokenize(regions)?;
        if tokenized.is_empty() {
            return Ok(vec![self.special_tokens.unk.clone()]);
        }
        Ok(tokenized.into_iter().map(|t| t.value).collect())
    }

    pub fn convert_token_to_id(&self, token: &str) -> Option<u32> {
        self.core.token_to_id(token)
    }

    pub fn convert_id_to_token(&self, id: u32) -> Option<String> {
        self.core.id_to_token(id)
    }

    pub fn get_vocab_size(&self) -> usize {
        self.core.get_vocab_size()
    }

    pub fn get_vocab(&self) -> HashMap<String, u32> {
        self.core.get_vocab()
    }

    pub fn get_unk_token(&self) -> String {
        self.special_tokens.unk.clone()
    }

    pub fn get_pad_token(&self) -> String {
        self.special_tokens.pad.clone()
    }

    pub fn get_mask_token(&self) -> String {
        self.special_tokens.mask.clone()
    }

    pub fn get_cls_token(&self) -> String {
        self.special_tokens.cls.clone()
    }

    pub fn get_eos_token(&self) -> String {
        self.special_tokens.eos.clone()
    }

    pub fn get_bos_token(&self) -> String {
        self.special_tokens.bos.clone()
    }

    pub fn get_sep_token(&self) -> String {
        self.special_tokens.sep.clone()
    }

    pub fn get_special_tokens(&self) -> &SpecialTokens {
        &self.special_tokens
    }

}

#[cfg(test)]
mod tokenizer_tests {
    use super::*;

    use pretty_assertions::assert_eq;
    use rstest::*;

    #[fixture]
    fn anndata_path() -> String {
        "tests/data/tokenizers/pbmc_hg38.h5ad".to_string()
    }

    #[rstest]
    fn test_tokenizer_creation_from_config() {
        let cfg_path = "tests/data/tokenizers/tokenizer.toml";
        let tokenizer =
            Tokenizer::from_config(cfg_path).expect("Failed to create tokenizer from config.");
        assert_eq!(tokenizer.get_vocab_size(), 32); // 25 regions + 7 special tokens
    }

    #[rstest]
    fn test_tokenizer_creation_from_bed() {
        let bed_path = "tests/data/tokenizers/peaks.bed";
        let tokenizer =
            Tokenizer::from_bed(bed_path).expect("Failed to create tokenizer from config.");
        assert_eq!(tokenizer.get_vocab_size(), 32); // 25 regions + 7 special tokens
    }

    #[rstest]
    fn test_tokenizer_creation_from_bed_gz() {
        let bed_path = "tests/data/tokenizers/peaks.bed.gz";
        let tokenizer =
            Tokenizer::from_bed(bed_path).expect("Failed to create tokenizer from config.");
        assert_eq!(tokenizer.get_vocab_size(), 32); // 25 regions + 7 special tokens
    }

    #[rstest]
    fn test_tokenizer_creation_auto_all() {
        let bed_path = "tests/data/tokenizers/peaks.bed";
        let tokenizer =
            Tokenizer::from_auto(bed_path).expect("Failed to create tokenizer from config.");
        assert_eq!(tokenizer.get_vocab_size(), 32); // 25 regions + 7 special tokens

        let cfg_path = "tests/data/tokenizers/tokenizer.toml";
        let tokenizer =
            Tokenizer::from_auto(cfg_path).expect("Failed to create tokenizer from config.");
        assert_eq!(tokenizer.get_vocab_size(), 32); // 25 regions + 7 special tokens

        let bed_path = "tests/data/tokenizers/peaks.bed.gz";
        let tokenizer =
            Tokenizer::from_auto(bed_path).expect("Failed to create tokenizer from config.");
        assert_eq!(tokenizer.get_vocab_size(), 32); // 25 regions + 7 special tokens
    }

    #[rstest]
    fn test_tokenizer_bad_tokenizer_type() {
        let cfg_path = "tests/data/tokenizers/tokenizer_bad_ttype.toml";
        let tokenizer = Tokenizer::from_config(cfg_path);
        assert_eq!(tokenizer.is_err(), true);
    }

    #[rstest]
    fn test_tokenizer_custom_special_tokens() {
        let cfg_path = "tests/data/tokenizers/tokenizer_custom_specials.toml";
        let tokenizer =
            Tokenizer::from_config(cfg_path).expect("Failed to create tokenizer from config.");

        assert_eq!(tokenizer.get_vocab_size(), 32); // 25 regions + 7 special tokens

        // check that unk was overridden
        assert_eq!(tokenizer.get_unk_token(), "<unk>");

        // check that pad didnt change
        assert_eq!(tokenizer.get_pad_token(), "<pad>");
    }

    #[rstest]
    fn test_tokenize_single_region_not_overlapping() {
        let cfg_path = "tests/data/tokenizers/tokenizer.toml";
        let tokenizer =
            Tokenizer::from_config(cfg_path).expect("Failed to create tokenizer from config.");

        let regions = vec![Region {
            chr: "chr1".to_string(),
            start: 50,
            end: 150,
            rest: None,
        }];

        let tokenized = tokenizer.tokenize(&regions);
        assert!(tokenized.is_ok());
        let tokenized = tokenized.unwrap();
        assert_eq!(tokenized.len(), 1);
        assert_eq!(tokenized[0], "<unk>");
    }

    #[rstest]
    fn test_tokenize_unk_chrom() {
        let cfg_path = "tests/data/tokenizers/tokenizer.toml";
        let tokenizer =
            Tokenizer::from_config(cfg_path).expect("Failed to create tokenizer from config.");

        let regions = vec![Region {
            chr: "chr999".to_string(),
            start: 50,
            end: 150,
            rest: None,
        }];

        let tokenized = tokenizer.tokenize(&regions);
        assert!(tokenized.is_ok());
        let tokenized = tokenized.unwrap();

        assert_eq!(tokenized.len(), 1);
    }

    #[rstest]
    fn test_tokenize_on_two_crhoms() {
        let cfg_path = "tests/data/tokenizers/tokenizer.toml";
        let tokenizer =
            Tokenizer::from_config(cfg_path).expect("Failed to create tokenizer from config.");

        let regions = vec![
            Region {
                chr: "chr1".to_string(),
                start: 151399441,
                end: 151399547,
                rest: None,
            },
            Region {
                chr: "chr2".to_string(),
                start: 203871220,
                end: 203871381,
                rest: None,
            },
        ];

        let tokenized = tokenizer.tokenize(&regions);
        assert!(tokenized.is_ok());

        let tokenized = tokenized.unwrap();
        assert_eq!(tokenized.len(), 2);

        // chr1:151399432-151399527 -- CONFIRMED IN IGV
        assert_eq!(tokenized[0], "chr1:151399431-151399527");
        assert_eq!(tokenizer.convert_token_to_id(&tokenized[0]), Some(6));

        // chr2:203871201-203871375 -- CONFIRMED IN IGV
        assert_eq!(tokenized[1], "chr2:203871200-203871375");
        assert_eq!(tokenizer.convert_token_to_id(&tokenized[1]), Some(7));
    }

    #[rstest]
    fn test_tokenize_with_multi_overlap() {
        let cfg_path = "tests/data/tokenizers/tokenizer.toml";
        let tokenizer =
            Tokenizer::from_config(cfg_path).expect("Failed to create tokenizer from config.");

        let regions = vec![Region {
            chr: "chr2".to_string(),
            start: 203871346,
            end: 203871616,
            rest: None,
        }];

        let tokenized = tokenizer.tokenize(&regions);
        assert!(tokenized.is_ok());

        let tokenized = tokenized.unwrap();
        assert_eq!(tokenized.len(), 2);

        // chr2:203871201-203871375 -- CONFIRMED IN IGV
        assert_eq!(tokenized[0],"chr2:203871200-203871375");
        assert_eq!(tokenizer.convert_token_to_id(&tokenized[0]), Some(7));

        // chr2:203871388-203871588 -- CONFIRMED IN IGV
        assert_eq!(tokenized[1],"chr2:203871387-203871588");
        assert_eq!(tokenizer.convert_token_to_id(&tokenized[1]), Some(8));
    }

    #[rstest]
    fn test_tokenize_with_order() {
        let cfg_path = "tests/data/tokenizers/peaks.scored.bed";
        let tokenizer =
            Tokenizer::from_auto(cfg_path).expect("Failed to create tokenizer from config.");

        let regions = vec![
            // overlaps chr9:3526184-3526269 (CONFIRMED IN IGV) (SCORE: 1.2233)
            Region {
                chr: "chr9".to_string(),
                start: 3_526_178,
                end: 3_526_249,
                rest: None,
            },
            // overlaps chr9:3526072-3526165 (CONFIRMED IN IGV) (SCORE: 23.8)
            Region {
                chr: "chr9".to_string(),
                start: 3_526_051,
                end: 3_526_145,
                rest: None,
            },
        ];

        let tokenized = tokenizer.tokenize(&regions);

        assert!(tokenized.is_ok());
        let tokenized = tokenized.unwrap();
        assert_eq!(tokenized.len(), 2);
        

        // chr9:3526071-3526165 -- CONFIRMED IN IGV
        assert_eq!(tokenized[0], "chr9:3526071-3526165");
        let first_id = tokenizer.convert_token_to_id(&tokenized[0]);
        assert_eq!(first_id, Some(11));

        // chr9:3526183-3526269 -- CONFIRMED IN IGV
        assert_eq!(tokenized[1], "chr9:3526183-3526269");
        let second_id = tokenizer.convert_token_to_id(&tokenized[1]);
        assert_eq!(second_id, Some(10));
    }
}

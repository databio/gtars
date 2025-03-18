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
    #[error("Universe error: {0}")]
    UniverseError(#[from] crate::tokenizers::universe::UniverseError),
    #[error(transparent)]
    Anyhow(#[from] anyhow::Error),
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

    ///
    /// Tokenize the given set of regions into tokens. This returns a `TokenizedRegionSet` which contains
    /// the IDs of the tokens corresponding to the input regions and a pointer to the universe.
    /// 
    pub fn tokenize(&self, regions: &[Region]) -> Result<TokenizedRegionSet, TokenizerError> {
        let mut res = self.core.tokenize(regions)?;
        if res.is_empty() {
            res.ids.push(self.core.token_to_id(&self.special_tokens.unk).unwrap());
        }
        Ok(res)
    }

    pub fn token_to_id(&self, token: &Region) -> Option<u32> {
        self.core.token_to_id(token)
    }

    pub fn id_to_token(&self, id: u32) -> Option<Region> {
        self.core.id_to_token(id)
    }

    pub fn get_vocab_size(&self) -> usize {
        self.core.get_vocab_size()
    }

    pub fn get_universe(&self) -> &Universe {
        self.core.get_universe()
    }

    pub fn get_special_tokens(&self) -> &SpecialTokens {
        &self.special_tokens
    }

    pub fn get_unk_token(&self) -> &Region {
        &self.special_tokens.unk
    }

    pub fn get_pad_token(&self) -> &Region {
        &self.special_tokens.pad
    }

    pub fn get_mask_token(&self) -> &Region {
        &self.special_tokens.mask
    }

    pub fn get_sep_token(&self) -> &Region {
        &self.special_tokens.sep
    }

    pub fn get_bos_token(&self) -> &Region {
        &self.special_tokens.bos
    }

    pub fn get_eos_token(&self) -> &Region {
        &self.special_tokens.eos
    }

    pub fn get_unk_token_id(&self) -> u32 {
        self.token_to_id(&self.special_tokens.unk).unwrap()
    }

    pub fn get_pad_token_id(&self) -> u32 {
        self.token_to_id(&self.special_tokens.pad).unwrap()
    }

    pub fn get_mask_token_id(&self) -> u32 {
        self.token_to_id(&self.special_tokens.mask).unwrap()
    }

    pub fn get_sep_token_id(&self) -> u32 {
        self.token_to_id(&self.special_tokens.sep).unwrap()
    }

    pub fn get_bos_token_id(&self) -> u32 {
        self.token_to_id(&self.special_tokens.bos).unwrap()
    }

    pub fn get_eos_token_id(&self) -> u32 {
        self.token_to_id(&self.special_tokens.eos).unwrap()
    }
}

#[cfg(test)]
mod tokenizer_tests {
    use super::*;

    use rstest::*;
    use pretty_assertions::assert_eq;

    #[fixture]
    fn anndata_path() -> String {
        "tests/data/tokenizers/pbmc_hg38.h5ad".to_string()
    }

    #[rstest]
    fn test_tokenizer_creation_from_config() {
        let cfg_path = "tests/data/tokenizers/tokenizer.toml";
        let tokenizer = Tokenizer::from_config(cfg_path)
            .expect("Failed to create tokenizer from config.");
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
        let tokenizer = Tokenizer::from_config(cfg_path)
            .expect("Failed to create tokenizer from config.");

        assert_eq!(tokenizer.get_vocab_size(), 32); // 25 regions + 7 special tokens

        // check that unk was overridden
        assert_eq!(tokenizer.get_unk_token().chr, "chrUNKNOWN");
        assert_eq!(tokenizer.get_unk_token().start, 100);
        assert_eq!(tokenizer.get_unk_token().end, 200);

        // check that pad didnt change
        assert_eq!(tokenizer.get_pad_token().chr, "chrPAD");
        assert_eq!(tokenizer.get_pad_token().start, 0);
        assert_eq!(tokenizer.get_pad_token().end, 0);
    }

    #[rstest]
    fn test_tokenize_single_region_not_overlapping() {
        let cfg_path = "tests/data/tokenizers/tokenizer.toml";
        let tokenizer = Tokenizer::from_config(cfg_path)
            .expect("Failed to create tokenizer from config.");

        let regions = vec![Region {
            chr: "chr1".to_string(),
            start: 50,
            end: 150,
            rest: None,
        }];

        let tokenized = tokenizer.tokenize(&regions);
        assert!(tokenized.is_ok());
        let tokenized = tokenized.unwrap();
        assert_eq!(tokenized.ids.len(), 1);
        assert_eq!(tokenized.ids[0], tokenizer.get_unk_token_id());
        assert_eq!(tokenizer.id_to_token(tokenized.ids[0]).unwrap().chr, "chrUNK");
    }

    #[rstest]
    fn test_tokenize_unk_chrom() {
        let cfg_path = "tests/data/tokenizers/tokenizer.toml";
        let tokenizer = Tokenizer::from_config(cfg_path)
            .expect("Failed to create tokenizer from config.");

        let regions = vec![Region {
            chr: "chr999".to_string(),
            start: 50,
            end: 150,
            rest: None,
        }];

        let tokenized = tokenizer.tokenize(&regions);
        assert!(tokenized.is_ok());
        let tokenized = tokenized.unwrap();

        assert_eq!(tokenized.ids.len(), 1);
    }

    #[rstest]
    fn test_tokenize_on_two_crhoms() {
        let cfg_path = "tests/data/tokenizers/tokenizer.toml";
        let tokenizer = Tokenizer::from_config(cfg_path)
            .expect("Failed to create tokenizer from config.");

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

        let tokenized = tokenized.into_region_vec();

        // chr1:151399432-151399527 -- CONFIRMED IN IGV
        assert_eq!(tokenized[0].chr, "chr1");
        assert_eq!(tokenized[0].start, 151399431); // igv shows 151399432 (but we are 0-based)
        assert_eq!(tokenized[0].end, 151399527);
        assert_eq!(tokenizer.token_to_id(&tokenized[0]), Some(6));

        // chr2:203871201-203871375 -- CONFIRMED IN IGV
        assert_eq!(tokenized[1].chr, "chr2");
        assert_eq!(tokenized[1].start, 203871200); // igv shows 203871201 (but we are 0-based)
        assert_eq!(tokenized[1].end, 203871375);
        assert_eq!(tokenizer.token_to_id(&tokenized[1]), Some(7));
    }

    #[rstest]
    fn test_tokenize_with_multi_overlap() {
        let cfg_path = "tests/data/tokenizers/tokenizer.toml";
        let tokenizer = Tokenizer::from_config(cfg_path)
            .expect("Failed to create tokenizer from config.");

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

        let tokenized = tokenized.into_region_vec();

        // chr2:203871201-203871375 -- CONFIRMED IN IGV
        assert_eq!(tokenized[0].chr, "chr2");
        assert_eq!(tokenized[0].start, 203871200); // igv shows 203871201 (but we are 0-based)
        assert_eq!(tokenized[0].end, 203871375);
        assert_eq!(tokenizer.token_to_id(&tokenized[0]), Some(7));

        // chr2:203871388-203871588 -- CONFIRMED IN IGV
        assert_eq!(tokenized[1].chr, "chr2");
        assert_eq!(tokenized[1].start, 203871387); // igv shows 203871388 (but we are 0-based)
        assert_eq!(tokenized[1].end, 203871588);
        assert_eq!(tokenizer.token_to_id(&tokenized[1]), Some(8));
    }
}
use std::collections::HashMap as StdHashMap;
use std::path::{Path, PathBuf};

use fxhash::FxHashMap as HashMap;

use gtars_core::models::Region;
use gtars_overlaprs::Overlapper;

use crate::utils::create_tokenize_core_from_universe;

use super::config::{TokenizerConfig, TokenizerInputFileType, TokenizerType};
use super::error::TokenizerError;
use super::universe::Universe;
use super::utils::prepare_universe_and_special_tokens;
use super::utils::special_tokens::SpecialTokens;

#[cfg(feature = "huggingface")]
use hf_hub::api::sync::Api;

pub const DEFAULT_UNIVERSE_FILENAME: &str = "universe.bed.gz";

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

pub struct Tokenizer {
    core: HashMap<String, Box<dyn Overlapper<u32, u32>>>,
    universe: Universe,
    special_tokens: SpecialTokens,
}

impl Tokenizer {
    ///
    /// Create a new tokenizer
    ///
    pub fn new(
        core: HashMap<String, Box<dyn Overlapper<u32, u32>>>,
        universe: Universe,
        special_tokens: SpecialTokens,
    ) -> Self {
        Self {
            core,
            universe,
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

        let core = create_tokenize_core_from_universe(
            &universe,
            config.tokenizer_type.unwrap_or(TokenizerType::Bits),
        );

        Ok(Tokenizer {
            core,
            special_tokens,
            universe,
        })
    }

    ///
    /// Create a new tokenizer from a bed file
    ///
    pub fn from_bed<P: AsRef<Path>>(bed_path: P) -> Result<Self, TokenizerError> {
        let special_tokens = SpecialTokens::default();
        let (universe, special_tokens) =
            prepare_universe_and_special_tokens(bed_path.as_ref(), special_tokens)?;
        let core = create_tokenize_core_from_universe(&universe, TokenizerType::Bits);

        Ok(Tokenizer {
            core,
            special_tokens,
            universe,
        })
    }

    ///
    /// Create a new tokenizer from a pre-trained model
    ///
    #[cfg(feature = "huggingface")]
    pub fn from_pretrained<P: AsRef<Path>>(path: P) -> Result<Self, TokenizerError> {
        // if local
        let universe_file_path: PathBuf = if path.as_ref().exists() {
            path.as_ref().join(DEFAULT_UNIVERSE_FILENAME)
        } else {
            let api = Api::new().unwrap();
            let repo = api.model(
                path.as_ref()
                    .to_str()
                    .expect("Path is not valid UTF-8")
                    .to_string(),
            );
            repo.get(DEFAULT_UNIVERSE_FILENAME).unwrap()
        };
        let file_type = TokenizerInputFileType::from_path(universe_file_path.as_path())?;
        match file_type {
            TokenizerInputFileType::Toml => Tokenizer::from_config(universe_file_path),
            TokenizerInputFileType::Bed => Tokenizer::from_bed(universe_file_path),
            TokenizerInputFileType::BedGz => Tokenizer::from_bed(universe_file_path),
        }
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
        let tokens = regions
            .iter()
            .filter_map(|region| {
                self.core.get(&region.chr).map(|tree| {
                    tree.find(region.start, region.end)
                        .iter()
                        .map(|i| i.val)
                        .collect::<Vec<u32>>()
                })
            })
            .flatten()
            .map(|id| Token {
                value: self.universe.convert_id_to_token(id).unwrap(),
                id,
            })
            .collect::<Vec<Token>>();

        if tokens.is_empty() {
            return Ok(vec![self.special_tokens.unk.clone()]);
        }

        Ok(tokens.into_iter().map(|t| t.value).collect())
    }

    pub fn encode(&self, regions: &[Region]) -> Result<Vec<u32>, TokenizerError> {
        let tokenized = self.tokenize(regions)?;
        Ok(tokenized
            .into_iter()
            .map(|t| self.convert_token_to_id(&t).unwrap())
            .collect())
    }

    pub fn decode(&self, ids: &[u32]) -> Result<Vec<String>, TokenizerError> {
        let decoded: Vec<String> = ids
            .iter()
            .map(|id| {
                self.universe
                    .convert_id_to_token(*id)
                    .unwrap_or(self.special_tokens.unk.clone())
            })
            .collect();
        Ok(decoded)
    }

    pub fn convert_token_to_id(&self, token: &str) -> Option<u32> {
        self.universe.convert_token_to_id(token)
    }

    pub fn convert_id_to_token(&self, id: u32) -> Option<String> {
        self.universe.convert_id_to_token(id)
    }

    pub fn get_vocab_size(&self) -> usize {
        self.universe.len()
    }

    pub fn get_vocab(&self) -> StdHashMap<String, u32> {
        self.universe.region_to_id.clone()
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

    // ids
    pub fn get_unk_token_id(&self) -> u32 {
        self.convert_token_to_id(&self.special_tokens.unk).unwrap()
    }

    pub fn get_pad_token_id(&self) -> u32 {
        self.convert_token_to_id(&self.special_tokens.pad).unwrap()
    }

    pub fn get_mask_token_id(&self) -> u32 {
        self.convert_token_to_id(&self.special_tokens.mask).unwrap()
    }

    pub fn get_cls_token_id(&self) -> u32 {
        self.convert_token_to_id(&self.special_tokens.cls).unwrap()
    }

    pub fn get_eos_token_id(&self) -> u32 {
        self.convert_token_to_id(&self.special_tokens.eos).unwrap()
    }

    pub fn get_bos_token_id(&self) -> u32 {
        self.convert_token_to_id(&self.special_tokens.bos).unwrap()
    }

    pub fn get_sep_token_id(&self) -> u32 {
        self.convert_token_to_id(&self.special_tokens.sep).unwrap()
    }

    pub fn get_special_tokens_mask(&self, tokens: &[String]) -> Vec<bool> {
        tokens
            .iter()
            .map(|token| {
                token == &self.special_tokens.unk
                    || token == &self.special_tokens.pad
                    || token == &self.special_tokens.mask
                    || token == &self.special_tokens.cls
                    || token == &self.special_tokens.eos
                    || token == &self.special_tokens.bos
                    || token == &self.special_tokens.sep
            })
            .collect()
    }

    pub fn get_special_tokens(&self) -> &SpecialTokens {
        &self.special_tokens
    }

    pub fn get_universe(&self) -> &Universe {
        &self.universe
    }
}

#[cfg(test)]
mod tokenizer_tests {
    use super::*;

    use pretty_assertions::assert_eq;
    use rstest::*;

    #[fixture]
    fn anndata_path() -> String {
        "../tests/data/tokenizers/pbmc_hg38.h5ad".to_string()
    }

    #[rstest]
    fn test_tokenizer_creation_from_config() {
        let cfg_path = "../tests/data/tokenizers/tokenizer.toml";
        let tokenizer =
            Tokenizer::from_config(cfg_path).expect("Failed to create tokenizer from config.");
        assert_eq!(tokenizer.get_vocab_size(), 32); // 25 regions + 7 special tokens
    }

    #[rstest]
    fn test_tokenizer_creation_from_bed() {
        let bed_path = "../tests/data/tokenizers/peaks.bed";
        let tokenizer =
            Tokenizer::from_bed(bed_path).expect("Failed to create tokenizer from config.");
        assert_eq!(tokenizer.get_vocab_size(), 32); // 25 regions + 7 special tokens
    }

    #[rstest]
    fn test_tokenizer_creation_from_bed_gz() {
        let bed_path = "../tests/data/tokenizers/peaks.bed.gz";
        let tokenizer =
            Tokenizer::from_bed(bed_path).expect("Failed to create tokenizer from config.");
        assert_eq!(tokenizer.get_vocab_size(), 32); // 25 regions + 7 special tokens
    }

    #[rstest]
    fn test_tokenizer_creation_auto_all() {
        let bed_path = "../tests/data/tokenizers/peaks.bed";
        let tokenizer =
            Tokenizer::from_auto(bed_path).expect("Failed to create tokenizer from config.");
        assert_eq!(tokenizer.get_vocab_size(), 32); // 25 regions + 7 special tokens

        let cfg_path = "../tests/data/tokenizers/tokenizer.toml";
        let tokenizer =
            Tokenizer::from_auto(cfg_path).expect("Failed to create tokenizer from config.");
        assert_eq!(tokenizer.get_vocab_size(), 32); // 25 regions + 7 special tokens

        let bed_path = "../tests/data/tokenizers/peaks.bed.gz";
        let tokenizer =
            Tokenizer::from_auto(bed_path).expect("Failed to create tokenizer from config.");
        assert_eq!(tokenizer.get_vocab_size(), 32); // 25 regions + 7 special tokens
    }

    #[rstest]
    fn test_tokenizer_bad_tokenizer_type() {
        let cfg_path = "../tests/data/tokenizers/tokenizer_bad_ttype.toml";
        let tokenizer = Tokenizer::from_config(cfg_path);
        assert_eq!(tokenizer.is_err(), true);
    }

    #[rstest]
    fn test_tokenizer_custom_special_tokens() {
        let cfg_path = "../tests/data/tokenizers/tokenizer_custom_specials.toml";
        let tokenizer =
            Tokenizer::from_config(cfg_path).expect("Failed to create tokenizer from config.");

        assert_eq!(tokenizer.get_vocab_size(), 32); // 25 regions + 7 special tokens

        // check that unk was overridden
        assert_eq!(tokenizer.get_unk_token(), "<UNKNOWN>");

        // check that pad didnt change
        assert_eq!(tokenizer.get_pad_token(), "<pad>");
    }

    #[rstest]
    fn test_tokenize_single_region_not_overlapping() {
        let cfg_path = "../tests/data/tokenizers/tokenizer.toml";
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
        let cfg_path = "../tests/data/tokenizers/tokenizer.toml";
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
        let cfg_path = "../tests/data/tokenizers/tokenizer.toml";
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
    fn test_tokenize_on_two_chroms_ailist() {
        let cfg_path = "../tests/data/tokenizers/tokenizer_ailist.toml";
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
        let cfg_path = "../tests/data/tokenizers/tokenizer.toml";
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
        assert_eq!(tokenized[0], "chr2:203871200-203871375");
        assert_eq!(tokenizer.convert_token_to_id(&tokenized[0]), Some(7));

        // chr2:203871388-203871588 -- CONFIRMED IN IGV
        assert_eq!(tokenized[1], "chr2:203871387-203871588");
        assert_eq!(tokenizer.convert_token_to_id(&tokenized[1]), Some(8));
    }
}

use std::collections::HashMap;
use std::path::Path;

use rust_lapper::Lapper;

use super::SpecialTokens;
use super::Tokenizer;
use super::TokenizerError;
use crate::common::models::Region;
use crate::tokenizers::config::{SpecialToken, TokenizerConfig};
use crate::tokenizers::tokens::TokenizedRegionSet;
use crate::tokenizers::universe::Universe;
use crate::tokenizers::utils::create_interval_tree_from_universe;

pub struct TreeTokenizer {
    /// The core interval tree. Actually, its **many** interval trees. The hash-map will map chrom names
    /// to an interval tree for querying. The hash-map lookup should be constant time (O(1)), while
    /// the interval tree is [reported to be NlogN](https://academic.oup.com/bioinformatics/article/29/1/1/273289?login=false)
    tree: HashMap<String, Lapper<u32, u32>>,
    universe: Universe,
    special_tokens: SpecialTokens,
}

impl TryFrom<TokenizerConfig> for TreeTokenizer {
    type Error = TokenizerError;

    fn try_from(value: TokenizerConfig) -> Result<Self, Self::Error> {
        let universe = Path::new(&value.universe);
        let mut universe = Universe::try_from(universe)?;
        let tree = create_interval_tree_from_universe(&universe);

        // we start with the default, then will replace as they exist in the config
        let mut special_tokens = SpecialTokens::default();
        if let Some(config_special_tokens) = value.special_tokens {
            for token in config_special_tokens {
                let token_type = token.name;
                let token_value = token.token;

                let parts = token_value.split(':').collect::<Vec<&str>>();
                if parts.len() != 2 {
                    return Err(TokenizerError::InvalidSpecialTokenConfig);
                }

                let chr = parts[0].to_string();
                let start_end = parts[1].split('-').collect::<Vec<&str>>();
                if start_end.len() != 2 {
                    return Err(TokenizerError::InvalidSpecialTokenConfig);
                }
                let start = start_end[0]
                    .parse::<u32>()
                    .map_err(|_| TokenizerError::InvalidSpecialTokenConfig)?;
                let end = start_end[1]
                    .parse::<u32>()
                    .map_err(|_| TokenizerError::InvalidSpecialTokenConfig)?;

                let token_value = Region {
                    chr,
                    start,
                    end,
                    rest: None,
                };

                match token_type {
                    SpecialToken::Unk => special_tokens.unk = token_value,
                    SpecialToken::Pad => special_tokens.pad = token_value,
                    SpecialToken::Mask => special_tokens.mask = token_value,
                    SpecialToken::Cls => special_tokens.cls = token_value,
                    SpecialToken::Sep => special_tokens.sep = token_value,
                    SpecialToken::Eos => special_tokens.eos = token_value,
                    SpecialToken::Bos => special_tokens.bos = token_value,
                }
            }
        }

        // insert all new special tokens into the universe
        let s_tokens: Vec<Region> = special_tokens.clone().into();
        for s_token in s_tokens.iter() {
            universe.insert_token(s_token);
        }

        Ok(TreeTokenizer {
            tree,
            universe,
            special_tokens,
        })
    }
}

impl Tokenizer for TreeTokenizer {
    fn tokenize<T: Into<Vec<Region>>>(
        &self,
        regions: T,
    ) -> Result<TokenizedRegionSet, TokenizerError> {
        let regions: Vec<Region> = regions.into();
        let mut tokenized_regions = Vec::new();
        for region in regions.iter() {
            if let Some(tree) = self.tree.get(&region.chr) {
                let overlapping_intervals = tree.find(region.start, region.end).map(|interval| interval.val).collect::<Vec<u32>>();
                if overlapping_intervals.is_empty() {
                    tokenized_regions.push(self.token_to_id(&self.special_tokens.unk).unwrap());
                } else {
                    // Assuming we take the first overlapping interval for simplicity
                    tokenized_regions.extend(overlapping_intervals);
                }
            } else {
                tokenized_regions.push(self.token_to_id(&self.special_tokens.unk).unwrap());
            }
        }

        Ok(TokenizedRegionSet {
            ids: tokenized_regions,
            universe: &self.universe,
        })
    }

    fn token_to_id(&self, token: &Region) -> Option<u32> {
        self.universe.convert_region_to_id(token)
    }

    fn id_to_token(&self, id: u32) -> Option<Region> {
        self.universe.convert_id_to_region(id)
    }

    fn get_vocab_size(&self) -> usize {
        self.universe.len()
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    use rstest::*;
    use pretty_assertions::assert_eq;

    #[fixture]
    fn tokenizer_config() -> TokenizerConfig {
        TokenizerConfig {
            universe: "tests/data/peaks.bed".to_string(),
            special_tokens: None
        }
    }

    #[rstest]
    fn test_tree_tokenizer_creation(tokenizer_config: TokenizerConfig) {
        let tokenizer = TreeTokenizer::try_from(tokenizer_config);
        assert_eq!(tokenizer.is_ok(), true);
    }

    #[rstest]
    fn test_universe_size(tokenizer_config: TokenizerConfig) {
        let tokenizer = TreeTokenizer::try_from(tokenizer_config);
        assert_eq!(tokenizer.is_ok(), true);

        let tokenizer = tokenizer.unwrap();
        assert_eq!(tokenizer.get_vocab_size(), 32); // 25 regions + 7 special tokens
    }

    #[rstest]
    fn test_special_tokens(tokenizer_config: TokenizerConfig) {
        let tokenizer = TreeTokenizer::try_from(tokenizer_config);
        assert_eq!(tokenizer.is_ok(), true);

        let tokenizer = tokenizer.unwrap();

        // confirm the special tokens
        assert_eq!(tokenizer.special_tokens.unk.chr, "chrUNK");
        assert_eq!(tokenizer.special_tokens.pad.chr, "chrPAD");
        assert_eq!(tokenizer.special_tokens.mask.chr, "chrMASK");
        assert_eq!(tokenizer.special_tokens.cls.chr, "chrCLS");
        assert_eq!(tokenizer.special_tokens.eos.chr, "chrEOS");
        assert_eq!(tokenizer.special_tokens.bos.chr, "chrBOS");
        assert_eq!(tokenizer.special_tokens.sep.chr, "chrSEP");

        assert_eq!(tokenizer.special_tokens.unk.start, 0);
        assert_eq!(tokenizer.special_tokens.pad.start, 0);
        assert_eq!(tokenizer.special_tokens.mask.start, 0);
        assert_eq!(tokenizer.special_tokens.cls.start, 0);
        assert_eq!(tokenizer.special_tokens.eos.start, 0);
        assert_eq!(tokenizer.special_tokens.bos.start, 0);
        assert_eq!(tokenizer.special_tokens.sep.start, 0);


        assert_eq!(tokenizer.special_tokens.unk.end, 0);
        assert_eq!(tokenizer.special_tokens.pad.end, 0);
        assert_eq!(tokenizer.special_tokens.mask.end, 0);
        assert_eq!(tokenizer.special_tokens.cls.end, 0);
        assert_eq!(tokenizer.special_tokens.eos.end, 0);
        assert_eq!(tokenizer.special_tokens.bos.end, 0);
        assert_eq!(tokenizer.special_tokens.sep.end, 0);
        
        // confirm id values
        assert_eq!(tokenizer.token_to_id(&tokenizer.special_tokens.unk).unwrap(), 25);
        assert_eq!(tokenizer.token_to_id(&tokenizer.special_tokens.pad).unwrap(), 26);
        assert_eq!(tokenizer.token_to_id(&tokenizer.special_tokens.mask).unwrap(), 27);
        assert_eq!(tokenizer.token_to_id(&tokenizer.special_tokens.cls).unwrap(), 28);
        assert_eq!(tokenizer.token_to_id(&tokenizer.special_tokens.eos).unwrap(), 29);
        assert_eq!(tokenizer.token_to_id(&tokenizer.special_tokens.bos).unwrap(), 30);
        assert_eq!(tokenizer.token_to_id(&tokenizer.special_tokens.sep).unwrap(), 31);
    }

    #[rstest]
    fn test_tokenize_single_region(tokenizer_config: TokenizerConfig) {
        let tokenizer = TreeTokenizer::try_from(tokenizer_config).unwrap();

        let regions = vec![
            Region {
                chr: "chr1".to_string(),
                start: 50,
                end: 150,
                rest: None,
            },
        ];

        let tokenized = tokenizer.tokenize(regions);
        assert!(tokenized.is_ok());
        let tokenized = tokenized.unwrap();
        assert_eq!(tokenized.ids.len(), 1);
    }

    #[rstest]
    fn test_tokenize_double_region(tokenizer_config: TokenizerConfig) {
        let tokenizer = TreeTokenizer::try_from(tokenizer_config).unwrap();

        let regions = vec![
            Region {
                chr: "chr1".to_string(),
                start: 50,
                end: 150,
                rest: None,
            },
            Region {
                chr: "chr2".to_string(),
                start: 50,
                end: 150,
                rest: None,
            },
        ];

        let tokenized = tokenizer.tokenize(regions);
        assert!(tokenized.is_ok());
        let tokenized = tokenized.unwrap();
        assert_eq!(tokenized.ids.len(), 2);
    }

    #[rstest]
    fn test_tokenize_unk_chrom(tokenizer_config: TokenizerConfig) {
        let tokenizer = TreeTokenizer::try_from(tokenizer_config).unwrap();

        let regions = vec![
            Region {
                chr: "chr999".to_string(),
                start: 50,
                end: 150,
                rest: None,
            },
        ];

        let tokenized = tokenizer.tokenize(regions);
        assert!(tokenized.is_ok());
        let tokenized = tokenized.unwrap();

        assert_eq!(tokenized.ids.len(), 1);
        assert_eq!(tokenized.ids[0], tokenizer.token_to_id(&tokenizer.special_tokens.unk).unwrap());
    }

    #[rstest]
    fn test_tokenize_on_two_crhoms(tokenizer_config: TokenizerConfig) {
        let tokenizer = TreeTokenizer::try_from(tokenizer_config).unwrap();

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

        let tokenized = tokenizer.tokenize(regions);
        assert!(tokenized.is_ok());
        
        let tokenized = tokenized.unwrap();
        assert_eq!(tokenized.len(), 2);
        
        let tokenized = tokenized.into_region_vec();

        // chr1:151399432-151399527 -- CONFIRMED IN IGV
        assert_eq!(tokenized[0].chr, "chr1");
        assert_eq!(tokenized[0].start, 151399431);
        assert_eq!(tokenized[0].end, 151399527);
        assert_eq!(tokenizer.token_to_id(&tokenized[0]), Some(6));

        // chr2:203871201-203871375 -- CONFIRMED IN IGV
        assert_eq!(tokenized[1].chr, "chr2");
        assert_eq!(tokenized[1].start, 203871200);
        assert_eq!(tokenized[1].end, 203871375);
        assert_eq!(tokenizer.token_to_id(&tokenized[1]), Some(7));


    }

    #[rstest]
    fn test_tokenize_with_multi_overlap(tokenizer_config: TokenizerConfig) {
        let tokenizer = TreeTokenizer::try_from(tokenizer_config).unwrap();

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

        let tokenized = tokenizer.tokenize(regions);
        assert!(tokenized.is_ok());
        
        let tokenized = tokenized.unwrap();
        assert_eq!(tokenized.len(), 2);
        
        let tokenized = tokenized.into_region_vec();

        // chr1:151399432-151399527 -- CONFIRMED IN IGV
        assert_eq!(tokenized[0].chr, "chr1");
        assert_eq!(tokenized[0].start, 151399431);
        assert_eq!(tokenized[0].end, 151399527);

        // chr2:203871201-203871375 -- CONFIRMED IN IGV
        assert_eq!(tokenized[1].chr, "chr2");
        assert_eq!(tokenized[1].start, 203871200);
        assert_eq!(tokenized[1].end, 203871375);


    }

    #[rstest]
    fn test_token_to_id_and_id_to_token(tokenizer_config: TokenizerConfig) {
        let tokenizer = TreeTokenizer::try_from(tokenizer_config).unwrap();

        let region = Region {
            chr: "chr1".to_string(),
            start: 0,
            end: 100,
            rest: None,
        };

        let id = tokenizer.token_to_id(&region);
        assert!(id.is_some());

        let converted_region = tokenizer.id_to_token(id.unwrap());
        assert!(converted_region.is_some());
        assert_eq!(converted_region.unwrap(), region);
    }

    #[rstest]
    fn test_get_vocab_size(tokenizer_config: TokenizerConfig) {
        let tokenizer = TreeTokenizer::try_from(tokenizer_config).unwrap();

        let vocab_size = tokenizer.get_vocab_size();
        assert!(vocab_size > 0);
    }
}

pub mod utils;

use std::collections::HashMap;

use rayon::prelude::*;
use rust_lapper::Lapper;

use self::utils::create_interval_tree_from_universe;
use super::GTokenize;
use super::Token;
use super::TokenizerError;
use crate::common::models::Region;
use crate::tokenizers::universe::Universe;

#[derive(Clone, Debug)]
pub struct BitsTree {
    /// The core interval tree. Actually, its **many** interval trees. The hash-map will map chrom names
    /// to an interval tree for querying. The hash-map lookup should be constant time (O(1)), while
    /// the interval tree is [reported to be NlogN](https://academic.oup.com/bioinformatics/article/29/1/1/273289?login=false)
    tree: HashMap<String, Lapper<u32, u32>>,
    universe: Universe,
}

impl From<Universe> for BitsTree {
    fn from(universe: Universe) -> Self {
        let tree = create_interval_tree_from_universe(&universe);

        BitsTree { tree, universe }
    }
}

impl GTokenize for BitsTree {
    fn tokenize(&self, regions: &[Region]) -> Result<Vec<Token>, TokenizerError> {
        let mut tokens = regions
            .par_iter()
            .filter_map(|region| {
                self.tree.get(&region.chr).map(|tree| {
                    tree.find(region.start, region.end)
                        .map(|i| i.val)
                        .collect::<Vec<u32>>()
                })
            })
            .flatten()
            .map(|id| Token {
                value: self.universe.convert_id_to_token(id).unwrap(),
                id
            })
            .collect::<Vec<Token>>();

        if let Some(scores) = &self.universe.scores {
            tokens.sort_by(|a, b| {
                let score_a = scores.get(&a.value).unwrap_or(&0.0);
                let score_b = scores.get(&b.value).unwrap_or(&0.0);
                score_b.partial_cmp(score_a).unwrap()
            });
            Ok(tokens)
        } else {
            Ok(tokens)
        }
    }

    fn token_to_id(&self, token: &str) -> Option<u32> {
        self.universe.convert_token_to_id(token)
    }

    fn id_to_token(&self, id: u32) -> Option<String> {
        self.universe.convert_id_to_token(id)
    }

    fn get_vocab_size(&self) -> usize {
        self.universe.len()
    }

    fn get_universe(&self) -> &Universe {
        &self.universe
    }

    fn get_vocab(&self) -> HashMap<String, u32> {
        self.universe.region_to_id.clone()
    }

}

#[cfg(test)]
mod tests {
    use crate::tokenizers::utils::{
        prepare_universe_and_special_tokens, special_tokens::SpecialTokens,
    };

    use super::*;

    use pretty_assertions::assert_eq;

    use rstest::*;

    #[fixture]
    fn universe() -> Universe {
        let universe_file = "tests/data/tokenizers/peaks.bed";
        let special_tokens = SpecialTokens::default();

        let (universe, _) = prepare_universe_and_special_tokens(universe_file, special_tokens)
            .expect("Failed to prepare universe and special tokens");

        universe
    }

    #[rstest]
    fn test_tree_tokenizer_creation(universe: Universe) {
        let _tokenizer = BitsTree::from(universe);
    }

    #[rstest]
    fn test_universe_size(universe: Universe) {
        let tokenizer = BitsTree::from(universe);
        assert_eq!(tokenizer.get_vocab_size(), 32); // 25 regions + 7 special tokens
    }

    #[rstest]
    fn test_special_tokens(universe: Universe) {
        let tokenizer = BitsTree::from(universe);
        let special_tokens = SpecialTokens::default();

        // confirm the special tokens
        assert_eq!(special_tokens.unk, "<unk>");
        assert_eq!(special_tokens.pad, "<pad>");
        assert_eq!(special_tokens.mask, "<mask>");
        assert_eq!(special_tokens.cls, "<cls>");
        assert_eq!(special_tokens.eos, "<eos>");
        assert_eq!(special_tokens.bos, "<bos>");
        assert_eq!(special_tokens.sep, "<sep>");

        // confirm id values
        assert_eq!(tokenizer.token_to_id(&special_tokens.unk).unwrap(), 25);
        assert_eq!(tokenizer.token_to_id(&special_tokens.pad).unwrap(), 26);
        assert_eq!(tokenizer.token_to_id(&special_tokens.mask).unwrap(), 27);
        assert_eq!(tokenizer.token_to_id(&special_tokens.cls).unwrap(), 28);
        assert_eq!(tokenizer.token_to_id(&special_tokens.eos).unwrap(), 29);
        assert_eq!(tokenizer.token_to_id(&special_tokens.bos).unwrap(), 30);
        assert_eq!(tokenizer.token_to_id(&special_tokens.sep).unwrap(), 31);
    }

    #[rstest]
    fn test_tokenize_single_region_not_overlapping(universe: Universe) {
        let tokenizer = BitsTree::from(universe);

        let regions = vec![Region {
            chr: "chr1".to_string(),
            start: 50,
            end: 150,
            rest: None,
        }];

        let tokenized = tokenizer.tokenize(&regions);
        assert!(tokenized.is_ok());
        let tokenized = tokenized.unwrap();
        assert_eq!(tokenized.len(), 0);
    }

    #[rstest]
    fn test_tokenize_unk_chrom(universe: Universe) {
        let tokenizer = BitsTree::from(universe);

        let regions = vec![Region {
            chr: "chr999".to_string(),
            start: 50,
            end: 150,
            rest: None,
        }];

        let tokenized = tokenizer.tokenize(&regions);
        assert!(tokenized.is_ok());
        let tokenized = tokenized.unwrap();

        assert_eq!(tokenized.len(), 0);
    }

    #[rstest]
    fn test_tokenize_on_two_crhoms(universe: Universe) {
        let tokenizer = BitsTree::from(universe);

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
        assert_eq!(tokenized[0].value, "chr1:151399431-151399527");  // igv shows 151399432 (but we are 0-based)
        assert_eq!(tokenizer.token_to_id(&tokenized[0].value), Some(6));

        // chr2:203871201-203871375 -- CONFIRMED IN IGV
        assert_eq!(tokenized[1].value, "chr2:203871200-203871375");  // igv shows 203871201 (but we are 0-based)
        assert_eq!(tokenizer.token_to_id(&tokenized[1].value), Some(7));
    }

    #[rstest]
    fn test_tokenize_with_multi_overlap(universe: Universe) {
        let tokenizer = BitsTree::from(universe);

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
        assert_eq!(tokenized[0].value, "chr2:203871200-203871375"); // igv shows 203871201 (but we are 0-based)
        assert_eq!(tokenizer.token_to_id(&tokenized[0].value), Some(7));

        // chr2:203871388-203871588 -- CONFIRMED IN IGV
        assert_eq!(tokenized[1].value, "chr2:203871387-203871588"); // igv shows 203871388 (but we are 0-based)
        assert_eq!(tokenizer.token_to_id(&tokenized[1].value), Some(8));
    }
}

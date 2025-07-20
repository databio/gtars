// #[cfg(test)]
// mod tests {
//     use crate::tokenizers::utils::{
//         prepare_universe_and_special_tokens, special_tokens::SpecialTokens,
//     };

//     use super::*;

//     use pretty_assertions::assert_eq;

//     use rstest::*;

//     #[fixture]
//     fn universe() -> Universe {
//         let universe_file = "tests/data/tokenizers/peaks.bed";
//         let special_tokens = SpecialTokens::default();

//         let (universe, _) = prepare_universe_and_special_tokens(universe_file, special_tokens)
//             .expect("Failed to prepare universe and special tokens");

//         universe
//     }

//     #[rstest]
//     fn test_tree_tokenizer_creation(universe: Universe) {
//         let _tokenizer = BitsTree::from(universe);
//     }

//     #[rstest]
//     fn test_universe_size(universe: Universe) {
//         let tokenizer = BitsTree::from(universe);
//         assert_eq!(tokenizer.get_vocab_size(), 32); // 25 regions + 7 special tokens
//     }

//     #[rstest]
//     fn test_special_tokens(universe: Universe) {
//         let tokenizer = BitsTree::from(universe);
//         let special_tokens = SpecialTokens::default();

//         // confirm the special tokens
//         assert_eq!(special_tokens.unk, "<unk>");
//         assert_eq!(special_tokens.pad, "<pad>");
//         assert_eq!(special_tokens.mask, "<mask>");
//         assert_eq!(special_tokens.cls, "<cls>");
//         assert_eq!(special_tokens.eos, "<eos>");
//         assert_eq!(special_tokens.bos, "<bos>");
//         assert_eq!(special_tokens.sep, "<sep>");

//         // confirm id values
//         assert_eq!(tokenizer.token_to_id(&special_tokens.unk).unwrap(), 25);
//         assert_eq!(tokenizer.token_to_id(&special_tokens.pad).unwrap(), 26);
//         assert_eq!(tokenizer.token_to_id(&special_tokens.mask).unwrap(), 27);
//         assert_eq!(tokenizer.token_to_id(&special_tokens.cls).unwrap(), 28);
//         assert_eq!(tokenizer.token_to_id(&special_tokens.eos).unwrap(), 29);
//         assert_eq!(tokenizer.token_to_id(&special_tokens.bos).unwrap(), 30);
//         assert_eq!(tokenizer.token_to_id(&special_tokens.sep).unwrap(), 31);
//     }

//     #[rstest]
//     fn test_tokenize_single_region_not_overlapping(universe: Universe) {
//         let tokenizer = BitsTree::from(universe);

//         let regions = vec![Region {
//             chr: "chr1".to_string(),
//             start: 50,
//             end: 150,
//             rest: None,
//         }];

//         let tokenized = tokenizer.tokenize(&regions);
//         assert!(tokenized.is_ok());
//         let tokenized = tokenized.unwrap();
//         assert_eq!(tokenized.len(), 0);
//     }

//     #[rstest]
//     fn test_tokenize_unk_chrom(universe: Universe) {
//         let tokenizer = BitsTree::from(universe);

//         let regions = vec![Region {
//             chr: "chr999".to_string(),
//             start: 50,
//             end: 150,
//             rest: None,
//         }];

//         let tokenized = tokenizer.tokenize(&regions);
//         assert!(tokenized.is_ok());
//         let tokenized = tokenized.unwrap();

//         assert_eq!(tokenized.len(), 0);
//     }

//     #[rstest]
//     fn test_tokenize_on_two_crhoms(universe: Universe) {
//         let tokenizer = BitsTree::from(universe);

//         let regions = vec![
//             Region {
//                 chr: "chr1".to_string(),
//                 start: 151399441,
//                 end: 151399547,
//                 rest: None,
//             },
//             Region {
//                 chr: "chr2".to_string(),
//                 start: 203871220,
//                 end: 203871381,
//                 rest: None,
//             },
//         ];

//         let tokenized = tokenizer.tokenize(&regions);
//         assert!(tokenized.is_ok());

//         let tokenized = tokenized.unwrap();
//         assert_eq!(tokenized.len(), 2);

//         // chr1:151399432-151399527 -- CONFIRMED IN IGV
//         assert_eq!(tokenized[0].value, "chr1:151399431-151399527"); // igv shows 151399432 (but we are 0-based)
//         assert_eq!(tokenizer.token_to_id(&tokenized[0].value), Some(6));

//         // chr2:203871201-203871375 -- CONFIRMED IN IGV
//         assert_eq!(tokenized[1].value, "chr2:203871200-203871375"); // igv shows 203871201 (but we are 0-based)
//         assert_eq!(tokenizer.token_to_id(&tokenized[1].value), Some(7));
//     }

//     #[rstest]
//     fn test_tokenize_with_multi_overlap(universe: Universe) {
//         let tokenizer = BitsTree::from(universe);

//         let regions = vec![Region {
//             chr: "chr2".to_string(),
//             start: 203871346,
//             end: 203871616,
//             rest: None,
//         }];

//         let tokenized = tokenizer.tokenize(&regions);
//         assert!(tokenized.is_ok());

//         let tokenized = tokenized.unwrap();
//         assert_eq!(tokenized.len(), 2);

//         // chr2:203871201-203871375 -- CONFIRMED IN IGV
//         assert_eq!(tokenized[0].value, "chr2:203871200-203871375"); // igv shows 203871201 (but we are 0-based)
//         assert_eq!(tokenizer.token_to_id(&tokenized[0].value), Some(7));

//         // chr2:203871388-203871588 -- CONFIRMED IN IGV
//         assert_eq!(tokenized[1].value, "chr2:203871387-203871588"); // igv shows 203871388 (but we are 0-based)
//         assert_eq!(tokenizer.token_to_id(&tokenized[1].value), Some(8));
//     }
// }

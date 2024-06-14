//! # Genomic data tokenizers and pre-processors to prepare interval data for machine learning pipelines.
//!
//! The tokenizers module is the most comprehensive module in `gtars`. It houses all tokenizers that implement
//! tokenization of genomic data into a known vocabulary. This is especially useful for genomic data machine
//! learning models that are based on NLP-models like tranformers.
//! 
//! ## Example
//! ### Create a tokenizer and tokenize a bed file
//! ```rust
//! use std::path::Path;
//! 
//! use gtars::tokenizers::{Tokenizer, TreeTokenizer};
//! use gtars::common::models::RegionSet;
//! 
//! let path_to_bed_file = "tests/data/peaks.bed.gz";
//! let tokenizer = TreeTokenizer::try_from(Path::new(path_to_bed_file)).unwrap();
//! 
//! let path_to_tokenize_bed_file = "tests/data/to_tokenize.bed";
//! let rs = RegionSet::try_from(Path::new(path_to_tokenize_bed_file)).unwrap();
//! 
//! let tokenized_regions = tokenizer.tokenize_region_set(&rs);
//! println!("{:?}", tokenized_regions.ids);
//! ```
pub mod cli;
pub mod config;
pub mod fragment_tokenizer;
pub mod soft_tokenizer;
pub mod special_tokens;
pub mod traits;
pub mod tree_tokenizer;
pub mod meta_tokenizer;

/// constants for the tokenizer module.
pub mod consts {
    /// command for the `gtars` cli
    pub const TOKENIZE_CMD: &str = "tokenize";
    pub const UNIVERSE_FILE_NAME: &str = "universe.bed";
}

// expose the TreeTokenizer struct to users of this crate
pub use config::TokenizerConfig;
pub use fragment_tokenizer::FragmentTokenizer;
pub use traits::{SingleCellTokenizer, Tokenizer};
pub use tree_tokenizer::TreeTokenizer;

#[cfg(test)]
mod tests {

    use crate::common::models::{Region, RegionSet};
    use crate::tokenizers::traits::SpecialTokens;
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

    #[fixture]
    fn path_to_bad_config_file() -> &'static str {
        "tests/data/tokenizer_bad.toml"
    }

    #[fixture]
    fn path_to_tokenize_bed_file() -> &'static str {
        "tests/data/to_tokenize.bed"
    }

    #[rstest]
    fn test_create_tokenizer_from_bed(path_to_bed_file: &str) {
        let tokenizer = TreeTokenizer::try_from(Path::new(path_to_bed_file)).unwrap();
        assert_eq!(tokenizer.vocab_size(), 32); // 25 regions + 7 special tokens
    }

    #[rstest]
    fn test_create_tokenizer_from_config(path_to_config_file: &str) {
        let tokenizer = TreeTokenizer::try_from(Path::new(path_to_config_file)).unwrap();
        assert_eq!(tokenizer.vocab_size(), 56); // 25 regions in main universe + 24 in hierarchical + 7 special tokens
    }

    #[rstest]
    #[should_panic]
    fn test_bad_config_file(path_to_bad_config_file: &str) {
        let tokenizer = TreeTokenizer::try_from(Path::new(path_to_bad_config_file));
        let _tokenizer = tokenizer.unwrap();
    }

    #[rstest]
    fn test_get_special_token_ids(path_to_bed_file: &str) {
        let tokenizer = TreeTokenizer::try_from(Path::new(path_to_bed_file)).unwrap();
        let unk_id = tokenizer.unknown_token_id();
        let pad_id = tokenizer.padding_token_id();
        let mask_id = tokenizer.mask_token_id();
        let eos_id = tokenizer.eos_token_id();
        let bos_id = tokenizer.bos_token_id();
        let cls_id = tokenizer.cls_token_id();
        let sep_id = tokenizer.sep_token_id();

        assert_eq!(unk_id, 25);
        assert_eq!(pad_id, 26);
        assert_eq!(mask_id, 27);
        assert_eq!(eos_id, 28);
        assert_eq!(bos_id, 29);
        assert_eq!(cls_id, 30);
        assert_eq!(sep_id, 31);
    }

    #[rstest]
    fn test_tokenize_bed_file(path_to_bed_file: &str, path_to_tokenize_bed_file: &str) {
        let tokenizer = TreeTokenizer::try_from(Path::new(path_to_bed_file)).unwrap();
        let rs = RegionSet::try_from(Path::new(path_to_tokenize_bed_file)).unwrap();
        let tokenized_regions = tokenizer.tokenize_region_set(&rs);

        println!("{}", tokenized_regions.len());
        assert_eq!(tokenized_regions.len(), 4);

        // last should be the unknown token
        let unknown_token = tokenizer
            .universe
            .convert_id_to_region(tokenized_regions[3])
            .unwrap();
        assert!(unknown_token.chr == "chrUNK");
    }

    #[rstest]
    fn test_hierarchical_universe_hit(path_to_config_file: &str) {
        let tokenizer = TreeTokenizer::try_from(Path::new(path_to_config_file)).unwrap();
        let res = tokenizer.tokenize_region(&Region {
            chr: "chr1".to_string(),
            start: 100,
            end: 200,
        });
        assert_eq!(res.len(), 1);

        // check the id, it should be len(primary_universe) + 1 (since its chr1)
        assert_eq!(res.ids, vec![25]);

        let res = res.into_region_vec();
        let region = &res[0];

        assert_eq!(region.chr, "chr1");
        assert_eq!(region.start, 0);
        assert_eq!(region.end, 248_956_422);
    }

    #[rstest]
    fn test_hierarchical_universe_no_hit(path_to_config_file: &str) {
        let tokenizer = TreeTokenizer::try_from(Path::new(path_to_config_file)).unwrap();
        let res = tokenizer.tokenize_region(&Region {
            chr: "chrFOO".to_string(),
            start: 100,
            end: 200,
        });
        assert_eq!(res.len(), 1);

        // check the id, it should be the id of the UNK token
        assert_eq!(res.ids, vec![49]);

        let res = res.into_region_vec();
        let region = &res[0];

        assert_eq!(region.chr, "chrUNK");
        assert_eq!(region.start, 0);
        assert_eq!(region.end, 0);
    }
}

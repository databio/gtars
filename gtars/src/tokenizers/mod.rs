//! # Tokenizers - tokenize new genomic intervals into a known universe for machine-learning pipelines
//!
//! There is currently only one tokenizer - the `TreeTokenizer`
pub mod cli;
pub mod config;
pub mod fragment_tokenizer;
pub mod soft_tokenizer;
pub mod special_tokens;
pub mod traits;
pub mod tree_tokenizer;

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

    use crate::common::models::RegionSet;
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
}

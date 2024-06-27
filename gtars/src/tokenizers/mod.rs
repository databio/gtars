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
pub mod builder;
pub mod cli;
pub mod config;
pub mod fragment_tokenizer;
pub mod meta_tokenizer;
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
pub use builder::TokenizerBuilder;
pub use config::TokenizerConfig;
pub use fragment_tokenizer::FragmentTokenizer;
pub use meta_tokenizer::MetaTokenizer;
pub use traits::{SingleCellTokenizer, Tokenizer};
pub use tree_tokenizer::TreeTokenizer;

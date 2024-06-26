//! # Tokenizers - tokenize new genomic intervals into a known universe for machine-learning pipelines
//!
//! There is currently only one tokenizer - the `TreeTokenizer`
pub mod cli;
pub mod fragment_tokenizer;
pub mod soft_tokenizer;
pub mod special_tokens;
pub mod traits;
pub mod tree_tokenizer;
pub mod config;

/// constants for the tokenizer module.
pub mod consts {
    /// command for the `gtars` cli
    pub const TOKENIZE_CMD: &str = "tokenize";
    pub const UNIVERSE_FILE_NAME: &str = "universe.bed";
}

// expose the TreeTokenizer struct to users of this crate
pub use fragment_tokenizer::FragmentTokenizer;
pub use traits::{SingleCellTokenizer, Tokenizer};
pub use tree_tokenizer::TreeTokenizer;
pub use config::TokenizerConfig;
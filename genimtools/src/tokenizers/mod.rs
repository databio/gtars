//! # Tokenizers - tokenize new genomic intervals into a known universe for machine-learning pipelines
//! 
//! There is currently only one tokenizer - the `TreeTokenizer`
pub mod traits;
pub mod tree_tokenizer;
pub mod cli;

/// constants for the tokenizer module.
pub mod consts {
    /// command for the `genimtools` cli
    pub const TOKENIZE_CMD: &str = "tokenize";
}

// expose the TreeTokenizer struct to users of this crate
pub use traits::Tokenizer;
pub use tree_tokenizer::TreeTokenizer;

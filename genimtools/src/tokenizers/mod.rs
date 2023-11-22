pub mod traits;
pub mod tree_tokenizer;
pub mod cli;

pub mod consts {
    pub const TOKENIZE_CMD: &str = "tokenize";
}

// expose the TreeTokenizer struct to users of this crate
pub use traits::Tokenizer;
pub use tree_tokenizer::TreeTokenizer;

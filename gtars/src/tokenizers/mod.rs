//! # Genomic data tokenizers and pre-processors to prepare interval data for machine learning pipelines.
//!
//! The tokenizers module is the most comprehensive module in `gtars`. It houses all tokenizers that implement
//! tokenization of genomic data into a known vocabulary. This is especially useful for genomic data machine
//! learning models that are based on NLP-models like tranformers.
//! 
pub mod config;
pub mod tokenizer_impl;
pub mod tokens;
pub mod universe;
pub mod utils;

// re-export things
pub use tokenizer_impl::bits_tree::*;
pub use tokenizer_impl::*;
pub use tokens::*;
pub use universe::*;
pub use utils::*;

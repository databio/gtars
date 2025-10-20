//! # Genomic data tokenizers and pre-processors to prepare interval data for machine learning pipelines.
//!
//! The tokenizers module is the most comprehensive module in `gtars`. It houses all tokenizers that implement
//! tokenization of genomic data into a known vocabulary. This is especially useful for genomic data machine
//! learning models that are based on NLP-models like tranformers.
//!
//! # Example
//! ```rust
//! use std::path::Path;
//! use std::error::Error; // Required for Box<dyn Error> in the example
//!
//! use gtars_tokenizers::{Tokenizer, TokenizerError}; // TokenizerError will be available via re-export
//! use gtars_core::models::Region;
//!
//! // The example has been updated to use Result propagation (`?` operator)
//! // instead of `unwrap()`, demonstrating the new error handling approach.
//! fn main() -> Result<(), Box<dyn Error>> {
//!     let tokenizer = Tokenizer::from_bed(Path::new("../tests/data/tokenizers/peaks.bed"))?;
//!
//!     let regions = vec![Region {
//!         chr: "chr1".to_string(),
//!         start: 100,  
//!         end: 200,
//!         rest: None,
//!     }];
//!     let tokens = tokenizer.tokenize(&regions);
//!
//!     // Example must return Ok(()) to compile as a doctest
//!     Ok(())
//! }
//! ```
//!
pub mod config;
pub mod encoding;
pub mod error;
pub mod tokenizer;
pub mod universe;
pub mod utils;

// re-export things
pub use encoding::*;
pub use error::*;
pub use tokenizer::*;
pub use universe::*;
pub use utils::*;

// contants
pub mod consts {
    pub const TOKENIZERS_CMD: &str = "tokenize";
}
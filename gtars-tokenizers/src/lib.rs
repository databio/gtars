//! # Genomic data tokenizers and pre-processors to prepare interval data for machine learning pipelines.
//!
//! The tokenizers module is the most comprehensive module in `gtars`. It houses all tokenizers that implement
//! tokenization of genomic data into a known vocabulary. This is especially useful for genomic data machine
//! learning models that are based on NLP-models like tranformers.
//!
//! # Example
//! ```rust
//! use std::path::Path;
//!
//! use gtars_tokenizers::Tokenizer;
//! use gtars_core::models::Region;
//!
//! let tokenizer = Tokenizer::from_bed(Path::new("../tests/data/tokenizers/peaks.bed")).unwrap();
//!
//! let regions = vec![Region {
//!     chr: "chr1".to_string(),
//!     start: 100,  
//!     end: 200,
//!     rest: None,
//! }];
//! let tokens = tokenizer.tokenize(&regions);
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

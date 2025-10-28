//! # gtars-tokenizers
//!
//! Wrapper around gtars-overlaprs for producing tokens for machine learning models.
//!
//! ## Purpose
//!
//! This module wraps the core overlap infrastructure from gtars-overlaprs to convert
//! genomic regions into vocabulary tokens for machine learning pipelines. It is
//! specifically designed for ML applications that need to represent genomic intervals
//! as discrete tokens.
//!
//! ## Design Philosophy
//!
//! All overlap computation is delegated to gtars-overlaprs. This module focuses on:
//! - Token vocabulary management
//! - Encoding/decoding strategies
//! - Integration with ML frameworks (HuggingFace, etc.)
//!
//! ## Use Cases
//!
//! - **Transformer Models**: Convert genomic regions to token sequences
//! - **Feature Extraction**: Represent intervals as discrete features for ML
//! - **Language Model Input**: Prepare genomic data for NLP-based models
//!
//! ## Main Components
//!
//! - **`Tokenizer`**: Maps regions to vocabulary tokens using overlap detection
//! - **`Universe`**: Vocabulary of genomic regions (peaks/intervals)
//!
//! ## Example
//!
//! ```rust
//! use std::path::Path;
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

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
pub mod config;
pub mod tokenizers;
pub mod tokens;
pub mod utils;
pub mod universe;
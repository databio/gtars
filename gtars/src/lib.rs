//! # gtars: *<small>Performance-critical tools to manipulate, analyze, and process genomic interval data. </small>*
//!
//! `gtars` is a rust crate that provides a set of tools for working with genomic interval data. Its primary goal is to provide
//! processors for our python package, [`geniml`](https:github.com/databio/geniml), a library for machine learning on genomic intervals.
//! However, it can be used as a standalone library for working with genomic intervals as well.
//!
//! There are several modules in this crate. The most comprehensive is the [tokenizers] modules which houses genomic region tokenizers
//! for use as pre-processors to machine learning pipelines.
//!
//! ## Examples
//! ### Create a tokenizer and tokenize a bed file
//! ```rust
//! use std::path::Path;
//!
//! use gtars::tokenizers::{Tokenizer, TreeTokenizer};
//! use gtars::common::models::RegionSet;
//!
//! let path_to_bed_file = "tests/data/peaks.bed";
//! let tokenizer = TreeTokenizer::try_from(Path::new(path_to_bed_file)).unwrap();
//!
//! let path_to_tokenize_bed_file = "tests/data/to_tokenize.bed";
//! let rs = RegionSet::try_from(Path::new(path_to_tokenize_bed_file)).unwrap();
//!
//! let tokenized_regions = tokenizer.tokenize_region_set(&rs);
//! println!("{:?}", tokenized_regions.ids);
//! ```
//!
//! You can save the result of this tokenization to a file for later use in machine learning model training:
//! ### Write tokens to `gtok` file for later use:
//! ```rust
//! use gtars::io::write_tokens_to_gtok;
//!
//! let ids = vec![42, 101, 999];
//! write_tokens_to_gtok("tests/data/out/tokens.gtok", &ids);
//! ```
pub mod ailist;
pub mod common;
pub mod digests;
pub mod fragsplit;
pub mod igd;
pub mod io;
pub mod scoring;
pub mod tokenizers;
pub mod uniwig;

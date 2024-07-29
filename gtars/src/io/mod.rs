//! # Input/Output utilities for genomic data.
//!
//! This small module provides some small, but convenient utility functions for writing and reading
//! genomic data to and from disk. Most importantly, it contains functions for saving and reading
//! `.gtok` files to disk - special files that store pre-tokenized genomic data for use in machine
//! learning pipelines.
//!
//! ## Examples
//! ### Save tokens to disk
//! ```rust
//! use gtars::io::write_tokens_to_gtok;
//!
//! let ids = vec![42, 101, 999];
//! write_tokens_to_gtok("tests/data/out/tokens.gtok", &ids);
//! ```
//! ### Read tokens from disk
//! ```rust
//! use gtars::io::read_tokens_from_gtok;
//! let ids = read_tokens_from_gtok("tests/data/out/tokens.gtok").unwrap();
//!
//! println!("{:?}", ids); // [42, 101, 999]
//! ```
pub mod consts;
pub mod gtok;

// re-expose core functions
pub use consts::*;
pub use gtok::*;

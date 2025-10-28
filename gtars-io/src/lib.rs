//! # Input/Output utilities for genomic data.
//!
//! This small module provides some small, but convenient utility functions for writing and reading
//! genomic data to and from disk. Most importantly, it contains functions for saving and reading
//! `.gtok` files to disk - special files that store pre-tokenized genomic data for use in machine
//! learning pipelines.
//!
pub mod consts;
pub mod error;
pub mod gtok;
pub mod bed;
#[cfg(feature = "bigbed")]
pub mod bigbed;

// re-expose core functions
pub use consts::*;
pub use error::*;
pub use gtok::*;
pub use bed::*;
#[cfg(feature = "bigbed")]
pub use bigbed::*;
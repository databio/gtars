//!
//! BM25 sparse embedding implementation for genomic intervals and information retrieval.
//!
//! This crate enables powerful hybrid search when paired with a dense embedding model like Atacformer.
//!
pub mod bm25;
pub mod sparse_vector;

pub use bm25::{Bm25, Bm25Builder};
pub use sparse_vector::SparseVector;
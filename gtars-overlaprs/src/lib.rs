//! # gtars-overlaprs
//!
//! High-performance genomic interval overlap operations in Rust.
//!
//! This crate provides efficient data structures and algorithms for finding overlapping intervals
//! in genomic data. It is part of the [gtars](https://github.com/databio/gtars) project, which
//! provides tools for working with genomic interval data in Rust, Python, and R.
//!
//! ## Features
//!
//! - **Fast overlap queries**: Efficiently find all intervals that overlap with a query interval
//! - **Iterator-based API**: Memory-efficient iteration over overlapping intervals
//! - **Thread-safe**: All data structures implement `Send` and `Sync` for concurrent access
//!
//! ## Quick Start
//!
//! ```rust
//! use gtars_overlaprs::{AiList, Overlapper, Interval};
//!
//! // create some genomic intervals (e.g., ChIP-seq peaks)
//! let intervals = vec![
//!     Interval { start: 100u32, end: 200, val: "gene1" },
//!     Interval { start: 150, end: 300, val: "gene2" },
//!     Interval { start: 400, end: 500, val: "gene3" },
//! ];
//!
//! // build the AIList data structure
//! let ailist = AiList::build(intervals);
//!
//! // query for overlapping intervals
//! let overlaps = ailist.find(180, 250);
//! assert_eq!(overlaps.len(), 2); // gene1 and gene2 overlap
//!
//! // or use an iterator for memory-efficient processing
//! for interval in ailist.find_iter(180, 250) {
//!     println!("Found overlap: {:?}", interval);
//! }
//! ```
//!
//! ## Performance
//!
//! The [`AiList`] data structure is optimized for queries on genomic-scale datasets and provides
//! excellent performance for typical genomic interval overlap operations. It uses a decomposition
//! strategy to handle intervals efficiently, particularly when dealing with high-coverage regions
//! common in genomic data.
//!
//! ## Examples
//!
//! ### Finding all genes that overlap a query region
//!
//! ```rust
//! use gtars_overlaprs::{AiList, Overlapper, Interval};
//!
//! let genes = vec![
//!     Interval { start: 1000u32, end: 2000, val: "BRCA1" },
//!     Interval { start: 3000, end: 4000, val: "TP53" },
//!     Interval { start: 5000, end: 6000, val: "EGFR" },
//! ];
//!
//! let gene_index = AiList::build(genes);
//!
//! // query a specific region (e.g., chr17:1500-3500)
//! let overlapping_genes: Vec<&str> = gene_index
//!     .find_iter(1500, 3500)
//!     .map(|interval| interval.val)
//!     .collect();
//!
//! println!("Genes in region: {:?}", overlapping_genes);
//! ```
pub mod ailist;
pub mod bits;
pub mod traits;

// re-exports
pub use self::ailist::AiList;
pub use self::bits::Bits;
pub use self::traits::{Interval, Overlapper};
// pub use self::iitree::IITree;

pub mod consts {
    pub const OVERLAP_CMD: &str = "overlap";
}

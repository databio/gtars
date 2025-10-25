//! # gtars-overlaprs
//!
//! Core infrastructure for all genomic interval overlap computations.
//!
//! ## Purpose
//!
//! This module provides the foundational overlap computation primitives that all other
//! modules build upon. It is the single source of truth for efficient interval operations
//! in the gtars ecosystem.
//!
//! ## Design Philosophy
//!
//! All overlap computation logic should live here. Higher-level modules (scoring, tokenizers)
//! wrap this functionality for their specific use cases but should not reimplement overlap
//! algorithms.
//!
//! ## Main Components
//!
//! - **`AIList`**: Augmented Interval List for fast interval queries
//! - **`Bits`**: Bit vector implementation for coverage calculations
//! - **`Overlapper`**: Common trait providing a consistent interface for all overlap operations
//!
//! ## Example
//!
//! ```rust,ignore
//! use gtars_overlaprs::{Bits, Overlapper};
//!
//! // Create an interval tree
//! let intervals = vec![
//!     Interval { start: 100, end: 200, val: 1 },
//!     Interval { start: 150, end: 250, val: 2 },
//! ];
//! let bits = Bits::build(intervals);
//!
//! // Query overlaps
//! let overlaps = bits.find(120, 180);
//! ```

pub mod ailist;
pub mod bits;
pub mod traits;

// re-exports
pub use self::ailist::AiList;
pub use self::bits::Bits;
pub use self::traits::Overlapper;
// pub use self::iitree::IITree;

pub mod consts {
    pub const OVERLAP_CMD: &str = "overlap";
}

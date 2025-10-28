//! # gtars-scoring
//!
//! Wrapper around gtars-overlaprs for generating X-by-peak count matrices.
//!
//! ## Purpose
//!
//! This module wraps the core overlap infrastructure from gtars-overlaprs to produce
//! count matrices for genomic analysis. It handles all use cases that need X-by-peak
//! matrices where X can be cells, pseudobulk aggregations, or samples.
//!
//! ## Design Philosophy
//!
//! All overlap computation is delegated to gtars-overlaprs. This module focuses on:
//! - Matrix data structures and I/O
//! - Aggregation strategies
//! - Output formats (dense/sparse)
//!
//! ## Use Cases
//!
//! - **Cell-by-peak** matrices for single-cell ATAC-seq
//! - **Sample-by-peak** matrices for bulk ATAC-seq/ChIP-seq
//! - **Pseudobulk-by-peak** matrices for aggregated single-cell data
//!
//! ## Main Components
//!
//! - **`ConsensusSet`**: Wraps consensus peaks using `gtars-overlaprs::Bits` for overlap detection
//! - **`CountMatrix`**: Dense matrix storage for count data
//! - **`region_scoring_from_fragments`**: Main scoring function for file-based input
//!
//! ## CLI Usage
//!
//! ```bash
//! # Score multiple fragment files against consensus peaks
//! gtars fscoring "fragments/*.bed.gz" peaks.bed --output matrix.csv.gz
//! ```
//!
//! ## Example
//!
//! ```rust,ignore
//! use gtars_scoring::{region_scoring_from_fragments, ConsensusSet, FragmentFileGlob, ScoringMode};
//!
//! // Set up inputs
//! let mut fragments = FragmentFileGlob::new("fragments/*.bed.gz")?;
//! let consensus = ConsensusSet::new("peaks.bed".into())?;
//!
//! // Generate count matrix
//! let count_matrix = region_scoring_from_fragments(
//!     &mut fragments,
//!     &consensus,
//!     ScoringMode::Atac
//! )?;
//!
//! // Write to file
//! count_matrix.write_to_file("output.csv.gz")?;
//! ```

pub mod consts;
pub mod counts;
pub mod files;
pub mod fragment_scoring;
pub mod scoring_modes;

// re-exports
pub use counts::*;
pub use files::*;
pub use fragment_scoring::*;
pub use scoring_modes::*;

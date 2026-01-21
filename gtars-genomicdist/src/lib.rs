//! Genomic distribution and statistics for region sets.
//!
//! This crate provides tools for analyzing the distribution of genomic regions
//! across chromosomes, including:
//!
//! - Computing summary statistics (min, max, mean, median) for region lengths per chromosome
//! - Binning genomes into fixed-size windows and counting region overlaps
//! - Analyzing genomic coverage patterns
//!
//! # Example
//!
//! ```no_run
//! use gtars_genomicdist::GenomicIntervalSetStatistics;
//! use gtars_core::models::RegionSet;
//!
//! let regions = RegionSet::try_from("input.bed").unwrap();
//!
//! // Get statistics per chromosome
//! let stats = regions.chromosome_statistics();
//!
//! // Get region distribution across 10 bins
//! let distribution = regions.region_distribution();
//! ```

pub mod bed_classifier;
pub mod errors;
pub mod models;
pub mod statistics;
pub mod utils;

// re-exports
#[cfg(feature = "bedclassifier")]
pub use bed_classifier::classify_bed;
pub use statistics::GenomicIntervalSetStatistics;
pub use statistics::{calc_dinucl_freq, calc_gc_content};

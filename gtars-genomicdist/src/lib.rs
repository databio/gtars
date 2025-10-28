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
//! use gtars_gd::GenomicIntervalSetStatistics;
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

pub mod models;
pub mod statistics;
pub mod utils;

// re-exports
pub use statistics::GenomicIntervalSetStatistics;

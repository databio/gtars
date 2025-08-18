//! # Scatrs: Single-cell ATAC-seq Simulation Module
//! 
//! A high-performance module for simulating realistic single-cell ATAC-seq data 
//! from bulk ATAC-seq experiments, including biological artifacts like doublets 
//! and ambient DNA contamination.
//! 
//! ## Overview
//! 
//! Scatrs transforms bulk ATAC-seq data into simulated single-cell data that closely
//! mimics real scATAC-seq experiments. It's designed for:
//! - Benchmarking computational methods (doublet detection, clustering, etc.)
//! - Understanding scATAC-seq data characteristics
//! - Testing analysis pipelines with ground truth labels
//! 
//! ## Key Features
//! 
//! - **Weighted cell type sampling** - Control cell type proportions
//! - **Doublet simulation** - Create multiplets with realistic mixing ratios
//! - **Ambient contamination** - Add background DNA from lysed cells
//! - **Peak-based or fragment-only** - Works with or without peak calls
//! - **Pipeline agnostic** - Compatible with any peak caller
//! 
//! ## Example Usage
//! 
//! ```rust,ignore
//! use gtars::scatrs::{BamProcessor, WeightedSampler, UniformSampler, DoubletSimulator};
//! use gtars::scatrs::{FragmentDistribution, ScatrsConfig};
//! use std::path::Path;
//! 
//! // Process BAM file to extract fragments
//! let fragments = BamProcessor::filter_fragments(Path::new("sample.bam"))?;
//! 
//! // Create samplers
//! let mut weighted = WeightedSampler::new(Some(42)); // seed for reproducibility
//! let uniform = UniformSampler::new(Some(42));
//! 
//! // Sample cells (simplified - normally you'd have peaks and weights)
//! let cell_regions = weighted.sample_regions_for_cells(
//!     None, None, None, None, 0.8, 
//!     &FragmentDistribution::Uniform { count: 10000 },
//!     100 // number of cells
//! )?;
//! 
//! // Create doublets (10% rate)
//! let doublet_sim = DoubletSimulator::new(Some(42));
//! let mut cells = vec![]; // your simulated cells
//! let config = ScatrsConfig::default();
//! let annotations = doublet_sim.simulate_doublets(&mut cells, 0.1, &config)?;
//! ```
//! 
//! ## Module Structure
//! 
//! - [`models`] - Core data structures (Fragment, Region, Cell, Config)
//! - [`io`] - File I/O for various formats (BED, narrowPeak, fragments)
//! - [`staging`] - BAM processing and peak operations
//! - [`sampling`] - Weighted and uniform sampling algorithms
//! - [`cli`] - Command-line interface implementation

pub mod cli;
pub mod consts;
pub mod models;
pub mod staging;
pub mod sampling;
pub mod io;

#[cfg(test)]
mod tests;

// Re-export commonly used types
pub use models::{ScatrsFragment, ScatrsRegion, Peak, SimulatedCell, ScatrsConfig, FragmentDistribution, 
                 StageConfig, StageParameters, StageStats};
pub use staging::{BamProcessor, PeakMerger, BackgroundGenerator};
pub use sampling::{WeightedSampler, UniformSampler, SimulationStats, DoubletSimulator, AmbientSimulator};
pub use io::{BedWriter, FragmentWriter, NarrowPeakReader, ManifestReader};
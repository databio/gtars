use crate::scatrs::models::{ScatrsRegion, ScatrsFragment, SimulatedCell, FragmentDistribution, ScatrsConfig};
use rand::prelude::*;
use rand::distributions::{Distribution, Uniform};
use rand_distr::{Normal, Beta};
use anyhow::{Result, Context};
use rayon::prelude::*;
use indicatif::{ProgressBar, ProgressStyle, MultiProgress};
use std::sync::Arc;
use std::collections::HashMap;

pub mod wrswr_skip;
use wrswr_skip::OptimizedWRSWRSkip;

// ============================================================================
// Weighted Sampling
// ============================================================================

/// Weighted sampling of genomic regions based on fragment counts
/// 
/// Uses weighted reservoir sampling to select regions proportional to their
/// fragment density, ensuring realistic peak enrichment in simulated cells.
pub struct WeightedSampler {}

#[derive(Debug, Clone)]
pub struct SampledRegion {
    pub region: ScatrsRegion,
    pub is_peak: bool,
    pub weight: f64,
    pub fragment_count: u32,
}

impl WeightedSampler {
    pub fn new(_seed: Option<u64>) -> Self {
        Self {}
    }
    
    /// Sample genomic regions for multiple cells
    /// 
    /// Assigns regions to cells based on weights and signal-to-noise ratio.
    /// Returns a vector of region assignments, one per cell.
    /// 
    /// # Arguments
    /// * `peaks` - Optional peak regions to sample from
    /// * `peak_weights` - Weights for peak sampling (e.g., from fragment counts)
    /// * `background` - Optional background regions
    /// * `background_weights` - Weights for background sampling
    /// * `signal_to_noise` - Ratio of peak to background fragments (0-1)
    /// * `fragment_dist` - Distribution for fragment counts per cell
    /// * `num_cells` - Number of cells to generate
    pub fn sample_regions_for_cells(
        &mut self,
        peaks: Option<&[ScatrsRegion]>,
        peak_weights: Option<&[f64]>,
        background: Option<&[ScatrsRegion]>,
        background_weights: Option<&[f64]>,
        signal_to_noise: f64,
        fragment_dist: &FragmentDistribution,
        num_cells: usize,
    ) -> Result<Vec<Vec<SampledRegion>>> {
        let multi_progress = Arc::new(MultiProgress::new());
        let main_pb = multi_progress.add(ProgressBar::new(num_cells as u64));
        main_pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} cells sampled")
                .unwrap()
        );
        
        // Check if we have peaks to sample from
        if peaks.is_none() || background.is_none() {
            // No peaks mode - return empty regions for uniform fragment sampling
            println!("No peaks provided - fragments will be sampled uniformly");
            let empty_cells = vec![Vec::new(); num_cells];
            main_pb.finish_with_message("Fragment-only mode (no peak regions)");
            return Ok(empty_cells);
        }
        
        let peaks = peaks.unwrap();
        let peak_weights = peak_weights.unwrap();
        let background = background.unwrap();
        let background_weights = background_weights.unwrap();
        
        // Pre-normalize weights
        let peak_sum: f64 = peak_weights.iter().sum();
        let bg_sum: f64 = background_weights.iter().sum();
        
        let norm_peak_weights: Vec<f64> = if peak_sum > 0.0 {
            peak_weights.iter().map(|w| w / peak_sum).collect()
        } else {
            vec![1.0 / peaks.len() as f64; peaks.len()]
        };
        
        let norm_bg_weights: Vec<f64> = if bg_sum > 0.0 {
            background_weights.iter().map(|w| w / bg_sum).collect()
        } else {
            vec![1.0 / background.len() as f64; background.len()]
        };
        
        // Generate cells in parallel batches
        let batch_size = 10;
        let mut all_cells = Vec::with_capacity(num_cells);
        
        for batch_start in (0..num_cells).step_by(batch_size) {
            let batch_end = (batch_start + batch_size).min(num_cells);
            let batch_results: Vec<Vec<SampledRegion>> = (batch_start..batch_end)
                .into_par_iter()
                .map(|_| {
                    let mut local_rng = StdRng::from_entropy();
                    self.sample_regions_for_cell_optimized(
                        peaks,
                        &norm_peak_weights,
                        background,
                        &norm_bg_weights,
                        signal_to_noise,
                        fragment_dist,
                        &mut local_rng,
                    ).unwrap_or_else(|e| {
                        eprintln!("Warning: Failed to sample cell: {}", e);
                        Vec::new()
                    })
                })
                .collect();
            
            all_cells.extend(batch_results);
            main_pb.inc((batch_end - batch_start) as u64);
        }
        
        main_pb.finish_with_message("Cell sampling complete");
        
        Ok(all_cells)
    }
    
    /// Sample genomic regions for a single cell using WRSWR-SKIP algorithm
    /// 
    /// This optimized version uses the WRSWR-SKIP algorithm for O(n + k log k)
    /// time complexity and without-replacement sampling for biological accuracy.
    pub fn sample_regions_for_cell_optimized(
        &self,
        peaks: &[ScatrsRegion],
        peak_weights: &[f64],
        background: &[ScatrsRegion],
        background_weights: &[f64],
        signal_to_noise: f64,
        fragment_dist: &FragmentDistribution,
        rng: &mut StdRng,
    ) -> Result<Vec<SampledRegion>> {
        // Determine number of fragments for this cell
        let total_fragments = match fragment_dist {
            FragmentDistribution::Uniform { count } => *count,
            FragmentDistribution::Gaussian { mean, variance, min, max } => {
                let normal = Normal::new(*mean, variance.sqrt())
                    .context("Invalid Gaussian parameters")?;
                let sampled = normal.sample(rng) as u32;
                sampled.clamp(*min, *max)
            }
        };
        
        // Split fragments between peaks and background based on signal-to-noise ratio
        let peak_fragments = (total_fragments as f64 * signal_to_noise) as usize;
        let background_fragments = (total_fragments - peak_fragments as u32) as usize;
        
        let mut sampled_regions = Vec::new();
        
        // Sample peaks using WRSWR-SKIP
        if peak_fragments > 0 && !peaks.is_empty() {
            let mut peak_sampler = OptimizedWRSWRSkip::new(
                peak_fragments.min(peaks.len()),
                Some(rng.gen())
            );
            
            let peak_iter = peaks.iter().zip(peak_weights.iter())
                .map(|(region, &weight)| (region.clone(), weight));
            
            let sampled_peaks = peak_sampler.sample_stream(
                peak_iter,
                |(_, weight)| *weight
            );
            
            // Convert to SampledRegion format
            for (region, weight) in sampled_peaks {
                sampled_regions.push(SampledRegion {
                    region,
                    is_peak: true,
                    weight,
                    fragment_count: 1, // Will be adjusted based on actual sampling
                });
            }
        }
        
        // Sample background using WRSWR-SKIP
        if background_fragments > 0 && !background.is_empty() {
            let mut bg_sampler = OptimizedWRSWRSkip::new(
                background_fragments.min(background.len()),
                Some(rng.gen())
            );
            
            let bg_iter = background.iter().zip(background_weights.iter())
                .map(|(region, &weight)| (region.clone(), weight));
            
            let sampled_bg = bg_sampler.sample_stream(
                bg_iter,
                |(_, weight)| *weight
            );
            
            for (region, weight) in sampled_bg {
                sampled_regions.push(SampledRegion {
                    region,
                    is_peak: false,
                    weight,
                    fragment_count: 1,
                });
            }
        }
        
        Ok(sampled_regions)
    }
    
}

// ============================================================================
// Uniform Sampling
// ============================================================================

pub struct UniformSampler {
    seed: Option<u64>,
}

impl UniformSampler {
    pub fn new(seed: Option<u64>) -> Self {
        Self { seed }
    }
    
    pub fn sample_fragments_for_cells(
        &self,
        cell_regions: Vec<Vec<SampledRegion>>,
        region_fragments: &HashMap<String, Vec<ScatrsFragment>>,
        cell_types: &[String],
        cell_type_weights: Option<&[f64]>,         // Weights for cell type proportions
        all_fragments: Option<&[ScatrsFragment]>,  // For no-peaks mode
        fragments_per_cell: Option<u32>,           // For no-peaks mode
    ) -> Result<Vec<SimulatedCell>> {
        let num_cells = cell_regions.len();
        let pb = ProgressBar::new(num_cells as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} cells processed")
                .unwrap()
        );
        pb.set_message("Sampling fragments for cells");
        
        // Pre-compute cell type assignments based on weights
        let cell_type_assignments = if let Some(weights) = cell_type_weights {
            self.assign_cell_types_weighted(num_cells, cell_types, weights)?
        } else {
            // Default: round-robin assignment
            (0..num_cells)
                .map(|i| cell_types[i % cell_types.len()].clone())
                .collect()
        };
        
        // Process cells in parallel
        let cells: Vec<SimulatedCell> = cell_regions
            .into_par_iter()
            .enumerate()
            .map(|(cell_idx, regions)| {
                let mut rng = match self.seed {
                    Some(s) => StdRng::seed_from_u64(s + cell_idx as u64),
                    None => StdRng::from_entropy(),
                };
                
                // Use pre-assigned cell type
                let cell_type = &cell_type_assignments[cell_idx];
                let cell_id = format!("cell_{:06}", cell_idx + 1);
                
                let mut cell = SimulatedCell::new(cell_id, cell_type.clone());
                
                // Check if this is no-peaks mode (empty regions)
                if regions.is_empty() && all_fragments.is_some() {
                    // No-peaks mode: sample uniformly from all fragments
                    let fragments = all_fragments.unwrap();
                    let num_to_sample = fragments_per_cell.unwrap_or(10000) as usize;
                    
                    if !fragments.is_empty() {
                        let sampled = self.uniform_sample_fragments(
                            fragments,
                            num_to_sample,
                            &cell.cell_id,
                            &mut rng,
                        );
                        
                        for frag in sampled {
                            cell.add_fragment(frag);
                        }
                    }
                } else {
                    // Normal mode: sample from regions
                    for sampled_region in regions {
                    let region_key = format!(
                        "{}:{}-{}",
                        sampled_region.region.chrom,
                        sampled_region.region.start,
                        sampled_region.region.end
                    );
                    
                    // Get fragments that overlap this region
                    let available_fragments = region_fragments.get(&region_key)
                        .or_else(|| {
                            // Try alternative key format
                            region_fragments.get(&format!(
                                "{}_{}_{}", 
                                sampled_region.region.chrom,
                                sampled_region.region.start,
                                sampled_region.region.end
                            ))
                        });
                    
                    if let Some(fragments) = available_fragments {
                        if !fragments.is_empty() {
                            // Uniformly sample fragments from this region
                            let sampled = self.uniform_sample_fragments(
                                fragments,
                                sampled_region.fragment_count as usize,
                                &cell.cell_id,
                                &mut rng,
                            );
                            
                            for frag in sampled {
                                if sampled_region.is_peak {
                                    cell.add_peak_fragment(frag);
                                } else {
                                    cell.add_background_fragment(frag);
                                }
                            }
                        }
                    } else {
                        // If no real fragments available, generate synthetic ones
                        let synthetic = self.generate_synthetic_fragments(
                            &sampled_region.region,
                            sampled_region.fragment_count as usize,
                            &cell.cell_id,
                            &mut rng,
                        );
                        
                        for frag in synthetic {
                            if sampled_region.is_peak {
                                cell.add_peak_fragment(frag);
                            } else {
                                cell.add_background_fragment(frag);
                            }
                        }
                    }
                    }  // Close the else block for normal mode
                }
                
                // Calculate FRiP score
                cell.calculate_frip();
                
                pb.inc(1);
                cell
            })
            .collect();
        
        pb.finish_with_message(format!("Generated {} cells", cells.len()));
        
        Ok(cells)
    }
    
    fn assign_cell_types_weighted(
        &self,
        num_cells: usize,
        cell_types: &[String],
        weights: &[f64],
    ) -> Result<Vec<String>> {
        if cell_types.len() != weights.len() {
            anyhow::bail!("Cell types and weights must have the same length");
        }
        
        // Normalize weights
        let total_weight: f64 = weights.iter().sum();
        let normalized_weights: Vec<f64> = weights.iter()
            .map(|w| w / total_weight)
            .collect();
        
        // Calculate number of cells per type
        let mut cell_counts: Vec<usize> = normalized_weights.iter()
            .map(|w| (w * num_cells as f64).round() as usize)
            .collect();
        
        // Adjust for rounding errors
        let total_assigned: usize = cell_counts.iter().sum();
        if total_assigned < num_cells {
            // Add remaining cells to the type with highest weight
            let max_idx = weights.iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                .map(|(i, _)| i)
                .unwrap_or(0);
            cell_counts[max_idx] += num_cells - total_assigned;
        } else if total_assigned > num_cells {
            // Remove excess cells from the type with most cells
            let max_idx = cell_counts.iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.cmp(b))
                .map(|(i, _)| i)
                .unwrap_or(0);
            cell_counts[max_idx] -= total_assigned - num_cells;
        }
        
        // Create cell type assignments
        let mut assignments = Vec::with_capacity(num_cells);
        for (cell_type, count) in cell_types.iter().zip(cell_counts.iter()) {
            for _ in 0..*count {
                assignments.push(cell_type.clone());
            }
        }
        
        // Shuffle to randomize positions
        let mut rng = match self.seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };
        assignments.shuffle(&mut rng);
        
        // Print distribution info
        println!("Cell type distribution:");
        for (cell_type, count) in cell_types.iter().zip(cell_counts.iter()) {
            let percentage = (*count as f64 / num_cells as f64) * 100.0;
            println!("  {}: {} cells ({:.1}%)", cell_type, count, percentage);
        }
        
        Ok(assignments)
    }
    
    fn uniform_sample_fragments(
        &self,
        available: &[ScatrsFragment],
        count: usize,
        cell_id: &str,
        rng: &mut StdRng,
    ) -> Vec<ScatrsFragment> {
        if available.is_empty() || count == 0 {
            return Vec::new();
        }
        
        let mut sampled = Vec::with_capacity(count);
        let uniform_dist = Uniform::new(0, available.len());
        
        // Sample with replacement (diploid genome)
        for _ in 0..count {
            let idx = uniform_dist.sample(rng);
            let mut fragment = available[idx].clone();
            fragment.cell_id = cell_id.to_string();
            sampled.push(fragment);
        }
        
        sampled
    }
    
    fn generate_synthetic_fragments(
        &self,
        region: &ScatrsRegion,
        count: usize,
        cell_id: &str,
        rng: &mut StdRng,
    ) -> Vec<ScatrsFragment> {
        let mut fragments = Vec::with_capacity(count);
        
        // Typical ATAC-seq fragment size distribution
        let size_dist = Uniform::new(50, 500);
        let pos_dist = Uniform::new(
            region.start,
            region.end.saturating_sub(500),
        );
        
        for _ in 0..count {
            let start = pos_dist.sample(rng);
            let size = size_dist.sample(rng) as u64;
            let end = (start + size).min(region.end);
            
            let fragment = ScatrsFragment::new(
                region.chrom.clone(),
                start,
                end,
            )
            .with_cell_id(cell_id.to_string())
            .with_quality(30); // Default quality score
            
            fragments.push(fragment);
        }
        
        fragments
    }
    
    pub fn build_fragment_index(
        fragments: &[ScatrsFragment],
        regions: &[ScatrsRegion],
    ) -> HashMap<String, Vec<ScatrsFragment>> {
        let mut index = HashMap::new();
        
        for region in regions {
            let region_key = format!(
                "{}:{}-{}",
                region.chrom,
                region.start,
                region.end
            );
            
            let overlapping: Vec<ScatrsFragment> = fragments
                .iter()
                .filter(|f| {
                    f.chrom == region.chrom
                        && f.start < region.end
                        && f.end > region.start
                })
                .cloned()
                .collect();
            
            if !overlapping.is_empty() {
                index.insert(region_key, overlapping);
            }
        }
        
        index
    }
}

// ============================================================================
// Statistics
// ============================================================================

pub struct SimulationStats {
    pub total_cells: usize,
    pub total_fragments: u64,
    pub mean_fragments_per_cell: f64,
    pub median_fragments_per_cell: f64,
    pub mean_frip: f64,
    pub median_frip: f64,
    pub min_fragments: u32,
    pub max_fragments: u32,
}

impl SimulationStats {
    pub fn from_cells(cells: &[SimulatedCell]) -> Self {
        if cells.is_empty() {
            return Self {
                total_cells: 0,
                total_fragments: 0,
                mean_fragments_per_cell: 0.0,
                median_fragments_per_cell: 0.0,
                mean_frip: 0.0,
                median_frip: 0.0,
                min_fragments: 0,
                max_fragments: 0,
            };
        }
        
        let total_cells = cells.len();
        
        // Fragment counts
        let mut fragment_counts: Vec<u32> = cells
            .iter()
            .map(|c| c.metadata.fragment_count)
            .collect();
        fragment_counts.sort_unstable();
        
        let total_fragments: u64 = fragment_counts.iter().map(|&c| c as u64).sum();
        let mean_fragments = total_fragments as f64 / total_cells as f64;
        let median_fragments = if total_cells % 2 == 0 {
            let mid = total_cells / 2;
            (fragment_counts[mid - 1] + fragment_counts[mid]) as f64 / 2.0
        } else {
            fragment_counts[total_cells / 2] as f64
        };
        
        // FRiP scores
        let mut frip_scores: Vec<f64> = cells
            .iter()
            .filter_map(|c| c.metadata.frip_score)
            .collect();
        frip_scores.sort_by(|a, b| a.partial_cmp(b).unwrap());
        
        let mean_frip = if !frip_scores.is_empty() {
            frip_scores.iter().sum::<f64>() / frip_scores.len() as f64
        } else {
            0.0
        };
        
        let median_frip = if !frip_scores.is_empty() {
            if frip_scores.len() % 2 == 0 {
                let mid = frip_scores.len() / 2;
                (frip_scores[mid - 1] + frip_scores[mid]) / 2.0
            } else {
                frip_scores[frip_scores.len() / 2]
            }
        } else {
            0.0
        };
        
        Self {
            total_cells,
            total_fragments,
            mean_fragments_per_cell: mean_fragments,
            median_fragments_per_cell: median_fragments,
            mean_frip,
            median_frip,
            min_fragments: *fragment_counts.first().unwrap_or(&0),
            max_fragments: *fragment_counts.last().unwrap_or(&0),
        }
    }
    
    pub fn print_summary(&self) {
        println!("=== Simulation Statistics ===");
        println!("Total cells: {}", self.total_cells);
        println!("Total fragments: {}", self.total_fragments);
        println!("Fragments per cell:");
        println!("  Mean: {:.2}", self.mean_fragments_per_cell);
        println!("  Median: {:.2}", self.median_fragments_per_cell);
        println!("  Min: {}", self.min_fragments);
        println!("  Max: {}", self.max_fragments);
        println!("FRiP scores:");
        println!("  Mean: {:.4}", self.mean_frip);
        println!("  Median: {:.4}", self.median_frip);
    }
}

pub fn calculate_weight_from_counts(counts: &[u32]) -> Vec<f64> {
    let total: u32 = counts.iter().sum();
    if total == 0 {
        vec![1.0 / counts.len() as f64; counts.len()]
    } else {
        counts.iter().map(|&c| c as f64 / total as f64).collect()
    }
}

pub fn normalize_weights(weights: &[f64]) -> Vec<f64> {
    let sum: f64 = weights.iter().sum();
    if sum > 0.0 {
        weights.iter().map(|w| w / sum).collect()
    } else {
        vec![1.0 / weights.len() as f64; weights.len()]
    }
}

// ============================================================================
// Doublet Simulation
// ============================================================================

/// Simulates doublet cells by merging fragments from multiple cells
/// 
/// Creates realistic multiplets with weighted mixing ratios using a Beta
/// distribution, mimicking the biological process where two cells are
/// captured in the same droplet.
/// 
/// # Example
/// ```no_run
/// # use anyhow::Result;
/// # fn main() -> Result<()> {
/// use gtars::scatrs::{DoubletSimulator, ScatrsConfig};
/// 
/// let simulator = DoubletSimulator::new(Some(42));
/// let mut cells = vec![]; // your cells
/// let config = ScatrsConfig::default();
/// 
/// // Create 10% doublets
/// let annotations = simulator.simulate_doublets(&mut cells, 0.1, &config)?;
/// # Ok(())
/// # }
/// ```
pub struct DoubletSimulator {
    seed: Option<u64>,
}

impl DoubletSimulator {
    pub fn new(seed: Option<u64>) -> Self {
        Self { seed }
    }
    
    pub fn simulate_doublets(
        &self,
        cells: &mut Vec<SimulatedCell>,
        doublet_rate: f64,
        _config: &ScatrsConfig,
    ) -> Result<Vec<(String, Vec<String>)>> {  // Returns doublet annotations
        if doublet_rate <= 0.0 || cells.is_empty() {
            return Ok(Vec::new());
        }
        
        let mut rng = match self.seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };
        
        let num_doublets = (cells.len() as f64 * doublet_rate) as usize;
        let mut doublet_annotations = Vec::new();
        
        println!("Creating {} doublets ({:.1}% of cells)", num_doublets, doublet_rate * 100.0);
        
        // Beta distribution for weighted doublet mixing (alpha=2, beta=2 gives a bell curve centered at 0.5)
        let beta_dist = Beta::new(2.0, 2.0).context("Invalid beta distribution parameters")?;
        
        // Randomly select cells to convert to doublets
        let mut indices: Vec<usize> = (0..cells.len()).collect();
        indices.shuffle(&mut rng);
        
        for i in 0..num_doublets.min(cells.len() / 2) {
            let idx1 = indices[i * 2];
            let idx2 = indices[i * 2 + 1];
            
            // Determine mixing proportion (what fraction comes from cell 1 vs cell 2)
            let mix_ratio = beta_dist.sample(&mut rng);
            
            // Create doublet by merging fragments
            let cell1_type = cells[idx1].cell_type.clone();
            let cell2_type = cells[idx2].cell_type.clone();
            
            // Sample fragments from both cells based on mixing ratio
            let total_fragments = cells[idx1].fragments.len() + cells[idx2].fragments.len();
            let from_cell1 = (total_fragments as f64 * mix_ratio) as usize;
            let from_cell2 = total_fragments - from_cell1;
            
            let mut combined_fragments = Vec::new();
            
            // Sample from cell 1
            if !cells[idx1].fragments.is_empty() && from_cell1 > 0 {
                let sampled_indices = self.sample_indices(cells[idx1].fragments.len(), from_cell1, &mut rng);
                for idx in sampled_indices {
                    combined_fragments.push(cells[idx1].fragments[idx].clone());
                }
            }
            
            // Sample from cell 2
            if !cells[idx2].fragments.is_empty() && from_cell2 > 0 {
                let sampled_indices = self.sample_indices(cells[idx2].fragments.len(), from_cell2, &mut rng);
                for idx in sampled_indices {
                    combined_fragments.push(cells[idx2].fragments[idx].clone());
                }
            }
            
            // Update the first cell to be the doublet
            cells[idx1].fragments = combined_fragments;
            cells[idx1].cell_type = format!("{}__{}", cell1_type, cell2_type);
            cells[idx1].metadata.is_doublet = true;
            cells[idx1].metadata.doublet_types = Some(vec![
                (cell1_type.clone(), mix_ratio),
                (cell2_type.clone(), 1.0 - mix_ratio),
            ]);
            cells[idx1].metadata.fragment_count = cells[idx1].fragments.len() as u32;
            
            // Clear the second cell (it was merged into the doublet)
            cells[idx2].fragments.clear();
            cells[idx2].metadata.fragment_count = 0;
            
            // Record annotation
            doublet_annotations.push((
                cells[idx1].cell_id.clone(),
                vec![cell1_type, cell2_type],
            ));
        }
        
        // Remove empty cells (those that were merged into doublets)
        cells.retain(|c| c.metadata.fragment_count > 0);
        
        Ok(doublet_annotations)
    }
    
    fn sample_indices(&self, total: usize, count: usize, rng: &mut StdRng) -> Vec<usize> {
        if count >= total {
            return (0..total).collect();
        }
        
        let mut indices: Vec<usize> = (0..total).collect();
        indices.shuffle(rng);
        indices.truncate(count);
        indices
    }
}

// ============================================================================
// Ambient DNA Contamination
// ============================================================================

/// Adds ambient DNA contamination to cells
/// 
/// Simulates background DNA from lysed cells that contaminates droplets,
/// a common artifact in single-cell experiments. Fragments are sampled
/// uniformly from the total pool.
pub struct AmbientSimulator {
    seed: Option<u64>,
}

impl AmbientSimulator {
    pub fn new(seed: Option<u64>) -> Self {
        Self { seed }
    }
    
    pub fn add_ambient_contamination(
        &self,
        cells: &mut Vec<SimulatedCell>,
        contamination_rate: f64,
        all_fragments: &[ScatrsFragment],  // Pool of all fragments from staging
    ) -> Result<()> {
        if contamination_rate <= 0.0 || all_fragments.is_empty() {
            return Ok(());
        }
        
        println!("Adding {:.1}% ambient DNA contamination", contamination_rate * 100.0);
        
        let pb = ProgressBar::new(cells.len() as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} cells")
                .unwrap()
        );
        
        // Process each cell
        for (idx, cell) in cells.iter_mut().enumerate() {
            let mut rng = match self.seed {
                Some(s) => StdRng::seed_from_u64(s + idx as u64),
                None => StdRng::from_entropy(),
            };
            
            let original_count = cell.fragments.len();
            if original_count == 0 {
                continue;
            }
            
            // Calculate number of ambient fragments to add
            let ambient_count = (original_count as f64 * contamination_rate / (1.0 - contamination_rate)) as usize;
            
            if ambient_count > 0 {
                // Randomly sample from the ambient pool
                let uniform = Uniform::new(0, all_fragments.len());
                for _ in 0..ambient_count {
                    let idx = uniform.sample(&mut rng);
                    let mut ambient_frag = all_fragments[idx].clone();
                    ambient_frag.cell_id = cell.cell_id.clone();
                    cell.fragments.push(ambient_frag);
                    cell.metadata.ambient_fragments += 1;
                }
                
                cell.metadata.fragment_count = cell.fragments.len() as u32;
            }
            
            pb.inc(1);
        }
        
        pb.finish_with_message("Ambient contamination added");
        Ok(())
    }
}
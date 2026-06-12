use crate::models::{ScatrsRegion, ScatrsFragment, SimulatedCell, FragmentDistribution, SignalToNoiseDistribution, ScatrsConfig};
use rand::prelude::*;
use rand::distributions::{Distribution, Uniform};
use rand_distr::{Normal, Beta, Gamma, Poisson};
use anyhow::{Result, Context};
use rayon::prelude::*;
use indicatif::{ProgressBar, ProgressStyle, MultiProgress};
use std::sync::Arc;
use std::collections::HashMap;

pub mod a_res;
pub mod a_expj;
use a_expj::AExpJSampler as VectorAResSampler;  // Alias for compatibility

#[cfg(test)]
mod tests;

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

    /// Sample fragment count from a distribution
    fn sample_fragment_count(
        fragment_dist: &FragmentDistribution,
        rng: &mut StdRng,
    ) -> Result<u32> {
        match fragment_dist {
            FragmentDistribution::Fixed { fragments_per_cell } => {
                Ok(*fragments_per_cell as u32)
            },
            FragmentDistribution::Uniform { min, max } => {
                let uniform_dist = Uniform::new_inclusive(*min, *max);
                Ok(uniform_dist.sample(rng) as u32)
            },
            FragmentDistribution::Gaussian { mean_fragments_per_cell, sd, min, max } => {
                let normal = Normal::new(*mean_fragments_per_cell, *sd)
                    .context("Invalid Gaussian parameters")?;
                let sampled = normal.sample(rng) as u32;
                let min_val = min.map(|v| v as u32).unwrap_or(0);
                match max {
                    Some(max_val) => Ok(sampled.clamp(min_val, *max_val as u32)),
                    None => Ok(sampled.max(min_val))
                }
            },
            FragmentDistribution::NegativeBinomial { mean, sd, min, max } => {
                let variance = sd * sd;
                if variance <= *mean {
                    return Err(anyhow::anyhow!(
                        "Invalid NegativeBinomial parameters: variance (sd²={}) must be greater than mean ({}) for overdispersion",
                        variance, mean
                    ));
                }

                let r = (mean * mean) / (variance - mean);
                let p = mean / variance;
                let gamma_shape = r;
                let gamma_scale = (1.0 - p) / p;

                let gamma = Gamma::new(gamma_shape as f64, gamma_scale)
                    .context("Invalid Gamma parameters for NegativeBinomial simulation")?;
                let lambda = gamma.sample(rng);

                let poisson = Poisson::new(lambda)
                    .context("Invalid Poisson parameter")?;
                let sampled = poisson.sample(rng) as u32;

                let min_val = min.map(|v| v as u32).unwrap_or(0);
                match max {
                    Some(max_val) => Ok(sampled.clamp(min_val, *max_val as u32)),
                    None => Ok(sampled.max(min_val))
                }
            }
        }
    }

    /// Sample genomic regions for multiple cells with cell-type-specific weights
    /// 
    /// Assigns regions to cells based on per-sample weights and signal-to-noise ratio.
    /// Returns a vector of region assignments, one per cell.
    /// 
    /// # Arguments
    /// * `peaks` - Optional peak regions to sample from
    /// * `sample_peak_weights` - Per-sample weights for peak sampling
    /// * `background` - Optional background regions
    /// * `sample_bg_weights` - Per-sample weights for background sampling
    /// * `signal_to_noise` - Ratio of peak to background fragments (0-1)
    /// * `fragment_dist` - Distribution for fragment counts per cell
    /// * `cell_type_assignments` - Pre-assigned cell types for each cell
    pub fn sample_regions_for_cells_by_type(
        &mut self,
        peaks: Option<&[ScatrsRegion]>,
        sample_peak_weights: Option<&HashMap<String, Vec<f64>>>,
        background: Option<&[ScatrsRegion]>,
        sample_bg_weights: Option<&HashMap<String, Vec<f64>>>,
        signal_to_noise_dist: &SignalToNoiseDistribution,
        fragment_dist: &FragmentDistribution,
        cell_type_assignments: &[String],
        ambient_contamination: f64,
    ) -> Result<Vec<Vec<SampledRegion>>> {
        let num_cells = cell_type_assignments.len();
        let multi_progress = Arc::new(MultiProgress::new());
        let main_pb = multi_progress.add(ProgressBar::new(num_cells as u64));
        main_pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} cells sampled")
                .unwrap()
        );
        
        // Check if we have peaks to sample from
        if peaks.is_none() || sample_peak_weights.is_none() {
            // No peaks mode - return empty regions for uniform fragment sampling
            println!("No peaks provided - fragments will be sampled uniformly");
            let empty_cells = vec![Vec::new(); num_cells];
            main_pb.finish_with_message("Fragment-only mode (no peak regions)");
            return Ok(empty_cells);
        }
        
        let peaks = peaks.unwrap();
        let sample_peak_weights = sample_peak_weights.unwrap();
        
        // Background regions are optional - many datasets may have no background
        let has_background = background.is_some() && sample_bg_weights.is_some() && 
                           background.map_or(false, |b| !b.is_empty());
        
        // Pre-normalize weights for each sample type
        let mut normalized_peak_weights: HashMap<String, Vec<f64>> = HashMap::new();
        let mut normalized_bg_weights: HashMap<String, Vec<f64>> = HashMap::new();
        
        for (sample_name, weights) in sample_peak_weights {
            let sum: f64 = weights.iter().sum();
            let normalized = if sum > 0.0 {
                weights.iter().map(|w| w / sum).collect()
            } else {
                vec![1.0 / peaks.len() as f64; peaks.len()]
            };
            normalized_peak_weights.insert(sample_name.clone(), normalized);
        }
        
        // Only process background weights if we have background regions
        if has_background {
            let background = background.unwrap();
            let sample_bg_weights = sample_bg_weights.unwrap();
            for (sample_name, weights) in sample_bg_weights {
                let sum: f64 = weights.iter().sum();
                let normalized = if sum > 0.0 {
                    weights.iter().map(|w| w / sum).collect()
                } else {
                    vec![1.0 / background.len() as f64; background.len()]
                };
                normalized_bg_weights.insert(sample_name.clone(), normalized);
            }
        }
        
        // Generate cells in parallel batches
        let batch_size = 10;
        let mut all_cells = Vec::with_capacity(num_cells);
        
        for batch_start in (0..num_cells).step_by(batch_size) {
            let batch_end = (batch_start + batch_size).min(num_cells);
            let batch_results: Vec<Vec<SampledRegion>> = (batch_start..batch_end)
                .into_par_iter()
                .map(|cell_idx| {
                    let cell_type = &cell_type_assignments[cell_idx];

                    // Get weights for this cell's type
                    let peak_weights = normalized_peak_weights.get(cell_type)
                        .ok_or_else(|| anyhow::anyhow!("No peak weights for cell type: {}", cell_type))
                        .unwrap();

                    let mut local_rng = StdRng::from_entropy();

                    // Sample signal-to-noise for this specific cell
                    let cell_signal_to_noise = sample_signal_to_noise(signal_to_noise_dist, &mut local_rng)
                        .unwrap_or(0.8);  // Fallback to default if sampling fails

                    if has_background {
                        let bg_weights = normalized_bg_weights.get(cell_type)
                            .ok_or_else(|| anyhow::anyhow!("No background weights for cell type: {}", cell_type))
                            .unwrap();
                        self.sample_regions_for_cell_optimized(
                            peaks,
                            peak_weights,
                            background.unwrap(),
                            bg_weights,
                            cell_signal_to_noise,
                            fragment_dist,
                            &mut local_rng,
                            ambient_contamination,
                        ).unwrap_or_else(|e| {
                            eprintln!("Warning: Failed to sample cell: {}", e);
                            Vec::new()
                        })
                    } else {
                        // No background regions - all fragments come from peaks
                        self.sample_regions_for_cell_peaks_only(
                            peaks,
                            peak_weights,
                            fragment_dist,
                            &mut local_rng,
                            ambient_contamination,
                        ).unwrap_or_else(|e| {
                            eprintln!("Warning: Failed to sample cell: {}", e);
                            Vec::new()
                        })
                    }
                })
                .collect();
            
            all_cells.extend(batch_results);
            main_pb.inc((batch_end - batch_start) as u64);
        }
        
        main_pb.finish_with_message("Cell sampling complete");
        
        Ok(all_cells)
    }
    
    /// Sample genomic regions for a single cell using A-Res algorithm
    /// 
    /// This version uses the basic Algorithm A for O(n + k log k)
    /// time complexity and without-replacement sampling for biological accuracy.
    /// This is the fundamental reservoir algorithm without skip optimizations.
    pub fn sample_regions_for_cell_optimized(
        &self,
        peaks: &[ScatrsRegion],
        peak_weights: &[f64],
        background: &[ScatrsRegion],
        background_weights: &[f64],
        signal_to_noise: f64,
        fragment_dist: &FragmentDistribution,
        rng: &mut StdRng,
        ambient_contamination: f64,
    ) -> Result<Vec<SampledRegion>> {
        // Determine total number of fragments for this cell
        let total_fragments = Self::sample_fragment_count(fragment_dist, rng)?;

        // Calculate source fragments (excluding ambient contamination)
        // This ensures the total fragment count stays as configured
        let source_fragments = ((total_fragments as f64) * (1.0 - ambient_contamination)) as u32;
        
        // Split source fragments between peaks and background based on signal-to-noise ratio
        let peak_fragments_uncapped = (source_fragments as f64 * signal_to_noise) as usize;
        let background_fragments_uncapped = (source_fragments - peak_fragments_uncapped as u32) as usize;

        // Cap at available regions (biological constraint: max 1 fragment per region per cell)
        let peak_fragments = peak_fragments_uncapped.min(peaks.len());
        let background_fragments = background_fragments_uncapped.min(background.len());

        // Always show fragment/region statistics (useful for understanding simulation parameters)
        if std::env::var("SCATRS_DEBUG").is_ok() {
            eprintln!("Cell sampling stats: {} peak fragments requested, {} peaks available (ratio: {:.2}%)",
                peak_fragments_uncapped, peaks.len(),
                (peak_fragments_uncapped as f64 / peaks.len() as f64) * 100.0);
            if !background.is_empty() {
                eprintln!("                    {} background fragments requested, {} background regions available (ratio: {:.2}%)",
                    background_fragments_uncapped, background.len(),
                    (background_fragments_uncapped as f64 / background.len() as f64) * 100.0);
            }
        }

        // Warn if we had to cap (this indicates unrealistic simulation parameters)
        if peak_fragments < peak_fragments_uncapped {
            eprintln!("WARNING: Requested {} peak fragments but only {} peaks available. Capped at {} fragments (ratio: {:.2}%).",
                peak_fragments_uncapped, peaks.len(), peak_fragments,
                (peak_fragments as f64 / peaks.len() as f64) * 100.0);
            eprintln!("         This is biologically unrealistic (hitting same peak multiple times in one cell).");
            eprintln!("         Consider reducing fragments_per_cell or increasing the peak set size.");
        }
        if background_fragments < background_fragments_uncapped {
            eprintln!("WARNING: Requested {} background fragments but only {} background regions available. Capped at {} fragments (ratio: {:.2}%).",
                background_fragments_uncapped, background.len(), background_fragments,
                (background_fragments as f64 / background.len() as f64) * 100.0);
            eprintln!("         This is biologically unrealistic (hitting same region multiple times in one cell).");
            eprintln!("         Consider reducing fragments_per_cell or increasing the background region set size.");
        }
        
        let mut sampled_regions = Vec::new();
        
        // Sample peaks using A-Res
        if peak_fragments > 0 && !peaks.is_empty() {
            let mut peak_sampler = VectorAResSampler::new(
                peak_fragments,  // Already capped at peaks.len() above
                Some(rng.gen())
            );

            let peak_iter = peaks.iter().zip(peak_weights.iter())
                .map(|(region, &weight)| (region.clone(), weight));

            // Debug: Check input iterator for duplicates
            if std::env::var("SCATRS_DEBUG").is_ok() {
                let iter_copy = peaks.iter().zip(peak_weights.iter())
                    .map(|(region, &weight)| (region.clone(), weight));
                let regions_vec: Vec<_> = iter_copy.collect();
                let mut region_set = std::collections::HashSet::new();
                for (region, _) in &regions_vec {
                    let key = format!("{}:{}-{}", region.chrom, region.start, region.end);
                    if !region_set.insert(key.clone()) {
                        eprintln!("WARNING: Duplicate region in A-Res input: {}", key);
                    }
                }
            }

            let sampled_peaks = peak_sampler.sample_stream(
                peak_iter,
                |(_, weight)| *weight
            );

            // Debug diagnostics (if enabled)
            if std::env::var("SCATRS_DEBUG").is_ok() && !sampled_peaks.is_empty() {
                let mut chr_counts = std::collections::HashMap::new();
                let mut region_set = std::collections::HashSet::new();
                for (region, _) in &sampled_peaks {
                    *chr_counts.entry(region.chrom.clone()).or_insert(0) += 1;
                    let key = format!("{}:{}-{}", region.chrom, region.start, region.end);
                    if !region_set.insert(key.clone()) {
                        eprintln!("ERROR: A-Res returned duplicate region: {}", key);
                    }
                }
                eprintln!("Sampled peaks by chromosome: {:?}", chr_counts);
            }

            // Each sampled region contributes exactly 1 fragment
            let peak_regions: Vec<SampledRegion> = sampled_peaks
                .into_iter()
                .map(|(region, weight)| SampledRegion {
                    region,
                    is_peak: true,
                    weight,
                    fragment_count: 1,
                })
                .collect();
            sampled_regions.extend(peak_regions);
        }
        
        // Sample background using A-Res
        if background_fragments > 0 && !background.is_empty() {
            let mut bg_sampler = VectorAResSampler::new(
                background_fragments,  // Already capped at background.len() above
                Some(rng.gen())
            );

            let bg_iter = background.iter().zip(background_weights.iter())
                .map(|(region, &weight)| (region.clone(), weight));

            let sampled_bg = bg_sampler.sample_stream(
                bg_iter,
                |(_, weight)| *weight
            );

            // Each sampled region contributes exactly 1 fragment
            let bg_regions: Vec<SampledRegion> = sampled_bg
                .into_iter()
                .map(|(region, weight)| SampledRegion {
                    region,
                    is_peak: false,
                    weight,
                    fragment_count: 1,
                })
                .collect();
            sampled_regions.extend(bg_regions);
        }
        
        Ok(sampled_regions)
    }
    
    /// Sample regions for a cell when only peaks are available (no background)
    /// 
    /// This is common for small test datasets or highly focused analyses
    /// where all fragments should come from peak regions.
    pub fn sample_regions_for_cell_peaks_only(
        &self,
        peaks: &[ScatrsRegion],
        peak_weights: &[f64],
        fragment_dist: &FragmentDistribution,
        rng: &mut StdRng,
        ambient_contamination: f64,
    ) -> Result<Vec<SampledRegion>> {
        // Determine total number of fragments for this cell
        let total_fragments = Self::sample_fragment_count(fragment_dist, rng)?;

        // Calculate source fragments (excluding ambient contamination)
        // This ensures the total fragment count stays as configured
        let source_fragments = ((total_fragments as f64) * (1.0 - ambient_contamination)) as u32;
        
        let mut sampled_regions = Vec::new();
        
        // All fragments come from peaks when no background is available
        if source_fragments > 0 && !peaks.is_empty() {
            // Cap at available regions (biological constraint: max 1 fragment per region per cell)
            let fragments_to_sample = (source_fragments as usize).min(peaks.len());

            let mut peak_sampler = VectorAResSampler::new(
                fragments_to_sample,
                Some(rng.gen())
            );

            let peak_iter = peaks.iter().zip(peak_weights.iter())
                .map(|(region, &weight)| (region.clone(), weight));

            let sampled_peaks = peak_sampler.sample_stream(
                peak_iter,
                |(_, weight)| *weight
            );

            // Debug diagnostics (if enabled)
            if std::env::var("SCATRS_DEBUG").is_ok() && !sampled_peaks.is_empty() {
                let mut chr_counts = std::collections::HashMap::new();
                for (region, _) in &sampled_peaks {
                    *chr_counts.entry(region.chrom.clone()).or_insert(0) += 1;
                }
                eprintln!("Sampled peaks by chromosome (peaks-only mode): {:?}", chr_counts);
            }

            // Each sampled region contributes exactly 1 fragment
            sampled_regions = sampled_peaks
                .into_iter()
                .map(|(region, weight)| SampledRegion {
                    region,
                    is_peak: true,
                    weight,
                    fragment_count: 1,
                })
                .collect();
        }
        
        Ok(sampled_regions)
    }

}

/// Sample a signal-to-noise value for a single cell
fn sample_signal_to_noise(
    distribution: &SignalToNoiseDistribution,
    rng: &mut StdRng,
) -> Result<f64> {
    let value = match distribution {
        SignalToNoiseDistribution::Fixed { signal_to_noise } => *signal_to_noise,
        SignalToNoiseDistribution::Beta { mean, sd } => {
            // Convert mean/sd to alpha/beta parameters for Beta distribution
            // mean = alpha / (alpha + beta)
            // variance = (alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1))
            let variance = sd * sd;

            // Solve for alpha and beta from mean and variance
            // alpha = mean * ((mean * (1 - mean) / variance) - 1)
            // beta = (1 - mean) * ((mean * (1 - mean) / variance) - 1)
            let mean_variance_ratio = (mean * (1.0 - mean)) / variance;
            if mean_variance_ratio <= 1.0 {
                return Err(anyhow::anyhow!(
                    "Invalid Beta parameters: variance too large for mean={}, sd={}. Max sd for this mean is {:.4}",
                    mean, sd, (mean * (1.0 - mean)).sqrt() * 0.99
                ));
            }

            let alpha = mean * (mean_variance_ratio - 1.0);
            let beta_param = (1.0 - mean) * (mean_variance_ratio - 1.0);

            let beta_dist = Beta::new(alpha, beta_param)
                .context("Invalid Beta distribution parameters")?;
            beta_dist.sample(rng)
        },
        SignalToNoiseDistribution::Gaussian { mean, sd, min, max } => {
            let normal = Normal::new(*mean, *sd)
                .context("Invalid Gaussian parameters")?;
            let sampled = normal.sample(rng);
            let min_val = min.unwrap_or(0.0);
            let max_val = max.unwrap_or(1.0);
            sampled.clamp(min_val, max_val)
        },
    };

    // Ensure value is in valid range [0, 1]
    Ok(value.clamp(0.0, 1.0))
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
    
    /// Sample fragments for cells with pre-assigned cell types
    /// 
    /// This method accepts pre-assigned cell types instead of computing them,
    /// ensuring consistency with the region sampling phase.
    pub fn sample_fragments_for_cells_with_types(
        &self,
        cell_regions: Vec<Vec<SampledRegion>>,
        region_fragments: &HashMap<String, Vec<ScatrsFragment>>,
        cell_type_assignments: &[String],  // Pre-assigned cell types
        _all_fragments: Option<&[ScatrsFragment]>,  // Kept for compatibility, not used
        _fragments_per_cell: Option<u32>,           // Kept for compatibility, not used
    ) -> Result<Vec<SimulatedCell>> {
        let num_cells = cell_regions.len();
        let pb = ProgressBar::new(num_cells as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} cells processed")
                .unwrap()
        );
        pb.set_message("Sampling fragments for cells");
        
        // Process cells in parallel
        let cells_result: Result<Vec<SimulatedCell>> = cell_regions
            .into_par_iter()
            .enumerate()
            .map(|(cell_idx, regions)| -> Result<SimulatedCell> {
                let mut rng = match self.seed {
                    Some(s) => StdRng::seed_from_u64(s + cell_idx as u64),
                    None => StdRng::from_entropy(),
                };
                
                // Use pre-assigned cell type
                let cell_type = &cell_type_assignments[cell_idx];
                let cell_id = format!("cell_{:06}", cell_idx + 1);
                
                let mut cell = SimulatedCell::new(cell_id, cell_type.clone());
                
                // Only redistribute fragments from real regions with data
                if !regions.is_empty() {
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
                            // Region has no fragments - this is expected for low-coverage regions
                            // Simply skip this region (no synthetic fragments will be generated)
                            // This is normal behavior - many regions have zero coverage
                        }
                    }
                }
                
                // Calculate FRiP score
                cell.calculate_frip();
                
                pb.inc(1);
                Ok(cell)
            })
            .collect();
        
        let cells = cells_result?;
        pb.finish_with_message(format!("Generated {} cells", cells.len()));
        
        Ok(cells)
    }
    
    /// Sample fragments for cells using type-specific fragment pools
    /// 
    /// This method ensures fragments are only sampled from the correct cell type's pool,
    /// preventing cross-contamination between cell types.
    pub fn sample_fragments_for_cells_with_type_pools(
        &self,
        cell_regions: Vec<Vec<SampledRegion>>,
        fragments_by_type: &HashMap<String, HashMap<String, Vec<ScatrsFragment>>>,
        cell_type_assignments: &[String],
        _cells_by_type: &HashMap<String, Vec<usize>>,
    ) -> Result<Vec<SimulatedCell>> {
        let num_cells = cell_regions.len();
        let pb = ProgressBar::new(num_cells as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} cells processed")
                .unwrap()
        );
        pb.set_message("Sampling fragments for cells from type-specific pools");
        
        // Process cells in parallel
        let cells_result: Result<Vec<SimulatedCell>> = cell_regions
            .into_par_iter()
            .enumerate()
            .map(|(cell_idx, regions)| -> Result<SimulatedCell> {
                let mut rng = match self.seed {
                    Some(s) => StdRng::seed_from_u64(s + cell_idx as u64),
                    None => StdRng::from_entropy(),
                };
                
                // Use pre-assigned cell type
                let cell_type = &cell_type_assignments[cell_idx];
                let cell_id = format!("cell_{:06}", cell_idx + 1);
                
                // Get the fragment pool for this cell's type
                // CRITICAL: Use type-specific fragment pool to prevent cross-contamination
                // Each cell must only sample fragments from its assigned cell type's BAM file
                // This preserves biological signal and prevents mixing of cell-type-specific patterns
                let type_fragments = fragments_by_type.get(cell_type)
                    .ok_or_else(|| anyhow::anyhow!("No fragments for cell type: {}", cell_type))?;
                
                let mut cell = SimulatedCell::new(cell_id, cell_type.clone());
                
                // Only redistribute fragments from real regions with data
                if !regions.is_empty() {
                    // Normal mode: sample from regions using type-specific fragments
                    for sampled_region in regions {
                        let region_key = format!(
                            "{}:{}-{}",
                            sampled_region.region.chrom,
                            sampled_region.region.start,
                            sampled_region.region.end
                        );
                        
                        // Get fragments for this region FROM THE CORRECT TYPE'S POOL
                        let available_fragments = type_fragments.get(&region_key)
                            .or_else(|| {
                                // Try alternative key format
                                type_fragments.get(&format!(
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
                            // Region has no fragments - this is expected for low-coverage regions
                            // Simply skip this region (no synthetic fragments will be generated)
                            // This is normal behavior - many regions have zero coverage
                        }
                    }
                }
                
                // Calculate FRiP score
                cell.calculate_frip();
                
                pb.inc(1);
                Ok(cell)
            })
            .collect();
        
        let cells = cells_result?;
        pb.finish_with_message(format!("Generated {} cells with type-specific fragments", cells.len()));
        
        Ok(cells)
    }
    
    pub fn sample_fragments_for_cells(
        &self,
        cell_regions: Vec<Vec<SampledRegion>>,
        region_fragments: &HashMap<String, Vec<ScatrsFragment>>,
        cell_types: &[String],
        cell_type_weights: Option<&[f64]>,         // Weights for cell type proportions
        _all_fragments: Option<&[ScatrsFragment]>,  // Kept for compatibility, not used
        _fragments_per_cell: Option<u32>,           // Kept for compatibility, not used
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
        let cells_result: Result<Vec<SimulatedCell>> = cell_regions
            .into_par_iter()
            .enumerate()
            .map(|(cell_idx, regions)| -> Result<SimulatedCell> {
                let mut rng = match self.seed {
                    Some(s) => StdRng::seed_from_u64(s + cell_idx as u64),
                    None => StdRng::from_entropy(),
                };
                
                // Use pre-assigned cell type
                let cell_type = &cell_type_assignments[cell_idx];
                let cell_id = format!("cell_{:06}", cell_idx + 1);
                
                let mut cell = SimulatedCell::new(cell_id, cell_type.clone());
                
                // Only redistribute fragments from real regions with data
                if !regions.is_empty() {
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
                        // Region has no fragments - this is expected for low-coverage regions
                        // Simply skip this region (no synthetic fragments will be generated)
                        // This is normal behavior - many regions have zero coverage
                    }
                    }  // Close the else block for normal mode
                }
                
                // Calculate FRiP score
                cell.calculate_frip();
                
                pb.inc(1);
                Ok(cell)
            })
            .collect();
        
        let cells = cells_result?;
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
/// use gtars_scatrs::sampling::DoubletSimulator;
/// use gtars_scatrs::models::ScatrsConfig;
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
    
    /// Add ambient contamination to simulated cells
    /// 
    /// This method adds ambient fragments to cells to simulate contamination from lysed cells.
    /// The number of ambient fragments is calculated to achieve the desired contamination rate
    /// while maintaining the total fragment count as configured.
    /// 
    /// # Arguments
    /// * `cells` - The cells to add contamination to
    /// * `contamination_rate` - Fraction of total fragments that should be ambient (0-1)
    /// * `all_fragments` - Pool of all fragments from which to sample ambient contamination
    /// * `target_fragments_per_cell` - The configured total fragments per cell
    /// 
    /// # Note
    /// This method assumes cells already have source fragments sampled, and adds ambient
    /// fragments to reach the target total. The source fragment sampling should have
    /// already accounted for the contamination rate by reducing the number of source fragments.
    pub fn add_ambient_contamination(
        &self,
        cells: &mut Vec<SimulatedCell>,
        contamination_rate: f64,
        all_fragments: &[ScatrsFragment],  // Pool of all fragments from staging
        target_fragments_per_cell: u32,  // The configured fragments_per_cell
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
            
            // Calculate number of ambient fragments based on the total fragment target
            // The cell currently has source fragments, we need to add ambient to reach total
            let current_count = cell.fragments.len();
            let ambient_count = (target_fragments_per_cell as f64 * contamination_rate) as usize;
            
            // Only add ambient if we haven't already reached the target
            let fragments_to_add = ambient_count.min(target_fragments_per_cell as usize - current_count);
            
            if fragments_to_add > 0 && !all_fragments.is_empty() {
                // Randomly sample from the ambient pool
                let uniform = Uniform::new(0, all_fragments.len());
                for _ in 0..fragments_to_add {
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
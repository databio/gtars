//! GenomicDistributions functions and extensions for RegionSet module
//!
//! This file includes popular statistics calculated on RegionSets and functions involving
//! TSS information and Reference Genome
//!

use std::collections::HashMap;

use gtars_core::models::{Region, RegionSet};

use crate::errors::GtarsGenomicDistError;
use crate::interval_ranges::IntervalRanges;
use crate::models::{
    ChromosomeStatistics, ClusterStats, DensityHomogeneity, DensityVector, Dinucleotide,
    RegionBin, SequenceAccess, SpacingStats,
};
use crate::utils::chrom_karyotype_key;

/// Trait for computing statistics and distributions of genomic intervals.
pub trait GenomicIntervalSetStatistics {
    /// Calculate basic statistics for regions on each chromosome.
    ///
    /// Returns a map from chromosome name to statistics including:
    /// - Region counts
    /// - Chromosome bounds (min start, max end)
    /// - Region length statistics (min, max, mean, median)
    fn chromosome_statistics(&self) -> HashMap<String, ChromosomeStatistics>;

    /// Compute the distribution of regions across chromosome bins.
    ///
    /// The genome is partitioned into `n_bins` fixed-size windows, where bin size
    /// is determined by the longest chromosome. Each region is assigned to the bin
    /// containing its midpoint (matching R GenomicDistributions behavior), so each
    /// region is counted exactly once regardless of width.
    fn region_distribution_with_bins(&self, n_bins: u32) -> HashMap<String, RegionBin>;

    /// Compute the distribution of regions across per-chromosome bins.
    ///
    /// Like `region_distribution_with_bins` but uses actual chromosome sizes
    /// to create bins per-chromosome, matching R's `getGenomeBins(chromSizes)`.
    /// Each chromosome gets `n_bins` bins sized to that chromosome's length.
    ///
    /// Regions on chromosomes not present in `chrom_sizes` are skipped.
    /// Regions whose midpoint falls beyond the stated chromosome size are also
    /// skipped (common with assembly mismatches, e.g. an hg19 BED paired with
    /// hg38 chrom_sizes). The total bin count may therefore be lower than the
    /// input region count; callers who need to detect mismatches can compare
    /// `sum(bin.n)` against their input region count.
    fn region_distribution_with_chrom_sizes(
        &self,
        n_bins: u32,
        chrom_sizes: &HashMap<String, u32>,
    ) -> HashMap<String, RegionBin>;

    /// Compute distances between consecutive regions on each chromosome.
    ///
    /// For each pair of adjacent regions on the same chromosome, returns the
    /// signed gap: `next.start - current.end`. Positive values indicate a gap,
    /// negative values indicate overlapping regions, zero means adjacent/abutting.
    ///
    /// Regions on chromosomes with fewer than 2 regions are skipped (they have
    /// no neighbors to measure against). Output length equals the total number
    /// of *gaps* across all multi-region chromosomes — generally shorter than
    /// the input region count. Output is not aligned 1:1 with input regions.
    /// No sentinel values are emitted.
    fn calc_neighbor_distances(&self) -> Result<Vec<i64>, GtarsGenomicDistError>;

    /// Compute the distance from each region to its nearest neighbor.
    ///
    /// For each region, takes the minimum absolute distance to its upstream
    /// and downstream neighbors. First and last regions on each chromosome
    /// use their only neighbor's distance. Overlapping neighbors have distance 0.
    ///
    /// Regions on chromosomes with only one region are skipped (they have no
    /// neighbors). Output length equals the total number of regions across all
    /// multi-region chromosomes — generally shorter than the input region
    /// count, and NOT aligned 1:1 with input regions. No sentinel values are
    /// emitted. Callers who need 1:1 alignment must filter their input to
    /// multi-region chromosomes first.
    ///
    /// Port of R GenomicDistributions `calcNearestNeighbors()`.
    fn calc_nearest_neighbors(&self) -> Result<Vec<u32>, GtarsGenomicDistError>;

    /// Compute region widths (end - start for each region).
    ///
    /// Port of R GenomicDistributions `calcWidth()`.
    fn calc_widths(&self) -> Vec<u32>;

    /// Compute summary statistics over the distribution of inter-region spacings.
    ///
    /// Wraps `calc_neighbor_distances()`: computes mean, median, standard
    /// deviation, IQR, and log-space mean/std over the positive inter-region
    /// gaps. Overlapping and abutting neighbors are excluded (matching
    /// `calc_neighbor_distances`); cross-chromosome pairs are never counted.
    ///
    /// Returns a `SpacingStats` with NaN float fields for inputs with zero
    /// positive gaps (empty, singleton per chromosome, or all-overlapping).
    ///
    /// Useful as a descriptive scalar summary of how evenly spaced peaks
    /// are — regular arrays (CTCF-style) have small IQR and low log_std,
    /// clustered pileups have large IQR and high log_std.
    fn calc_inter_peak_spacing(&self) -> SpacingStats;

    /// Compute cluster-level summary statistics at a given stitching radius.
    ///
    /// Wraps `IntervalRanges::cluster(radius_bp)`: groups regions into
    /// single-linkage connected components where two regions link if the
    /// bp gap between `prev.end` and `next.start` is at most `radius_bp`.
    /// Chromosome-scoped (no cross-chromosome linking).
    ///
    /// Returns a `ClusterStats` summarizing the cluster size distribution.
    /// For empty inputs `mean_cluster_size` and `fraction_clustered` are NaN.
    ///
    /// Useful for characterizing enhancer clusters, super-enhancer stitching
    /// distances, or any ChIP-seq-style "how clustered are my peaks?" question.
    fn calc_peak_clusters(&self, radius_bp: u32) -> ClusterStats;

    /// Compute the dense per-window peak count vector across the genome.
    ///
    /// Unlike `region_distribution_with_chrom_sizes`, which returns only
    /// bins containing ≥1 region, this returns the full zero-padded vector
    /// with one entry per window on every chromosome in `chrom_sizes`.
    /// Output is ordered by karyotypic chromosome order and then by bin
    /// index within chromosome, so feature vectors are stable across files
    /// (suitable for ML input).
    ///
    /// Regions on chromosomes not present in `chrom_sizes`, or whose midpoint
    /// falls past the stated chromosome size, are skipped (matching
    /// `region_distribution_with_chrom_sizes`).
    ///
    /// Returns an empty `DensityVector` if `chrom_sizes` is empty or
    /// `n_bins == 0`.
    ///
    /// # `n_bins` semantics
    ///
    /// `n_bins` is the **target bin count for the longest chromosome in
    /// `chrom_sizes`**, not the total number of bins returned. Bin width is
    /// derived as `max(chrom_sizes.values()) / n_bins` (floored, minimum 1
    /// bp), and every chromosome is then tiled with windows of that width.
    /// Shorter chromosomes therefore receive proportionally fewer bins —
    /// specifically `ceil(chrom_size / bin_width)` — so the total length of
    /// `counts` is
    ///
    /// ```text
    /// sum(ceil(chrom_size / bin_width)) for chrom in chrom_sizes
    /// ```
    ///
    /// which can substantially exceed `n_bins` when `chrom_sizes` contains
    /// many chromosomes. To target a specific bin width in bp, set
    /// `n_bins` to `max_chrom_len / desired_bin_width_bp`.
    ///
    /// # Per-chromosome bin width
    ///
    /// The last bin on each chromosome is narrower than `bin_width` whenever
    /// `chrom_size` is not an exact multiple of `bin_width`. Chromosomes
    /// shorter than `bin_width` (common with UCSC alt / random / unplaced
    /// contigs like `chrUn_*`, `*_random`, `*_alt`) reduce to a single bin
    /// whose effective width equals the chromosome size rather than
    /// `bin_width`. `counts[i]` is therefore a count per bin, not a count
    /// per `bin_width` bp — bins of different effective widths are not
    /// directly comparable as densities, and consumers computing per-bp
    /// rates should be aware of this when `chrom_sizes` includes contigs
    /// significantly shorter than `bin_width`.
    fn calc_density_vector(
        &self,
        chrom_sizes: &HashMap<String, u32>,
        n_bins: u32,
    ) -> DensityVector;

    /// Compute summary statistics over the dense per-window count vector.
    ///
    /// Builds the `DensityVector` internally (see `calc_density_vector`),
    /// then computes mean, population variance, coefficient of variation,
    /// Gini coefficient, and the count of nonzero windows. Useful as a
    /// scalar measure of peak spatial evenness.
    ///
    /// See `calc_density_vector` for the definition of `n_bins` (it is the
    /// target bin count for the longest chromosome, not the total window
    /// count) and for the treatment of chromosomes shorter than the derived
    /// `bin_width`. Both affect the interpretation of the statistics below —
    /// short contigs each contribute a narrow single-bin entry which dilutes
    /// `mean_count`, inflates `n_windows`, and raises `gini`.
    ///
    /// **Gini bias note:** the Gini coefficient is biased toward high
    /// values for sparse count distributions (many zero-count windows).
    /// Callers characterizing very sparse peak sets should also consult
    /// `n_nonzero_windows` to judge whether the Gini is meaningful.
    fn calc_density_homogeneity(
        &self,
        chrom_sizes: &HashMap<String, u32>,
        n_bins: u32,
    ) -> DensityHomogeneity;
}

impl GenomicIntervalSetStatistics for RegionSet {
    fn chromosome_statistics(&self) -> HashMap<String, ChromosomeStatistics> {
        let mut widths_by_chr: HashMap<&String, Vec<u32>> = HashMap::new();
        let mut bounds_by_chr: HashMap<&String, (u32, u32)> = HashMap::new();

        // single pass iterator: collect widths and track chromosome bounds
        for region in &self.regions {
            let width = region.width();
            widths_by_chr.entry(&region.chr).or_default().push(width);

            bounds_by_chr
                .entry(&region.chr)
                .and_modify(|(min_start, max_end)| {
                    *min_start = (*min_start).min(region.start);
                    *max_end = (*max_end).max(region.end);
                })
                .or_insert((region.start, region.end));
        }

        // compute statistics from sorted widths
        widths_by_chr
            .into_iter()
            .map(|(chr, mut widths)| {
                let count = widths.len() as u32;
                widths.sort_unstable();

                let minimum = widths[0];
                let maximum = widths[widths.len() - 1];
                let sum: u64 = widths.iter().map(|&w| w as u64).sum();
                let mean = sum as f64 / count as f64;

                let median = if count % 2 == 0 {
                    (widths[(count / 2 - 1) as usize] + widths[(count / 2) as usize]) as f64 / 2.0
                } else {
                    widths[(count / 2) as usize] as f64
                };

                let (start, end) = bounds_by_chr[&chr];

                (
                    chr.clone(),
                    ChromosomeStatistics {
                        chromosome: chr.clone(),
                        number_of_regions: count,
                        start_nucleotide_position: start,
                        end_nucleotide_position: end,
                        minimum_region_length: minimum,
                        maximum_region_length: maximum,
                        mean_region_length: mean,
                        median_region_length: median,
                    },
                )
            })
            .collect()
    }

    fn region_distribution_with_bins(&self, n_bins: u32) -> HashMap<String, RegionBin> {
        if self.regions.is_empty() {
            return HashMap::new();
        }

        // Use midpoint of each region for bin assignment, matching R GenomicDistributions.
        // This ensures each region is counted in exactly one bin regardless of width,
        // answering "where are the regions located?" rather than "how much coverage?"
        let chrom_maxes = self.get_max_end_per_chr();
        let chrom_max_length = match chrom_maxes.values().max() {
            Some(&v) => v,
            None => return HashMap::new(),
        };
        let bin_size = if n_bins == 0 {
            chrom_max_length.max(1)
        } else {
            (chrom_max_length / n_bins).max(1)
        };

        let mut plot_results: HashMap<String, RegionBin> = HashMap::new();

        for region in &self.regions {
            let mid = region.mid_point();
            let rid = mid / bin_size;
            let bin_start = rid * bin_size;
            let chrom_end = chrom_maxes.get(&region.chr).copied().unwrap_or(0);
            let bin_end = (bin_start + bin_size).min(chrom_end);

            let key = format!("{}-{}-{}", region.chr, bin_start, bin_end);
            if let Some(bin) = plot_results.get_mut(&key) {
                bin.n += 1;
            } else {
                plot_results.insert(
                    key,
                    RegionBin {
                        chr: region.chr.clone(),
                        start: bin_start,
                        end: bin_end,
                        n: 1,
                        rid,
                    },
                );
            }
        }

        plot_results
    }

    fn region_distribution_with_chrom_sizes(
        &self,
        n_bins: u32,
        chrom_sizes: &HashMap<String, u32>,
    ) -> HashMap<String, RegionBin> {
        if self.regions.is_empty() || n_bins == 0 {
            return HashMap::new();
        }

        // Proportional binning: the longest chromosome gets n_bins bins,
        // shorter chromosomes get proportionally fewer. This produces a
        // uniform bin width (in bp) across all chromosomes so that
        // positional heatmaps are reference-aligned.
        let max_chrom_len = chrom_sizes.values().copied().max().unwrap_or(1) as u64;
        let bin_width = (max_chrom_len / n_bins as u64).max(1);

        let mut plot_results: HashMap<String, RegionBin> = HashMap::new();

        for region in &self.regions {
            let chrom_size = match chrom_sizes.get(&region.chr) {
                Some(&s) => s,
                None => continue, // skip regions on chromosomes not in chrom_sizes
            };
            let bin_size = bin_width as u32;

            let mid = region.mid_point();
            // Skip regions whose midpoint falls beyond the stated chromosome size
            // (e.g. BED file assembled against a different reference than the one
            // supplied by chrom_sizes). Without this, rid would exceed n_bins and
            // produce bins with end < start.
            if mid >= chrom_size {
                continue;
            }
            let rid = mid / bin_size;
            let bin_start = rid * bin_size;
            let bin_end = (bin_start + bin_size).min(chrom_size);

            let key = format!("{}-{}-{}", region.chr, bin_start, bin_end);
            if let Some(bin) = plot_results.get_mut(&key) {
                bin.n += 1;
            } else {
                plot_results.insert(
                    key,
                    RegionBin {
                        chr: region.chr.clone(),
                        start: bin_start,
                        end: bin_end,
                        n: 1,
                        rid,
                    },
                );
            }
        }

        plot_results
    }

    fn calc_neighbor_distances(&self) -> Result<Vec<i64>, GtarsGenomicDistError> {
        let mut distances: Vec<i64> = vec![];

        for chr in self.iter_chroms() {
            let mut chr_regions: Vec<&Region> = self.iter_chr_regions(chr).collect();
            chr_regions.sort_by_key(|r| (r.start, r.end));

            if chr_regions.len() < 2 {
                continue;
            }

            for window in chr_regions.windows(2) {
                let distance = window[1].start as i64 - window[0].end as i64;
                // Only include positive distances (non-overlapping gaps), matching R
                if distance > 0 {
                    distances.push(distance);
                }
            }
        }

        Ok(distances)
    }

    fn calc_nearest_neighbors(&self) -> Result<Vec<u32>, GtarsGenomicDistError> {
        let mut nearest: Vec<u32> = vec![];

        for chr in self.iter_chroms() {
            let mut chr_regions: Vec<&Region> = self.iter_chr_regions(chr).collect();
            chr_regions.sort_by_key(|r| (r.start, r.end));

            if chr_regions.len() < 2 {
                // Single region on this chromosome — skip (R drops these)
                continue;
            }

            // compute absolute neighbor distances for this chromosome
            // overlapping regions get distance 0
            let distances: Vec<u32> = chr_regions
                .windows(2)
                .map(|w| {
                    let d = w[1].start as i64 - w[0].end as i64;
                    if d > 0 { d as u32 } else { 0 }
                })
                .collect();

            // first region: only has right neighbor
            nearest.push(distances[0]);

            // middle regions: min of left and right neighbor distances
            for pair in distances.windows(2) {
                nearest.push(pair[0].min(pair[1]));
            }

            // last region: only has left neighbor
            nearest.push(distances[distances.len() - 1]);
        }

        Ok(nearest)
    }

    fn calc_widths(&self) -> Vec<u32> {
        self.regions.iter().map(|r| r.width()).collect()
    }

    fn calc_inter_peak_spacing(&self) -> SpacingStats {
        // calc_neighbor_distances already per-chromosome and filters to positive gaps.
        let gaps: Vec<f64> = match self.calc_neighbor_distances() {
            Ok(v) => v.into_iter().map(|g| g as f64).collect(),
            Err(_) => Vec::new(),
        };

        if gaps.is_empty() {
            return SpacingStats {
                n_gaps: 0,
                mean: f64::NAN,
                median: f64::NAN,
                std: f64::NAN,
                iqr: f64::NAN,
                log_mean: f64::NAN,
                log_std: f64::NAN,
            };
        }

        let n = gaps.len();
        let mean = gaps.iter().sum::<f64>() / n as f64;
        let variance = gaps.iter().map(|g| (g - mean).powi(2)).sum::<f64>() / n as f64;
        let std = variance.sqrt();

        // Sorted copy for quantiles.
        let mut sorted = gaps.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let median = quantile_sorted(&sorted, 0.5);
        let q1 = quantile_sorted(&sorted, 0.25);
        let q3 = quantile_sorted(&sorted, 0.75);
        let iqr = q3 - q1;

        // log-space statistics: log10(gap + 1) handles zero and heavy tails.
        let log_gaps: Vec<f64> = gaps.iter().map(|g| (g + 1.0).log10()).collect();
        let log_mean = log_gaps.iter().sum::<f64>() / n as f64;
        let log_var =
            log_gaps.iter().map(|g| (g - log_mean).powi(2)).sum::<f64>() / n as f64;
        let log_std = log_var.sqrt();

        SpacingStats {
            n_gaps: n,
            mean,
            median,
            std,
            iqr,
            log_mean,
            log_std,
        }
    }

    fn calc_peak_clusters(&self, radius_bp: u32) -> ClusterStats {
        let ids = self.cluster(radius_bp);
        let total = ids.len();

        if total == 0 {
            return ClusterStats {
                radius_bp,
                n_clusters: 0,
                n_clustered_peaks: 0,
                mean_cluster_size: f64::NAN,
                max_cluster_size: 0,
                fraction_clustered: f64::NAN,
            };
        }

        // Count cluster sizes. Cluster IDs are dense 0..n_clusters from cluster().
        let max_id = *ids.iter().max().unwrap();
        let mut sizes: Vec<usize> = vec![0; (max_id as usize) + 1];
        for &id in &ids {
            sizes[id as usize] += 1;
        }

        let n_clusters = sizes.len();
        let max_cluster_size = *sizes.iter().max().unwrap_or(&0);

        // Clusters of size > 1 and peaks belonging to them.
        let nontrivial_sizes: Vec<usize> =
            sizes.iter().copied().filter(|&s| s > 1).collect();
        let n_clustered_peaks: usize = nontrivial_sizes.iter().sum();
        let mean_cluster_size = if nontrivial_sizes.is_empty() {
            f64::NAN
        } else {
            nontrivial_sizes.iter().sum::<usize>() as f64 / nontrivial_sizes.len() as f64
        };

        let fraction_clustered = n_clustered_peaks as f64 / total as f64;

        ClusterStats {
            radius_bp,
            n_clusters,
            n_clustered_peaks,
            mean_cluster_size,
            max_cluster_size,
            fraction_clustered,
        }
    }

    fn calc_density_vector(
        &self,
        chrom_sizes: &HashMap<String, u32>,
        n_bins: u32,
    ) -> DensityVector {
        if chrom_sizes.is_empty() || n_bins == 0 {
            return DensityVector {
                n_bins,
                bin_width: 0,
                counts: Vec::new(),
                chrom_offsets: Vec::new(),
            };
        }

        // Same bin-width formula as region_distribution_with_chrom_sizes:
        // longest chromosome gets n_bins bins; shorter chromosomes get
        // proportionally fewer. bin_width is clamped at 1 to avoid
        // zero-sized bins when n_bins > max_chrom_len.
        let max_chrom_len = *chrom_sizes.values().max().unwrap() as u64;
        let bin_width_u64 = (max_chrom_len / n_bins as u64).max(1);
        let bin_width = bin_width_u64 as u32;

        // Karyotypic chromosome order so feature vectors are stable across files.
        let mut chr_order: Vec<(&String, u32)> = chrom_sizes
            .iter()
            .map(|(c, s)| (c, *s))
            .collect();
        chr_order.sort_by(|a, b| {
            chrom_karyotype_key(a.0).cmp(&chrom_karyotype_key(b.0))
        });

        // Compute per-chromosome bin counts and global offsets.
        let mut chrom_offsets: Vec<(String, usize)> = Vec::with_capacity(chr_order.len());
        let mut per_chrom_bin_counts: HashMap<&str, (usize, usize)> = HashMap::new();
        // chr -> (offset, n_bins_in_chrom)
        let mut running: usize = 0;
        for (chr, size) in &chr_order {
            // Number of bins that cover [0, size): ceil(size / bin_width).
            // If size == 0, the chromosome contributes zero bins.
            let n_chrom_bins = if *size == 0 {
                0
            } else {
                ((*size as u64).div_ceil(bin_width_u64)) as usize
            };
            chrom_offsets.push(((*chr).clone(), running));
            per_chrom_bin_counts.insert(chr.as_str(), (running, n_chrom_bins));
            running += n_chrom_bins;
        }

        let total = running;
        let mut counts = vec![0u32; total];

        // Assign each region to its bin based on midpoint. Mirrors
        // region_distribution_with_chrom_sizes so the two calls agree on
        // which bin a peak lands in.
        for region in &self.regions {
            let (offset, n_chrom_bins) = match per_chrom_bin_counts.get(region.chr.as_str()) {
                Some(v) => *v,
                None => continue, // unknown chromosome
            };
            if n_chrom_bins == 0 {
                continue;
            }
            let chrom_size = chrom_sizes[&region.chr];
            let mid = region.mid_point();
            if mid >= chrom_size {
                continue;
            }
            let rid = (mid / bin_width) as usize;
            // rid can in principle reach n_chrom_bins - 1; safety clamp.
            if rid >= n_chrom_bins {
                continue;
            }
            counts[offset + rid] += 1;
        }

        DensityVector {
            n_bins,
            bin_width,
            counts,
            chrom_offsets,
        }
    }

    fn calc_density_homogeneity(
        &self,
        chrom_sizes: &HashMap<String, u32>,
        n_bins: u32,
    ) -> DensityHomogeneity {
        let dv = self.calc_density_vector(chrom_sizes, n_bins);
        let n_windows = dv.counts.len();

        if n_windows == 0 {
            return DensityHomogeneity {
                bin_width: dv.bin_width,
                n_windows: 0,
                n_nonzero_windows: 0,
                mean_count: f64::NAN,
                variance: f64::NAN,
                cv: f64::NAN,
                gini: f64::NAN,
            };
        }

        let n_nonzero_windows = dv.counts.iter().filter(|&&c| c > 0).count();
        let sum: u64 = dv.counts.iter().map(|&c| c as u64).sum();
        let mean_count = sum as f64 / n_windows as f64;

        // Population variance.
        let variance: f64 = if n_windows == 0 {
            0.0
        } else {
            dv.counts
                .iter()
                .map(|&c| (c as f64 - mean_count).powi(2))
                .sum::<f64>()
                / n_windows as f64
        };

        // Coefficient of variation. Undefined when the mean is zero
        // (e.g. empty RegionSet over populated chrom_sizes).
        let cv = if mean_count > 0.0 {
            variance.sqrt() / mean_count
        } else {
            f64::NAN
        };

        // Classical Gini coefficient on the sorted count vector:
        //   G = (2 * Σ i * x_i  −  (n + 1) * Σ x_i)  /  (n * Σ x_i)
        // where x is sorted ascending and i is 1-based.
        // For an all-zero vector the denominator is zero; by convention
        // a perfectly even (all-zero) distribution has Gini = 0.
        let gini = if sum == 0 {
            0.0
        } else {
            let mut sorted = dv.counts.clone();
            sorted.sort_unstable();
            let n = n_windows as f64;
            let cum: f64 = sorted
                .iter()
                .enumerate()
                .map(|(i, &c)| (i + 1) as f64 * c as f64)
                .sum();
            let total = sum as f64;
            (2.0 * cum - (n + 1.0) * total) / (n * total)
        };

        DensityHomogeneity {
            bin_width: dv.bin_width,
            n_windows,
            n_nonzero_windows,
            mean_count,
            variance,
            cv,
            gini,
        }
    }
}

/// Linear-interpolated quantile on an already-sorted slice.
///
/// Uses numpy's default ("linear") interpolation so results match
/// `np.quantile(..., interpolation='linear')`. Panics if the slice is
/// empty — callers must handle that upstream.
fn quantile_sorted(sorted: &[f64], q: f64) -> f64 {
    debug_assert!(!sorted.is_empty(), "quantile_sorted: empty slice");
    let n = sorted.len();
    if n == 1 {
        return sorted[0];
    }
    let pos = q * (n as f64 - 1.0);
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    if lo == hi {
        sorted[lo]
    } else {
        let frac = pos - lo as f64;
        sorted[lo] * (1.0 - frac) + sorted[hi] * frac
    }
}

///
///  Calculate GC content for bed file
///
/// Arguments:
/// - region_set: RegionSet object
/// - genome: GenomeAssembly object (reference genome)
/// - ignore_unk_chroms: bool to ignore unknown chromosomes for reference genome
///
pub fn calc_gc_content(
    region_set: &RegionSet,
    genome: &(impl SequenceAccess + ?Sized),
    ignore_unk_chroms: bool,
) -> Result<Vec<f64>, GtarsGenomicDistError> {
    // for region in region_set
    let mut gc_contents: Vec<f64> = vec![];
    for chr in region_set.iter_chroms() {
        // check if the chrom is even in genome
        if ignore_unk_chroms && !genome.contains_chr(chr) {
            continue;
        }

        for region in region_set.iter_chr_regions(chr) {
            let mut gc_count: u32 = 0;
            let mut total_count: u32 = 0;
            let seq = genome.get_sequence(region);

            match seq {
                Ok(seq) => {
                    for &base in &seq {
                        match base.to_ascii_lowercase() {
                            b'g' | b'c' => {
                                gc_count += 1;
                            }
                            _ => {}
                        }
                        total_count += 1;
                    }
                    if total_count > 0 {
                        gc_contents.push(gc_count as f64 / total_count as f64);
                    } else {
                        gc_contents.push(0.0);
                    }
                }
                Err(e) => {
                    if ignore_unk_chroms {
                        continue;
                    } else {
                        return Err(GtarsGenomicDistError::GCContentError(
                            region.chr.to_string(),
                            region.start,
                            region.end,
                            format!("{}", e),
                        ));
                    }
                }
            }
        }
    }

    Ok(gc_contents)
}

/// Canonical ordering of dinucleotides, matching GenomicDistributions' column order.
pub const DINUCL_ORDER: [Dinucleotide; 16] = [
    Dinucleotide::Aa, Dinucleotide::Ac, Dinucleotide::Ag, Dinucleotide::At,
    Dinucleotide::Ca, Dinucleotide::Cc, Dinucleotide::Cg, Dinucleotide::Ct,
    Dinucleotide::Ga, Dinucleotide::Gc, Dinucleotide::Gg, Dinucleotide::Gt,
    Dinucleotide::Ta, Dinucleotide::Tc, Dinucleotide::Tg, Dinucleotide::Tt,
];

/// Per-region dinucleotide frequencies.
///
/// Matches R GenomicDistributions `calcDinuclFreq`: one row per region,
/// 16 dinucleotide columns in [`DINUCL_ORDER`] order.
///
/// Arguments:
/// - `region_set`: RegionSet object
/// - `genome`: GenomeAssembly object (reference genome)
/// - `raw_counts`: if `true`, return raw integer-valued counts;
///   if `false`, return percentages (0–100) per row (matches R default)
/// - `ignore_unk_chroms`: if `true`, skip regions on chromosomes not in
///   the assembly; if `false`, error on unknown chromosomes
///
/// Returns a tuple `(region_labels, frequency_matrix)`:
/// - `region_labels`: `chr_start_end` for each region
/// - `frequency_matrix`: `Vec<[f64; 16]>` — one row per region
///
/// For pooled global counts across all regions, sum the columns of the
/// raw-counts matrix:
/// ```no_run
/// # use gtars_genomicdist::{calc_dinucl_freq, DINUCL_ORDER};
/// # use gtars_genomicdist::models::GenomeAssembly;
/// # use gtars_core::models::RegionSet;
/// # fn example(rs: &RegionSet, assembly: &GenomeAssembly) -> Result<(), Box<dyn std::error::Error>> {
/// let (_, matrix) = calc_dinucl_freq(rs, assembly, true, false)?;
/// let mut totals = [0.0f64; 16];
/// for row in &matrix {
///     for (i, &c) in row.iter().enumerate() {
///         totals[i] += c;
///     }
/// }
/// # Ok(()) }
/// ```
pub fn calc_dinucl_freq(
    region_set: &RegionSet,
    genome: &(impl SequenceAccess + ?Sized),
    raw_counts: bool,
    ignore_unk_chroms: bool,
) -> Result<(Vec<String>, Vec<[f64; 16]>), GtarsGenomicDistError> {
    let mut labels: Vec<String> = Vec::new();
    let mut matrix: Vec<[f64; 16]> = Vec::new();

    for chr in region_set.iter_chroms() {
        if ignore_unk_chroms && !genome.contains_chr(chr) {
            continue;
        }
        for region in region_set.iter_chr_regions(chr) {
            let seq = match genome.get_sequence(region) {
                Ok(s) => s,
                Err(e) => {
                    if ignore_unk_chroms {
                        continue;
                    }
                    return Err(e);
                }
            };
            let mut counts = [0u64; 16];
            let mut total: u64 = 0;

            for window in seq.windows(2) {
                if let Some(dinucl) = Dinucleotide::from_bytes(window) {
                    let idx = DINUCL_ORDER.iter().position(|d| *d == dinucl).unwrap();
                    counts[idx] += 1;
                    total += 1;
                }
            }

            let row: [f64; 16] = if raw_counts {
                let mut r = [0.0f64; 16];
                for (i, &c) in counts.iter().enumerate() {
                    r[i] = c as f64;
                }
                r
            } else if total > 0 {
                let mut r = [0.0f64; 16];
                for (i, &c) in counts.iter().enumerate() {
                    r[i] = (c as f64 / total as f64) * 100.0;
                }
                r
            } else {
                [0.0; 16]
            };

            labels.push(format!("{}_{}_{}",
                region.chr, region.start, region.end));
            matrix.push(row);
        }
    }

    Ok((labels, matrix))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::GenomeAssembly;

    use pretty_assertions::assert_eq;
    use rstest::*;
    use std::path::PathBuf;

    fn get_test_path(file_name: &str) -> Result<PathBuf, std::io::Error> {
        let file_path: PathBuf = std::env::current_dir()
            .unwrap()
            .join("../tests/data/regionset")
            .join(file_name);
        Ok(file_path)
    }

    #[rstest]
    fn test_statistics() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        let stats = region_set.chromosome_statistics();
        let chr1_stats = stats.get("chr1").unwrap();
        assert_eq!(chr1_stats.number_of_regions, 9);
        assert_eq!(chr1_stats.maximum_region_length, 9);
        assert_eq!(chr1_stats.minimum_region_length, 2);
        assert_eq!(chr1_stats.median_region_length, 3f64);
        assert_eq!(chr1_stats.start_nucleotide_position, 5);
        assert_eq!(chr1_stats.end_nucleotide_position, 36);
    }

    #[rstest]
    fn test_distribution_plot() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        let distribution = region_set.region_distribution_with_bins(5);
        assert_eq!(distribution.len(), 5);
        assert!((distribution.values().next().unwrap().rid as i32 > -1));
    }

    #[rstest]
    fn test_calculate_distances() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        let distances = region_set.calc_neighbor_distances().unwrap();
        // 9 regions on chr1 → 8 consecutive pairs, but only 4 have positive gaps.
        // Sorted: (5,7)(8,10)(11,13)(14,20)(16,18)(17,22)(25,28)(25,32)(27,36)
        // Gaps:    1     1     1    -4    -1     3    -3    -5
        assert_eq!(distances, vec![1, 1, 1, 3]);
    }

    #[rstest]
    fn test_calc_nearest_neighbors() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        let nearest = region_set.calc_nearest_neighbors().unwrap();
        // 9 regions on chr1, all-pairs absolute distances (neg→0):
        //   [1, 1, 1, 0, 0, 3, 0, 0]
        // Nearest = min(left, right) for interior; single neighbor for endpoints:
        //   first=1, min(1,1)=1, min(1,1)=1, min(1,0)=0, min(0,0)=0,
        //   min(0,3)=0, min(3,0)=0, min(0,0)=0, last=0
        assert_eq!(nearest, vec![1, 1, 1, 0, 0, 0, 0, 0, 0]);
    }

    #[rstest]
    fn test_calc_widths() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        let widths = region_set.calc_widths();
        assert_eq!(widths.len(), 9);
        assert_eq!(*widths.iter().min().unwrap(), 2);
        assert_eq!(*widths.iter().max().unwrap(), 9);
    }

    fn get_fasta_path(file_name: &str) -> PathBuf {
        std::env::current_dir()
            .unwrap()
            .join("../tests/data/fasta")
            .join(file_name)
    }

    // --- calc_nearest_neighbors bug regression ---

    #[rstest]
    fn test_nearest_neighbors_single_region_chrom() {
        // Single-region chromosomes are skipped (matching R behavior).
        let regions = vec![
            Region { chr: "chr1".into(), start: 10, end: 20, rest: None },
            Region { chr: "chr1".into(), start: 30, end: 40, rest: None },
            Region { chr: "chr2".into(), start: 100, end: 200, rest: None }, // lone region
        ];
        let rs = RegionSet::from(regions);
        let nearest = rs.calc_nearest_neighbors().unwrap();

        // Only chr1 regions contribute (2 values), chr2 lone region is skipped
        assert_eq!(nearest.len(), 2);
        assert_eq!(nearest, vec![10, 10]);
    }

    // --- GC content ---

    #[rstest]
    fn test_calc_gc_content() {
        // base.fa: chr1=GGAA (2G, 0C → 50%), chr2=GCGC (2G, 2C → 100%)
        let path = get_fasta_path("base.fa");
        let ga = GenomeAssembly::try_from(path.to_str().unwrap()).unwrap();

        let regions = vec![
            Region { chr: "chr1".into(), start: 0, end: 4, rest: None },
            Region { chr: "chr2".into(), start: 0, end: 4, rest: None },
        ];
        let rs = RegionSet::from(regions);
        let gc = calc_gc_content(&rs, &ga, false).unwrap();

        assert_eq!(gc.len(), 2);
        // iter_chroms preserves insertion order: chr1 first, chr2 second
        assert!((gc[0] - 0.5).abs() < 1e-10);  // chr1: GGAA → 50%
        assert!((gc[1] - 1.0).abs() < 1e-10);  // chr2: GCGC → 100%
    }

    #[rstest]
    fn test_calc_gc_content_ignore_unknown_chroms() {
        let path = get_fasta_path("base.fa");
        let ga = GenomeAssembly::try_from(path.to_str().unwrap()).unwrap();

        let regions = vec![
            Region { chr: "chr1".into(), start: 0, end: 4, rest: None },
            Region { chr: "chrUnknown".into(), start: 0, end: 10, rest: None },
        ];
        let rs = RegionSet::from(regions);

        // With ignore_unk_chroms=true, should skip unknown chromosome
        let gc = calc_gc_content(&rs, &ga, true).unwrap();
        assert_eq!(gc.len(), 1);

        // With ignore_unk_chroms=false, should error
        let gc_err = calc_gc_content(&rs, &ga, false);
        assert!(gc_err.is_err());
    }

    // --- Dinucleotide frequency ---

    #[rstest]
    fn test_calc_dinucl_freq_raw_counts() {
        // base.fa: chr1=GGAA → dinucleotides: GG, GA, AA (3 total)
        let path = get_fasta_path("base.fa");
        let ga = GenomeAssembly::try_from(path.to_str().unwrap()).unwrap();

        let regions = vec![
            Region { chr: "chr1".into(), start: 0, end: 4, rest: None },
        ];
        let rs = RegionSet::from(regions);
        let (labels, matrix) = calc_dinucl_freq(&rs, &ga, true, false).unwrap();

        assert_eq!(labels, vec!["chr1_0_4"]);
        assert_eq!(matrix.len(), 1);
        let row = &matrix[0];
        let gg_idx = DINUCL_ORDER.iter().position(|d| *d == Dinucleotide::Gg).unwrap();
        let ga_idx = DINUCL_ORDER.iter().position(|d| *d == Dinucleotide::Ga).unwrap();
        let aa_idx = DINUCL_ORDER.iter().position(|d| *d == Dinucleotide::Aa).unwrap();
        assert_eq!(row[gg_idx], 1.0);
        assert_eq!(row[ga_idx], 1.0);
        assert_eq!(row[aa_idx], 1.0);
        // total should be 3 (4 bases → 3 dinucleotides)
        let total: f64 = row.iter().sum();
        assert_eq!(total, 3.0);
    }

    #[rstest]
    fn test_calc_dinucl_freq_percentages() {
        // base.fa: chr2=GCGC → dinucleotides: GC, CG, GC (percentages)
        let path = get_fasta_path("base.fa");
        let ga = GenomeAssembly::try_from(path.to_str().unwrap()).unwrap();

        let regions = vec![
            Region { chr: "chr2".into(), start: 0, end: 4, rest: None },
        ];
        let rs = RegionSet::from(regions);
        let (labels, matrix) = calc_dinucl_freq(&rs, &ga, false, false).unwrap();

        assert_eq!(labels.len(), 1);
        assert_eq!(labels[0], "chr2_0_4");
        assert_eq!(matrix.len(), 1);

        // GC appears 2/3 times, CG appears 1/3 times
        let gc_idx = DINUCL_ORDER.iter().position(|d| *d == Dinucleotide::Gc).unwrap();
        let cg_idx = DINUCL_ORDER.iter().position(|d| *d == Dinucleotide::Cg).unwrap();

        let row = &matrix[0];
        assert!((row[gc_idx] - 200.0 / 3.0).abs() < 0.1); // ~66.67%
        assert!((row[cg_idx] - 100.0 / 3.0).abs() < 0.1);  // ~33.33%

        // percentages sum to 100
        let total: f64 = row.iter().sum();
        assert!((total - 100.0).abs() < 0.1);
    }

    #[rstest]
    fn test_calc_dinucl_freq_global_derivable() {
        // Global counts are derivable by column-summing the raw-counts matrix.
        // Two regions: chr1=GGAA, chr2=GCGC → pooled: GG×1, GA×1, AA×1, GC×2, CG×1 = 6 total
        let path = get_fasta_path("base.fa");
        let ga = GenomeAssembly::try_from(path.to_str().unwrap()).unwrap();
        let regions = vec![
            Region { chr: "chr1".into(), start: 0, end: 4, rest: None },
            Region { chr: "chr2".into(), start: 0, end: 4, rest: None },
        ];
        let rs = RegionSet::from(regions);
        let (_, matrix) = calc_dinucl_freq(&rs, &ga, true, false).unwrap();

        let mut totals = [0.0f64; 16];
        for row in &matrix {
            for (i, &c) in row.iter().enumerate() {
                totals[i] += c;
            }
        }
        let grand: f64 = totals.iter().sum();
        assert_eq!(grand, 6.0);
    }

    // --- Empty RegionSet edge cases ---

    #[rstest]
    fn test_empty_regionset_chromosome_statistics() {
        let rs = RegionSet::from(Vec::<Region>::new());
        let stats = rs.chromosome_statistics();
        assert!(stats.is_empty());
    }

    #[rstest]
    fn test_empty_regionset_region_distribution() {
        // Empty RegionSet should return empty distribution, not panic
        let rs = RegionSet::from(Vec::<Region>::new());
        let dist = rs.region_distribution_with_bins(10);
        assert!(dist.is_empty());
    }

    #[rstest]
    fn test_empty_regionset_calc_widths() {
        let rs = RegionSet::from(Vec::<Region>::new());
        assert!(rs.calc_widths().is_empty());
    }

    #[rstest]
    fn test_empty_regionset_neighbor_distances() {
        let rs = RegionSet::from(Vec::<Region>::new());
        let dists = rs.calc_neighbor_distances().unwrap();
        assert!(dists.is_empty());
    }

    #[rstest]
    fn test_empty_regionset_nearest_neighbors() {
        let rs = RegionSet::from(Vec::<Region>::new());
        let nearest = rs.calc_nearest_neighbors().unwrap();
        assert!(nearest.is_empty());
    }

    // --- GC content edge cases ---

    #[rstest]
    fn test_calc_gc_content_zero_length_region() {
        // A zero-length region (start == end) should return 0.0, not NaN
        let path = get_fasta_path("base.fa");
        let ga = GenomeAssembly::try_from(path.to_str().unwrap()).unwrap();

        let regions = vec![
            Region { chr: "chr1".into(), start: 2, end: 2, rest: None },
        ];
        let rs = RegionSet::from(regions);
        let gc = calc_gc_content(&rs, &ga, false).unwrap();
        assert_eq!(gc.len(), 1);
        assert!(!gc[0].is_nan());
        assert!((gc[0] - 0.0).abs() < 1e-10);
    }

    // --- Property-based tests inspired by R GenomicDistributions ---

    #[rstest]
    fn test_neighbor_distances_shift_invariant() {
        // R GenomicDistributions tests that shifting all regions by a constant
        // produces the same neighbor distances. This verifies the algorithm
        // depends on relative positions, not absolute coordinates.
        let regions = vec![
            Region { chr: "chr1".into(), start: 100, end: 200, rest: None },
            Region { chr: "chr1".into(), start: 300, end: 400, rest: None },
            Region { chr: "chr1".into(), start: 500, end: 700, rest: None },
        ];
        let rs1 = RegionSet::from(regions);

        let shifted = vec![
            Region { chr: "chr1".into(), start: 10100, end: 10200, rest: None },
            Region { chr: "chr1".into(), start: 10300, end: 10400, rest: None },
            Region { chr: "chr1".into(), start: 10500, end: 10700, rest: None },
        ];
        let rs2 = RegionSet::from(shifted);

        let d1 = rs1.calc_neighbor_distances().unwrap();
        let d2 = rs2.calc_neighbor_distances().unwrap();
        assert_eq!(d1, d2);

        let nn1 = rs1.calc_nearest_neighbors().unwrap();
        let nn2 = rs2.calc_nearest_neighbors().unwrap();
        assert_eq!(nn1, nn2);
    }

    #[rstest]
    fn test_overlapping_regions_neighbor_distance() {
        // Overlapping regions: negative distances are filtered out (matching R).
        let regions = vec![
            Region { chr: "chr1".into(), start: 100, end: 300, rest: None },
            Region { chr: "chr1".into(), start: 200, end: 400, rest: None },
            Region { chr: "chr1".into(), start: 500, end: 600, rest: None },
        ];
        let rs = RegionSet::from(regions);

        let dists = rs.calc_neighbor_distances().unwrap();
        // Only positive distances kept: [200,400) to [500,600) = 100
        assert_eq!(dists.len(), 1);
        assert_eq!(dists[0], 100);

        let nearest = rs.calc_nearest_neighbors().unwrap();
        // First region: only right neighbor, distance clamped to 0 (overlap)
        assert_eq!(nearest[0], 0);
        // Middle region: min(0, 100) = 0
        assert_eq!(nearest[1], 0);
        // Last region: only left neighbor, distance 100
        assert_eq!(nearest[2], 100);
    }

    #[rstest]
    fn test_region_distribution_total_count() {
        // Midpoint assignment: each region counted exactly once.
        // Matches R GenomicDistributions: sum(result$N) == length(query)
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let rs = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();
        let n_regions = rs.regions.len();

        let dist = rs.region_distribution_with_bins(5);
        let total_n: u32 = dist.values().map(|b| b.n).sum();
        assert_eq!(
            total_n as usize, n_regions,
            "total bin count ({}) should equal number of regions ({})",
            total_n, n_regions
        );
    }

    #[rstest]
    fn test_region_distribution_with_chrom_sizes_skips_out_of_range() {
        // Regions whose midpoint falls beyond the stated chromosome size are
        // silently skipped (assembly mismatch case). Previously this produced
        // bins with end < start or rid >= n_bins.
        let regions = vec![
            Region { chr: "chr1".into(), start: 100, end: 200, rest: None },  // mid=150, in range
            Region { chr: "chr1".into(), start: 800, end: 900, rest: None },  // mid=850, in range
            Region { chr: "chr1".into(), start: 1200, end: 1300, rest: None },// mid=1250, out of range (chrom_size=1000)
            Region { chr: "chr2".into(), start: 300, end: 400, rest: None },  // mid=350, in range
            Region { chr: "chr2".into(), start: 2000, end: 2100, rest: None },// mid=2050, out of range
            Region { chr: "chr3".into(), start: 0, end: 100, rest: None },    // chr3 not in chrom_sizes
        ];
        let rs = RegionSet::from(regions);
        let mut chrom_sizes = HashMap::new();
        chrom_sizes.insert("chr1".to_string(), 1000u32);
        chrom_sizes.insert("chr2".to_string(), 500u32);

        let bins = rs.region_distribution_with_chrom_sizes(10, &chrom_sizes);

        // Every bin should have end > start and rid < n_bins
        for bin in bins.values() {
            assert!(bin.end > bin.start, "bin {:?} has end <= start", bin);
            assert!(bin.rid < 10, "bin {:?} has rid >= n_bins", bin);
        }

        // Total counted regions: 2 (chr1) + 1 (chr2) = 3
        // (chr1's third region, chr2's second, and chr3's only region are all skipped)
        let total: u32 = bins.values().map(|b| b.n).sum();
        assert_eq!(total, 3, "expected 3 in-range regions counted");
    }

    #[rstest]
    fn test_gc_content_in_valid_range() {
        // R GenomicDistributions tests: all GC values should be in [0, 1]
        let path = get_fasta_path("base.fa");
        let ga = GenomeAssembly::try_from(path.to_str().unwrap()).unwrap();

        let regions = vec![
            Region { chr: "chr1".into(), start: 0, end: 4, rest: None },
            Region { chr: "chr2".into(), start: 0, end: 4, rest: None },
            Region { chr: "chrX".into(), start: 0, end: 8, rest: None },
        ];
        let rs = RegionSet::from(regions);
        let gc = calc_gc_content(&rs, &ga, false).unwrap();

        assert_eq!(gc.len(), 3);
        for &val in &gc {
            assert!(val >= 0.0 && val <= 1.0, "GC content out of range: {}", val);
            assert!(!val.is_nan(), "GC content should not be NaN");
        }
    }

    #[rstest]
    fn test_widths_match_region_coordinates() {
        // Width should equal end - start for each region
        let regions = vec![
            Region { chr: "chr1".into(), start: 10, end: 20, rest: None },
            Region { chr: "chr1".into(), start: 0, end: 0, rest: None },
            Region { chr: "chr2".into(), start: 100, end: 350, rest: None },
        ];
        let rs = RegionSet::from(regions);
        let widths = rs.calc_widths();
        assert_eq!(widths, vec![10, 0, 250]);
    }

    #[test]
    fn test_chromosome_statistics_large_widths() {
        // Regression: many wide regions whose total width exceeds u32::MAX
        // Use separate chromosomes to avoid coordinate overlap issues
        let regions: Vec<Region> = (0..5)
            .map(|i| Region {
                chr: format!("chr{}", i + 1),
                start: 0,
                end: 1_000_000_000,
                rest: None,
            })
            .collect();
        let rs = RegionSet::from(regions);
        let stats = rs.chromosome_statistics();
        assert_eq!(stats.len(), 5);
        let s = &stats["chr1"];
        // Mean width should be 1 billion
        assert!((s.mean_region_length - 1_000_000_000.0).abs() < 1.0);
    }

    // ── spatial-arrangement feature tests ───────────────────────────────

    fn make_rs(regions: &[(&str, u32, u32)]) -> RegionSet {
        let regs: Vec<Region> = regions
            .iter()
            .map(|(chr, s, e)| Region {
                chr: chr.to_string(),
                start: *s,
                end: *e,
                rest: None,
            })
            .collect();
        RegionSet::from(regs)
    }

    fn sizes(entries: &[(&str, u32)]) -> HashMap<String, u32> {
        entries.iter().map(|(c, s)| (c.to_string(), *s)).collect()
    }

    // --- calc_inter_peak_spacing ---

    #[test]
    fn test_inter_peak_spacing_regular_array() {
        // Regions spaced evenly at 10bp gaps → zero std, zero IQR.
        let rs = make_rs(&[
            ("chr1", 0, 10),
            ("chr1", 20, 30),
            ("chr1", 40, 50),
            ("chr1", 60, 70),
        ]);
        let s = rs.calc_inter_peak_spacing();
        assert_eq!(s.n_gaps, 3);
        assert!((s.mean - 10.0).abs() < 1e-10);
        assert!((s.median - 10.0).abs() < 1e-10);
        assert!(s.std.abs() < 1e-10);
        assert!(s.iqr.abs() < 1e-10);
    }

    #[test]
    fn test_inter_peak_spacing_variable_gaps() {
        // Gaps of 5, 10, 25 → mean 13.333, median 10.
        let rs = make_rs(&[
            ("chr1", 0, 10),
            ("chr1", 15, 25),  // gap 5
            ("chr1", 35, 45),  // gap 10
            ("chr1", 70, 80),  // gap 25
        ]);
        let s = rs.calc_inter_peak_spacing();
        assert_eq!(s.n_gaps, 3);
        assert!((s.mean - 40.0 / 3.0).abs() < 1e-9);
        assert!((s.median - 10.0).abs() < 1e-10);
        assert!(s.std > 0.0);
    }

    #[test]
    fn test_inter_peak_spacing_empty() {
        let rs = make_rs(&[]);
        let s = rs.calc_inter_peak_spacing();
        assert_eq!(s.n_gaps, 0);
        assert!(s.mean.is_nan());
        assert!(s.std.is_nan());
        assert!(s.iqr.is_nan());
    }

    #[test]
    fn test_inter_peak_spacing_singleton() {
        let rs = make_rs(&[("chr1", 10, 20)]);
        let s = rs.calc_inter_peak_spacing();
        assert_eq!(s.n_gaps, 0);
        assert!(s.mean.is_nan());
    }

    #[test]
    fn test_inter_peak_spacing_per_chromosome_only() {
        // Two regions each on two chromosomes; gaps must not span the boundary.
        let rs = make_rs(&[
            ("chr1", 0, 10),
            ("chr1", 50, 60),   // gap 40 on chr1
            ("chr2", 0, 10),
            ("chr2", 100, 110), // gap 90 on chr2
        ]);
        let s = rs.calc_inter_peak_spacing();
        assert_eq!(s.n_gaps, 2);
        assert!((s.mean - 65.0).abs() < 1e-10);
    }

    // --- calc_peak_clusters ---

    #[test]
    fn test_peak_clusters_all_singletons() {
        // Regions far apart, no clustering at radius 0.
        let rs = make_rs(&[("chr1", 0, 10), ("chr1", 100, 110), ("chr1", 200, 210)]);
        let c = rs.calc_peak_clusters(0);
        assert_eq!(c.n_clusters, 3);
        assert_eq!(c.n_clustered_peaks, 0);
        assert_eq!(c.max_cluster_size, 1);
        assert!(c.mean_cluster_size.is_nan());
        assert!((c.fraction_clustered - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_peak_clusters_one_big_cluster() {
        // Five regions all within 5bp of each other → single cluster of size 5.
        let rs = make_rs(&[
            ("chr1", 0, 10),
            ("chr1", 13, 20),
            ("chr1", 22, 30),
            ("chr1", 33, 40),
            ("chr1", 42, 50),
        ]);
        let c = rs.calc_peak_clusters(5);
        assert_eq!(c.n_clusters, 1);
        assert_eq!(c.n_clustered_peaks, 5);
        assert_eq!(c.max_cluster_size, 5);
        assert!((c.mean_cluster_size - 5.0).abs() < 1e-10);
        assert!((c.fraction_clustered - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_peak_clusters_mixed() {
        // Two clusters + one singleton at radius 5:
        //   cluster A: (0,10),(13,20)  — size 2
        //   cluster B: (100,110),(113,120),(122,130) — size 3
        //   singleton: (500, 510)
        let rs = make_rs(&[
            ("chr1", 0, 10),
            ("chr1", 13, 20),
            ("chr1", 100, 110),
            ("chr1", 113, 120),
            ("chr1", 122, 130),
            ("chr1", 500, 510),
        ]);
        let c = rs.calc_peak_clusters(5);
        assert_eq!(c.n_clusters, 3);
        assert_eq!(c.n_clustered_peaks, 5); // 2 + 3
        assert_eq!(c.max_cluster_size, 3);
        // Mean over nontrivial clusters only: (2 + 3) / 2 = 2.5
        assert!((c.mean_cluster_size - 2.5).abs() < 1e-10);
        assert!((c.fraction_clustered - 5.0 / 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_peak_clusters_cross_chromosome() {
        // Two regions on different chromosomes should never cluster.
        let rs = make_rs(&[("chr1", 0, 10), ("chr2", 0, 10)]);
        let c = rs.calc_peak_clusters(1_000_000);
        assert_eq!(c.n_clusters, 2);
        assert_eq!(c.n_clustered_peaks, 0);
    }

    #[test]
    fn test_peak_clusters_empty() {
        let rs = make_rs(&[]);
        let c = rs.calc_peak_clusters(100);
        assert_eq!(c.n_clusters, 0);
        assert_eq!(c.n_clustered_peaks, 0);
        assert_eq!(c.max_cluster_size, 0);
        assert!(c.mean_cluster_size.is_nan());
        assert!(c.fraction_clustered.is_nan());
        assert_eq!(c.radius_bp, 100);
    }

    // --- calc_density_vector ---

    #[test]
    fn test_density_vector_single_chrom() {
        // chr1 of size 100, 5 bins → bin_width = 20.
        // Peaks at 10,30,50 (midpoints) → bins 0,1,2.
        let rs = make_rs(&[("chr1", 5, 15), ("chr1", 25, 35), ("chr1", 45, 55)]);
        let cs = sizes(&[("chr1", 100)]);
        let dv = rs.calc_density_vector(&cs, 5);
        assert_eq!(dv.bin_width, 20);
        assert_eq!(dv.counts, vec![1, 1, 1, 0, 0]);
        assert_eq!(dv.chrom_offsets, vec![("chr1".to_string(), 0)]);
    }

    #[test]
    fn test_density_vector_zero_padding() {
        // Empty RegionSet, 5 bins on chr1 → all-zero vector of length 5.
        let rs = make_rs(&[]);
        let cs = sizes(&[("chr1", 100)]);
        let dv = rs.calc_density_vector(&cs, 5);
        assert_eq!(dv.counts, vec![0, 0, 0, 0, 0]);
    }

    #[test]
    fn test_density_vector_multi_chrom_karyotypic_order() {
        // chr1 size 100, chr2 size 60. bin_width = 100/5 = 20.
        // chr1 gets 5 bins, chr2 gets ceil(60/20) = 3 bins.
        // Order: chr1 first (offset 0), chr2 second (offset 5), total length 8.
        let rs = make_rs(&[
            ("chr2", 10, 20),  // chr2 mid=15 → bin 0 (offset 5)
            ("chr1", 50, 60),  // chr1 mid=55 → bin 2 (offset 0)
        ]);
        let cs = sizes(&[("chr1", 100), ("chr2", 60)]);
        let dv = rs.calc_density_vector(&cs, 5);
        assert_eq!(dv.counts.len(), 8);
        assert_eq!(dv.counts[2], 1); // chr1 bin 2
        assert_eq!(dv.counts[5], 1); // chr2 bin 0
        assert_eq!(
            dv.chrom_offsets,
            vec![("chr1".to_string(), 0), ("chr2".to_string(), 5)]
        );
    }

    #[test]
    fn test_density_vector_unknown_chrom_skipped() {
        // Regions on chr3 should be ignored when chrom_sizes only lists chr1.
        let rs = make_rs(&[("chr3", 10, 20)]);
        let cs = sizes(&[("chr1", 100)]);
        let dv = rs.calc_density_vector(&cs, 5);
        assert_eq!(dv.counts, vec![0, 0, 0, 0, 0]);
    }

    #[test]
    fn test_density_vector_past_chrom_end_skipped() {
        // A region whose midpoint is past chrom_size should be skipped.
        let rs = make_rs(&[("chr1", 95, 200)]); // mid = 147, > 100
        let cs = sizes(&[("chr1", 100)]);
        let dv = rs.calc_density_vector(&cs, 5);
        assert_eq!(dv.counts, vec![0, 0, 0, 0, 0]);
    }

    #[test]
    fn test_density_vector_empty_chrom_sizes() {
        let rs = make_rs(&[("chr1", 10, 20)]);
        let cs: HashMap<String, u32> = HashMap::new();
        let dv = rs.calc_density_vector(&cs, 5);
        assert!(dv.counts.is_empty());
        assert!(dv.chrom_offsets.is_empty());
    }

    #[test]
    fn test_density_vector_zero_bins() {
        let rs = make_rs(&[("chr1", 10, 20)]);
        let cs = sizes(&[("chr1", 100)]);
        let dv = rs.calc_density_vector(&cs, 0);
        assert!(dv.counts.is_empty());
    }

    // --- calc_density_homogeneity ---

    #[test]
    fn test_density_homogeneity_even_distribution() {
        // One peak in each of 5 bins → perfectly even. CV=0, Gini=0.
        let rs = make_rs(&[
            ("chr1", 5, 15),
            ("chr1", 25, 35),
            ("chr1", 45, 55),
            ("chr1", 65, 75),
            ("chr1", 85, 95),
        ]);
        let cs = sizes(&[("chr1", 100)]);
        let h = rs.calc_density_homogeneity(&cs, 5);
        assert_eq!(h.n_windows, 5);
        assert_eq!(h.n_nonzero_windows, 5);
        assert!((h.mean_count - 1.0).abs() < 1e-10);
        assert!(h.variance.abs() < 1e-10);
        assert!(h.cv.abs() < 1e-10);
        assert!(h.gini.abs() < 1e-10);
    }

    #[test]
    fn test_density_homogeneity_maximum_concentration() {
        // All peaks in the same bin → maximum concentration.
        // 5 peaks all in chr1 bin 0 (mid 5..10 → rid 0 with bin_width 20).
        let rs = make_rs(&[
            ("chr1", 1, 2),
            ("chr1", 3, 4),
            ("chr1", 5, 6),
            ("chr1", 7, 8),
            ("chr1", 9, 10),
        ]);
        let cs = sizes(&[("chr1", 100)]);
        let h = rs.calc_density_homogeneity(&cs, 5);
        assert_eq!(h.n_windows, 5);
        assert_eq!(h.n_nonzero_windows, 1);
        assert!(h.cv > 0.0);
        // Gini for [5,0,0,0,0] (sorted: [0,0,0,0,5]):
        //   2*(1*0 + 2*0 + 3*0 + 4*0 + 5*5) = 50
        //   (n+1)*sum = 6*5 = 30
        //   (50 - 30) / (5 * 5) = 20 / 25 = 0.8
        assert!((h.gini - 0.8).abs() < 1e-10);
    }

    #[test]
    fn test_density_homogeneity_empty_regionset() {
        // Empty RegionSet over populated chrom_sizes: all zeros, Gini=0, CV=NaN.
        let rs = make_rs(&[]);
        let cs = sizes(&[("chr1", 100)]);
        let h = rs.calc_density_homogeneity(&cs, 5);
        assert_eq!(h.n_windows, 5);
        assert_eq!(h.n_nonzero_windows, 0);
        assert_eq!(h.mean_count, 0.0);
        assert_eq!(h.variance, 0.0);
        assert!(h.cv.is_nan());
        assert_eq!(h.gini, 0.0);
    }

    #[test]
    fn test_density_homogeneity_empty_chrom_sizes() {
        let rs = make_rs(&[("chr1", 10, 20)]);
        let cs: HashMap<String, u32> = HashMap::new();
        let h = rs.calc_density_homogeneity(&cs, 5);
        assert_eq!(h.n_windows, 0);
        assert!(h.mean_count.is_nan());
        assert!(h.gini.is_nan());
    }
}

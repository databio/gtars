//! GenomicDistributions functions and extensions for RegionSet module
//!
//! This file includes popular statistics calculated on RegionSets and functions involving
//! TSS information and Reference Genome
//!

use std::collections::HashMap;

use gtars_core::models::{Region, RegionSet};

use crate::errors::GtarsGenomicDistError;
use crate::models::{ChromosomeStatistics, Dinucleotide, GenomeAssembly, RegionBin};

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
    /// Regions on chromosomes with fewer than 2 regions are skipped.
    fn calc_neighbor_distances(&self) -> Result<Vec<i64>, GtarsGenomicDistError>;

    /// Compute the distance from each region to its nearest neighbor.
    ///
    /// For each region, takes the minimum absolute distance to its upstream
    /// and downstream neighbors. First and last regions on each chromosome
    /// use their only neighbor's distance. Overlapping neighbors have distance 0.
    ///
    /// Port of R GenomicDistributions `calcNearestNeighbors()`.
    fn calc_nearest_neighbors(&self) -> Result<Vec<u32>, GtarsGenomicDistError>;

    /// Compute region widths (end - start for each region).
    ///
    /// Port of R GenomicDistributions `calcWidth()`.
    fn calc_widths(&self) -> Vec<u32>;
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

        let mut plot_results: HashMap<String, RegionBin> = HashMap::new();

        for region in &self.regions {
            let chrom_size = match chrom_sizes.get(&region.chr) {
                Some(&s) => s,
                None => continue, // skip regions on chromosomes not in chrom_sizes
            };
            let bin_size = (chrom_size / n_bins).max(1);

            let mid = region.mid_point();
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
    genome: &GenomeAssembly,
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
            let seq = genome.seq_from_region(region);

            match seq {
                Ok(seq) => {
                    for base in seq {
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

///
///  Calculate Dinucleotide frequencies
///
/// Arguments:
/// - region_set: RegionSet object
/// - genome: GenomeAssembly object (reference genome)
/// Return: Result of hashmap of dinucleotide and frequencies e.g. 'Aa: 13142'
pub fn calc_dinucl_freq(
    region_set: &RegionSet,
    genome: &GenomeAssembly,
) -> Result<HashMap<Dinucleotide, u64>, GtarsGenomicDistError> {
    let mut dinucl_freqs: HashMap<Dinucleotide, u64> = HashMap::new();

    for chr in region_set.iter_chroms() {
        for region in region_set.iter_chr_regions(chr) {
            let seq = genome.seq_from_region(region)?;
            for aas in seq.windows(2) {
                let dinucl = Dinucleotide::from_bytes(aas);
                match dinucl {
                    Some(dinucl) => {
                        let current_freq = dinucl_freqs.entry(dinucl).or_insert(0);
                        *current_freq += 1;
                    }
                    None => continue,
                }
            }
        }
    }

    Ok(dinucl_freqs)
}

/// Canonical ordering of dinucleotides, matching GenomicDistributions' column order.
pub const DINUCL_ORDER: [Dinucleotide; 16] = [
    Dinucleotide::Aa, Dinucleotide::Ac, Dinucleotide::Ag, Dinucleotide::At,
    Dinucleotide::Ca, Dinucleotide::Cc, Dinucleotide::Cg, Dinucleotide::Ct,
    Dinucleotide::Ga, Dinucleotide::Gc, Dinucleotide::Gg, Dinucleotide::Gt,
    Dinucleotide::Ta, Dinucleotide::Tc, Dinucleotide::Tg, Dinucleotide::Tt,
];

/// Per-region dinucleotide frequencies as percentages.
///
/// Returns a tuple of (region_labels, frequency_matrix) where:
/// - `region_labels` is `chr_start_end` for each region
/// - `frequency_matrix` is a `Vec<[f64; 16]>` — one row per region, columns
///   in [`DINUCL_ORDER`] order, values are percentages (0–100).
///
/// This matches the output format of R's GenomicDistributions::calcDinuclFreq.
pub fn calc_dinucl_freq_per_region(
    region_set: &RegionSet,
    genome: &GenomeAssembly,
) -> Result<(Vec<String>, Vec<[f64; 16]>), GtarsGenomicDistError> {
    let mut labels: Vec<String> = Vec::new();
    let mut matrix: Vec<[f64; 16]> = Vec::new();

    for chr in region_set.iter_chroms() {
        for region in region_set.iter_chr_regions(chr) {
            let seq = genome.seq_from_region(region)?;
            let mut counts = [0u64; 16];
            let mut total: u64 = 0;

            for window in seq.windows(2) {
                if let Some(dinucl) = Dinucleotide::from_bytes(window) {
                    let idx = DINUCL_ORDER.iter().position(|d| *d == dinucl).unwrap();
                    counts[idx] += 1;
                    total += 1;
                }
            }

            let freqs: [f64; 16] = if total > 0 {
                let mut f = [0.0f64; 16];
                for (i, &c) in counts.iter().enumerate() {
                    f[i] = (c as f64 / total as f64) * 100.0;
                }
                f
            } else {
                [0.0; 16]
            };

            labels.push(format!("{}_{}_{}",
                region.chr, region.start, region.end));
            matrix.push(freqs);
        }
    }

    Ok((labels, matrix))
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn test_calc_dinucl_freq() {
        // base.fa: chr1=GGAA → dinucleotides: GG, GA, AA
        let path = get_fasta_path("base.fa");
        let ga = GenomeAssembly::try_from(path.to_str().unwrap()).unwrap();

        let regions = vec![
            Region { chr: "chr1".into(), start: 0, end: 4, rest: None },
        ];
        let rs = RegionSet::from(regions);
        let freq = calc_dinucl_freq(&rs, &ga).unwrap();

        assert_eq!(*freq.get(&Dinucleotide::Gg).unwrap_or(&0), 1);
        assert_eq!(*freq.get(&Dinucleotide::Ga).unwrap_or(&0), 1);
        assert_eq!(*freq.get(&Dinucleotide::Aa).unwrap_or(&0), 1);
        // total should be 3 (4 bases → 3 dinucleotides)
        let total: u64 = freq.values().sum();
        assert_eq!(total, 3);
    }

    #[rstest]
    fn test_calc_dinucl_freq_per_region() {
        // base.fa: chr2=GCGC → dinucleotides: GC, CG, GC
        let path = get_fasta_path("base.fa");
        let ga = GenomeAssembly::try_from(path.to_str().unwrap()).unwrap();

        let regions = vec![
            Region { chr: "chr2".into(), start: 0, end: 4, rest: None },
        ];
        let rs = RegionSet::from(regions);
        let (labels, matrix) = calc_dinucl_freq_per_region(&rs, &ga).unwrap();

        assert_eq!(labels.len(), 1);
        assert_eq!(labels[0], "chr2_0_4");
        assert_eq!(matrix.len(), 1);

        // GC appears 2/3 times, CG appears 1/3 times
        // Find indices in DINUCL_ORDER
        let gc_idx = DINUCL_ORDER.iter().position(|d| *d == Dinucleotide::Gc).unwrap();
        let cg_idx = DINUCL_ORDER.iter().position(|d| *d == Dinucleotide::Cg).unwrap();

        let row = &matrix[0];
        assert!((row[gc_idx] - 200.0 / 3.0).abs() < 0.1); // ~66.67%
        assert!((row[cg_idx] - 100.0 / 3.0).abs() < 0.1);  // ~33.33%

        // all other dinucleotides should be 0
        let total: f64 = row.iter().sum();
        assert!((total - 100.0).abs() < 0.1);
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
}

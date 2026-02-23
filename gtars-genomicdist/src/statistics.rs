//! GenomicDistributions functions and extensions for RegionSet module
//!
//! This file includes popular statistics calculated on RegionSets and functions involving
//! TSS information and Reference Genome
//!

use std::collections::HashMap;

use gtars_core::models::{Region, RegionSet};
use gtars_overlaprs::multi_chrom_overlapper::IntoMultiChromOverlapper;
use gtars_overlaprs::OverlapperType;

use crate::errors::GtarsGenomicDistError;
use crate::models::{ChromosomeStatistics, Dinucleotide, GenomeAssembly, RegionBin};
use crate::utils::partition_genome_into_bins;

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
    /// is determined by the longest chromosome. Returns bins with overlap counts.
    fn region_distribution_with_bins(&self, n_bins: u32) -> HashMap<String, RegionBin>;

    /// Compute region distribution using 10 bins (default).
    ///
    /// Convenience method that calls `region_distribution_with_bins(10)`.
    fn region_distribution(&self) -> HashMap<String, RegionBin> {
        self.region_distribution_with_bins(250)
    }

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
                let sum: u32 = widths.iter().sum();
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
        let binned_genome = partition_genome_into_bins(&self.get_max_end_per_chr(), n_bins);
        let binned_genome_overlapper =
            binned_genome.into_multi_chrom_overlapper(OverlapperType::AIList);

        let region_hits = binned_genome_overlapper
            .find_overlaps_iter(self)
            .map(|(chr, iv)| Region {
                chr,
                start: iv.start,
                end: iv.end,
                rest: None,
            })
            .collect::<Vec<Region>>();

        let mut plot_results: HashMap<String, RegionBin> = HashMap::new();
        let region_length = region_hits.first().unwrap().end - region_hits.first().unwrap().start;

        for k in &region_hits {
            if let Some(region_bin) = plot_results.get_mut(&k.digest()) {
                region_bin.n += 1;
            } else {
                plot_results.insert(
                    k.digest().clone(),
                    RegionBin {
                        chr: k.chr.clone(),
                        start: k.start,
                        end: k.end,
                        n: 1,
                        rid: k.start / region_length,
                    },
                );
            }
        }
        plot_results
    }

    fn calc_neighbor_distances(&self) -> Result<Vec<i64>, GtarsGenomicDistError> {
        let mut distances: Vec<i64> = vec![];

        for chr in self.iter_chroms() {
            let chr_regions: Vec<&Region> = self.iter_chr_regions(chr).collect();

            if chr_regions.len() < 2 {
                continue;
            }

            for window in chr_regions.windows(2) {
                let distance = window[1].start as i64 - window[0].end as i64;
                distances.push(distance);
            }
        }

        Ok(distances)
    }

    fn calc_nearest_neighbors(&self) -> Result<Vec<u32>, GtarsGenomicDistError> {
        let mut nearest: Vec<u32> = vec![];

        for chr in self.iter_chroms() {
            let chr_regions: Vec<&Region> = self.iter_chr_regions(chr).collect();

            if chr_regions.len() < 2 {
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
                    gc_contents.push(gc_count as f64 / total_count as f64);
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

        let distribution = region_set.calc_neighbor_distances().unwrap();
        assert_eq!(distribution.len(), 8);
    }

    #[rstest]
    fn test_calc_nearest_neighbors() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        let nearest = region_set.calc_nearest_neighbors().unwrap();
        // 9 regions on chr1 → 9 nearest-neighbor distances (one per region)
        assert_eq!(nearest.len(), 9);
        // each nearest-neighbor distance should be ≤ the max of its neighbor distances
        // neighbor_dists are signed (i64), nearest are absolute (u32)
        let neighbor_dists = region_set.calc_neighbor_distances().unwrap();
        let abs_dists: Vec<u32> = neighbor_dists
            .iter()
            .map(|&d| if d > 0 { d as u32 } else { 0 })
            .collect();
        for i in 0..nearest.len() {
            if i == 0 {
                assert_eq!(nearest[i], abs_dists[0]);
            } else if i == nearest.len() - 1 {
                assert_eq!(nearest[i], abs_dists[abs_dists.len() - 1]);
            } else {
                assert_eq!(
                    nearest[i],
                    abs_dists[i - 1].min(abs_dists[i])
                );
            }
        }
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

    #[rstest]
    #[ignore] // only for local testing
    fn test_calc_dinucl_freq() {
        let ga = GenomeAssembly::try_from("/home/bnt4me/virginia/repos/bedboss/data/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/fasta/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.fa").unwrap();
        let rs =
            RegionSet::try_from("/home/bnt4me/Downloads/dcc005e8761ad5599545cc538f6a2a4d.bed.gz")
                .unwrap();

        let result = calc_dinucl_freq(&rs, &ga).unwrap();

        assert_ne!(result.len(), 0);
    }
}

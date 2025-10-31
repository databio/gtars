use std::collections::HashMap;

use gtars_core::models::{Region, RegionSet};
use gtars_overlaprs::multi_chrom_overlapper::IntoMultiChromOverlapper;
use gtars_overlaprs::OverlapperType;

use crate::models::{ChromosomeStatistics, RegionBin};
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
}

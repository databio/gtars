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

    /// Compute Neighbor_distances
    ///
    /// Returns a vector of vectors between the regions
    fn calc_neighbor_distances(&self) -> Result<Vec<u32>, GtarsGenomicDistError>;
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

    fn calc_neighbor_distances(&self) -> Result<Vec<u32>, GtarsGenomicDistError> {
        let mut distances: Vec<u32> = vec![];

        for chr in self.iter_chroms() {
            let chr_regions: Vec<&Region> = self.iter_chr_regions(chr).collect();

            // if there is only one region on the chromosome, skip it, can't calculate distance between one region
            if chr_regions.len() < 2 {
                continue;
            }

            for window in chr_regions.windows(2) {
                let distance: f32 = window[1].start as f32 - window[0].end as f32;
                let absolute_dist: u32 = if distance > 0f32 {
                    distance as u32
                } else {
                    0u32
                };
                distances.push(absolute_dist);
            }
        }

        Ok(distances)
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

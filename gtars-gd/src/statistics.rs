use crate::models::{ChromosomeStats, RegionBin};
use crate::utils::create_bin_regionset;
use gtars_core::models::{Region, RegionSet};
use gtars_overlaprs::genome_index::IntoGenomeIndex;
use gtars_overlaprs::OverlapperType;
use std::collections::HashMap;

pub trait Statistics {
    /// Calculate basic statistic for each region in the chromosome
    fn calculate_chr_statistics(&self) -> HashMap<String, ChromosomeStats>;

    /// Generate chrom distribution data based on regions used in the RegionSet
    fn generate_region_distribution(&self, n_bins: u32) -> HashMap<String, RegionBin>;
}

impl Statistics for RegionSet {
    fn calculate_chr_statistics(&self) -> HashMap<String, ChromosomeStats> {
        let mut stats: HashMap<String, ChromosomeStats> = HashMap::new();

        let mut regions_by_chr: HashMap<&String, Vec<&Region>> = HashMap::new();
        for region in &self.regions {
            regions_by_chr.entry(&region.chr).or_default().push(region);
        }

        for (chr, regions) in regions_by_chr {
            let count = regions.len() as u32;
            let widths: Vec<u32> = regions.iter().map(|r| r.width()).collect();
            let minimum = *widths.iter().min().unwrap_or(&0);
            let maximum = *widths.iter().max().unwrap_or(&0);

            let earliest_position = regions.iter().map(|r| r.start).min().unwrap_or(0);
            let end_position = regions.iter().map(|r| r.end).max().unwrap_or(0);

            let mean = widths.iter().sum::<u32>() as f64 / count as f64;

            let mut sorted_widths = widths.clone();
            sorted_widths.sort_unstable();
            let median = if count % 2 == 0 {
                (sorted_widths[(count / 2 - 1) as usize] + sorted_widths[(count / 2) as usize])
                    as f64
                    / 2.0
            } else {
                sorted_widths[(count / 2) as usize] as f64
            };

            stats.insert(
                chr.clone(),
                ChromosomeStats {
                    chromosome: chr.clone(),
                    number_of_regions: count,
                    start_nucleotide_position: earliest_position,
                    end_nucleotide_position: end_position,
                    minimum_region_length: minimum,
                    maximum_region_length: maximum,
                    mean_region_length: mean,
                    median_region_length: median,
                },
            );
        }

        stats
    }

    fn generate_region_distribution(&self, n_bins: u32) -> HashMap<String, RegionBin> {
        let universe: RegionSet = create_bin_regionset(self.get_max_end_per_chr(), n_bins);

        let universe_overlap = universe.into_genome_index(OverlapperType::Bits);

        let regions = universe_overlap.find_overlaps_to_rs(&self).unwrap();

        let mut plot_results: HashMap<String, RegionBin> = HashMap::new();

        let region_length: u32 = regions.get(0).unwrap().end - regions.get(0).unwrap().start;
        for k in &regions {
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

    #[rstest]
    fn test_chrom_bins() {
        let region_set = RegionSet::try_from(
            "/home/bnt4me/Downloads/test_chrombins/test_chrombins/GSM6732293_Con_liver-IP2.bed",
        )
        .unwrap();

        // let k = create_bin_regionset(region_set.get_max_end_per_chr());
        let k = region_set.generate_region_distribution();
        println!("{:?}", k);
    }
}

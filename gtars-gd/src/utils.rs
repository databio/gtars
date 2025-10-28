use std::collections::HashMap;

use gtars_core::models::region::Region;
use gtars_core::models::RegionSet;

/// Partitions a genome into bins of a fixed size.
///
/// The bin size is determined by dividing the length of the longest chromosome by `n_bins`.
/// Each chromosome is then tiled with regions of this calculated bin size. The final bin
/// on each chromosome may be smaller to accommodate the exact chromosome length.
pub fn partition_genome_into_bins(chrom_sizes: &HashMap<String, u32>, n_bins: u32) -> RegionSet {
    let mut regions = Vec::new();
    let chrom_max_length = *chrom_sizes.values().max().unwrap();
    let bin_size: u32 = chrom_max_length / n_bins;

    for (chr, size) in chrom_sizes {
        let mut start = 1;
        while start <= *size {
            let end = (start + bin_size - 1).min(*size);
            regions.push(Region {
                chr: chr.clone(),
                start,
                end,
                rest: None,
            });
            start = end + 1;
        }
    }

    RegionSet {
        regions,
        header: None,
        path: None,
    }
}

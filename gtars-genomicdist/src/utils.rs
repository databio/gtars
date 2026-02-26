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

    RegionSet::from(regions)
}

/// Returns a sort key that orders chromosome names karyotypically:
/// numeric (1, 2, …, 22) → X → Y → M/MT → everything else alphabetically.
pub fn chrom_karyotype_key(chr: &str) -> (u8, u32, String) {
    let bare = chr.strip_prefix("chr").unwrap_or(chr);
    match bare {
        "X" => (1, 0, String::new()),
        "Y" => (2, 0, String::new()),
        "M" | "MT" => (3, 0, String::new()),
        _ => match bare.parse::<u32>() {
            Ok(n) => (0, n, String::new()),
            Err(_) => (4, 0, bare.to_string()),
        },
    }
}

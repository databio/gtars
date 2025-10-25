use gtars_core::models::region::Region;
use gtars_core::models::RegionSet;
use std::collections::HashMap;
use std::hash::Hash;

pub fn create_bin_regionset(chrom_sizes: HashMap<String, u32>, n_bins: u32) -> RegionSet {
    let mut regions = Vec::new();
    let chrom_max_length = chrom_sizes.values().max().unwrap().clone();
    let bin_size: u32 = chrom_max_length / n_bins;

    for (chr, size) in chrom_sizes {
        let mut start = 1;
        while start <= size {
            let end = (start + bin_size - 1).min(size);
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

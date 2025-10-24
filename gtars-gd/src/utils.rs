use gtars_core::models::region::Region;
use gtars_core::models::RegionSet;
use std::collections::HashMap;
use std::hash::Hash;

pub fn create_bin_regionset(chrom_sizes: HashMap<String, u32>) -> RegionSet {
    let bin_count: u32 = 1000;

    let mut regions = Vec::new();

    for (chr, size) in chrom_sizes {
        let bin_size = (size as f64 / bin_count as f64).ceil() as u32;

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

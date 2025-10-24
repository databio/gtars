use gtars_core::models::region::Region;
use gtars_core::models::RegionSet;
use gtars_overlaprs::genome_index::IntoGenomeIndex;
use gtars_overlaprs::OverlapperType;
use std::collections::HashMap;
use std::hash::Hash;

#[derive(Debug, Clone)]
pub struct RegionBin {
    pub chr: String,
    pub start: u32,
    pub end: u32,
    pub n: u32,
    pub rid: u32,
}

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

pub fn generate_chrom_distribution(rs: RegionSet) -> HashMap<String, RegionBin> {
    let universe: RegionSet = create_bin_regionset(rs.get_max_end_per_chr());

    let universe_overlap = universe.into_genome_index(OverlapperType::Bits);

    let regions = universe_overlap.find_overlaps_to_rs(&rs).unwrap();

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
        let k = generate_chrom_distribution(region_set);
        println!("{:?}", k);
    }
}

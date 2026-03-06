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
    let chrom_max_length = match chrom_sizes.values().max() {
        Some(&v) => v,
        None => return RegionSet::from(regions),
    };
    let bin_size: u32 = if n_bins == 0 { chrom_max_length.max(1) } else { (chrom_max_length / n_bins).max(1) };

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

/// Compute the median of absolute values from a slice of signed distances,
/// filtering out sentinel values (`i64::MAX`) that indicate missing data.
///
/// Returns `None` if the slice is empty or contains only sentinels.
pub fn median_abs_distance(dists: &[i64]) -> Option<f64> {
    let mut sorted: Vec<f64> = dists
        .iter()
        .filter(|&&d| d != i64::MAX)
        .map(|&d| (d as f64).abs())
        .collect();
    if sorted.is_empty() {
        return None;
    }
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = sorted.len();
    if n % 2 == 0 {
        Some((sorted[n / 2 - 1] + sorted[n / 2]) / 2.0)
    } else {
        Some(sorted[n / 2])
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_partition_genome_into_bins() {
        let mut chrom_sizes = HashMap::new();
        chrom_sizes.insert("chr1".to_string(), 100);
        chrom_sizes.insert("chr2".to_string(), 50);

        let rs = partition_genome_into_bins(&chrom_sizes, 5);
        // bin_size = 100 / 5 = 20
        // chr1: 5 bins (1-20, 21-40, 41-60, 61-80, 81-100)
        // chr2: 3 bins (1-20, 21-40, 41-50)
        assert_eq!(rs.regions.len(), 8);

        let chr1_regions: Vec<_> = rs.regions.iter().filter(|r| r.chr == "chr1").collect();
        assert_eq!(chr1_regions.len(), 5);

        let chr2_regions: Vec<_> = rs.regions.iter().filter(|r| r.chr == "chr2").collect();
        assert_eq!(chr2_regions.len(), 3);
        // last bin on chr2 should be truncated
        let last = chr2_regions.iter().max_by_key(|r| r.start).unwrap();
        assert_eq!(last.end, 50);
    }

    #[test]
    fn test_chrom_karyotype_sort_order() {
        let mut chroms = vec!["chrM", "chrX", "chr2", "chr10", "chr1", "chrY", "chrUn_gl"];
        chroms.sort_by_key(|c| chrom_karyotype_key(c));
        assert_eq!(chroms, vec!["chr1", "chr2", "chr10", "chrX", "chrY", "chrM", "chrUn_gl"]);
    }

    #[test]
    fn test_partition_genome_bins_larger_than_chrom() {
        // Regression: when n_bins > chrom_max_length, bin_size was 0,
        // causing an infinite loop and OOM.
        let mut chrom_sizes = HashMap::new();
        chrom_sizes.insert("chr1".to_string(), 10);

        // 250 bins on a 10bp chromosome should not hang
        let rs = partition_genome_into_bins(&chrom_sizes, 250);
        // bin_size clamped to 1, so we get 10 bins
        assert_eq!(rs.regions.len(), 10);
    }

    #[test]
    fn test_median_abs_distance_filters_sentinels() {
        // i64::MAX sentinels should be excluded
        let dists = vec![10, -20, i64::MAX, 30];
        let median = median_abs_distance(&dists).unwrap();
        // abs values: [10, 20, 30] → median = 20
        assert!((median - 20.0).abs() < 1e-10);
    }

    #[test]
    fn test_median_abs_distance_all_sentinels() {
        let dists = vec![i64::MAX, i64::MAX];
        assert!(median_abs_distance(&dists).is_none());
    }

    #[test]
    fn test_median_abs_distance_empty() {
        assert!(median_abs_distance(&[]).is_none());
    }

    #[test]
    fn test_median_abs_distance_even_count() {
        let dists = vec![-10, 20, -30, 40];
        let median = median_abs_distance(&dists).unwrap();
        // abs sorted: [10, 20, 30, 40] → median = (20+30)/2 = 25
        assert!((median - 25.0).abs() < 1e-10);
    }

    #[test]
    fn test_chrom_karyotype_without_prefix() {
        // works without "chr" prefix
        let mut chroms = vec!["MT", "X", "2", "1", "Y"];
        chroms.sort_by_key(|c| chrom_karyotype_key(c));
        assert_eq!(chroms, vec!["1", "2", "X", "Y", "MT"]);
    }

    #[test]
    fn test_partition_genome_empty_chrom_sizes() {
        // Empty chrom_sizes should return empty RegionSet, not panic
        let chrom_sizes = HashMap::new();
        let rs = partition_genome_into_bins(&chrom_sizes, 10);
        assert!(rs.regions.is_empty());
    }

    #[test]
    fn test_partition_genome_zero_bins() {
        // n_bins=0 should not divide by zero; treat as 1 bin per chromosome
        let mut chrom_sizes = HashMap::new();
        chrom_sizes.insert("chr1".to_string(), 100);
        let rs = partition_genome_into_bins(&chrom_sizes, 0);
        // bin_size = max(100, 1) = 100, so chr1 gets 1 bin
        assert_eq!(rs.regions.len(), 1);
        assert_eq!(rs.regions[0].start, 1);
        assert_eq!(rs.regions[0].end, 100);
    }
}

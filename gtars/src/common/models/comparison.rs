use crate::common::models::{Region,RegionSet};
use std::collections::{HashSet};

#[derive(Debug, Default)]
pub struct ChromNameStats {
    pub xs: f64, // Extra Sequences (equivalent to sensitivity)
    pub q_and_m: u32,
    pub q_and_not_m: u32,
    pub not_q_and_m: u32,
    pub jaccard_index: f64,
    pub jaccard_index_binary: f64,
    pub passed_chrom_names: bool,
}

#[derive(Debug, Default)]
pub struct ChromLengthStats {
    pub oobr: f64, // Out of Bounds Range (equivalent to sensitivity)
    pub beyond_range: bool,
    pub num_of_chrom_beyond: u32,
    pub percentage_bed_chrom_beyond: f64,
    pub percentage_genome_chrom_beyond: f64,
}

#[derive(Debug, Default)]
pub struct SequenceFitStats {
    pub sequence_fit: f64,
}

#[derive(Debug, Default)]
pub struct CompatibilityStats {
    pub chrom_name_stats: ChromNameStats,
    pub chrom_length_stats: ChromLengthStats,
    pub chrom_sequence_fit_stats: SequenceFitStats,
}

pub fn compare_regionset_chrom_sizes(regionset_query: RegionSet, regionset_ref: RegionSet)-> CompatibilityStats {

    // Get the max end per chromosome for both regionsets.
    let query_map = regionset_query.get_max_end_per_chr();
    let ref_map = regionset_ref.get_max_end_per_chr();

    let query_keys: HashSet<_> = query_map.keys().cloned().collect();
    let ref_keys: HashSet<_> = ref_map.keys().cloned().collect();

    // Layer 1: Check chrom names and calculate stats.
    let mut chrom_name_stats = ChromNameStats::default();

    // Calculate intersections and differences using Rust's HashSet operations.
    let q_and_m_keys: HashSet<_> = query_keys.intersection(&ref_keys).cloned().collect();
    let q_and_not_m_keys: HashSet<_> = query_keys.difference(&ref_keys).cloned().collect();
    let not_q_and_m_keys: HashSet<_> = ref_keys.difference(&query_keys).cloned().collect();

    chrom_name_stats.q_and_m = q_and_m_keys.len() as u32;
    chrom_name_stats.q_and_not_m = q_and_not_m_keys.len() as u32;
    chrom_name_stats.not_q_and_m = not_q_and_m_keys.len() as u32;

    // Calculate Jaccard Index. Handle the case where the union is empty to avoid division by zero.
    let union_len = (query_keys.len() + ref_keys.len() - q_and_m_keys.len()) as f64;
    if union_len > 0.0 {
        chrom_name_stats.jaccard_index = q_and_m_keys.len() as f64 / union_len;
    }

    // Calculate binary Jaccard Index.
    let denominator_binary =
        chrom_name_stats.q_and_m + chrom_name_stats.not_q_and_m + chrom_name_stats.q_and_not_m;
    if denominator_binary > 0 {
        chrom_name_stats.jaccard_index_binary =
            chrom_name_stats.q_and_m as f64 / denominator_binary as f64;
    }

    // Determine if chrom names passed the first layer.
    chrom_name_stats.passed_chrom_names = chrom_name_stats.q_and_not_m == 0;

    // Calculate sensitivity.
    let sensitivity_denominator = chrom_name_stats.q_and_m + chrom_name_stats.q_and_not_m;
    if sensitivity_denominator > 0 {
        chrom_name_stats.xs = chrom_name_stats.q_and_m as f64 / sensitivity_denominator as f64;
    }

    // Layer 2: Check lengths, but only if layer 1 passed.
    let mut chrom_length_stats = ChromLengthStats::default();
    if chrom_name_stats.passed_chrom_names {
        let mut num_of_chrom_beyond = 0;
        let mut num_chrom_within_bounds = 0;

        for key in &q_and_m_keys {
            if let (Some(&query_len), Some(&ref_len)) = (query_map.get(key), ref_map.get(key)) {
                if query_len > ref_len {
                    num_of_chrom_beyond += 1;
                } else {
                    num_chrom_within_bounds += 1;
                }
            }
        }

        chrom_length_stats.beyond_range = num_of_chrom_beyond > 0;
        chrom_length_stats.num_of_chrom_beyond = num_of_chrom_beyond;

        let total_query_chroms = q_and_m_keys.len() as f64;
        if total_query_chroms > 0.0 {
            chrom_length_stats.percentage_bed_chrom_beyond =
                100.0 * num_of_chrom_beyond as f64 / total_query_chroms;
        }

        let total_ref_chroms = ref_keys.len() as f64;
        if total_ref_chroms > 0.0 {
            chrom_length_stats.percentage_genome_chrom_beyond =
                100.0 * num_of_chrom_beyond as f64 / total_ref_chroms;
        }

        let oobr_denominator = num_chrom_within_bounds + num_of_chrom_beyond;
        if oobr_denominator > 0 {
            chrom_length_stats.oobr = num_chrom_within_bounds as f64 / oobr_denominator as f64;
        }
    }

    // Layer 3: Calculate sequence fit if any query chrom names were present.
    let mut seq_fit_stats = SequenceFitStats::default();
    if !q_and_m_keys.is_empty() {
        let mut bed_sum: u64 = 0;
        let mut ref_genome_sum: u64 = 0;

        for q_chr in &q_and_m_keys {
            if let Some(&length) = ref_map.get(q_chr) {
                bed_sum += length as u64;
            }
        }

        for (_, &length) in &ref_map {
            ref_genome_sum += length as u64;
        }

        if ref_genome_sum > 0 {
            seq_fit_stats.sequence_fit = bed_sum as f64 / ref_genome_sum as f64;
        }
    }

    // Combine and return the final stats struct.
    CompatibilityStats {
        chrom_name_stats,
        chrom_length_stats,
        chrom_sequence_fit_stats: seq_fit_stats,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // A helper function to create a basic RegionSet for tests.
    fn create_regionset(regions: Vec<Region>) -> RegionSet {
        RegionSet {
            regions,
            header: None,
            path: None,
        }
    }

    #[test]
    fn test_compatible_regionsets() {
        // Query regions are fully contained within the reference regions.
        let query_set = create_regionset(vec![
            Region { chr: "chr1".to_string(), start: 1, end: 100, rest: None },
            Region { chr: "chr2".to_string(), start: 1, end: 150, rest: None },
        ]);
        let ref_set = create_regionset(vec![
            Region { chr: "chr1".to_string(), start: 1, end: 200, rest: None },
            Region { chr: "chr2".to_string(), start: 1, end: 300, rest: None },
            Region { chr: "chr3".to_string(), start: 1, end: 50, rest: None },
        ]);

        let stats = compare_regionset_chrom_sizes(query_set, ref_set);

        assert!(stats.chrom_name_stats.passed_chrom_names);
        assert_eq!(stats.chrom_name_stats.q_and_m, 2);
        assert_eq!(stats.chrom_name_stats.q_and_not_m, 0);
        assert_eq!(stats.chrom_name_stats.not_q_and_m, 1);
        assert!(!stats.chrom_length_stats.beyond_range);
        // The query's chr1 (100) and chr2 (150) are within the ref's lengths.
        assert_eq!(stats.chrom_length_stats.oobr, 1.0);
    }

    #[test]
    fn test_query_has_extra_chroms() {
        // The query contains a chromosome ("chr4") not present in the reference.
        let query_set = create_regionset(vec![
            Region { chr: "chr1".to_string(), start: 1, end: 100, rest: None },
            Region { chr: "chr4".to_string(), start: 1, end: 200, rest: None },
        ]);
        let ref_set = create_regionset(vec![
            Region { chr: "chr1".to_string(), start: 1, end: 200, rest: None },
            Region { chr: "chr2".to_string(), start: 1, end: 300, rest: None },
        ]);

        let stats = compare_regionset_chrom_sizes(query_set, ref_set);

        assert!(!stats.chrom_name_stats.passed_chrom_names);
        assert_eq!(stats.chrom_name_stats.q_and_m, 1);
        assert_eq!(stats.chrom_name_stats.q_and_not_m, 1);
        assert_eq!(stats.chrom_name_stats.not_q_and_m, 1);
        // Length stats should be default since chrom names didn't pass.
        assert_eq!(stats.chrom_length_stats.oobr, 0.0);
    }

    #[test]
    fn test_query_has_beyond_range_chrom() {
        // The query's "chr1" length is greater than the reference's "chr1" length.
        let query_set = create_regionset(vec![
            Region { chr: "chr1".to_string(), start: 1, end: 300, rest: None },
        ]);
        let ref_set = create_regionset(vec![
            Region { chr: "chr1".to_string(), start: 1, end: 200, rest: None },
            Region { chr: "chr2".to_string(), start: 1, end: 300, rest: None },
        ]);

        let stats = compare_regionset_chrom_sizes(query_set, ref_set);

        assert!(stats.chrom_name_stats.passed_chrom_names);
        assert!(stats.chrom_length_stats.beyond_range);
        assert_eq!(stats.chrom_length_stats.num_of_chrom_beyond, 1);
        assert_eq!(stats.chrom_length_stats.oobr, 0.0);
    }
    
}

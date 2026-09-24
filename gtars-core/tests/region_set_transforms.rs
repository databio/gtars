mod common;

use common::*;
use gtars_core::models::{Region, RegionSet};
use pretty_assertions::assert_eq;
use rstest::rstest;
use std::collections::HashMap;

#[rstest]
#[case::clamps_past_end(vec![("chr1", 90, 150)], vec![("chr1", 90, 100)])]
#[case::within_bounds(vec![("chr1", 10, 50)], vec![("chr1", 10, 50)])]
#[case::drops_unknown_chrom(vec![("chrX", 0, 50), ("chr1", 10, 20)], vec![("chr1", 10, 20)])]
#[case::keeps_zero_width_after_clamp(vec![("chr1", 100, 200)], vec![("chr1", 100, 100)])]
#[case::drops_inverted(vec![("chr1", 60, 40)], vec![])]
fn test_trim(#[case] input: Vec<Iv>, #[case] expected: Vec<Iv>) {
    let trimmed = make_regionset(input).trim(&chrom_sizes(&[("chr1", 100)]));
    assert_eq!(coords(&trimmed), expected);
}

#[rstest]
#[case::standard(("chr1", 1000, 2000), 500, 200, ("chr1", 500, 1200))]
#[case::saturates_at_zero(("chr1", 100, 500), 200, 50, ("chr1", 0, 150))]
#[case::at_origin(("chr1", 0, 100), 500, 200, ("chr1", 0, 200))]
fn test_promoters(#[case] input: Iv, #[case] up: u32, #[case] down: u32, #[case] expected: Iv) {
    let result = make_regionset(vec![input]).promoters(up, down);
    assert_eq!(coords(&result), vec![expected]);
}

#[rstest]
#[case::non_overlapping(
    vec![("chr1", 0, 5), ("chr1", 10, 15), ("chr1", 20, 25)],
    vec![("chr1", 0, 5), ("chr1", 10, 15), ("chr1", 20, 25)]
)]
// half-open touching intervals merge
#[case::adjacent(vec![("chr1", 0, 10), ("chr1", 10, 20)], vec![("chr1", 0, 20)])]
#[case::multi_chrom(
    vec![("chr1", 0, 10), ("chr1", 5, 15), ("chr2", 0, 10), ("chr2", 20, 30)],
    vec![("chr1", 0, 15), ("chr2", 0, 10), ("chr2", 20, 30)]
)]
// output chroms are sorted lexicographically, not karyotypically
#[case::lexicographic_chrom_order(
    vec![("chr10", 0, 10), ("chr2", 0, 10), ("chr1", 0, 10), ("chrX", 0, 10), ("chrM", 0, 10), ("chrY", 0, 10)],
    vec![("chr1", 0, 10), ("chr10", 0, 10), ("chr2", 0, 10), ("chrM", 0, 10), ("chrX", 0, 10), ("chrY", 0, 10)]
)]
#[case::unsorted_interleaved(
    vec![("chr10", 5, 15), ("chr2", 0, 10), ("chr10", 0, 8), ("chr2", 5, 20)],
    vec![("chr10", 0, 15), ("chr2", 0, 20)]
)]
#[case::empty(vec![], vec![])]
fn test_reduce(#[case] input: Vec<Iv>, #[case] expected: Vec<Iv>) {
    assert_eq!(coords(&make_regionset(input).reduce()), expected);
}

#[rstest]
fn test_shift_negative_saturates() {
    let result = make_regionset(vec![("chr1", 3, 10)]).shift(-5);
    assert_eq!(coords(&result), vec![("chr1", 0, 5)]);
}

// strand-unaware: "start" is always the lower coordinate
#[rstest]
#[case::upstream_of_start(true, false, ("chr1", 50, 100))]
#[case::downstream_of_end(false, false, ("chr1", 200, 250))]
#[case::both_around_start(true, true, ("chr1", 50, 150))]
fn test_flank(#[case] use_start: bool, #[case] both: bool, #[case] expected: Iv) {
    let result = make_regionset(vec![("chr1", 100, 200)]).flank(50, use_start, both);
    assert_eq!(coords(&result), vec![expected]);
}

#[rstest]
#[case::from_start(("chr1", 100, 200), 50, "start", ("chr1", 100, 150))]
#[case::from_end(("chr1", 100, 200), 50, "end", ("chr1", 150, 200))]
#[case::from_center(("chr1", 100, 200), 40, "center", ("chr1", 130, 170))]
// midpoint must not overflow u32
#[case::center_large_coords(
    ("chr1", 2_000_000_000, 2_100_000_000), 1000, "center",
    ("chr1", 2_049_999_500, 2_050_000_500)
)]
fn test_resize(#[case] input: Iv, #[case] width: u32, #[case] fix: &str, #[case] expected: Iv) {
    let result = make_regionset(vec![input]).resize(width, fix);
    assert_eq!(coords(&result), vec![expected]);
}

// narrow's start/end are 1-based and relative to the region
#[rstest]
#[case::start_and_end(Some(1), Some(50), None, ("chr1", 100, 150))]
#[case::start_and_width(Some(1), None, Some(30), ("chr1", 100, 130))]
// start=0 behaves like 1 instead of underflowing
#[case::start_zero(Some(0), Some(50), None, ("chr1", 100, 150))]
fn test_narrow(
    #[case] start: Option<u32>,
    #[case] end: Option<u32>,
    #[case] width: Option<u32>,
    #[case] expected: Iv,
) {
    let result = make_regionset(vec![("chr1", 100, 200)]).narrow(start, end, width);
    assert_eq!(coords(&result), vec![expected]);
}

#[rstest]
fn test_gaps_basic() {
    // Three peaks on chr1 with gaps between them; leading + trailing
    // gaps also present.
    let rs = make_regionset(vec![("chr1", 10, 20), ("chr1", 30, 40), ("chr1", 50, 60)]);
    let cs = chrom_sizes(&[("chr1", 100)]);
    let result = rs.gaps(&cs);
    assert_eq!(
        coords(&result),
        vec![
            ("chr1", 0, 10),   // leading
            ("chr1", 20, 30),  // between peak 1 and 2
            ("chr1", 40, 50),  // between peak 2 and 3
            ("chr1", 60, 100), // trailing
        ]
    );
}

#[rstest]
fn test_gaps_peak_at_origin_no_leading() {
    // First peak starts at 0 — no leading gap.
    let rs = make_regionset(vec![("chr1", 0, 10), ("chr1", 20, 30)]);
    let cs = chrom_sizes(&[("chr1", 100)]);
    let result = rs.gaps(&cs);
    let gaps: Vec<(u32, u32)> = result.regions.iter().map(|r| (r.start, r.end)).collect();
    assert_eq!(gaps, vec![(10, 20), (30, 100)]);
}

#[rstest]
fn test_gaps_peak_at_chrom_end_no_trailing() {
    // Last peak ends at chrom_size — no trailing gap.
    let rs = make_regionset(vec![("chr1", 10, 20), ("chr1", 80, 100)]);
    let cs = chrom_sizes(&[("chr1", 100)]);
    let result = rs.gaps(&cs);
    let gaps: Vec<(u32, u32)> = result.regions.iter().map(|r| (r.start, r.end)).collect();
    assert_eq!(gaps, vec![(0, 10), (20, 80)]);
}

#[rstest]
fn test_gaps_peak_past_chrom_end_clipped() {
    // Last peak extends past chrom_size — should be clipped, no trailing.
    let rs = make_regionset(vec![("chr1", 10, 20), ("chr1", 80, 150)]);
    let cs = chrom_sizes(&[("chr1", 100)]);
    let result = rs.gaps(&cs);
    let gaps: Vec<(u32, u32)> = result.regions.iter().map(|r| (r.start, r.end)).collect();
    assert_eq!(gaps, vec![(0, 10), (20, 80)]);
}

#[rstest]
fn test_gaps_empty_regionset_populated_chrom_sizes() {
    // No regions, but chrom_sizes has entries — emit whole-chrom gaps.
    let rs = RegionSet::from(Vec::<Region>::new());
    let cs = chrom_sizes(&[("chr1", 100), ("chr2", 50)]);
    let result = rs.gaps(&cs);
    let mut gaps = coords(&result);
    gaps.sort();
    assert_eq!(gaps, vec![("chr1", 0, 100), ("chr2", 0, 50)]);
}

#[rstest]
fn test_gaps_empty_regionset_empty_chrom_sizes() {
    let rs = RegionSet::from(Vec::<Region>::new());
    assert!(rs.gaps(&HashMap::new()).regions.is_empty());
}

#[rstest]
fn test_gaps_chromosome_not_in_chrom_sizes_skipped() {
    // Peak on chr2 with no chr2 entry in chrom_sizes — should be ignored.
    let rs = make_regionset(vec![("chr1", 10, 20), ("chr2", 5, 15)]);
    let cs = chrom_sizes(&[("chr1", 100)]);
    let result = rs.gaps(&cs);
    assert_eq!(coords(&result), vec![("chr1", 0, 10), ("chr1", 20, 100)]);
}

#[rstest]
fn test_gaps_full_chrom_gap_for_unrepresented_chrom() {
    // chrom_sizes has chr2 but input has no chr2 peaks — emit whole chr2.
    let rs = make_regionset(vec![("chr1", 10, 20)]);
    let cs = chrom_sizes(&[("chr1", 100), ("chr2", 200)]);
    let result = rs.gaps(&cs);
    let chr2_gaps: Vec<(u32, u32)> = result
        .regions
        .iter()
        .filter(|r| r.chr == "chr2")
        .map(|r| (r.start, r.end))
        .collect();
    assert_eq!(chr2_gaps, vec![(0, 200)]);
}

#[rstest]
fn test_gaps_overlapping_peaks_reduced() {
    // Overlapping peaks get merged by reduce() before gap computation.
    let rs = make_regionset(vec![
        ("chr1", 10, 30),
        ("chr1", 25, 40), // overlaps with previous
        ("chr1", 50, 60),
    ]);
    let cs = chrom_sizes(&[("chr1", 100)]);
    let result = rs.gaps(&cs);
    let gaps: Vec<(u32, u32)> = result.regions.iter().map(|r| (r.start, r.end)).collect();
    // After reduce: [10,40], [50,60] → gaps: [0,10], [40,50], [60,100]
    assert_eq!(gaps, vec![(0, 10), (40, 50), (60, 100)]);
}

#[rstest]
fn test_gaps_karyotypic_ordering() {
    // Output should be karyotypically ordered regardless of chrom_sizes insertion order.
    let rs = make_regionset(vec![("chr2", 10, 20), ("chr1", 10, 20), ("chr10", 10, 20)]);
    let cs = chrom_sizes(&[("chr10", 100), ("chr1", 100), ("chr2", 100)]);
    let result = rs.gaps(&cs);
    let mut order: Vec<&str> = result.regions.iter().map(|r| r.chr.as_str()).collect();
    order.dedup();
    assert_eq!(order, vec!["chr1", "chr2", "chr10"]);
}

#[rstest]
fn test_gaps_fully_covered_chrom_no_gaps() {
    // A single region spanning the whole chromosome yields zero gaps.
    let rs = make_regionset(vec![("chr1", 0, 100)]);
    let cs = chrom_sizes(&[("chr1", 100)]);
    assert!(rs.gaps(&cs).regions.is_empty());
}

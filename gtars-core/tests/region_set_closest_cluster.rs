mod common;

use common::*;
use pretty_assertions::assert_eq;
use rstest::rstest;

// Result rows are (query_idx, other_idx, distance); other_idx indexes the caller's `other`.
#[rstest]
#[case::overlapping(vec![("chr1", 100, 200)], vec![("chr1", 150, 250)], vec![(0, 0, 0)])]
#[case::upstream(vec![("chr1", 200, 300)], vec![("chr1", 100, 150)], vec![(0, 0, 50)])]
#[case::downstream(vec![("chr1", 100, 150)], vec![("chr1", 200, 300)], vec![(0, 0, 50)])]
#[case::choose_nearer(vec![("chr1", 100, 110)], vec![("chr1", 50, 60), ("chr1", 115, 120)], vec![(0, 1, 5)])]
#[case::unsorted_other(
    vec![("chr1", 100, 110)],
    vec![("chr1", 500, 510), ("chr1", 120, 130), ("chr1", 900, 910)],
    vec![(0, 1, 10)]
)]
// query chroms absent from `other` are omitted
#[case::multi_chrom_no_cross(vec![("chr1", 100, 200), ("chr2", 100, 200)], vec![("chr1", 300, 400)], vec![(0, 0, 100)])]
#[case::absent_chrom(vec![("chr2", 100, 200)], vec![("chr1", 100, 200)], vec![])]
#[case::nearest_far_left_of_insertion(
    vec![("chr1", 1000, 1010)],
    vec![("chr1", 0, 999), ("chr1", 1, 2), ("chr1", 3, 4), ("chr1", 5, 6), ("chr1", 1020, 2000)],
    vec![(0, 0, 1)]
)]
#[case::overlap_far_right_of_insertion(
    vec![("chr1", 500, 510)],
    vec![("chr1", 0, 490), ("chr1", 520, 521), ("chr1", 522, 523), ("chr1", 524, 525),
         ("chr1", 526, 527), ("chr1", 528, 529), ("chr1", 530, 531), ("chr1", 508, 509)],
    vec![(0, 7, 0)]
)]
#[case::wide_region_far_left(
    vec![("chr1", 1000, 1010)],
    vec![("chr1", 10, 998), ("chr1", 500, 501), ("chr1", 600, 601), ("chr1", 700, 701), ("chr1", 1050, 1060)],
    vec![(0, 0, 2)]
)]
// a long interval far to the left still overlaps the query
#[case::long_interval_overlaps(
    vec![("chr1", 500, 501)],
    vec![("chr1", 0, 10000), ("chr1", 100, 101), ("chr1", 200, 201), ("chr1", 300, 301), ("chr1", 400, 401)],
    vec![(0, 0, 0)]
)]
fn test_closest(
    #[case] query: Vec<Iv>,
    #[case] other: Vec<Iv>,
    #[case] expected: Vec<(usize, usize, i64)>,
) {
    let result = make_regionset(query).closest(&make_regionset(other));
    assert_eq!(result, expected);
}

#[rstest]
fn test_closest_variable_density() {
    let mut other: Vec<Iv> = (0..20).map(|i| ("chr1", i * 2, i * 2 + 1)).collect();
    other.push(("chr1", 4990, 4999));
    other.push(("chr1", 6000, 6001));
    let result = make_regionset(vec![("chr1", 5000, 5010)]).closest(&make_regionset(other));
    assert_eq!(result, vec![(0, 20, 1)]);
}

// Cluster ids follow sorted position but are returned in input order.
#[rstest]
#[case::separate(vec![("chr1", 0, 10), ("chr1", 20, 30), ("chr1", 40, 50)], 0, vec![0, 1, 2])]
#[case::overlapping(vec![("chr1", 0, 15), ("chr1", 10, 25)], 0, vec![0, 0])]
// half-open touching intervals have gap 0
#[case::adjacent(vec![("chr1", 0, 10), ("chr1", 10, 20)], 0, vec![0, 0])]
#[case::within_max_gap(vec![("chr1", 0, 10), ("chr1", 15, 25)], 5, vec![0, 0])]
#[case::multi_chrom(vec![("chr1", 0, 10), ("chr2", 0, 10)], 0, vec![0, 1])]
#[case::unsorted_input(vec![("chr1", 100, 110), ("chr1", 0, 10)], 0, vec![1, 0])]
#[case::single(vec![("chr1", 100, 200)], 0, vec![0])]
#[case::empty(vec![], 0, vec![])]
fn test_cluster(#[case] input: Vec<Iv>, #[case] max_gap: u32, #[case] expected: Vec<u32>) {
    assert_eq!(make_regionset(input).cluster(max_gap), expected);
}

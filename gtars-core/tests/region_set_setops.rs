mod common;

use common::*;
use gtars_core::models::{IntervalSetOps, Region, RegionSet};
use pretty_assertions::assert_eq;
use rstest::rstest;

#[rstest]
fn test_setops_on_bed_files() {
    // dummy.bed reduces to chr1:2-12; dummy_b.bed is chr1:3-5, chr1:8-10
    let a = load_bed("dummy.bed");
    let b = load_bed("dummy_b.bed");
    assert_eq!(coords(&a.reduce()), vec![("chr1", 2, 12)]);
    assert_eq!(coords(&a.union(&b)), vec![("chr1", 2, 12)]);
    assert_eq!(
        coords(&a.setdiff(&b)),
        vec![("chr1", 2, 3), ("chr1", 5, 8), ("chr1", 10, 12)]
    );
    assert_approx(a.jaccard(&b), 0.4);
}

#[rstest]
#[case::disjoint(vec![("chr1", 0, 10)], vec![("chr1", 20, 30)], vec![("chr1", 0, 10), ("chr1", 20, 30)])]
#[case::overlapping(vec![("chr1", 0, 15)], vec![("chr1", 10, 25)], vec![("chr1", 0, 25)])]
// half-open touching intervals merge
#[case::adjacent(vec![("chr1", 0, 10)], vec![("chr1", 10, 20)], vec![("chr1", 0, 20)])]
#[case::one_bp_gap(vec![("chr1", 0, 10)], vec![("chr1", 11, 20)], vec![("chr1", 0, 10), ("chr1", 11, 20)])]
#[case::contained(vec![("chr1", 0, 100)], vec![("chr1", 20, 50)], vec![("chr1", 0, 100)])]
#[case::identical(vec![("chr1", 10, 20), ("chr1", 30, 40)], vec![("chr1", 10, 20), ("chr1", 30, 40)], vec![("chr1", 10, 20), ("chr1", 30, 40)])]
#[case::multi_chrom(vec![("chr2", 0, 10)], vec![("chr1", 0, 10)], vec![("chr1", 0, 10), ("chr2", 0, 10)])]
#[case::one_empty(vec![("chr1", 0, 10), ("chr1", 20, 30)], vec![], vec![("chr1", 0, 10), ("chr1", 20, 30)])]
#[case::both_empty(vec![], vec![], vec![])]
fn test_union(#[case] a: Vec<Iv>, #[case] b: Vec<Iv>, #[case] expected: Vec<Iv>) {
    let result = make_regionset(a).union(&make_regionset(b));
    assert_eq!(coords(&result), expected);
}

#[rstest]
fn test_union_into_matches_union() {
    // Use overlapping regions on the same chr to exercise the merge in
    // `reduce`, plus a separate chr to exercise ordering.
    let a_regions = vec![make_region("chr1", 100, 200), make_region("chr2", 0, 50)];
    let b_regions = vec![make_region("chr1", 150, 250), make_region("chr3", 10, 20)];

    let a = RegionSet::from(a_regions.clone());
    let b = RegionSet::from(b_regions.clone());
    let borrowed = a.union(&b);

    let a2 = RegionSet::from(a_regions);
    let b2 = RegionSet::from(b_regions);
    let consumed = a2.union_into(b2);

    assert_eq!(borrowed.regions, consumed.regions);
}

#[rstest]
#[case::overlapping(vec![("chr1", 100, 200)], vec![("chr1", 150, 250)], vec![("chr1", 150, 200)])]
#[case::no_overlap(vec![("chr1", 100, 200)], vec![("chr1", 300, 400)], vec![])]
// half-open touching intervals share no bases
#[case::adjacent(vec![("chr1", 0, 10)], vec![("chr1", 10, 20)], vec![])]
#[case::contained(vec![("chr1", 0, 100)], vec![("chr1", 30, 70)], vec![("chr1", 30, 70)])]
#[case::one_vs_many(vec![("chr1", 100, 300)], vec![("chr1", 120, 150), ("chr1", 200, 250), ("chr1", 400, 500)], vec![("chr1", 120, 150), ("chr1", 200, 250)])]
#[case::many_vs_one(vec![("chr1", 0, 10), ("chr1", 20, 30)], vec![("chr1", 5, 25)], vec![("chr1", 5, 10), ("chr1", 20, 25)])]
#[case::multi_chrom(vec![("chr1", 0, 20), ("chr2", 0, 10)], vec![("chr1", 10, 30), ("chr2", 5, 15)], vec![("chr1", 10, 20), ("chr2", 5, 10)])]
#[case::empty_other(vec![("chr1", 100, 200)], vec![], vec![])]
#[case::empty_self(vec![], vec![("chr1", 100, 200)], vec![])]
fn test_intersect(#[case] a: Vec<Iv>, #[case] b: Vec<Iv>, #[case] expected: Vec<Iv>) {
    let result = make_regionset(a).intersect(&make_regionset(b));
    let mut got = coords(&result);
    got.sort();
    assert_eq!(got, expected);
}

#[rstest]
#[case::middle(vec![("chr1", 0, 10)], vec![("chr1", 3, 7)], vec![("chr1", 0, 3), ("chr1", 7, 10)])]
#[case::complete(vec![("chr1", 3, 7)], vec![("chr1", 0, 10)], vec![])]
#[case::no_overlap(vec![("chr1", 0, 5)], vec![("chr1", 10, 20)], vec![("chr1", 0, 5)])]
#[case::multi_chrom(vec![("chr1", 0, 10), ("chr2", 0, 10)], vec![("chr1", 5, 15)], vec![("chr1", 0, 5), ("chr2", 0, 10)])]
#[case::multiple_holes(
    vec![("chr1", 0, 20)],
    vec![("chr1", 2, 5), ("chr1", 8, 12), ("chr1", 15, 18)],
    vec![("chr1", 0, 2), ("chr1", 5, 8), ("chr1", 12, 15), ("chr1", 18, 20)]
)]
fn test_setdiff(#[case] a: Vec<Iv>, #[case] b: Vec<Iv>, #[case] expected: Vec<Iv>) {
    let result = make_regionset(a).setdiff(&make_regionset(b));
    assert_eq!(coords(&result), expected);
}

#[rstest]
#[case::overlapping(vec![("chr1", 0, 10)], vec![("chr1", 5, 15)], vec![("chr1", 5, 10)])]
#[case::contained(vec![("chr1", 0, 100)], vec![("chr1", 30, 70)], vec![("chr1", 30, 70)])]
// non-overlapping pairs yield a zero-width region at b's start
#[case::no_overlap(vec![("chr1", 0, 5)], vec![("chr1", 10, 20)], vec![("chr1", 10, 10)])]
// chrom mismatch yields zero-width at a's start
#[case::chrom_mismatch(vec![("chr1", 0, 10)], vec![("chr2", 0, 10)], vec![("chr1", 0, 0)])]
#[case::multiple_pairs(
    vec![("chr1", 0, 10), ("chr1", 20, 30), ("chr2", 0, 100)],
    vec![("chr1", 5, 15), ("chr1", 25, 35), ("chr2", 50, 60)],
    vec![("chr1", 5, 10), ("chr1", 25, 30), ("chr2", 50, 60)]
)]
fn test_pintersect(#[case] a: Vec<Iv>, #[case] b: Vec<Iv>, #[case] expected: Vec<Iv>) {
    let result = make_regionset(a).pintersect(&make_regionset(b));
    assert_eq!(coords(&result), expected);
}

// concat appends without sorting, merging, or deduplicating
#[rstest]
#[case::keeps_order_no_merge(
    vec![("chr1", 0, 10), ("chr1", 20, 30)],
    vec![("chr1", 5, 15), ("chr2", 50, 60)],
    vec![("chr1", 0, 10), ("chr1", 20, 30), ("chr1", 5, 15), ("chr2", 50, 60)]
)]
#[case::duplicates(vec![("chr1", 0, 10)], vec![("chr1", 0, 10)], vec![("chr1", 0, 10), ("chr1", 0, 10)])]
#[case::empty_self(vec![], vec![("chr1", 0, 10)], vec![("chr1", 0, 10)])]
#[case::empty_other(vec![("chr1", 0, 10)], vec![], vec![("chr1", 0, 10)])]
#[case::both_empty(vec![], vec![], vec![])]
fn test_concat(#[case] a: Vec<Iv>, #[case] b: Vec<Iv>, #[case] expected: Vec<Iv>) {
    let result = make_regionset(a).concat(&make_regionset(b));
    assert_eq!(coords(&result), expected);
}

#[rstest]
fn test_concat_into_matches_concat() {
    // Borrowing `concat` and consuming `concat_into` must produce identical
    // results, including order: self's regions first, then other's.
    let a_regions = vec![make_region("chr1", 100, 200), make_region("chr2", 50, 60)];
    let b_regions = vec![make_region("chr1", 150, 250), make_region("chr3", 10, 20)];

    let a = RegionSet::from(a_regions.clone());
    let b = RegionSet::from(b_regions.clone());
    let borrowed = a.concat(&b);

    let a2 = RegionSet::from(a_regions);
    let b2 = RegionSet::from(b_regions);
    let consumed = a2.concat_into(b2);

    assert_eq!(borrowed.regions, consumed.regions);
}

#[rstest]
fn test_concat_into_empty() {
    let non_empty = vec![make_region("chr1", 100, 200), make_region("chr2", 5, 15)];

    // empty + non-empty
    let empty = RegionSet::from(Vec::<Region>::new());
    let ne = RegionSet::from(non_empty.clone());
    let result = empty.concat_into(ne);
    assert_eq!(result.regions, non_empty);

    // non-empty + empty
    let ne2 = RegionSet::from(non_empty.clone());
    let empty2 = RegionSet::from(Vec::<Region>::new());
    let result2 = ne2.concat_into(empty2);
    assert_eq!(result2.regions, non_empty);
}

// Columns: jaccard (symmetric), coverage(a, b), overlap_coefficient (symmetric).
#[rstest]
#[case::identical(vec![("chr1", 0, 100)], vec![("chr1", 0, 100)], 1.0, 1.0, 1.0)]
#[case::disjoint(vec![("chr1", 0, 10)], vec![("chr1", 20, 30)], 0.0, 0.0, 0.0)]
// half-open touching intervals share no bases
#[case::adjacent(vec![("chr1", 0, 10)], vec![("chr1", 10, 20)], 0.0, 0.0, 0.0)]
#[case::different_chroms(vec![("chr1", 0, 100)], vec![("chr2", 0, 100)], 0.0, 0.0, 0.0)]
#[case::partial(vec![("chr1", 0, 10)], vec![("chr1", 5, 15)], 1.0 / 3.0, 0.5, 0.5)]
#[case::superset(vec![("chr1", 0, 100)], vec![("chr1", 20, 50)], 0.3, 0.3, 1.0)]
#[case::subset(vec![("chr1", 20, 50)], vec![("chr1", 0, 100)], 0.3, 1.0, 1.0)]
#[case::multi_chrom(vec![("chr1", 0, 20), ("chr2", 0, 10)], vec![("chr1", 10, 30), ("chr2", 5, 15)], 1.0 / 3.0, 0.5, 0.5)]
// a is reduced to [0,25) first
#[case::self_overlapping(vec![("chr1", 0, 15), ("chr1", 10, 25)], vec![("chr1", 5, 20)], 0.6, 0.6, 1.0)]
#[case::empty_self(vec![], vec![("chr1", 0, 100)], 0.0, 0.0, 0.0)]
#[case::empty_other(vec![("chr1", 0, 100)], vec![], 0.0, 0.0, 0.0)]
#[case::both_empty(vec![], vec![], 0.0, 0.0, 0.0)]
fn test_similarity_metrics(
    #[case] a: Vec<Iv>,
    #[case] b: Vec<Iv>,
    #[case] jaccard: f64,
    #[case] coverage: f64,
    #[case] overlap_coef: f64,
) {
    let (a, b) = (make_regionset(a), make_regionset(b));
    assert_approx(a.jaccard(&b), jaccard);
    assert_approx(b.jaccard(&a), jaccard);
    assert_approx(a.coverage(&b), coverage);
    assert_approx(a.overlap_coefficient(&b), overlap_coef);
    assert_approx(b.overlap_coefficient(&a), overlap_coef);
}

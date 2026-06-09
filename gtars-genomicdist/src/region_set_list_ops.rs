//! Batch pairwise Jaccard and indexed pair operations on collections of
//! region sets ([`RegionSetList`]).
//!
//! These let callers operate on pairs by index without cloning full RegionSets
//! across an FFI boundary (wasm, Python, R).

use gtars_core::models::IntervalSetOps;
use gtars_core::models::RegionSet;

/// Pairwise Jaccard similarity.
///
/// Pre-reduces each set once (N reduces total) and feeds the already-reduced
/// sets to each pairwise [`IntervalSetOps::jaccard`] call. Note `jaccard`
/// re-reduces its inputs internally, so this does NOT skip those reductions
/// entirely; the win is that `reduce()` on an already-reduced set is cheap and
/// the O(N^2) `union()` calls in the inner loop operate on reduced inputs.
/// Without the pre-reduction each raw set would be re-reduced ~N times across
/// the matrix instead of once. Returns a flat Vec<f64> of length N*N in
/// row-major order (symmetric matrix with 1.0 on the diagonal).
pub fn pairwise_jaccard(sets: &[RegionSet]) -> Vec<f64> {
    let n = sets.len();
    let mut matrix = vec![0.0f64; n * n];

    // Pre-reduce all sets once (see fn-level note on why this still helps even
    // though `jaccard` re-reduces internally).
    let reduced: Vec<RegionSet> = sets.iter().map(|s| s.reduce()).collect();

    // Diagonal
    for i in 0..n {
        matrix[i * n + i] = 1.0;
    }

    // Upper triangle
    for i in 0..n {
        for j in (i + 1)..n {
            let jac = reduced[i].jaccard(&reduced[j]);
            matrix[i * n + j] = jac;
            matrix[j * n + i] = jac;
        }
    }

    matrix
}

// --- Indexed operations on RegionSetList ---
//
// These let callers operate on pairs by index without cloning full RegionSets
// across an FFI boundary (wasm, Python, R).

use gtars_core::models::RegionSetList;

/// Indexed pair operations on a RegionSetList.
pub trait RegionSetListOps {
    fn pintersect_at(&self, i: usize, j: usize) -> Option<RegionSet>;
    fn pintersect_count(&self, i: usize, j: usize) -> Option<u32>;
    fn jaccard_at(&self, i: usize, j: usize) -> Option<f64>;
    fn union_at(&self, i: usize, j: usize) -> Option<RegionSet>;
    fn setdiff_at(&self, i: usize, j: usize) -> Option<RegionSet>;
    fn region_count(&self, i: usize) -> Option<u32>;
    fn union_except(&self, skip: usize) -> Option<RegionSet>;
    /// Compute union-of-all and all N union-except results in O(n) unions
    /// using prefix/suffix arrays. Returns (full_union, vec_of_union_except).
    fn bulk_union_except(&self) -> Option<(RegionSet, Vec<RegionSet>)>;
    /// Fold all sets into a single union.
    fn union_all(&self) -> Option<RegionSet>;
    /// Fold all sets into a single intersection.
    fn intersect_all(&self) -> Option<RegionSet>;
}

impl RegionSetListOps for RegionSetList {
    fn pintersect_at(&self, i: usize, j: usize) -> Option<RegionSet> {
        let a = self.get(i)?;
        let b = self.get(j)?;
        Some(a.pintersect(b))
    }

    fn pintersect_count(&self, i: usize, j: usize) -> Option<u32> {
        self.pintersect_at(i, j).map(|rs| rs.len() as u32)
    }

    fn jaccard_at(&self, i: usize, j: usize) -> Option<f64> {
        let a = self.get(i)?;
        let b = self.get(j)?;
        Some(a.jaccard(b))
    }

    fn union_at(&self, i: usize, j: usize) -> Option<RegionSet> {
        let a = self.get(i)?;
        let b = self.get(j)?;
        Some(a.union(b))
    }

    fn setdiff_at(&self, i: usize, j: usize) -> Option<RegionSet> {
        let a = self.get(i)?;
        let b = self.get(j)?;
        Some(a.setdiff(b))
    }

    fn region_count(&self, i: usize) -> Option<u32> {
        self.get(i).map(|rs| rs.len() as u32)
    }

    fn union_except(&self, skip: usize) -> Option<RegionSet> {
        let n = self.len();
        if n < 2 || skip >= n { return None; }
        let first = if skip == 0 { 1 } else { 0 };
        let mut acc = self.get(first)?.clone();
        for k in (first + 1)..n {
            if k == skip { continue; }
            if let Some(other) = self.get(k) {
                acc = acc.union(other);
            }
        }
        Some(acc)
    }

    fn bulk_union_except(&self) -> Option<(RegionSet, Vec<RegionSet>)> {
        let n = self.len();
        if n < 2 { return None; }

        // prefix[i] = union(set[0]..=set[i])
        let mut prefix = Vec::with_capacity(n);
        prefix.push(self.get(0)?.clone());
        for i in 1..n {
            let prev = &prefix[i - 1];
            prefix.push(prev.union(self.get(i)?));
        }

        // suffix[i] = union(set[i]..=set[n-1]), built incrementally from right
        let mut suffix = vec![None; n];
        suffix[n - 1] = Some(self.get(n - 1)?.clone());
        for i in (0..n - 1).rev() {
            // INVARIANT: the loop fills suffix from n-1 downward, so suffix[i+1]
            // was set on the previous iteration (or as the seed at n-1); it is
            // always Some here.
            suffix[i] = Some(self.get(i)?.union(suffix[i + 1].as_ref().unwrap()));
        }
        // INVARIANT: every index 0..n was assigned Some above (seed at n-1, the
        // loop covers 0..n-1), so each unwrap is sound.
        let suffix: Vec<RegionSet> = suffix.into_iter().map(|s| s.unwrap()).collect();

        let full_union = prefix[n - 1].clone();

        // union_except[i] = union(prefix[i-1], suffix[i+1])
        let mut results = Vec::with_capacity(n);
        for i in 0..n {
            let except = match (i > 0, i < n - 1) {
                (false, true) => suffix[1].clone(),
                (true, false) => prefix[i - 1].clone(),
                (true, true) => prefix[i - 1].union(&suffix[i + 1]),
                (false, false) => unreachable!(), // n >= 2
            };
            results.push(except);
        }

        Some((full_union, results))
    }

    fn union_all(&self) -> Option<RegionSet> {
        let n = self.len();
        if n == 0 { return None; }
        let mut acc = self.get(0)?.clone();
        for i in 1..n {
            if let Some(other) = self.get(i) {
                acc = acc.union(other);
            }
        }
        Some(acc)
    }

    fn intersect_all(&self) -> Option<RegionSet> {
        let n = self.len();
        if n == 0 { return None; }
        let mut acc = self.get(0)?.clone();
        for i in 1..n {
            if let Some(other) = self.get(i) {
                acc = acc.intersect(other);
            }
        }
        Some(acc)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gtars_core::models::IntervalSetOps;
    use gtars_core::models::{Region, RegionSetList};
    use pretty_assertions::assert_eq;
    use rstest::*;

    fn make_region(chr: &str, start: u32, end: u32) -> Region {
        Region {
            chr: chr.to_string(),
            start,
            end,
            rest: None,
        }
    }

    fn make_regionset(regions: Vec<(&str, u32, u32)>) -> RegionSet {
        let regions: Vec<Region> = regions
            .into_iter()
            .map(|(chr, start, end)| make_region(chr, start, end))
            .collect();
        RegionSet::from(regions)
    }

    // ── trim tests ──────────────────────────────────────────────────────


    #[rstest]
    fn test_pairwise_jaccard_matches_single() {
        // Pairwise matrix entries should match individual jaccard() calls
        let a = make_regionset(vec![("chr1", 0, 20), ("chr2", 0, 10)]);
        let b = make_regionset(vec![("chr1", 10, 30), ("chr2", 5, 15)]);
        let c = make_regionset(vec![("chr1", 50, 60)]);

        let matrix = super::pairwise_jaccard(&[a.clone(), b.clone(), c.clone()]);
        assert_eq!(matrix.len(), 9); // 3x3

        // Compare each pair against the single jaccard
        let eps = 1e-10;
        assert!((matrix[0 * 3 + 1] - a.jaccard(&b)).abs() < eps);
        assert!((matrix[0 * 3 + 2] - a.jaccard(&c)).abs() < eps);
        assert!((matrix[1 * 3 + 2] - b.jaccard(&c)).abs() < eps);
    }

    #[rstest]
    fn test_pairwise_jaccard_diagonal() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 50, 150)]);
        let matrix = super::pairwise_jaccard(&[a, b]);
        assert!((matrix[0] - 1.0).abs() < 1e-10); // [0,0]
        assert!((matrix[3] - 1.0).abs() < 1e-10); // [1,1]
    }

    #[rstest]
    fn test_pairwise_jaccard_symmetric() {
        let a = make_regionset(vec![("chr1", 0, 15), ("chr1", 30, 50)]);
        let b = make_regionset(vec![("chr1", 10, 40)]);
        let c = make_regionset(vec![("chr2", 0, 100)]);
        let matrix = super::pairwise_jaccard(&[a, b, c]);
        let n = 3;
        for i in 0..n {
            for j in 0..n {
                assert!((matrix[i * n + j] - matrix[j * n + i]).abs() < 1e-10);
            }
        }
    }

    #[rstest]
    fn test_pairwise_jaccard_empty_sets() {
        let a = RegionSet::from(Vec::<Region>::new());
        let b = RegionSet::from(Vec::<Region>::new());
        let matrix = super::pairwise_jaccard(&[a, b]);
        // Empty sets: intersection=0, union=0 -> jaccard=0
        // Diagonal: 1.0 (by convention)
        assert!((matrix[0] - 1.0).abs() < 1e-10);
        assert!((matrix[3] - 1.0).abs() < 1e-10);
        assert!(matrix[1].abs() < 1e-10);
        assert!(matrix[2].abs() < 1e-10);
    }

    #[rstest]
    fn test_pairwise_jaccard_single_set() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let matrix = super::pairwise_jaccard(&[a]);
        assert_eq!(matrix.len(), 1);
        assert!((matrix[0] - 1.0).abs() < 1e-10);
    }

    #[rstest]
    fn test_pairwise_jaccard_no_overlap_different_chroms() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr2", 0, 100)]);
        let matrix = super::pairwise_jaccard(&[a, b]);
        assert!(matrix[1].abs() < 1e-10); // [0,1]
        assert!(matrix[2].abs() < 1e-10); // [1,0]
    }

    #[rstest]
    fn test_pairwise_jaccard_identical_sets() {
        let a = make_regionset(vec![("chr1", 0, 50), ("chr2", 10, 30)]);
        let b = make_regionset(vec![("chr1", 0, 50), ("chr2", 10, 30)]);
        let matrix = super::pairwise_jaccard(&[a, b]);
        assert!((matrix[1] - 1.0).abs() < 1e-10); // identical -> jaccard=1
    }

    #[rstest]
    fn test_pairwise_jaccard_overlapping_input_regions() {
        // Input regions that need reducing (self-overlapping)
        let a = make_regionset(vec![("chr1", 0, 15), ("chr1", 10, 25)]); // reduces to [0,25)
        let b = make_regionset(vec![("chr1", 5, 20)]);
        let matrix = super::pairwise_jaccard(&[a.clone(), b.clone()]);
        let expected = a.jaccard(&b);
        assert!((matrix[1] - expected).abs() < 1e-10);
    }

    // ── RegionSetListOps tests ─────────────────────────────────────────


    fn make_rsl(sets: Vec<RegionSet>) -> RegionSetList {
        RegionSetList::from(sets)
    }


    #[rstest]
    fn test_rsl_pintersect_count() {
        let a = make_regionset(vec![("chr1", 0, 100), ("chr1", 200, 300)]);
        let b = make_regionset(vec![("chr1", 50, 150), ("chr1", 250, 350)]);
        let rsl = make_rsl(vec![a, b]);
        // a has 2 regions, b has 2 regions, both overlap each other
        let count = rsl.pintersect_count(0, 1).unwrap();
        assert_eq!(count, 2); // both pairs overlap
    }

    #[rstest]
    fn test_rsl_pintersect_count_no_overlap() {
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr1", 100, 200)]);
        let rsl = make_rsl(vec![a, b]);
        // paired by index: chr1:0-10 vs chr1:100-200 → no genomic overlap,
        // but pintersect produces a zero-width region (start=end) per pair
        let count = rsl.pintersect_count(0, 1).unwrap();
        assert_eq!(count, 1); // zero-width region still counted
    }

    #[rstest]
    fn test_rsl_pintersect_count_out_of_bounds() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let rsl = make_rsl(vec![a]);
        assert!(rsl.pintersect_count(0, 5).is_none());
    }

    #[rstest]
    fn test_rsl_jaccard_at() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 0, 100)]);
        let rsl = make_rsl(vec![a, b]);
        let j = rsl.jaccard_at(0, 1).unwrap();
        assert!((j - 1.0).abs() < 1e-9, "identical sets should have jaccard=1.0");
    }

    #[rstest]
    fn test_rsl_jaccard_at_disjoint() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 200, 300)]);
        let rsl = make_rsl(vec![a, b]);
        let j = rsl.jaccard_at(0, 1).unwrap();
        assert!((j - 0.0).abs() < 1e-9, "disjoint sets should have jaccard=0.0");
    }

    #[rstest]
    fn test_rsl_union_at() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 50, 150)]);
        let rsl = make_rsl(vec![a, b]);
        let u = rsl.union_at(0, 1).unwrap();
        assert_eq!(u.regions.len(), 1);
        assert_eq!(u.regions[0].start, 0);
        assert_eq!(u.regions[0].end, 150);
    }

    #[rstest]
    fn test_rsl_setdiff_at() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 50, 150)]);
        let rsl = make_rsl(vec![a, b]);
        let diff = rsl.setdiff_at(0, 1).unwrap();
        // a minus b: chr1:0-50
        assert_eq!(diff.regions.len(), 1);
        assert_eq!(diff.regions[0].start, 0);
        assert_eq!(diff.regions[0].end, 50);
    }

    #[rstest]
    fn test_rsl_region_count() {
        let a = make_regionset(vec![("chr1", 0, 100), ("chr1", 200, 300)]);
        let b = make_regionset(vec![("chr1", 50, 150)]);
        let rsl = make_rsl(vec![a, b]);
        assert_eq!(rsl.region_count(0).unwrap(), 2);
        assert_eq!(rsl.region_count(1).unwrap(), 1);
        assert!(rsl.region_count(5).is_none());
    }

    #[rstest]
    fn test_rsl_union_all() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 50, 200)]);
        let c = make_regionset(vec![("chr1", 150, 300)]);
        let rsl = make_rsl(vec![a, b, c]);
        let u = rsl.union_all().unwrap();
        assert_eq!(u.regions.len(), 1);
        assert_eq!(u.regions[0].start, 0);
        assert_eq!(u.regions[0].end, 300);
    }

    #[rstest]
    fn test_rsl_union_all_empty() {
        let rsl = make_rsl(vec![]);
        assert!(rsl.union_all().is_none());
    }

    #[rstest]
    fn test_rsl_union_all_single() {
        let a = make_regionset(vec![("chr1", 10, 50)]);
        let rsl = make_rsl(vec![a]);
        let u = rsl.union_all().unwrap();
        assert_eq!(u.regions.len(), 1);
        assert_eq!(u.regions[0].start, 10);
        assert_eq!(u.regions[0].end, 50);
    }

    #[rstest]
    fn test_rsl_intersect_all() {
        // Three overlapping sets — intersection is the region shared by all three
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 30, 200)]);
        let c = make_regionset(vec![("chr1", 60, 150)]);
        let rsl = make_rsl(vec![a, b, c]);
        let inter = rsl.intersect_all().unwrap();
        assert_eq!(inter.regions.len(), 1);
        assert_eq!(inter.regions[0].start, 60);
        assert_eq!(inter.regions[0].end, 100);
    }

    #[rstest]
    fn test_rsl_intersect_all_disjoint() {
        let a = make_regionset(vec![("chr1", 0, 50)]);
        let b = make_regionset(vec![("chr1", 100, 200)]);
        let rsl = make_rsl(vec![a, b]);
        let inter = rsl.intersect_all().unwrap();
        assert_eq!(inter.regions.len(), 0);
    }

    #[rstest]
    fn test_rsl_intersect_all_empty() {
        let rsl = make_rsl(vec![]);
        assert!(rsl.intersect_all().is_none());
    }

    #[rstest]
    fn test_rsl_intersect_all_different_sizes() {
        // This is the case where pintersect would give wrong results:
        // sets have different numbers of regions, but share genomic coverage
        let a = make_regionset(vec![("chr1", 0, 100), ("chr1", 200, 300)]);
        let b = make_regionset(vec![("chr1", 50, 250)]);
        let rsl = make_rsl(vec![a, b]);
        let inter = rsl.intersect_all().unwrap();
        // Shared coverage: [50,100) and [200,250)
        assert_eq!(inter.regions.len(), 2);
        assert_eq!(inter.regions[0].start, 50);
        assert_eq!(inter.regions[0].end, 100);
        assert_eq!(inter.regions[1].start, 200);
        assert_eq!(inter.regions[1].end, 250);
    }

    #[rstest]
    fn test_rsl_union_except() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 200, 300)]);
        let c = make_regionset(vec![("chr1", 400, 500)]);
        let rsl = make_rsl(vec![a, b, c]);
        // union_except(1) = union of sets 0 and 2 (skip set 1)
        let ue = rsl.union_except(1).unwrap();
        assert_eq!(ue.regions.len(), 2);
        assert_eq!(ue.regions[0].start, 0);
        assert_eq!(ue.regions[0].end, 100);
        assert_eq!(ue.regions[1].start, 400);
        assert_eq!(ue.regions[1].end, 500);
    }

    #[rstest]
    fn test_rsl_union_except_too_small() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let rsl = make_rsl(vec![a]);
        assert!(rsl.union_except(0).is_none());
    }

    #[rstest]
    fn test_rsl_bulk_union_except_n2() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 200, 300)]);
        let rsl = make_rsl(vec![a, b]);
        let (full_union, excepts) = rsl.bulk_union_except().unwrap();

        // Full union covers both regions
        assert_eq!(full_union.regions.len(), 2);

        // except[0] = union of everything except set 0 = set 1
        assert_eq!(excepts.len(), 2);
        assert_eq!(excepts[0].regions.len(), 1);
        assert_eq!(excepts[0].regions[0].start, 200);
        assert_eq!(excepts[0].regions[0].end, 300);

        // except[1] = union of everything except set 1 = set 0
        assert_eq!(excepts[1].regions.len(), 1);
        assert_eq!(excepts[1].regions[0].start, 0);
        assert_eq!(excepts[1].regions[0].end, 100);
    }

    #[rstest]
    fn test_rsl_bulk_union_except_n3() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 200, 300)]);
        let c = make_regionset(vec![("chr1", 400, 500)]);
        let rsl = make_rsl(vec![a, b, c]);
        let (full_union, excepts) = rsl.bulk_union_except().unwrap();

        assert_eq!(full_union.regions.len(), 3);
        assert_eq!(excepts.len(), 3);

        // except[0] = union(b, c)
        assert_eq!(excepts[0].regions.len(), 2);
        assert_eq!(excepts[0].regions[0].start, 200);
        assert_eq!(excepts[0].regions[1].start, 400);

        // except[1] = union(a, c)
        assert_eq!(excepts[1].regions.len(), 2);
        assert_eq!(excepts[1].regions[0].start, 0);
        assert_eq!(excepts[1].regions[1].start, 400);

        // except[2] = union(a, b)
        assert_eq!(excepts[2].regions.len(), 2);
        assert_eq!(excepts[2].regions[0].start, 0);
        assert_eq!(excepts[2].regions[1].start, 200);
    }

    #[rstest]
    fn test_rsl_bulk_union_except_too_small() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let rsl = make_rsl(vec![a]);
        assert!(rsl.bulk_union_except().is_none());

        let rsl_empty = make_rsl(vec![]);
        assert!(rsl_empty.bulk_union_except().is_none());
    }

    #[rstest]
    fn test_rsl_bulk_union_except_matches_union_except() {
        // Verify bulk algorithm produces same results as individual union_except calls
        let a = make_regionset(vec![("chr1", 0, 100), ("chr2", 50, 200)]);
        let b = make_regionset(vec![("chr1", 80, 180), ("chr2", 100, 300)]);
        let c = make_regionset(vec![("chr1", 150, 250)]);
        let d = make_regionset(vec![("chr2", 0, 150)]);
        let rsl = make_rsl(vec![a, b, c, d]);

        let (_, bulk_excepts) = rsl.bulk_union_except().unwrap();

        for i in 0..4 {
            let individual = rsl.union_except(i).unwrap();
            assert_eq!(
                bulk_excepts[i].regions.len(),
                individual.regions.len(),
                "region count mismatch at index {}",
                i
            );
            for (j, (bulk_r, indiv_r)) in bulk_excepts[i]
                .regions
                .iter()
                .zip(individual.regions.iter())
                .enumerate()
            {
                assert_eq!(
                    bulk_r.chr, indiv_r.chr,
                    "chr mismatch at except[{}][{}]", i, j
                );
                assert_eq!(
                    bulk_r.start, indiv_r.start,
                    "start mismatch at except[{}][{}]", i, j
                );
                assert_eq!(
                    bulk_r.end, indiv_r.end,
                    "end mismatch at except[{}][{}]", i, j
                );
            }
        }
    }
}

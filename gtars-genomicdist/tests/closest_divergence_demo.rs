//! TEMPORARY demonstration test — referenced from a PR review comment.
//! DELETE after review.
//!
//! Shows that the new gtars-core inherent `RegionSet::closest` and the existing
//! gtars-genomicdist `IntervalRanges::closest` return DIFFERENT `other` indices
//! for the same nearest region when `other` is not already sorted by start.
//!
//! Both find the correct nearest *distance*, but the core copy reports an index
//! into a sorted clone it builds internally, so the index maps to the wrong
//! region in the caller's `other`. The genomicdist version preserves the
//! original index.
//!
//! Run: cargo test -p gtars-genomicdist --test closest_divergence_demo -- --nocapture

use gtars_core::models::{Region, RegionSet};
use gtars_genomicdist::IntervalRanges;

fn r(chr: &str, start: u32, end: u32) -> Region {
    Region { chr: chr.into(), start, end, rest: None }
}

#[test]
fn closest_index_divergence() {
    // self: a single region near position 100
    let me = RegionSet::from(vec![r("chr1", 100, 110)]);

    // other, deliberately NOT in start order:
    //   original index 0 -> 500..510  (far)
    //   original index 1 -> 120..130  (NEAREST to self, gap = 10)
    //   original index 2 -> 900..910  (far)
    let other = RegionSet::from(vec![
        r("chr1", 500, 510), // idx 0
        r("chr1", 120, 130), // idx 1  <- the true nearest
        r("chr1", 900, 910), // idx 2
    ]);

    // gtars-core inherent method (inherent wins over the trait when both are in scope)
    let core = RegionSet::closest(&me, &other);
    // gtars-genomicdist trait method, called via UFCS to disambiguate
    let gd = IntervalRanges::closest(&me, &other);

    let (_, core_idx, core_dist) = core[0];
    let (_, gd_idx, gd_dist) = gd[0];

    let core_region = &other.regions[core_idx];
    let gd_region = &other.regions[gd_idx];

    println!("\n--- self region: {:?} ---", me.regions[0]);
    println!("core  -> other_idx={core_idx} dist={core_dist}  => other.regions[{core_idx}] = {:?}", core_region);
    println!("gdist -> other_idx={gd_idx} dist={gd_dist}  => other.regions[{gd_idx}] = {:?}", gd_region);

    // Both agree on the DISTANCE (10) ...
    assert_eq!(core_dist, 10, "core distance");
    assert_eq!(gd_dist, 10, "genomicdist distance");

    // ... but the returned indices map to DIFFERENT regions when fed back into `other`.
    // genomicdist returns the original index (1 -> 120..130, the true nearest).
    assert_eq!(gd_idx, 1, "genomicdist returns ORIGINAL index of the nearest");
    assert_eq!(*gd_region, r("chr1", 120, 130));

    // core returns the index into a SORTED copy (0 -> nearest is first after sort),
    // which maps to the WRONG region in the caller's `other`.
    assert_eq!(core_idx, 0, "core returns SORTED-copy index");
    assert_eq!(*core_region, r("chr1", 500, 510), "core index maps to the wrong region!");

    // Proof they diverge:
    assert_ne!(core_idx, gd_idx);
    assert_ne!(core_region, gd_region);
}

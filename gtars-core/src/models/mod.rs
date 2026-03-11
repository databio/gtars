pub mod bed_set;
pub mod fragments;
pub mod interval;
pub mod region;
pub mod region_set;

// re-export for cleaner imports
pub use self::fragments::Fragment;
pub use self::interval::Interval;
pub use self::region::Region;
pub use self::region_set::{
    sweep_intersect_chr, sweep_setdiff_chr, IntervalSetOps, RegionSet, SortedRegionSet,
};

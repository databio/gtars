pub mod bed_set;
pub mod fragments;
pub mod region;
pub mod region_set;
pub mod tokenized_region;
pub mod tokenized_regionset;
pub mod universe;
pub mod interval;

// re-export for cleaner imports
pub use self::fragments::Fragment;
pub use self::region::Region;
pub use self::region_set::RegionSet;
pub use self::interval::Interval;
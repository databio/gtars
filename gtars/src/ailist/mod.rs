//!
//! # Augmented Interval List: a novel data structure for efficient genomic interval search
//! This is a rust implementation of the Augmented Interval List (AIList): [https://academic.oup.com/bioinformatics/article/35/23/4907/5509521](https://academic.oup.com/bioinformatics/article/35/23/4907/5509521).
//!
//! The Augmented Interval List (AIList), enumerates intersections between a query interval q and an interval set R.
//!
//! It should be complete, but has not been rigorously tested.
//!
pub mod core;

// re-expose models
pub use core::{AIList, Interval};

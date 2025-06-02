//! This module provides an implementation of the Binary Interval Tree Search (BITS) algorithm.
//! 
//! [BITS](https://arxiv.org/pdf/1208.3407.pdf) algorithm
//! - Parallel friendly. Queries are on an immutable structure, even for seek
//! - Consumer / Adapter paradigm, Iterators are returned and serve as the main API for interacting
//!   with the lapper
//!
//! ## Details:
//!
//! ```text
//!       0  1  2  3  4  5  6  7  8  9  10 11
//! (0,10]X  X  X  X  X  X  X  X  X  X
//! (2,5]       X  X  X
//! (3,8]          X  X  X  X  X
//! (3,8]          X  X  X  X  X
//! (3,8]          X  X  X  X  X
//! (3,8]          X  X  X  X  X
//! (5,9]                X  X  X  X
//! (8,11]                        X  X  X
//!
//! Query: (8, 11]
//! Answer: ((0,10], (5,9], (8,11])
//! ```
//!
//! Most interaction with this crate will be through the [`Lapper`](struct.Lapper.html) struct
//! The main methods are [`find`](struct.Lapper.html#method.find),
//! [`seek`](struct.Lapper.html#method.seek), and [`count`](struct.Lapper.html#method.count)
//! where both `seek` and `count` are special cases allowing for very fast queries in certain scenarios.
//!
//! The overlap function for this assumes a zero based genomic coordinate system. So [start, end)
//! is not inclusive of the end position for neither the queries, nor the Intervals.
//!
//! Lapper does not use an interval tree, instead, it operates on the assumtion that most intervals are
//! of similar length; or, more exactly, that the longest interval in the set is not long compred to
//! the average distance between intervals.
//!
//! For cases where this holds true (as it often does with genomic data), we can sort by start and
//! use binary search on the starts, accounting for the length of the longest interval. The advantage
//! of this approach is simplicity of implementation and speed. In realistic tests queries returning
//! the overlapping intervals are 1000 times faster than brute force and queries that merely check
//! for the overlaps are > 5000 times faster.
//!
//! When this is not the case, if possible in your scenario, use merge_overlaps first, and then use
//! `find` or `seek`. The `count` method will be fast in all scenarios.
//!
//! # Examples
//!
//! ```rust
//!    use rust_lapper::{Interval, Lapper};
//!    use std::cmp;
//!    type Iv = Interval<usize, u32>;
//!
//!    // create some fake data
//!    let data: Vec<Iv> = (0..20).step_by(5).map(|x| Iv{start: x, end: x + 2, val: 0}).collect();
//!    println!("{:#?}", data);
//!
//!    // make lapper structure
//!    let laps = Lapper::new(data);
//!
//!    assert_eq!(laps.find(6, 11).next(), Some(&Iv{start: 5, end: 7, val: 0}));
//!
//!    // Demonstration of seek function. By passing in the &mut cursor, seek can have thread local
//!    // cursors going
//!    let mut sim: usize = 0;
//!    let mut cursor = 0;
//!    // Calculate the overlap between the query and the found intervals, sum total overlap
//!    for i in (0..10).step_by(3) {
//!        sim += laps
//!            .seek(i, i + 2, &mut cursor)
//!            .map(|iv| cmp::min(i + 2, iv.end) - cmp::max(i, iv.start))
//!            .sum::<usize>();
//!    }
//!    assert_eq!(sim, 4);
//! ```
pub mod core;

pub use core::*;
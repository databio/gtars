//! This module provides a simple data-structure for fast interval searches.
//! ## Features
//! - Extremely fast on most genomic datasets. (3-4x faster than other methods)
//! - Extremely fast on in order queries. (10x faster than other methods)
//! - Extremely fast intersections count method based on the
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
use num_traits::{
    identities::{one, zero},
    PrimInt, Unsigned,
};
use crate::tokenizers::intersect::intervals::Interval;
use std::collections::VecDeque;


/// Primary object of the library. The public intervals holds all the intervals and can be used for
/// iterating / pulling values out of the tree.
#[derive(Debug, Clone)]
pub struct Lapper<I, T>
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    /// List of intervals
    pub intervals: Vec<Interval<I, T>>,
    /// Sorted list of start positions,
    starts: Vec<I>,
    /// Sorted list of end positions,
    ends: Vec<I>,
    /// The length of the longest interval
    max_len: I,
    /// The calculated number of positions covered by the intervals
    cov: Option<I>,
    /// Whether or not overlaps have been merged
    pub overlaps_merged: bool,
}



impl<I, T> Lapper<I, T>
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    /// Create a new instance of Lapper by passing in a vector of Intervals. This vector will
    /// immediately be sorted by start order.
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let data = (0..20).step_by(5)
    ///                   .map(|x| Interval{start: x, end: x + 10, val: true})
    ///                   .collect::<Vec<Interval<usize, bool>>>();
    /// let lapper = Lapper::new(data);
    /// ```
    pub fn new(mut intervals: Vec<Interval<I, T>>) -> Self {
        intervals.sort();
        let (mut starts, mut ends): (Vec<_>, Vec<_>) =
            intervals.iter().map(|x| (x.start, x.end)).unzip();
        starts.sort();
        ends.sort();
        let mut max_len = zero::<I>();
        for interval in intervals.iter() {
            let i_len = interval
                .end
                .checked_sub(&interval.start)
                .unwrap_or_else(zero::<I>);
            if i_len > max_len {
                max_len = i_len;
            }
        }
        Lapper {
            intervals,
            starts,
            ends,
            max_len,
            cov: None,
            overlaps_merged: false,
        }
    }

    /// Insert a new interval after the Lapper has been created. This is very
    /// inefficient and should be avoided if possible.
    ///
    /// SIDE EFFECTS: This clears cov() and overlaps_merged
    /// meaning that those will have to be recomputed after a insert
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let data : Vec<Interval<usize, usize>>= vec!{
    ///     Interval{start:0,  end:5,  val:1},
    ///     Interval{start:6,  end:10, val:2},
    /// };
    /// let mut lapper = Lapper::new(data);
    /// lapper.insert(Interval{start:0, end:20, val:5});
    /// assert_eq!(lapper.len(), 3);
    /// assert_eq!(lapper.find(1,3).collect::<Vec<&Interval<usize,usize>>>(),
    ///     vec![
    ///         &Interval{start:0, end:5, val:1},
    ///         &Interval{start:0, end:20, val:5},
    ///     ]
    /// );
    ///
    /// ```
    pub fn insert(&mut self, elem: Interval<I, T>) {
        let starts_insert_index = Self::bsearch_seq(elem.start, &self.starts);
        let ends_insert_index = Self::bsearch_seq(elem.end, &self.ends);
        let intervals_insert_index = Self::bsearch_seq_ref(&elem, &self.intervals);
        let i_len = elem.end.checked_sub(&elem.start).unwrap_or_else(zero::<I>);
        if i_len > self.max_len {
            self.max_len = i_len;
        }
        self.starts.insert(starts_insert_index, elem.start);
        self.ends.insert(ends_insert_index, elem.end);
        self.intervals.insert(intervals_insert_index, elem);
        self.cov = None;
        self.overlaps_merged = false;
    }

    /// Get the number over intervals in Lapper
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let data = (0..20).step_by(5)
    ///                   .map(|x| Interval{start: x, end: x + 10, val: true})
    ///                   .collect::<Vec<Interval<usize, bool>>>();
    /// let lapper = Lapper::new(data);
    /// assert_eq!(lapper.len(), 4);
    /// ```
    #[inline]
    pub fn len(&self) -> usize {
        self.intervals.len()
    }

    /// Check if lapper is empty
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let data: Vec<Interval<usize, bool>> = vec![];
    /// let lapper = Lapper::new(data);
    /// assert_eq!(lapper.is_empty(), true);
    /// ```
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.intervals.is_empty()
    }

    /// Get the number of positions covered by the intervals in Lapper. This provides immutable
    /// access if it has already been set, or on the fly calculation.
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let data = (0..20).step_by(5)
    ///                   .map(|x| Interval{start: x, end: x + 10, val: true})
    ///                   .collect::<Vec<Interval<usize, bool>>>();
    /// let lapper = Lapper::new(data);
    /// assert_eq!(lapper.cov(), 25);
    #[inline]
    pub fn cov(&self) -> I {
        match self.cov {
            None => self.calculate_coverage(),
            Some(cov) => cov,
        }
    }

    /// Get the number fo positions covered by the intervals in Lapper and store it. If you are
    /// going to be using the coverage, you should set it to avoid calculating it over and over.
    pub fn set_cov(&mut self) -> I {
        let cov = self.calculate_coverage();
        self.cov = Some(cov);
        cov
    }

    /// Calculate the actual coverage behind the scenes.
    fn calculate_coverage(&self) -> I {
        let mut moving_interval = Interval {
            start: zero::<I>(),
            end: zero::<I>(),
            val: zero::<I>(),
        };
        let mut cov = zero::<I>();

        for interval in self.intervals.iter() {
            // If it overlaps, embrace, extend, extinguish
            if moving_interval.overlap(interval.start, interval.end) {
                moving_interval.start = std::cmp::min(moving_interval.start, interval.start);
                moving_interval.end = std::cmp::max(moving_interval.end, interval.end);
            } else {
                // add the set and move on
                cov = cov + (moving_interval.end - moving_interval.start);
                moving_interval.start = interval.start;
                moving_interval.end = interval.end;
            }
        }
        // add in the last bit
        cov = cov + (moving_interval.end - moving_interval.start);
        cov
    }

    /// Return an iterator over the intervals in Lapper
    #[inline]
    pub fn iter(&self) -> IterLapper<I, T> {
        IterLapper {
            inner: self,
            pos: 0,
        }
    }

    /// Merge any intervals that overlap with eachother within the Lapper. This is an easy way to
    /// speed up queries.
    pub fn merge_overlaps(&mut self) {
        let mut stack: VecDeque<&mut Interval<I, T>> = VecDeque::new();
        let mut ivs = self.intervals.iter_mut();
        if let Some(first) = ivs.next() {
            stack.push_back(first);
            for interval in ivs {
                let top = stack.pop_back().unwrap();
                if top.end < interval.start {
                    stack.push_back(top);
                    stack.push_back(interval);
                } else if top.end < interval.end {
                    top.end = interval.end;
                    //stack.pop_back();
                    stack.push_back(top);
                } else {
                    // they were equal
                    stack.push_back(top);
                }
            }
            self.overlaps_merged = true;
            self.intervals = stack
                .into_iter()
                .map(|x| Interval {
                    start: x.start,
                    end: x.end,
                    val: x.val.clone(),
                })
                .collect();
        }
        // Fix the starts and ends used by counts
        let (mut starts, mut ends): (Vec<_>, Vec<_>) =
            self.intervals.iter().map(|x| (x.start, x.end)).unzip();
        starts.sort();
        ends.sort();
        self.starts = starts;
        self.ends = ends;
        self.max_len = self
            .intervals
            .iter()
            .map(|x| x.end.checked_sub(&x.start).unwrap_or_else(zero::<I>))
            .max()
            .unwrap_or_else(zero::<I>);
    }

    /// Determine the first index that we should start checking for overlaps for via a binary
    /// search.
    /// Assumes that the maximum interval length in `intervals` has been subtracted from
    /// `start`, otherwise the result is undefined
    #[inline]
    pub fn lower_bound(start: I, intervals: &[Interval<I, T>]) -> usize {
        let mut size = intervals.len();
        let mut low = 0;

        while size > 0 {
            let half = size / 2;
            let other_half = size - half;
            let probe = low + half;
            let other_low = low + other_half;
            let v = &intervals[probe];
            size = half;
            low = if v.start < start { other_low } else { low }
        }
        low
    }

    #[inline]
    pub fn bsearch_seq<K>(key: K, elems: &[K]) -> usize
    where
        K: PartialEq + PartialOrd,
    {
        Self::bsearch_seq_ref(&key, elems)
    }

    #[inline]
    pub fn bsearch_seq_ref<K>(key: &K, elems: &[K]) -> usize
    where
        K: PartialEq + PartialOrd,
    {
        if elems.is_empty() {
            return 0;
        }
        if elems[0] > *key {
            return 0;
        }
        let mut high = elems.len();
        let mut low = 0;

        while high - low > 1 {
            let mid = (high + low) / 2;
            if elems[mid] < *key {
                low = mid;
            } else {
                high = mid;
            }
        }
        high
    }

    /// Find the union and the intersect of two lapper objects.
    /// Union: The set of positions found in both lappers
    /// Intersect: The number of positions where both lappers intersect. Note that a position only
    /// counts one time, multiple Intervals covering the same position don't add up.
    /// ``` rust
    /// use rust_lapper::{Lapper, Interval};
    /// type Iv = Interval<u32, u32>;
    /// let data1: Vec<Iv> = vec![
    ///     Iv{start: 70, end: 120, val: 0}, // max_len = 50
    ///     Iv{start: 10, end: 15, val: 0}, // exact overlap
    ///     Iv{start: 12, end: 15, val: 0}, // inner overlap
    ///     Iv{start: 14, end: 16, val: 0}, // overlap end
    ///     Iv{start: 68, end: 71, val: 0}, // overlap start
    /// ];
    /// let data2: Vec<Iv> = vec![
    ///
    ///     Iv{start: 10, end: 15, val: 0},
    ///     Iv{start: 40, end: 45, val: 0},
    ///     Iv{start: 50, end: 55, val: 0},
    ///     Iv{start: 60, end: 65, val: 0},
    ///     Iv{start: 70, end: 75, val: 0},
    /// ];
    ///
    /// let (mut lapper1, mut lapper2) = (Lapper::new(data1), Lapper::new(data2)) ;
    /// // Should be the same either way it's calculated
    /// let (union, intersect) = lapper1.union_and_intersect(&lapper2);
    /// assert_eq!(intersect, 10);
    /// assert_eq!(union, 73);
    /// let (union, intersect) = lapper2.union_and_intersect(&lapper1);
    /// assert_eq!(intersect, 10);
    /// assert_eq!(union, 73);
    /// lapper1.merge_overlaps();
    /// lapper1.set_cov();
    /// lapper2.merge_overlaps();
    /// lapper2.set_cov();
    ///
    /// // Should be the same either way it's calculated
    /// let (union, intersect) = lapper1.union_and_intersect(&lapper2);
    /// assert_eq!(intersect, 10);
    /// assert_eq!(union, 73);
    /// let (union, intersect) = lapper2.union_and_intersect(&lapper1);
    /// assert_eq!(intersect, 10);
    /// assert_eq!(union, 73);
    /// ```
    #[inline]
    pub fn union_and_intersect(&self, other: &Self) -> (I, I) {
        let mut cursor: usize = 0;

        if !self.overlaps_merged || !other.overlaps_merged {
            let mut intersections: Vec<Interval<I, bool>> = vec![];
            for self_iv in self.iter() {
                for other_iv in other.seek(self_iv.start, self_iv.end, &mut cursor) {
                    let start = std::cmp::max(self_iv.start, other_iv.start);
                    let end = std::cmp::min(self_iv.end, other_iv.end);
                    intersections.push(Interval {
                        start,
                        end,
                        val: true,
                    });
                }
            }
            let mut temp_lapper = Lapper::new(intersections);
            temp_lapper.merge_overlaps();
            temp_lapper.set_cov();
            let union = self.cov() + other.cov() - temp_lapper.cov();
            (union, temp_lapper.cov())
        } else {
            let mut intersect = zero::<I>();
            for c1_iv in self.iter() {
                for c2_iv in other.seek(c1_iv.start, c1_iv.end, &mut cursor) {
                    let local_intersect = c1_iv.intersect(c2_iv);
                    intersect = intersect + local_intersect;
                }
            }
            let union = self.cov() + other.cov() - intersect;
            (union, intersect)
        }
    }

    /// Find the intersect of two lapper objects.
    /// Intersect: The number of positions where both lappers intersect. Note that a position only
    /// counts one time, multiple Intervals covering the same position don't add up
    #[inline]
    pub fn intersect(&self, other: &Self) -> I {
        self.union_and_intersect(other).1
    }

    /// Find the union of two lapper objects.
    #[inline]
    pub fn union(&self, other: &Self) -> I {
        self.union_and_intersect(other).0
    }

    /// Return the contiguous intervals of coverage, `val` represents the number of intervals
    /// covering the returned interval.
    ///
    /// # Examples
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let data = (0..20).step_by(5)
    ///                   .map(|x| Interval{start: x, end: x + 10, val: true})
    ///                   .collect::<Vec<Interval<usize, bool>>>();
    /// let lapper = Lapper::new(data);
    /// assert_eq!(lapper.depth().collect::<Vec<Interval<usize, usize>>>(), vec![
    ///             Interval { start: 0, end: 5, val: 1 },
    ///             Interval { start: 5, end: 20, val: 2 },
    ///             Interval { start: 20, end: 25, val: 1 }]);
    /// ```
    #[inline]
    pub fn depth(&self) -> IterDepth<I, T> {
        let mut merged_lapper = Lapper::new(
            self.intervals
                .iter()
                .map(|i| Interval {
                    start: i.start,
                    end: i.end,
                    val: true,
                })
                .collect::<Vec<Interval<I, bool>>>(),
        );
        merged_lapper.merge_overlaps();
        let merged_len = merged_lapper.intervals.len();
        IterDepth {
            inner: self,
            merged: merged_lapper,
            curr_merged_pos: zero::<I>(),
            curr_pos: 0,
            cursor: 0,
            end: merged_len,
        }
    }

    /// Count all intervals that overlap start .. end. This performs two binary search in order to
    /// find all the excluded elements, and then deduces the intersection from there. See
    /// [BITS](https://arxiv.org/pdf/1208.3407.pdf) for more details.
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let lapper = Lapper::new((0..100).step_by(5)
    ///                                 .map(|x| Interval{start: x, end: x+2 , val: true})
    ///                                 .collect::<Vec<Interval<usize, bool>>>());
    /// assert_eq!(lapper.count(5, 11), 2);
    /// ```
    #[inline]
    pub fn count(&self, start: I, end: I) -> usize {
        let len = self.intervals.len();
        let mut first = Self::bsearch_seq(start, &self.ends);
        let last = Self::bsearch_seq(end, &self.starts);
        while first < len && self.ends[first] == start {
            first += 1;
        }
        let num_cant_after = len - last;
        len - first - num_cant_after
    }

    /// Find all intervals that overlap start .. end
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let lapper = Lapper::new((0..100).step_by(5)
    ///                                 .map(|x| Interval{start: x, end: x+2 , val: true})
    ///                                 .collect::<Vec<Interval<usize, bool>>>());
    /// assert_eq!(lapper.find(5, 11).count(), 2);
    /// ```
    #[inline]
    pub fn find(&self, start: I, end: I) -> IterFind<I, T> {
        IterFind {
            inner: self,
            off: Self::lower_bound(
                start.checked_sub(&self.max_len).unwrap_or_else(zero::<I>),
                &self.intervals,
            ),
            start,
            end,
        }
    }

    /// Find all intevals that overlap start .. end. This method will work when queries
    /// to this lapper are in sorted (start) order. It uses a linear search from the last query
    /// instead of a binary search. A reference to a cursor must be passed in. This reference will
    /// be modified and should be reused in the next query. This allows seek to not need to make
    /// the lapper object mutable, and thus use the same lapper accross threads.
    /// ```
    /// use rust_lapper::{Lapper, Interval};
    /// let lapper = Lapper::new((0..100).step_by(5)
    ///                                 .map(|x| Interval{start: x, end: x+2 , val: true})
    ///                                 .collect::<Vec<Interval<usize, bool>>>());
    /// let mut cursor = 0;
    /// for i in lapper.iter() {
    ///    assert_eq!(lapper.seek(i.start, i.end, &mut cursor).count(), 1);
    /// }
    /// ```
    #[inline]
    pub fn seek<'a>(&'a self, start: I, end: I, cursor: &mut usize) -> IterFind<'a, I, T> {
        if *cursor == 0 || (*cursor < self.intervals.len() && self.intervals[*cursor].start > start)
        {
            *cursor = Self::lower_bound(
                start.checked_sub(&self.max_len).unwrap_or_else(zero::<I>),
                &self.intervals,
            );
        }

        while *cursor + 1 < self.intervals.len()
            && self.intervals[*cursor + 1].start
                < start.checked_sub(&self.max_len).unwrap_or_else(zero::<I>)
        {
            *cursor += 1;
        }

        IterFind {
            inner: self,
            off: *cursor,
            start,
            end,
        }
    }
}

/// Find Iterator
#[derive(Debug)]
pub struct IterFind<'a, I, T>
where
    T: Eq + Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
{
    inner: &'a Lapper<I, T>,
    off: usize,
    start: I,
    end: I,
}

impl<'a, I, T> Iterator for IterFind<'a, I, T>
where
    T: Eq + Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
{
    type Item = &'a Interval<I, T>;

    #[inline]
    // interval.start < end && interval.end > start
    fn next(&mut self) -> Option<Self::Item> {
        while self.off < self.inner.intervals.len() {
            //let mut generator = self.inner.intervals[self.off..].iter();
            //while let Some(interval) = generator.next() {
            let interval = &self.inner.intervals[self.off];
            self.off += 1;
            if interval.overlap(self.start, self.end) {
                return Some(interval);
            } else if interval.start >= self.end {
                break;
            }
        }
        None
    }
}

/// Depth Iterator
#[derive(Debug)]
pub struct IterDepth<'a, I, T>
where
    T: Eq + Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
{
    inner: &'a Lapper<I, T>,
    merged: Lapper<I, bool>, // A lapper that is the merged_lapper of inner
    curr_merged_pos: I,      // Current start position in current interval
    curr_pos: usize,         // In merged list of non-overlapping intervals
    cursor: usize,           // cursor for seek over inner lapper
    end: usize,              // len of merged
}

impl<'a, I, T> Iterator for IterDepth<'a, I, T>
where
    T: Eq + Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
{
    type Item = Interval<I, I>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let mut interval: &Interval<I, bool> = &self.merged.intervals[self.curr_pos];
        if self.curr_merged_pos == zero::<I>() {
            self.curr_merged_pos = interval.start;
        }
        if interval.end == self.curr_merged_pos {
            if self.curr_pos + 1 != self.end {
                self.curr_pos += 1;
                interval = &self.merged.intervals[self.curr_pos];
                self.curr_merged_pos = interval.start;
            } else {
                return None;
            }
        }
        let start = self.curr_merged_pos;
        let depth_at_point = self
            .inner
            .seek(
                self.curr_merged_pos,
                self.curr_merged_pos + one::<I>(),
                &mut self.cursor,
            )
            .count();
        let mut new_depth_at_point = depth_at_point;
        while new_depth_at_point == depth_at_point && self.curr_merged_pos < interval.end {
            self.curr_merged_pos = self.curr_merged_pos + one::<I>();
            new_depth_at_point = self
                .inner
                .seek(
                    self.curr_merged_pos,
                    self.curr_merged_pos + one::<I>(),
                    &mut self.cursor,
                )
                .count();
        }
        Some(Interval {
            start,
            end: self.curr_merged_pos,
            val: I::from(depth_at_point).unwrap(), // from usize should always work
        })
    }
}
/// Lapper Iterator
pub struct IterLapper<'a, I, T>
where
    T: Eq + Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
{
    inner: &'a Lapper<I, T>,
    pos: usize,
}

impl<'a, I, T> Iterator for IterLapper<'a, I, T>
where
    T: Eq + Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
{
    type Item = &'a Interval<I, T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.inner.intervals.len() {
            None
        } else {
            self.pos += 1;
            self.inner.intervals.get(self.pos - 1)
        }
    }
}

impl<I, T> IntoIterator for Lapper<I, T>
where
    T: Eq + Clone + Send + Sync,
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
{
    type Item = Interval<I, T>;
    type IntoIter = ::std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.intervals.into_iter()
    }
}

impl<'a, I, T> IntoIterator for &'a Lapper<I, T>
where
    T: Eq + Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
{
    type Item = &'a Interval<I, T>;
    type IntoIter = std::slice::Iter<'a, Interval<I, T>>;

    fn into_iter(self) -> std::slice::Iter<'a, Interval<I, T>> {
        self.intervals.iter()
    }
}

impl<'a, I, T> IntoIterator for &'a mut Lapper<I, T>
where
    T: Eq + Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
{
    type Item = &'a mut Interval<I, T>;
    type IntoIter = std::slice::IterMut<'a, Interval<I, T>>;

    fn into_iter(self) -> std::slice::IterMut<'a, Interval<I, T>> {
        self.intervals.iter_mut()
    }
}

#[cfg(test)]
#[rustfmt::skip]
mod tests {
    use super::*;

    type Iv = Interval<usize, u32>;
    fn setup_nonoverlapping() -> Lapper<usize, u32> {
        let data: Vec<Iv> = (0..100)
            .step_by(20)
            .map(|x| Iv {
                start: x,
                end: x + 10,
                val: 0,
            })
            .collect();
        Lapper::new(data)
    }

    fn setup_overlapping() -> Lapper<usize, u32> {
        let data: Vec<Iv> = (0..100)
            .step_by(10)
            .map(|x| Iv {
                start: x,
                end: x + 15,
                val: 0,
            })
            .collect();
        Lapper::new(data)
    }

    fn setup_badlapper() -> Lapper<usize, u32> {
        let data: Vec<Iv> = vec![
            Iv{start: 70, end: 120, val: 0}, // max_len = 50
            Iv{start: 10, end: 15, val: 0},
            Iv{start: 10, end: 15, val: 0}, // exact overlap
            Iv{start: 12, end: 15, val: 0}, // inner overlap
            Iv{start: 14, end: 16, val: 0}, // overlap end
            Iv{start: 40, end: 45, val: 0},
            Iv{start: 50, end: 55, val: 0},
            Iv{start: 60, end: 65, val: 0},
            Iv{start: 68, end: 71, val: 0}, // overlap start
            Iv{start: 70, end: 75, val: 0},
        ];
        Lapper::new(data)
    }

    fn setup_single() -> Lapper<usize, u32> {
        let data: Vec<Iv> = vec![Iv {
            start: 10,
            end: 35,
            val: 0,
        }];
        Lapper::new(data)
    }

    // Test that inserting data ends up with the same lapper (nonoverlapping)
    #[test]
    fn insert_equality_nonoverlapping() {
        let data: Vec<Iv> = (0..100)
            .step_by(20)
            .map(|x| Iv {
                start: x,
                end: x + 10,
                val: 0,
            })
            .collect();
        let new_lapper = Lapper::new(data.clone());
        let mut insert_lapper = Lapper::new(vec![]);
        for elem in data {
            insert_lapper.insert(elem);
        }
        assert_eq!(new_lapper.starts, insert_lapper.starts);
        assert_eq!(new_lapper.ends, insert_lapper.ends);
        assert_eq!(new_lapper.intervals, insert_lapper.intervals);
        assert_eq!(new_lapper.max_len, insert_lapper.max_len);
    }

    // Test that inserting data ends up with the same lapper (overlapping)
    #[test]
    fn insert_equality_overlapping() {
        let data: Vec<Iv> = (0..100)
            .step_by(10)
            .map(|x| Iv {
                start: x,
                end: x + 15,
                val: 0,
            })
            .collect();
        let new_lapper = Lapper::new(data.clone());
        let mut insert_lapper = Lapper::new(vec![]);
        for elem in data {
            insert_lapper.insert(elem);
        }
        assert_eq!(new_lapper.starts, insert_lapper.starts);
        assert_eq!(new_lapper.ends, insert_lapper.ends);
        assert_eq!(new_lapper.intervals, insert_lapper.intervals);
        assert_eq!(new_lapper.max_len, insert_lapper.max_len);
    }

    // Test that inserting data half with new and half with insert
    // ends up with the same lapper
    #[test]
    fn insert_equality_half_and_half() {
        let data: Vec<Iv> = (0..100)
            .step_by(1)
            .map(|x| Iv {
                start: x,
                end: x + 15,
                val: 0,
            })
            .collect();
        let new_lapper = Lapper::new(data.clone());
        let (new_data, insert_data) = data.split_at(50);
        let mut insert_lapper = Lapper::new(new_data.to_vec());
        let mut insert_data = insert_data.to_vec();
        insert_data.reverse();
        for elem in insert_data {
            insert_lapper.insert(elem);
        }
        assert_eq!(new_lapper.starts, insert_lapper.starts);
        assert_eq!(new_lapper.ends, insert_lapper.ends);
        assert_eq!(new_lapper.intervals, insert_lapper.intervals);
        assert_eq!(new_lapper.max_len, insert_lapper.max_len);
    }

    // Test that inserting data ends up with the same lapper (badlapper)
    #[test]
    fn insert_equality_badlapper() {
        let data: Vec<Iv> = vec![
            Iv{start: 70, end: 120, val: 0}, // max_len = 50
            Iv{start: 10, end: 15, val: 0},
            Iv{start: 10, end: 15, val: 0}, // exact overlap
            Iv{start: 12, end: 15, val: 0}, // inner overlap
            Iv{start: 14, end: 16, val: 0}, // overlap end
            Iv{start: 40, end: 45, val: 0},
            Iv{start: 50, end: 55, val: 0},
            Iv{start: 60, end: 65, val: 0},
            Iv{start: 68, end: 71, val: 0}, // overlap start
            Iv{start: 70, end: 75, val: 0},
        ];
        let new_lapper = Lapper::new(data.clone());
        let mut insert_lapper = Lapper::new(vec![]);
        for elem in data {
            insert_lapper.insert(elem);
        }
        assert_eq!(new_lapper.starts, insert_lapper.starts);
        assert_eq!(new_lapper.ends, insert_lapper.ends);
        assert_eq!(new_lapper.intervals, insert_lapper.intervals);
        assert_eq!(new_lapper.max_len, insert_lapper.max_len);
    }

    // Test that inserting data ends up with the same lapper (single)
    #[test]
    fn insert_equality_single() {
        let data: Vec<Iv> = vec![Iv {
            start: 10,
            end: 35,
            val: 0,
        }];
        let new_lapper = Lapper::new(data.clone());
        let mut insert_lapper = Lapper::new(vec![]);
        for elem in data {
            insert_lapper.insert(elem);
        }
        assert_eq!(new_lapper.starts, insert_lapper.starts);
        assert_eq!(new_lapper.ends, insert_lapper.ends);
        assert_eq!(new_lapper.intervals, insert_lapper.intervals);
        assert_eq!(new_lapper.max_len, insert_lapper.max_len);
    }

    // Test that a query end that hits an interval start returns no interval
    #[test]
    fn test_query_end_interval_start() {
        let lapper = setup_nonoverlapping();
        let mut cursor = 0;
        assert_eq!(None, lapper.find(15, 20).next());
        assert_eq!(None, lapper.seek(15, 20, &mut cursor).next());
        assert_eq!(lapper.find(15, 20).count(), lapper.count(15, 20));
    }

    // Test that a query start that hits an interval end returns no interval
    #[test]
    fn test_query_start_interval_end() {
        let lapper = setup_nonoverlapping();
        let mut cursor = 0;
        assert_eq!(None, lapper.find(30, 35).next());
        assert_eq!(None, lapper.seek(30, 35, &mut cursor).next());
        assert_eq!(lapper.find(30, 35).count(), lapper.count(30, 35));
    }

    // Test that a query that overlaps the start of an interval returns that interval
    #[test]
    fn test_query_overlaps_interval_start() {
        let lapper = setup_nonoverlapping();
        let mut cursor = 0;
        let expected = Iv {
            start: 20,
            end: 30,
            val: 0,
        };
        assert_eq!(Some(&expected), lapper.find(15, 25).next());
        assert_eq!(Some(&expected), lapper.seek(15, 25, &mut cursor).next());
        assert_eq!(lapper.find(15, 25).count(), lapper.count(15, 25));
    }

    // Test that a query that overlaps the end of an interval returns that interval
    #[test]
    fn test_query_overlaps_interval_end() {
        let lapper = setup_nonoverlapping();
        let mut cursor = 0;
        let expected = Iv {
            start: 20,
            end: 30,
            val: 0,
        };
        assert_eq!(Some(&expected), lapper.find(25, 35).next());
        assert_eq!(Some(&expected), lapper.seek(25, 35, &mut cursor).next());
        assert_eq!(lapper.find(25, 35).count(), lapper.count(25, 35));
    }

    // Test that a query that is enveloped by interval returns interval
    #[test]
    fn test_interval_envelops_query() {
        let lapper = setup_nonoverlapping();
        let mut cursor = 0;
        let expected = Iv {
            start: 20,
            end: 30,
            val: 0,
        };
        assert_eq!(Some(&expected), lapper.find(22, 27).next());
        assert_eq!(Some(&expected), lapper.seek(22, 27, &mut cursor).next());
        assert_eq!(lapper.find(22, 27).count(), lapper.count(22, 27));
    }

    // Test that a query that envolops an interval returns that interval
    #[test]
    fn test_query_envolops_interval() {
        let lapper = setup_nonoverlapping();
        let mut cursor = 0;
        let expected = Iv {
            start: 20,
            end: 30,
            val: 0,
        };
        assert_eq!(Some(&expected), lapper.find(15, 35).next());
        assert_eq!(Some(&expected), lapper.seek(15, 35, &mut cursor).next());
        assert_eq!(lapper.find(15, 35).count(), lapper.count(15, 35));
    }

    #[test]
    fn test_overlapping_intervals() {
        let lapper = setup_overlapping();
        let mut cursor = 0;
        let e1 = Iv {
            start: 0,
            end: 15,
            val: 0,
        };
        let e2 = Iv {
            start: 10,
            end: 25,
            val: 0,
        };
        assert_eq!(vec![&e1, &e2], lapper.find(8, 20).collect::<Vec<&Iv>>());
        assert_eq!(
            vec![&e1, &e2],
            lapper.seek(8, 20, &mut cursor).collect::<Vec<&Iv>>()
        );
        assert_eq!(lapper.count(8, 20), 2);
    }

    #[test]
    fn test_merge_overlaps() {
        let mut lapper = setup_badlapper();
        let expected: Vec<&Iv> = vec![
            &Iv{start: 10, end: 16, val: 0},
            &Iv{start: 40, end: 45, val: 0},
            &Iv{start: 50, end: 55, val: 0},
            &Iv{start: 60, end: 65, val: 0},
            &Iv{start: 68, end: 120, val: 0}, // max_len = 50
        ];
        assert_eq!(lapper.intervals.len(), lapper.starts.len());
        lapper.merge_overlaps();
        assert_eq!(expected, lapper.iter().collect::<Vec<&Iv>>());
        assert_eq!(lapper.intervals.len(), lapper.starts.len())

    }

    // This test was added because this breakage was found in a library user's code, where after
    // calling merge_overlaps(), the find() call returned an empty iterator.
    #[test]
    fn test_merge_overlaps_find() {
        let data = vec![
                Iv{start: 2, end: 3, val: 0},
                Iv{start: 3, end: 4, val: 0},
                Iv{start: 4, end: 6, val: 0},
                Iv{start: 6, end: 7, val: 0},
                Iv{start: 7, end: 8, val: 0},
        ];
        let mut lapper = Lapper::new(data);

        let found = lapper.find(7, 9).collect::<Vec<&Interval<_,_>>>();
        assert_eq!(found, vec![
            &Iv{start:7, end: 8, val: 0},
        ]);

        // merge_overlaps should merge all intervals to one, which should be returned in the find call.
        lapper.merge_overlaps();

        let found = lapper.find(7, 9).collect::<Vec<&Interval<_,_>>>();
        assert_eq!(found, vec![
            &Iv{start:2, end: 8, val: 0},
        ]);
    }

    #[test]
    fn test_lapper_cov() {
        let mut lapper = setup_badlapper();
        let before = lapper.cov();
        lapper.merge_overlaps();
        let after = lapper.cov();
        assert_eq!(before, after);

        let mut lapper = setup_nonoverlapping();
        lapper.set_cov();
        assert_eq!(lapper.cov(), 50);
    }

    #[test]
    fn test_interval_intersects() {
        let i1 = Iv{start: 70, end: 120, val: 0}; // max_len = 50
        let i2 = Iv{start: 10, end: 15, val: 0};
        let i3 = Iv{start: 10, end: 15, val: 0}; // exact overlap
        let i4 = Iv{start: 12, end: 15, val: 0}; // inner overlap
        let i5 = Iv{start: 14, end: 16, val: 0}; // overlap end
        let i6 = Iv{start: 40, end: 50, val: 0};
        let i7 = Iv{start: 50, end: 55, val: 0};
        let i_8 = Iv{start: 60, end: 65, val: 0};
        let i9 = Iv{start: 68, end: 71, val: 0}; // overlap start
        let i10 = Iv{start: 70, end: 75, val: 0};

        assert_eq!(i2.intersect(&i3), 5); // exact match
        assert_eq!(i2.intersect(&i4), 3); // inner intersect
        assert_eq!(i2.intersect(&i5), 1); // end intersect
        assert_eq!(i9.intersect(&i10), 1); // start intersect
        assert_eq!(i7.intersect(&i_8), 0); // no intersect
        assert_eq!(i6.intersect(&i7), 0); // no intersect end = start
        assert_eq!(i1.intersect(&i10), 5); // inner intersect at start
    }

    #[test]
    fn test_union_and_intersect() {
        let data1: Vec<Iv> = vec![
            Iv{start: 70, end: 120, val: 0}, // max_len = 50
            Iv{start: 10, end: 15, val: 0}, // exact overlap
            Iv{start: 12, end: 15, val: 0}, // inner overlap
            Iv{start: 14, end: 16, val: 0}, // overlap end
            Iv{start: 68, end: 71, val: 0}, // overlap start
        ];
        let data2: Vec<Iv> = vec![

            Iv{start: 10, end: 15, val: 0},
            Iv{start: 40, end: 45, val: 0},
            Iv{start: 50, end: 55, val: 0},
            Iv{start: 60, end: 65, val: 0},
            Iv{start: 70, end: 75, val: 0},
        ];

        let (mut lapper1, mut lapper2) = (Lapper::new(data1), Lapper::new(data2)) ;
        // Should be the same either way it's calculated
        let (union, intersect) = lapper1.union_and_intersect(&lapper2);
        assert_eq!(intersect, 10);
        assert_eq!(union, 73);
        let (union, intersect) = lapper2.union_and_intersect(&lapper1);
        assert_eq!(intersect, 10);
        assert_eq!(union, 73);
        lapper1.merge_overlaps();
        lapper1.set_cov();
        lapper2.merge_overlaps();
        lapper2.set_cov();

        // Should be the same either way it's calculated
        let (union, intersect) = lapper1.union_and_intersect(&lapper2);
        assert_eq!(intersect, 10);
        assert_eq!(union, 73);
        let (union, intersect) = lapper2.union_and_intersect(&lapper1);
        assert_eq!(intersect, 10);
        assert_eq!(union, 73);
    }

    #[test]
    fn test_find_overlaps_in_large_intervals() {
        let data1: Vec<Iv> = vec![
            Iv{start: 0, end: 8, val: 0},
            Iv{start: 1, end: 10, val: 0},
            Iv{start: 2, end: 5, val: 0},
            Iv{start: 3, end: 8, val: 0},
            Iv{start: 4, end: 7, val: 0},
            Iv{start: 5, end: 8, val: 0},
            Iv{start: 8, end: 8, val: 0},
            Iv{start: 9, end: 11, val: 0},
            Iv{start: 10, end: 13, val: 0},
            Iv{start: 100, end: 200, val: 0},
            Iv{start: 110, end: 120, val: 0},
            Iv{start: 110, end: 124, val: 0},
            Iv{start: 111, end: 160, val: 0},
            Iv{start: 150, end: 200, val: 0},
        ];
        let lapper = Lapper::new(data1);
        let found = lapper.find(8, 11).collect::<Vec<&Iv>>();
        assert_eq!(found, vec![
            &Iv{start: 1, end: 10, val: 0},
            &Iv{start: 9, end: 11, val: 0},
            &Iv{start: 10, end: 13, val: 0},
        ]);
        assert_eq!(lapper.count(8, 11), 3);
        let found = lapper.find(145, 151).collect::<Vec<&Iv>>();
        assert_eq!(found, vec![
            &Iv{start: 100, end: 200, val: 0},
            &Iv{start: 111, end: 160, val: 0},
            &Iv{start: 150, end: 200, val: 0},
        ]);

        assert_eq!(lapper.count(145, 151), 3);
    }

    #[test]
    fn test_depth_sanity() {
        let data1: Vec<Iv> = vec![
            Iv{start: 0, end: 10, val: 0},
            Iv{start: 5, end: 10, val: 0}
        ];
        let lapper = Lapper::new(data1);
        let found = lapper.depth().collect::<Vec<Interval<usize, usize>>>();
        assert_eq!(found, vec![
                   Interval{start: 0, end: 5, val: 1},
                   Interval{start: 5, end: 10, val: 2}
        ]);
    }

    #[test]
    fn test_depth_hard() {
        let data1: Vec<Iv> = vec![
            Iv{start: 1, end: 10, val: 0},
            Iv{start: 2, end: 5, val: 0},
            Iv{start: 3, end: 8, val: 0},
            Iv{start: 3, end: 8, val: 0},
            Iv{start: 3, end: 8, val: 0},
            Iv{start: 5, end: 8, val: 0},
            Iv{start: 9, end: 11, val: 0},
        ];
        let lapper = Lapper::new(data1);
        let found = lapper.depth().collect::<Vec<Interval<usize, usize>>>();
        assert_eq!(found, vec![
                   Interval{start: 1, end: 2, val: 1},
                   Interval{start: 2, end: 3, val: 2},
                   Interval{start: 3, end: 8, val: 5},
                   Interval{start: 8, end: 9, val: 1},
                   Interval{start: 9, end: 10, val: 2},
                   Interval{start: 10, end: 11, val: 1},
        ]);
    }
    #[test]
    fn test_depth_harder() {
        let data1: Vec<Iv> = vec![
            Iv{start: 1, end: 10, val: 0},
            Iv{start: 2, end: 5, val: 0},
            Iv{start: 3, end: 8, val: 0},
            Iv{start: 3, end: 8, val: 0},
            Iv{start: 3, end: 8, val: 0},
            Iv{start: 5, end: 8, val: 0},
            Iv{start: 9, end: 11, val: 0},
            Iv{start: 15, end: 20, val: 0},
        ];
        let lapper = Lapper::new(data1);
        let found = lapper.depth().collect::<Vec<Interval<usize, usize>>>();
        assert_eq!(found, vec![
                   Interval{start: 1, end: 2, val: 1},
                   Interval{start: 2, end: 3, val: 2},
                   Interval{start: 3, end: 8, val: 5},
                   Interval{start: 8, end: 9, val: 1},
                   Interval{start: 9, end: 10, val: 2},
                   Interval{start: 10, end: 11, val: 1},
                   Interval{start: 15, end: 20, val: 1},
        ]);
    }
    // BUG TESTS - these are tests that came from real life

    // Test that it's not possible to induce index out of bounds by pushing the cursor past the end
    // of the lapper.
    #[test]
    fn test_seek_over_len() {
        let lapper = setup_nonoverlapping();
        let single = setup_single();
        let mut cursor: usize = 0;

        for interval in lapper.iter() {
            for o_interval in single.seek(interval.start, interval.end, &mut cursor) {
                println!("{:#?}", o_interval);
            }
        }
    }

    // Test that if lower_bound puts us before the first match, we still return a match
    #[test]
    fn test_find_over_behind_first_match() {
        let lapper = setup_badlapper();
        let e1 = Iv {start: 50, end: 55, val: 0};
        let found = lapper.find(50, 55).next();
        assert_eq!(found, Some(&e1));
        assert_eq!(lapper.find(50, 55).count(), lapper.count(50,55));
    }

    // When there is a very long interval that spans many little intervals, test that the little
    // intevals still get returne properly
    #[test]
    fn test_bad_skips() {
        let data = vec![
            Iv{start:25264912, end: 25264986, val: 0},
            Iv{start:27273024, end: 27273065	, val: 0},
            Iv{start:27440273, end: 27440318	, val: 0},
            Iv{start:27488033, end: 27488125	, val: 0},
            Iv{start:27938410, end: 27938470	, val: 0},
            Iv{start:27959118, end: 27959171	, val: 0},
            Iv{start:28866309, end: 33141404	, val: 0},
        ];
        let lapper = Lapper::new(data);

        let found = lapper.find(28974798, 33141355).collect::<Vec<&Iv>>();
        assert_eq!(found, vec![
            &Iv{start:28866309, end: 33141404	, val: 0},
        ]);
        assert_eq!(lapper.count(28974798, 33141355), 1);
    }

}

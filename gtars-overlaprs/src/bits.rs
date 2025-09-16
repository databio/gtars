use num_traits::{
    PrimInt, Unsigned,
    identities::{one, zero},
};

use super::Overlapper;
use gtars_core::models::Interval;

#[derive(Debug, Clone)]
pub struct Bits<I, T>
where
    I: PrimInt + Unsigned + Send + Sync,
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

impl<I, T> Overlapper<I, T> for Bits<I, T>
where
    I: PrimInt + Unsigned + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    /// Create a new instance of Bits by passing in a vector of Intervals. This vector will
    /// immediately be sorted by start order.
    /// ```
    /// use gtars_overlaprs::{Bits, Overlapper};
    /// use gtars_core::models::Interval;
    ///
    /// let data = (0..20).step_by(5)
    ///                   .map(|x| Interval{start: x, end: x + 10, val: true})
    ///                   .collect::<Vec<Interval<usize, bool>>>();
    /// let bits = Bits::build(data);
    /// ```
    fn build(mut intervals: Vec<Interval<I, T>>) -> Self
    where
        Self: Sized,
    {
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
        Bits {
            intervals,
            starts,
            ends,
            max_len,
            cov: None,
            overlaps_merged: false,
        }
    }

    /// Find all intervals that overlap start .. stop
    /// ```
    /// use gtars_overlaprs::{Bits, Overlapper};
    /// use gtars_core::models::Interval;
    ///
    /// let bits = Bits::build((0..100).step_by(5)
    ///                                 .map(|x| Interval{start: x, end: x+2 , val: true})
    ///                                 .collect::<Vec<Interval<usize, bool>>>());
    /// assert_eq!(bits.find_iter(5, 11).count(), 2);
    /// ```
    #[inline]
    fn find(&self, start: I, stop: I) -> Vec<Interval<I, T>> {
        let finder = IterFind {
            inner: self,
            off: Self::lower_bound(
                start.checked_sub(&self.max_len).unwrap_or_else(zero::<I>),
                &self.intervals,
            ),
            start,
            stop,
        };

        // TODO: we are literally just collection the iterator, which feels like smell
        // but how do we write the trait in such a way that we return an actual
        // iterator instead of a vector of intervals?
        finder.into_iter().cloned().collect()
    }

    fn find_iter<'a>(
        &'a self,
        start: I,
        stop: I,
    ) -> Box<dyn Iterator<Item = &'a Interval<I, T>> + 'a> {
        let finder = IterFind {
            inner: self,
            off: Self::lower_bound(
                start.checked_sub(&self.max_len).unwrap_or_else(zero::<I>),
                &self.intervals,
            ),
            start,
            stop,
        };
        Box::new(finder)
    }
}

impl<I, T> Bits<I, T>
where
    I: PrimInt + Unsigned + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    /// Insert a new interval after the BITS has been created. This is very
    /// inefficient and should be avoided if possible.
    ///
    /// SIDE EFFECTS: This clears cov() and overlaps_merged
    /// meaning that those will have to be recomputed after a insert
    /// ```
    /// use gtars_overlaprs::{Bits, Overlapper};
    /// use gtars_core::models::Interval;
    ///
    /// let data : Vec<Interval<usize, usize>>= vec!{
    ///     Interval{start:0,  end:5,  val:1},
    ///     Interval{start:6,  end:10, val:2},
    /// };
    /// let mut bits = Bits::build(data);
    /// bits.insert(Interval{start:0, end:20, val:5});
    /// assert_eq!(bits.len(), 3);
    /// assert_eq!(bits.find_iter(1,3).collect::<Vec<&Interval<usize,usize>>>(),
    ///     vec![
    ///         &Interval{start:0, end:5, val:1},
    ///         &Interval{start:0, end:20, val:5},
    ///     ]
    /// );
    ///
    /// ```
    pub fn insert(&mut self, elem: Interval<I, T>) {
        let starts_insert_index = Self::bsearch_seq(elem.start, &self.starts);
        let stops_insert_index = Self::bsearch_seq(elem.end, &self.ends);
        let intervals_insert_index = Self::bsearch_seq_ref(&elem, &self.intervals);
        let i_len = elem.end.checked_sub(&elem.start).unwrap_or_else(zero::<I>);
        if i_len > self.max_len {
            self.max_len = i_len;
        }
        self.starts.insert(starts_insert_index, elem.start);
        self.ends.insert(stops_insert_index, elem.end);
        self.intervals.insert(intervals_insert_index, elem);
        self.cov = None;
        self.overlaps_merged = false;
    }

    /// Get the number over intervals in Bits
    /// ```
    /// use gtars_overlaprs::{Bits, Overlapper};
    /// use gtars_core::models::Interval;
    ///
    /// let data = (0..20).step_by(5)
    ///                   .map(|x| Interval{start: x, end: x + 10, val: true})
    ///                   .collect::<Vec<Interval<usize, bool>>>();
    ///
    /// let bits = Bits::build(data);
    /// assert_eq!(bits.len(), 4);
    /// ```
    #[inline]
    pub fn len(&self) -> usize {
        self.intervals.len()
    }

    /// Check if BITS is empty
    /// ```
    /// use gtars_overlaprs::{Bits, Overlapper};
    /// use gtars_core::models::Interval;
    ///
    /// let data: Vec<Interval<usize, bool>> = vec![];
    /// let bits = Bits::build(data);
    /// assert_eq!(bits.is_empty(), true);
    /// ```
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.intervals.is_empty()
    }

    /// Return an iterator over the intervals in Bits
    #[inline]
    pub fn iter(&'_ self) -> IterBits<'_, I, T> {
        IterBits {
            inner: self,
            pos: 0,
        }
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
        if elems.is_empty() || elems[0] >= *key {
            return 0;
        } else if elems[elems.len() - 1] < *key {
            return elems.len();
        }

        let mut cursor = 0;
        let mut length = elems.len();
        while length > 1 {
            let half = length >> 1;
            length -= half;
            cursor += (usize::from(elems[cursor + half - 1] < *key)) * half;
        }
        cursor
    }

    /// Count all intervals that overlap start .. stop. This performs two binary search in order to
    /// find all the excluded elements, and then deduces the intersection from there. See
    /// [BITS](https://arxiv.org/pdf/1208.3407.pdf) for more details.
    /// ```
    /// use gtars_overlaprs::{Bits, Overlapper};
    /// use gtars_core::models::Interval;
    ///
    /// let bits = Bits::build((0..100).step_by(5)
    ///                                 .map(|x| Interval{start: x, end: x+2 , val: true})
    ///                                 .collect::<Vec<Interval<usize, bool>>>());
    /// assert_eq!(bits.count(5, 11), 2);
    /// ```
    #[inline]
    pub fn count(&self, start: I, stop: I) -> usize {
        let len = self.intervals.len();
        // Plus one to account for half-openness of bits intervals compared to BITS paper
        let first = Self::bsearch_seq(start + one::<I>(), &self.ends);
        let last = Self::bsearch_seq(stop, &self.starts);
        let num_cant_after = len - last;
        len - first - num_cant_after
    }

    /// Find all intevals that overlap start .. stop. This method will work when queries
    /// to this Bits are in sorted (start) order. It uses a linear search from the last query
    /// instead of a binary search. A reference to a cursor must be passed in. This reference will
    /// be modified and should be reused in the next query. This allows seek to not need to make
    /// the Bits object mutable, and thus use the same Bits accross threads.
    /// ```
    /// use gtars_overlaprs::{Bits, Overlapper};
    /// use gtars_core::models::Interval;
    ///
    /// let bits = Bits::build((0..100).step_by(5)
    ///                                 .map(|x| Interval{start: x, end: x+2 , val: true})
    ///                                 .collect::<Vec<Interval<usize, bool>>>());
    /// let mut cursor = 0;
    /// for i in bits.iter() {
    ///    assert_eq!(bits.seek(i.start, i.end, &mut cursor).count(), 1);
    /// }
    /// ```
    #[inline]
    pub fn seek<'a>(&'a self, start: I, stop: I, cursor: &mut usize) -> IterFind<'a, I, T> {
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
            stop,
        }
    }
}

/// Find Iterator
#[derive(Debug)]
pub struct IterFind<'a, I, T>
where
    T: Eq + Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Send + Sync,
{
    inner: &'a Bits<I, T>,
    off: usize,
    start: I,
    stop: I,
}

impl<'a, I, T> Iterator for IterFind<'a, I, T>
where
    T: Eq + Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Send + Sync,
{
    type Item = &'a Interval<I, T>;

    #[inline]
    // interval.start < stop && interval.end > start
    fn next(&mut self) -> Option<Self::Item> {
        while self.off < self.inner.intervals.len() {
            //let mut generator = self.inner.intervals[self.off..].iter();
            //while let Some(interval) = generator.next() {
            let interval = &self.inner.intervals[self.off];
            self.off += 1;
            if interval.overlap(self.start, self.stop) {
                return Some(interval);
            } else if interval.start >= self.stop {
                break;
            }
        }
        None
    }
}

/// Bits Iterator
pub struct IterBits<'a, I, T>
where
    T: Eq + Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Send + Sync,
{
    inner: &'a Bits<I, T>,
    pos: usize,
}

impl<'a, I, T> Iterator for IterBits<'a, I, T>
where
    T: Eq + Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Send + Sync,
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

impl<I, T> IntoIterator for Bits<I, T>
where
    T: Eq + Clone + Send + Sync,
    I: PrimInt + Unsigned + Send + Sync,
{
    type Item = Interval<I, T>;
    type IntoIter = ::std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.intervals.into_iter()
    }
}

impl<'a, I, T> IntoIterator for &'a Bits<I, T>
where
    T: Eq + Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Send + Sync,
{
    type Item = &'a Interval<I, T>;
    type IntoIter = std::slice::Iter<'a, Interval<I, T>>;

    fn into_iter(self) -> std::slice::Iter<'a, Interval<I, T>> {
        self.intervals.iter()
    }
}

impl<'a, I, T> IntoIterator for &'a mut Bits<I, T>
where
    T: Eq + Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Send + Sync,
{
    type Item = &'a mut Interval<I, T>;
    type IntoIter = std::slice::IterMut<'a, Interval<I, T>>;

    fn into_iter(self) -> std::slice::IterMut<'a, Interval<I, T>> {
        self.intervals.iter_mut()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use pretty_assertions::{assert_eq, assert_ne};
    use rstest::{fixture, rstest};

    #[fixture]
    fn intervals() -> Vec<Interval<u32, &'static str>> {
        vec![
            Interval {
                start: 1,
                end: 5,
                val: "a",
            },
            Interval {
                start: 3,
                end: 7,
                val: "b",
            },
            Interval {
                start: 6,
                end: 10,
                val: "c",
            },
            Interval {
                start: 8,
                end: 12,
                val: "d",
            },
        ]
    }

    #[rstest]
    fn test_build_and_len(intervals: Vec<Interval<u32, &'static str>>) {
        let ailist = Bits::build(intervals.clone());
        assert_eq!(ailist.len(), intervals.len());
        assert_ne!(ailist.is_empty(), true);
    }

    #[rstest]
    fn test_find_overlapping_intervals(intervals: Vec<Interval<u32, &'static str>>) {
        let ailist = Bits::build(intervals);

        // Query that overlaps with "a" and "b"
        let results = ailist.find(2, 4);
        let vals: Vec<&str> = results.iter().map(|i| i.val).collect();
        assert_eq!(vals.contains(&"a"), true);
        assert_eq!(vals.contains(&"b"), true);
        assert_eq!(vals.contains(&"c"), false);

        // Query that overlaps with "c" and "d"
        let results = ailist.find(9, 11);
        let vals: Vec<&str> = results.iter().map(|i| i.val).collect();
        assert_eq!(vals.contains(&"c"), true);
        assert_eq!(vals.contains(&"d"), true);
        assert_eq!(vals.contains(&"a"), false);
    }

    #[rstest]
    fn test_find_no_overlap(intervals: Vec<Interval<u32, &'static str>>) {
        let ailist = Bits::build(intervals);

        // Query outside all intervals
        let results = ailist.find(13, 15);
        assert_eq!(results.is_empty(), true);
    }

    #[rstest]
    fn test_empty_ailist() {
        let ailist: Bits<u32, &str> = Bits::build(vec![]);

        assert_eq!(ailist.len(), 0);
        assert_eq!(ailist.is_empty(), true);

        let results = ailist.find(1, 2);

        assert_eq!(results.is_empty(), true);
    }
}

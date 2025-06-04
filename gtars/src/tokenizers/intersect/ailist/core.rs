use num_traits::{PrimInt, Unsigned};

use crate::tokenizers::intersect::{Intersect, Interval, IntervalCount};

///
/// The Augmented Interval List (AIList), enumerates intersections between a query interval q and an interval set R.
///
pub struct AIList<I, T>
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    pub intervals: Vec<Interval<I, T>>,
    starts: Vec<I>,
    ends: Vec<I>,
    max_ends: Vec<I>,
    header_list: Vec<usize>,
}

impl<I, T> AIList<I, T>
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    ///
    /// Create a new AIList struct
    ///
    /// # Arguments
    /// - intervals: list of intervals to create from
    ///
    /// # Returns
    /// - AIList struct
    pub fn new(intervals: Vec<Interval<I, T>>, minimum_coverage_length: usize) -> AIList<I, T> {
        // in the future, clone and sort...
        let mut intervals = intervals;
        intervals.sort_by_key(|key| key.start);

        let mut starts = Vec::new();
        let mut ends = Vec::new();
        let mut max_ends = Vec::new();
        let mut header_list = vec![0];

        loop {
            let mut results = Self::decompose(&intervals, minimum_coverage_length);

            starts.append(&mut results.0);
            ends.append(&mut results.1);
            max_ends.append(&mut results.2);
            intervals = results.3;

            if intervals.is_empty() {
                break;
            } else {
                header_list.push(starts.len());
            }
        }

        AIList {
            starts,
            ends,
            max_ends,
            header_list,
            intervals,
        }
    }

    #[allow(clippy::type_complexity)]
    fn decompose(
        intervals: &[Interval<I, T>],
        minimum_coverage_length: usize,
    ) -> (Vec<I>, Vec<I>, Vec<I>, Vec<Interval<I, T>>) {
        // look at the next minL*2 intervals
        let mut starts = Vec::new();
        let mut ends = Vec::new();
        let mut max_ends = Vec::new();
        let mut l2 = Vec::new();

        for (index, interval) in intervals.iter().enumerate() {
            let mut count = 0;
            for i in 1..(minimum_coverage_length * 2) {
                match intervals.get(index + i) {
                    Some(interval2) => {
                        if interval.end > interval2.end {
                            count += 1;
                        }
                    }
                    None => break,
                }
            }
            if count >= minimum_coverage_length {
                l2.push(Interval {
                    start: interval.start,
                    end: interval.end,
                    val: interval.val.clone(),
                });
            } else {
                starts.push(interval.start);
                ends.push(interval.end);
            }
        }

        let mut max: I = I::zero();

        for end in ends.iter() {
            max = if max > *end { max } else { *end };
            max_ends.push(max);
        }

        (starts, ends, max_ends, l2)
    }

    pub fn find_overlaps(&self, start: I, end: I) -> IterFind<I, T> {
        IterFind {
            inner: self,
            start,
            end,
            idx:  0,
            slice_pos: 0
        }
    }

    #[inline]
    fn slice_bounds(&self, slice: usize) -> (usize, usize) {
        let lo = self.header_list[slice];
        let hi = *self.header_list.get(slice + 1).unwrap_or(&self.starts.len());
        (lo, hi)
    }
}

pub struct IterFind<'a, I, T>
where
    T: Eq + Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
{
    inner: &'a AIList<I, T>,
    start: I,
    end: I,
    idx: usize,
    slice_pos: usize
}

impl<'a, I, T> Iterator for IterFind<'a, I, T>
where
    T: Eq + Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
{
    type Item = &'a Interval<I, T>;

    fn next(&mut self) -> Option<Self::Item> {
        // outer loop
        while self.idx < self.inner.header_list.len() - 1 {
            // inner loop
            self.slice_pos = self.inner.starts.partition_point(|&x| x < self.end);

            // get slice starts and ends
            let starts = &self.inner.starts[self.inner.header_list[self.idx]..self.inner.header_list[self.idx + 1]];
            let ends = &self.inner.ends[self.inner.header_list[self.idx]..self.inner.header_list[self.idx + 1]];
            let max_ends = &self.inner.max_ends[self.inner.header_list[self.idx]..self.inner.header_list[self.idx + 1]];

            while self.slice_pos > 0 {
                self.slice_pos -= 1;
                if starts[self.slice_pos] > ends[self.slice_pos] {
                    if self.start > max_ends[self.slice_pos] {
                        break
                    }
                } else {
                    // TODO: this needs to be addressed, this is not right...
                    // we need to return a reference to an interval.
                    return Some(&self.inner.intervals[self.slice_pos])
                }
            }
            self.idx += 1;
        }

        // final pass on decomposed list
        if self.idx == self.inner.header_list.len() - 1 {
            self.slice_pos = self.idx;

            // get slice starts and ends
            let starts = &self.inner.starts[self.inner.header_list[self.idx]..];
            let ends = &self.inner.ends[self.inner.header_list[self.idx]..];
            let max_ends = &self.inner.max_ends[self.inner.header_list[self.idx]..];
            
            while self.slice_pos > 0 {
                self.slice_pos -= 1;
                if starts[self.slice_pos] > ends[self.slice_pos] {
                    if self.start > max_ends[self.slice_pos] {
                        break
                    }
                } else {
                    // TODO: this needs to be addressed, this is not right...
                    // we need to return a reference to an interval.
                    return Some(&self.inner.intervals[self.slice_pos])
                }
            }
        }

        None
    }
}

impl<I, T> IntervalCount for AIList<I, T>
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    fn len(&self) -> usize {
        self.starts.len()
    }

    fn is_empty(&self) -> bool {
        self.starts.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;

    #[fixture]
    fn consensus_set() -> &'static str {
        "tests/data/consensus/consensus1.bed"
    }

    #[rstest]
    fn test_create_ailist() {
        let intervals = vec![
            Interval {
                start: 0_u32,
                end: 5_u32,
                val: 0_u32,
            },
            Interval {
                start: 5_u32,
                end: 20_u32,
                val: 0_u32,
            },
            Interval {
                start: 20_u32,
                end: 25_u32,
                val: 0_u32,
            },
        ];

        let _ailist = AIList::new(intervals, 10);
    }

    #[rstest]
    fn test_query_ailist() {
        let universe_intervals = vec![
            Interval {
                start: 0_u32,
                end: 5_u32,
                val: 0_u32,
            },
            Interval {
                start: 5_u32,
                end: 20_u32,
                val: 0_u32,
            },
            Interval {
                start: 20_u32,
                end: 25_u32,
                val: 0_u32,
            },
        ];
        let ailist = AIList::new(universe_intervals, 10);

        let query_interval = (6, 11);

        let res = ailist.find_overlaps(query_interval.0, query_interval.1);
        let res = res.collect::<Vec<&Interval<u32, u32>>>();

        assert_eq!(
            *res[0],
            Interval {
                start: 5,
                end: 20,
                val: 0
            }
        );
    }
}

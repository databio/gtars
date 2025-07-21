use num_traits::{PrimInt, Unsigned};

use crate::common::models::Interval;
use super::Overlapper;

///
/// The Augmented Interval List (AiList), enumerates intersections between a query interval q and an interval set R.
///
pub struct AiList<I, T>
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    starts: Vec<I>,
    ends: Vec<I>,
    max_ends: Vec<I>,
    header_list: Vec<usize>,
    data_list: Vec<T>,
}

type DecomposeResult<I, T> = (Vec<I>, Vec<I>, Vec<I>, Vec<T>, Vec<Interval<I, T>>);

impl<I, T> Overlapper<I, T> for AiList<I, T>
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
    fn build(intervals: Vec<Interval<I, T>>) -> Self  where Self: Sized {
        // in the future, clone and sort...
        let mut intervals = intervals;
        intervals.sort_by_key(|key| key.start);

        let mut starts = Vec::new();
        let mut ends = Vec::new();
        let mut max_ends = Vec::new();
        let mut data_list = Vec::new();
        let mut header_list = vec![0];

        loop {
            let mut results = Self::decompose(&intervals, 10);

            starts.append(&mut results.0);
            ends.append(&mut results.1);
            max_ends.append(&mut results.2);
            data_list.append(&mut results.3);
            intervals = results.4;

            if intervals.is_empty() {
                break;
            } else {
                header_list.push(starts.len());
            }
        }

        AiList {
            starts,
            ends,
            max_ends,
            header_list,
            data_list,
        }
    }

    fn find(&self, start: I, end: I) -> Vec<Interval<I,T>> {
        let mut results_list = Vec::new();

        for i in 0..(self.header_list.len() - 1) {
            results_list.append(&mut Self::query_slice(
                start,
                end,
                &self.starts[self.header_list[i]..self.header_list[i + 1]],
                &self.ends[self.header_list[i]..self.header_list[i + 1]],
                &self.max_ends[self.header_list[i]..self.header_list[i + 1]],
                &self.data_list[self.header_list[i]..self.header_list[i + 1]],
            ));
        }
        // now do the last decomposed AiList
        let i = self.header_list.len() - 1;
        results_list.extend(Self::query_slice(
            start,
            end,
            &self.starts[self.header_list[i]..],
            &self.ends[self.header_list[i]..],
            &self.max_ends[self.header_list[i]..],
            &self.data_list[self.header_list[i]..],
        ));

        results_list
    }
}

impl<I, T> AiList<I, T>
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    T: Eq + Clone + Send + Sync,
{

    fn decompose(
        intervals: &[Interval<I, T>],
        minimum_coverage_length: usize,
    ) -> DecomposeResult<I, T> {
        // look at the next minL*2 intervals
        let mut starts = Vec::new();
        let mut ends = Vec::new();
        let mut max_ends = Vec::new();
        let mut data_list = Vec::new();
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
                data_list.push(interval.val.clone());
            }
        }

        let mut max: I = I::zero();

        for end in ends.iter() {
            max = if max > *end { max } else { *end };
            max_ends.push(max);
        }

        (starts, ends, max_ends, data_list, l2)
    }

    fn query_slice(
        start: I,
        end: I,
        starts: &[I],
        ends: &[I],
        max_ends: &[I],
        data_list: &[T],
    ) -> Vec<Interval<I, T>> {
        let mut results_list = Vec::new();
        let mut i = starts.partition_point(|&x| x < end);

        while i > 0 {
            i -= 1;
            if start > ends[i] {
                //this means that there is no intersection
                if start > max_ends[i] {
                    //there is no further intersection
                    return results_list;
                }
            } else {
                results_list.push(Interval {
                    start: starts[i],
                    end: ends[i],
                    val: data_list[i].clone(),
                })
            }
        }
        results_list
    }

    pub fn len(&self) -> usize {
        self.starts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.starts.is_empty()
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    use pretty_assertions::{assert_eq, assert_ne};
    use rstest::{rstest, fixture};

    #[fixture]
    fn intervals() -> Vec<Interval<u32, &'static str>> {
        vec![
            Interval { start: 1, end: 5, val: "a" },
            Interval { start: 3, end: 7, val: "b" },
            Interval { start: 6, end: 10, val: "c" },
            Interval { start: 8, end: 12, val: "d" },
        ]
    }

    #[rstest]
    fn test_build_and_len(intervals: Vec<Interval<u32, &'static str>>) {
        let ailist = AiList::build(intervals.clone());
        assert_eq!(ailist.len(), intervals.len());
        assert_ne!(ailist.is_empty(), true);
    }

    #[rstest]
    fn test_find_overlapping_intervals(intervals: Vec<Interval<u32, &'static str>>) {

        let ailist = AiList::build(intervals);

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
        
        let ailist = AiList::build(intervals);

        // Query outside all intervals
        let results = ailist.find(13, 15);
        assert_eq!(results.is_empty(), true);
    }

    #[rstest]
    fn test_empty_ailist() {

        let ailist: AiList<u32, &str> = AiList::build(vec![]);

        assert_eq!(ailist.len(), 0);
        assert_eq!(ailist.is_empty(), true);

        let results = ailist.find(1, 2);

        assert_eq!(results.is_empty(), true);
    }
}

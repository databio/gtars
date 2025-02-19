use num_traits::{PrimInt, Unsigned};
use std::fmt;

#[derive(Eq, Debug, Clone, PartialEq)]
pub struct Interval<I, T>
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    pub start: I,
    pub end: I,
    pub data: T,
}

///
/// The Augmented Interval List (AIList), enumerates intersections between a query interval q and an interval set R.
///
pub struct AIList<I, T>
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
        let mut data_list = Vec::new();
        let mut header_list = vec![0];

        loop {
            let mut results = Self::decompose(&intervals, minimum_coverage_length);

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

        AIList {
            starts,
            ends,
            max_ends,
            header_list,
            data_list,
        }
    }

    fn decompose(
        intervals: &[Interval<I, T>],
        minimum_coverage_length: usize,
    ) -> (Vec<I>, Vec<I>, Vec<I>, Vec<T>, Vec<Interval<I, T>>) {
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
                    data: interval.data.clone(),
                });
            } else {
                starts.push(interval.start);
                ends.push(interval.end);
                data_list.push(interval.data.clone());
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
                    data: data_list[i].clone(),
                })
            }
        }
        results_list
    }

    pub fn query(&self, start: I, end: I) -> Vec<Interval<I, T>> {
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
        // now do the last decomposed ailist
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
    use crate::common::utils::extract_regions_from_bed_file;
    use std::path::Path;
    // use pretty_assertions::assert_eq;
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
                data: 0_u32,
            },
            Interval {
                start: 5_u32,
                end: 20_u32,
                data: 0_u32,
            },
            Interval {
                start: 20_u32,
                end: 25_u32,
                data: 0_u32,
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
                data: 0_u32,
            },
            Interval {
                start: 5_u32,
                end: 20_u32,
                data: 0_u32,
            },
            Interval {
                start: 20_u32,
                end: 25_u32,
                data: 0_u32,
            },
        ];
        let ailist = AIList::new(universe_intervals, 10);

        let query_interval = (6, 11);

        let res = ailist.query(query_interval.0, query_interval.1);

        assert_eq!(
            res.first(),
            Some(&Interval {
                start: 5,
                end: 20,
                data: 0
            })
        );
    }

    #[rstest]
    fn test_ailist_tokenizer() {
        let bed_file_path = Path::new("tests/data/peaks.bed.gz");
    }
}

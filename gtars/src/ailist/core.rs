use std::fmt;

#[derive(Debug, PartialEq)]
pub struct Interval<T: Clone> {
    pub start: u32,
    pub end: u32,
    pub data: T
}

impl<T: Clone> fmt::Display for Interval<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {})", self.start, self.end)
    }
}

///
/// The Augmented Interval List (AIList), enumerates intersections between a query interval q and an interval set R.
///
pub struct AIList<T: Clone> {
    starts: Vec<u32>,
    ends: Vec<u32>,
    max_ends: Vec<u32>,
    header_list: Vec<usize>,
    data_list: Vec<T>
}

impl<T: Clone> AIList<T> {
    ///
    /// Create a new AIList struct
    ///
    /// # Arguments
    /// - intervals: list of intervals to create from
    ///
    /// # Returns
    /// - AIList struct
    pub fn new(intervals: Vec<Interval<T>>, minimum_coverage_length: usize) -> AIList<T> {
        // in the future, clone and sort...
        let mut intervals = intervals;
        intervals.sort_by_key(|key| key.start);

        let mut starts = Vec::new();
        let mut ends = Vec::new();
        let mut max_ends= Vec::new();
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
            data_list
        }
    }

    fn decompose(
        intervals: &[Interval<T>],
        minimum_coverage_length: usize,
    ) -> (Vec<u32>, Vec<u32>, Vec<u32>, Vec<T>, Vec<Interval<T>>) {
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
                    data: interval.data.clone()
                });
            } else {
                starts.push(interval.start);
                ends.push(interval.end);
                data_list.push(interval.data.clone());
            }
        }

        let mut max: u32 = 0;

        for end in ends.iter() {
            max = if max > *end { max } else { *end };
            max_ends.push(max);
        }

        (starts, ends, max_ends, data_list, l2)
    }

    
    fn query_slice(
        start: u32,
        end: u32,
        starts: &[u32],
        ends: &[u32],
        max_ends: &[u32],
        data_list: &[T]
    ) -> Vec<Interval<T>> {
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
                    data: data_list[i].clone()
                })
            }
        }
        results_list
    }

    pub fn query(&self, start: u32, end: u32) -> Vec<Interval<T>> {
        let mut results_list = Vec::new();

        for i in 0..(self.header_list.len() - 1) {
            results_list.append(&mut Self::query_slice(
                start, end,
                &self.starts[self.header_list[i]..self.header_list[i + 1]],
                &self.ends[self.header_list[i]..self.header_list[i + 1]],
                &self.max_ends[self.header_list[i]..self.header_list[i + 1]],
                &self.data_list[self.header_list[i]..self.header_list[i + 1]],
            ));
        }
        // now do the last decomposed ailist
        let i = self.header_list.len() - 1;
        results_list.extend(Self::query_slice(
            start, end,
            &self.starts[self.header_list[i]..],
            &self.ends[self.header_list[i]..],
            &self.max_ends[self.header_list[i]..],
            &self.data_list[self.header_list[i]..],
        ));

        results_list
    }
}

impl<T: Clone> fmt::Display for AIList<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut string = String::new();
        string.push('\n');
        for element in self.starts.iter() {
            string.push_str(format!("{element}").as_str());
        }
        string.push('\n');
        for element in self.ends.iter() {
            string.push_str(format!("{element}").as_str());
        }
        string.push('\n');
        for element in self.max_ends.iter() {
            string.push_str(format!("{element}").as_str());
        }
        string.push('\n');
        for element in self.header_list.iter() {
            string.push_str(format!("{element}").as_str());
        }
        string.push('\n');
        write!(f, "{string}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // use pretty_assertions::assert_eq;
    use rstest::*;

    #[fixture]
    fn path_to_fragment_files() -> &'static str {
        "tests/data/fragments/region_scoring/*.bed.gz"
    }

    #[fixture]
    fn consensus_set() -> &'static str {
        "tests/data/consensus/consensus1.bed"
    }

    #[fixture]
    fn output_file() -> &'static str {
        "tests/data/out/region_scoring_count.csv.gz"
    }

    #[rstest]
    fn test_create_ailist() {
        let intervals = vec![
            Interval { start: 0, end: 5, data: 0 },
            Interval { start: 5, end: 20, data: 0 },
            Interval { start: 20, end: 25, data: 0 }
        ];

        let _ailist = AIList::new(intervals, 10);
    }

    #[rstest]
    fn test_query_ailist() {
        let universe_intervals = vec![
            Interval { start: 0, end: 5, data: 0 },
            Interval { start: 5, end: 20, data: 0 },
            Interval { start: 20, end: 25, data: 0 }
        ];
        let ailist = AIList::new(universe_intervals, 10);

        let query_interval = (6, 11);

        let res = ailist.query(query_interval.0, query_interval.1);

        assert_eq!(res.first(), Some(&Interval { start: 5, end: 20, data: 0 }));
    }
}
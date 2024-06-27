use std::fmt;

pub struct Interval {
    pub start: u32,
    pub end: u32,
}

impl fmt::Display for Interval {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {})", self.start, self.end)
    }
}

///
/// The Augmented Interval List (AIList), enumerates intersections between a query interval q and an interval set R.
///
pub struct AIList {
    starts: Vec<u32>,
    ends: Vec<u32>,
    max_ends: Vec<u32>,
    header_list: Vec<usize>,
}

impl AIList {
    ///
    /// Create a new AIList struct
    ///
    /// # Arguments
    /// - intervals: list of intervals to create from
    ///
    /// # Returns
    /// - AIList struct
    pub fn new(intervals: &mut Vec<Interval>, minimum_coverage_length: usize) -> AIList {
        // in the future, clone and sort...
        intervals.sort_by_key(|key| key.start);

        let mut starts: Vec<u32> = Vec::new();
        let mut ends: Vec<u32> = Vec::new();
        let mut max_ends: Vec<u32> = Vec::new();
        let mut header_list: Vec<usize> = vec![0];

        loop {
            let mut results = Self::decompose(intervals, minimum_coverage_length);

            starts.append(&mut results.0);
            ends.append(&mut results.1);
            max_ends.append(&mut results.2);

            *intervals = results.3;

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
        }
    }

    fn decompose(
        intervals: &mut [Interval],
        minimum_coverage_length: usize,
    ) -> (Vec<u32>, Vec<u32>, Vec<u32>, Vec<Interval>) {
        // look at the next minL*2 intervals
        let mut starts: Vec<u32> = Vec::new();
        let mut ends: Vec<u32> = Vec::new();
        let mut max_ends: Vec<u32> = Vec::new();
        let mut l2: Vec<Interval> = Vec::new();

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
                });
            } else {
                starts.push(interval.start);
                ends.push(interval.end)
            }
        }

        let mut max: u32 = 0;

        for end in ends.iter() {
            max = if max > *end { max } else { *end };
            max_ends.push(max);
        }

        (starts, ends, max_ends, l2)
    }

    fn query_slice(
        interval: &Interval,
        starts: &[u32],
        ends: &[u32],
        max_ends: &[u32],
    ) -> Vec<Interval> {
        let mut results_list: Vec<Interval> = Vec::new();
        let mut i = starts.partition_point(|&x| x < interval.end);

        while i > 0 {
            i -= 1;
            if interval.start > ends[i] {
                //this means that there is no intersection
                if interval.start > max_ends[i] {
                    //there is no further intersection
                    return results_list;
                }
            } else {
                results_list.push(Interval {
                    start: starts[i],
                    end: ends[i],
                })
            }
        }
        results_list
    }

    pub fn query(&self, interval: &Interval) -> Vec<Interval> {
        let mut results_list: Vec<Interval> = Vec::new();

        for i in 0..(self.header_list.len() - 1) {
            results_list.append(&mut Self::query_slice(
                interval,
                &self.starts[self.header_list[i]..self.header_list[i + 1]],
                &self.ends[self.header_list[i]..self.header_list[i + 1]],
                &self.max_ends[self.header_list[i]..self.header_list[i + 1]],
            ));
        }
        // now do the last decomposed ailist
        let i = self.header_list.len() - 1;
        results_list.extend(Self::query_slice(
            interval,
            &self.starts[self.header_list[i]..],
            &self.ends[self.header_list[i]..],
            &self.max_ends[self.header_list[i]..],
        ));

        results_list
    }
}

impl fmt::Display for AIList {
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

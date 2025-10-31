use std::mem::swap;

use num_traits::{PrimInt, Unsigned};

use super::Overlapper;
use gtars_core::models::Interval;

/// An Augmented Interval List for efficient genomic interval overlap queries.
///
/// From the following article: <https://academic.oup.com/bioinformatics/article/35/23/4907/5509521>
///
/// The Augmented Interval List (AIList) is a data structure optimized for finding overlaps
/// between a query interval and a large collection of genomic intervals. It is particularly
/// efficient for datasets with high-coverage regions, which are common in genomic data such
/// as ChIP-seq peaks, gene annotations, or aligned reads.
///
/// # Examples
///
/// ```
/// use gtars_overlaprs::{AIList, Overlapper, Interval};
///
/// // Create intervals for genomic features
/// let genes = vec![
///     Interval { start: 1000u32, end: 2000, val: "GENE1" },
///     Interval { start: 1500, end: 2500, val: "GENE2" },
///     Interval { start: 5000, end: 6000, val: "GENE3" },
/// ];
///
/// let ailist = AIList::build(genes);
///
/// // Query for genes overlapping position 1800-2200
/// let overlaps = ailist.find(1800, 2200);
/// assert_eq!(overlaps.len(), 2); // GENE1 and GENE2
#[derive(Debug, Clone)]
pub struct AIList<I, T>
where
    I: PrimInt + Unsigned + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    starts: Vec<I>,
    ends: Vec<I>,
    max_ends: Vec<I>,
    header_list: Vec<usize>,
    stored_intervals: Vec<Interval<I, T>>,
}

/// Storage for the intermediate results from [`AIList::decompose`].
#[derive(Debug, Default)]
struct DecomposeResult<I, T>
where
    I: PrimInt + Unsigned + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    /// Start positions.
    starts: Vec<I>,
    /// End positions.
    ends: Vec<I>,
    /// The max end position seen up to this index.
    max_ends: Vec<I>,
    /// The associated Interval.
    stored_intervals: Vec<Interval<I, T>>,
    /// The remaining intervals to be decomposed.
    l2: Vec<Interval<I, T>>,
}

impl<I, T> DecomposeResult<I, T>
where
    I: PrimInt + Unsigned + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    /// Clear the contents of the [`DecomposeResult`], maintaining capacity.
    fn clear(&mut self) {
        self.starts.clear();
        self.ends.clear();
        self.max_ends.clear();
        self.stored_intervals.clear();
        self.l2.clear();
    }

    /// Create an empty [`DecomposeResult`] with the given `cap` capacity.
    fn with_capacity(cap: usize) -> Self {
        Self {
            starts: Vec::with_capacity(cap),
            ends: Vec::with_capacity(cap),
            max_ends: Vec::with_capacity(cap),
            stored_intervals: Vec::with_capacity(cap),
            l2: Vec::with_capacity(cap),
        }
    }
}

impl<I, T> Overlapper<I, T> for AIList<I, T>
where
    I: PrimInt + Unsigned + Send + Sync,
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
    fn build(intervals: Vec<Interval<I, T>>) -> Self
    where
        Self: Sized,
    {
        // in the future, clone and sort...
        let mut intervals = intervals;
        intervals.sort_by_key(|key| key.start);

        let mut starts = Vec::with_capacity(intervals.len());
        let mut ends = Vec::with_capacity(intervals.len());
        let mut max_ends = Vec::with_capacity(intervals.len());
        let mut stored_intervals = Vec::with_capacity(intervals.len());

        // Scratch space for construction
        // The scratch vecs will get drained into the above final vectors, but the capacity
        // in the scratch space will remain on each call, avoiding any additional allocations.
        //
        // Creating with cap of `intervals.len()` is overkill, but it means there will only every be the
        // on allocation and prevents any possible re-allocs if a given sub-list is larger than the previous.
        let mut results = DecomposeResult::with_capacity(intervals.len());

        let mut header_list = vec![0];

        loop {
            Self::decompose(&intervals, 10, &mut results);

            starts.append(&mut results.starts);
            ends.append(&mut results.ends);
            max_ends.append(&mut results.max_ends);
            stored_intervals.append(&mut results.stored_intervals);
            swap(&mut intervals, &mut results.l2);

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
            stored_intervals,
        }
    }

    fn find(&self, start: I, end: I) -> Vec<Interval<I, T>> {
        let mut results_list = Vec::new();

        for i in 0..(self.header_list.len() - 1) {
            results_list.append(&mut Self::query_slice(
                start,
                end,
                &self.starts[self.header_list[i]..self.header_list[i + 1]],
                &self.ends[self.header_list[i]..self.header_list[i + 1]],
                &self.max_ends[self.header_list[i]..self.header_list[i + 1]],
                &self.stored_intervals[self.header_list[i]..self.header_list[i + 1]],
            ));
        }
        // now do the last decomposed AIList
        let i = self.header_list.len() - 1;
        results_list.extend(Self::query_slice(
            start,
            end,
            &self.starts[self.header_list[i]..],
            &self.ends[self.header_list[i]..],
            &self.max_ends[self.header_list[i]..],
            &self.stored_intervals[self.header_list[i]..],
        ));

        results_list
    }

    fn find_iter<'a>(
        &'a self,
        start: I,
        stop: I,
    ) -> Box<dyn Iterator<Item = &'a Interval<I, T>> + 'a> {
        Box::new(IterFind::new(self, start, stop))
    }
}

impl<I, T> AIList<I, T>
where
    I: PrimInt + Unsigned + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    fn decompose(
        intervals: &[Interval<I, T>],
        minimum_coverage_length: usize,
        scratch: &mut DecomposeResult<I, T>,
    ) {
        scratch.clear();

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
                scratch.l2.push(Interval {
                    start: interval.start,
                    end: interval.end,
                    val: interval.val.clone(),
                });
            } else {
                scratch.starts.push(interval.start);
                scratch.ends.push(interval.end);
                scratch.stored_intervals.push(interval.clone());
            }
        }

        let mut max: I = I::zero();

        for end in scratch.ends.iter() {
            max = if max > *end { max } else { *end };
            scratch.max_ends.push(max);
        }
    }

    fn query_slice(
        start: I,
        end: I,
        starts: &[I],
        ends: &[I],
        max_ends: &[I],
        stored_intervals: &[Interval<I, T>],
    ) -> Vec<Interval<I, T>> {
        let mut results_list = Vec::new();
        let mut i = starts.partition_point(|&x| x < end);

        while i > 0 {
            i -= 1;
            // maintain start inclusive, end exclusive
            if start >= ends[i] {
                //this means that there is no intersection
                if start > max_ends[i] {
                    //there is no further intersection
                    return results_list;
                }
            } else {
                results_list.push(stored_intervals[i].clone())
            }
        }
        results_list
    }

    /// Returns the number of intervals in the AIList.
    pub fn len(&self) -> usize {
        self.starts.len()
    }

    /// Returns `true` if the AIList contains no intervals.
    pub fn is_empty(&self) -> bool {
        self.starts.is_empty()
    }
}

/// An iterator over intervals in an [`AIList`] that overlap with a query range.
///
/// This struct is created by the [`find_iter`](Overlapper::find_iter) method on [`AIList`].
/// It lazily yields references to intervals that overlap with the specified query range.
///
/// The iterator maintains state to traverse the decomposed sublists within the `AIList`
/// efficiently, yielding overlapping intervals one at a time without allocating a vector.
///
/// # Examples
///
/// ```
/// use gtars_overlaprs::{AIList, Overlapper, Interval};
///
/// let intervals = vec![
///     Interval { start: 10u32, end: 20, val: "a" },
///     Interval { start: 15, end: 25, val: "b" },
/// ];
///
/// let ailist = AIList::build(intervals);
///
/// // The iterator is created by find_iter
/// for interval in ailist.find_iter(12, 18) {
///     println!("Found: {}", interval.val);
/// }
/// ```
#[derive(Debug)]
pub struct IterFind<'a, I, T>
where
    T: Eq + Clone + Send + Sync + 'a,
    I: PrimInt + Unsigned + Send + Sync,
{
    inner: &'a AIList<I, T>,
    header_list_idx: usize,
    list_idx: Option<usize>,
    start: I,
    stop: I,
}

impl<'a, I, T> IterFind<'a, I, T>
where
    I: PrimInt + Unsigned + Send + Sync + 'a,
    T: Eq + Clone + Send + Sync,
{
    fn new(ailist: &'a AIList<I, T>, start: I, stop: I) -> Self {
        Self {
            inner: ailist,
            header_list_idx: 0,
            list_idx: None,
            start,
            stop,
        }
    }
}

impl<'a, I, T> Iterator for IterFind<'a, I, T>
where
    I: PrimInt + Unsigned + Send + Sync,
    T: Eq + Clone + Send + Sync + 'a,
{
    type Item = &'a Interval<I, T>;

    fn next(&mut self) -> Option<Self::Item> {
        while self.header_list_idx < self.inner.header_list.len() {
            let range = if self.header_list_idx == self.inner.header_list.len() - 1 {
                // is last
                self.inner.header_list[self.header_list_idx]..self.inner.starts.len()
            } else {
                self.inner.header_list[self.header_list_idx]
                    ..self.inner.header_list[self.header_list_idx + 1]
            };
            let starts = &self.inner.starts[range.clone()];
            let ends = &self.inner.ends[range.clone()];
            let max_ends = &self.inner.max_ends[range.clone()];
            let stored_intervalss = &self.inner.stored_intervals[range.clone()];

            let i = if let Some(list_idx) = self.list_idx.as_mut() {
                list_idx
            } else {
                self.list_idx = Some(starts.partition_point(|&x| x < self.stop));
                self.list_idx.as_mut().unwrap()
            };

            while *i > 0 {
                *i -= 1;
                // maintain start inclusive, end exclusive
                if self.start >= ends[*i] {
                    // there is no further intersection, try the next header_list_idx
                    if self.start > max_ends[*i] {
                        break;
                    }
                } else {
                    return Some(&stored_intervalss[*i]);
                }
            }
            self.list_idx = None;
            self.header_list_idx += 1;
        }
        None
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
        let ailist = AIList::build(intervals.clone());
        assert_eq!(ailist.len(), intervals.len());
        assert_ne!(ailist.is_empty(), true);
    }

    #[rstest]
    fn test_find_overlapping_intervals(intervals: Vec<Interval<u32, &'static str>>) {
        let ailist = AIList::build(intervals);

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
        let ailist = AIList::build(intervals);

        // Query outside all intervals
        let results = ailist.find(13, 15);
        assert_eq!(results.is_empty(), true);
    }

    #[rstest]
    fn test_empty_ailist() {
        let ailist: AIList<u32, &str> = AIList::build(vec![]);

        assert_eq!(ailist.len(), 0);
        assert_eq!(ailist.is_empty(), true);

        let results = ailist.find(1, 2);

        assert_eq!(results.is_empty(), true);
    }

    #[rstest]
    fn test_find_iter_overlapping_intervals(intervals: Vec<Interval<u32, &'static str>>) {
        let ailist = AIList::build(intervals);

        // Query that overlaps with "a" and "b"
        let results: Vec<&Interval<u32, &str>> = ailist.find_iter(2, 4).collect();
        let vals: Vec<&str> = results.iter().map(|i| i.val).collect();
        assert_eq!(vals.contains(&"a"), true);
        assert_eq!(vals.contains(&"b"), true);
        assert_eq!(vals.contains(&"c"), false);
        assert_eq!(results.len(), 2);

        // Query that overlaps with "c" and "d"
        let results: Vec<&Interval<u32, &str>> = ailist.find_iter(9, 11).collect();
        let vals: Vec<&str> = results.iter().map(|i| i.val).collect();
        assert_eq!(vals.contains(&"c"), true);
        assert_eq!(vals.contains(&"d"), true);
        assert_eq!(vals.contains(&"a"), false);
        assert_eq!(results.len(), 2);
    }

    #[rstest]
    fn test_find_iter_no_overlap(intervals: Vec<Interval<u32, &'static str>>) {
        let ailist = AIList::build(intervals);

        // Query outside all intervals
        let results: Vec<&Interval<u32, &str>> = ailist.find_iter(13, 15).collect();
        assert_eq!(results.is_empty(), true);

        // Query before all intervals
        let results: Vec<&Interval<u32, &str>> = ailist.find_iter(0, 1).collect();
        assert_eq!(results.is_empty(), true);
    }

    #[rstest]
    fn test_find_iter_empty_ailist() {
        let ailist: AIList<u32, &str> = AIList::build(vec![]);

        let results: Vec<&Interval<u32, &str>> = ailist.find_iter(1, 2).collect();
        assert_eq!(results.is_empty(), true);
    }

    #[rstest]
    fn test_find_iter_matches_find(intervals: Vec<Interval<u32, &'static str>>) {
        let ailist = AIList::build(intervals);

        // Test multiple queries to ensure find_iter produces same results as find
        let test_queries = vec![(2, 4), (5, 8), (9, 11), (0, 15), (7, 9)];

        for (start, end) in test_queries {
            let find_results = ailist.find(start, end);
            let find_iter_results: Vec<&Interval<u32, &str>> =
                ailist.find_iter(start, end).collect();

            // Check that both methods return the same number of intervals
            assert_eq!(
                find_results.len(),
                find_iter_results.len(),
                "Mismatch in number of results for query ({}, {})",
                start,
                end
            );

            // Check that all intervals from find are present in find_iter results
            for interval in &find_results {
                assert!(
                    find_iter_results.contains(&interval),
                    "Interval {interval:?} from find() not found in find_iter() results",
                );
            }
        }
    }

    #[rstest]
    fn test_find_iter_single_interval() {
        let intervals = vec![Interval {
            start: 5,
            end: 10,
            val: "single",
        }];
        let ailist = AIList::build(intervals);

        // Overlapping query
        let results: Vec<&Interval<u32, &str>> = ailist.find_iter(6, 8).collect();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].val, "single");

        // Non-overlapping query
        let results: Vec<&Interval<u32, &str>> = ailist.find_iter(11, 15).collect();
        assert_eq!(results.is_empty(), true);
    }

    #[rstest]
    fn test_find_iter_complex_interval() {
        let iv = |start: usize, end: usize| -> Interval<usize, ()> {
            Interval {
                start,
                end,
                val: (),
            }
        };
        let intervals = vec![
            iv(0, 30), // Span many ivs
            iv(0, 10),
            iv(0, 10),
            iv(5, 15),
            iv(5, 15),
            iv(10, 20),
            iv(10, 20),
            iv(15, 25),
            iv(15, 25),
            iv(21, 22),
            iv(22, 23),
            iv(20, 30),
            iv(20, 30),
            iv(25, 100), // Span many ivs
            iv(26, 27),
            iv(27, 28),
            iv(29, 30),
            iv(30, 31),
            iv(32, 33),
            iv(50, 51),
            iv(51, 52),
            iv(52, 53),
            iv(53, 54),
            iv(55, 56),
            iv(60, 61),
            iv(70, 71),
        ];

        let ailist = AIList::build(intervals);
        // Confirm we are at least iterating over header list values a bit.
        assert_eq!(ailist.header_list.len(), 2);

        // Overlapping query
        let results: Vec<&Interval<usize, ()>> = ailist.find_iter(6, 8).collect();
        assert_eq!(results.len(), 5);

        let results: Vec<&Interval<usize, ()>> = ailist.find_iter(30, 35).collect();
        assert_eq!(results.len(), 3);

        // Non-overlapping query
        let results: Vec<&Interval<usize, ()>> = ailist.find_iter(101, 150).collect();
        assert_eq!(results.is_empty(), true);
    }
}

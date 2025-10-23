use std::{collections::HashMap, fmt::Debug};

use gtars_core::models::{Interval, RegionSet};
use num_traits::{PrimInt, Unsigned};
use thiserror::Error;

use crate::{AiList, Bits, Overlapper, OverlapperType};

#[derive(Debug, Error)]
pub enum GenomeIndexError {
    #[error("Error parsing region: {0}")]
    RegionParsingError(String),
    #[error("Error converting interval coordinates to u32: start={0}, end={1}")]
    CoordinateConversionError(String, String),
}

pub struct GenomeIndex<I, T> {
    index_maps: HashMap<String, Box<dyn Overlapper<I, T>>>,
    #[allow(dead_code)]
    overlapper_type: OverlapperType,
}

pub struct IterFindOverlaps<'a, 'b, I, T>
where
    I: PrimInt + Unsigned + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    inner: &'a HashMap<String, Box<dyn Overlapper<I, T>>>,
    rs: &'b RegionSet,
    region_idx: usize,
    current_chr: Option<String>,
    current_iter: Option<Box<dyn Iterator<Item = &'a Interval<I, T>> + 'a>>,
}

impl<'a, 'b, I, T> Iterator for IterFindOverlaps<'a, 'b, I, T>
where
    I: PrimInt + Unsigned + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    type Item = (String, &'a Interval<I, T>);

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // first, try to get next item from current iterator
            #[allow(clippy::collapsible_if)]
            if let Some(ref mut iter) = self.current_iter {
                if let Some(interval) = iter.next() {
                    return Some((self.current_chr.as_ref().unwrap().clone(), interval));
                }
            }

            // current iterator exhausted or doesn't exist, move to next region
            if self.region_idx >= self.rs.regions.len() {
                // we are done afte this!
                // this is the terminal point
                return None;
            }

            let region = &self.rs.regions[self.region_idx];
            self.region_idx += 1;

            // try to get overlapper for this chromosome
            if let Some(lapper) = self.inner.get(&region.chr) {
                // convert coordinates
                if let (Some(start), Some(end)) = (I::from(region.start), I::from(region.end)) {
                    self.current_chr = Some(region.chr.clone());
                    self.current_iter = Some(lapper.find_iter(start, end));
                    // continue loop to get first item from new iterator
                } else {
                    // This is a programming error: the GenomeIndex type I cannot represent
                    // the Region's u32 coordinates. This should never happen in practice since
                    // genomic coordinates are u32 and the index should be GenomeIndex<u32, T>.
                    panic!(
                        "Type conversion error: cannot convert Region coordinates to index type. \
                         Region: {}:{}-{}, expected type: {}",
                        region.chr,
                        region.start,
                        region.end,
                        std::any::type_name::<I>()
                    );
                }
            } else {
                // no overlapper for this chromosome, skip to next region
                continue;
            }
        }
    }
}

impl<I, T> GenomeIndex<I, T>
where
    I: PrimInt + Unsigned + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    /// Returns an iterator over all overlapping intervals for the query regions.
    ///
    /// Each item is a tuple of (chromosome, interval reference).
    /// Invalid regions (coordinates that can't convert to type I) are silently skipped.
    pub fn find_overlaps_iter<'a, 'b>(
        &'a self,
        rs: &'b RegionSet,
    ) -> IterFindOverlaps<'a, 'b, I, T> {
        IterFindOverlaps {
            inner: &self.index_maps,
            rs,
            region_idx: 0,
            current_chr: None,
            current_iter: None,
        }
    }

    /// Collect all overlaps into a Vec for convenience.
    ///
    /// This is a helper method that collects the iterator results.
    pub fn find_overlaps(&self, rs: &RegionSet) -> Vec<(String, Interval<I, T>)> {
        self.find_overlaps_iter(rs)
            .map(|(chr, interval)| (chr, interval.clone()))
            .collect()
    }
}

pub trait IntoGenomeIndex<I, T>
where
    I: PrimInt + Unsigned + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    fn into_genome_index(self, overlapper_type: OverlapperType) -> GenomeIndex<I, T>;
}

impl IntoGenomeIndex<u32, Option<String>> for RegionSet {
    fn into_genome_index(
        self,
        overlapper_type: OverlapperType,
    ) -> GenomeIndex<u32, Option<String>> {
        // instantiate the tree and list of intervals
        let mut core: HashMap<String, Box<dyn Overlapper<u32, Option<String>>>> =
            HashMap::default();
        let mut intervals: HashMap<String, Vec<Interval<u32, Option<String>>>> = HashMap::default();

        // STEP 1: filter/organize/sort regions into vectors, one for each chrom
        for region in self.regions.into_iter() {
            // create interval
            let interval = Interval {
                start: region.start,
                end: region.end,
                val: region.rest,
            };

            // use chr to get the vector of intervals
            let chr_intervals = intervals.entry(region.chr.clone()).or_default();

            // push interval to vector
            chr_intervals.push(interval);
        }

        //STEP 2: take each vector (one for each chrom) and build the overlapper
        for (chr, chr_intervals) in intervals.into_iter() {
            let lapper: Box<dyn Overlapper<u32, Option<String>>> = match overlapper_type {
                OverlapperType::Bits => Box::new(Bits::build(chr_intervals)),
                OverlapperType::AiList => Box::new(AiList::build(chr_intervals)),
            };
            core.insert(chr.to_string(), lapper);
        }

        GenomeIndex {
            index_maps: core,
            overlapper_type,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gtars_core::models::Region;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[rstest]
    #[case(OverlapperType::AiList)]
    #[case(OverlapperType::Bits)]
    fn test_basic_overlaps(#[case] overlapper_type: OverlapperType) {
        let regions = vec![
            Region {
                chr: "chr1".to_string(),
                start: 100,
                end: 200,
                rest: None,
            },
            Region {
                chr: "chr1".to_string(),
                start: 300,
                end: 400,
                rest: None,
            },
            Region {
                chr: "chr1".to_string(),
                start: 600,
                end: 800,
                rest: None,
            },
        ];
        let rs = RegionSet::from(regions);
        let gi = rs.into_genome_index(overlapper_type);

        let query = RegionSet::from(vec![Region {
            chr: "chr1".to_string(),
            start: 110,
            end: 210,
            rest: None,
        }]);

        let hits = gi.find_overlaps(&query);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].0, "chr1");
        assert_eq!(hits[0].1.start, 100);
        assert_eq!(hits[0].1.end, 200);
    }

    #[rstest]
    #[case(OverlapperType::AiList)]
    #[case(OverlapperType::Bits)]
    fn test_multiple_overlaps_single_query(#[case] overlapper_type: OverlapperType) {
        let regions = vec![
            Region {
                chr: "chr1".to_string(),
                start: 100,
                end: 200,
                rest: None,
            },
            Region {
                chr: "chr1".to_string(),
                start: 150,
                end: 250,
                rest: None,
            },
            Region {
                chr: "chr1".to_string(),
                start: 180,
                end: 300,
                rest: None,
            },
        ];
        let rs = RegionSet::from(regions);
        let gi = rs.into_genome_index(overlapper_type);

        let query = RegionSet::from(vec![Region {
            chr: "chr1".to_string(),
            start: 160,
            end: 190,
            rest: None,
        }]);

        let hits = gi.find_overlaps(&query);
        assert_eq!(hits.len(), 3);
    }

    #[rstest]
    #[case(OverlapperType::AiList)]
    #[case(OverlapperType::Bits)]
    fn test_no_overlaps(#[case] overlapper_type: OverlapperType) {
        let regions = vec![
            Region {
                chr: "chr1".to_string(),
                start: 100,
                end: 200,
                rest: None,
            },
            Region {
                chr: "chr1".to_string(),
                start: 300,
                end: 400,
                rest: None,
            },
        ];
        let rs = RegionSet::from(regions);
        let gi = rs.into_genome_index(overlapper_type);

        let query = RegionSet::from(vec![Region {
            chr: "chr1".to_string(),
            start: 500,
            end: 600,
            rest: None,
        }]);

        let hits = gi.find_overlaps(&query);
        assert_eq!(hits.len(), 0);
    }

    #[rstest]
    #[case(OverlapperType::AiList)]
    #[case(OverlapperType::Bits)]
    fn test_multiple_chromosomes(#[case] overlapper_type: OverlapperType) {
        let regions = vec![
            Region {
                chr: "chr1".to_string(),
                start: 100,
                end: 200,
                rest: None,
            },
            Region {
                chr: "chr2".to_string(),
                start: 300,
                end: 400,
                rest: None,
            },
            Region {
                chr: "chr3".to_string(),
                start: 500,
                end: 600,
                rest: None,
            },
        ];
        let rs = RegionSet::from(regions);
        let gi = rs.into_genome_index(overlapper_type);

        let query = RegionSet::from(vec![
            Region {
                chr: "chr1".to_string(),
                start: 150,
                end: 250,
                rest: None,
            },
            Region {
                chr: "chr2".to_string(),
                start: 350,
                end: 450,
                rest: None,
            },
        ]);

        let hits = gi.find_overlaps(&query);
        assert_eq!(hits.len(), 2);

        let chr1_hits: Vec<_> = hits.iter().filter(|(chr, _)| chr == "chr1").collect();
        let chr2_hits: Vec<_> = hits.iter().filter(|(chr, _)| chr == "chr2").collect();

        assert_eq!(chr1_hits.len(), 1);
        assert_eq!(chr2_hits.len(), 1);
    }

    #[rstest]
    #[case(OverlapperType::AiList)]
    #[case(OverlapperType::Bits)]
    fn test_exact_boundary_overlaps(#[case] overlapper_type: OverlapperType) {
        let regions = vec![Region {
            chr: "chr1".to_string(),
            start: 100,
            end: 200,
            rest: None,
        }];
        let rs = RegionSet::from(regions);
        let gi = rs.into_genome_index(overlapper_type);

        // Query starts exactly at region end
        let query = RegionSet::from(vec![Region {
            chr: "chr1".to_string(),
            start: 200,
            end: 300,
            rest: None,
        }]);

        let hits = gi.find_overlaps(&query);
        // Typically intervals are half-open [start, end), so start=200 shouldn't overlap with end=200
        assert_eq!(hits.len(), 0);
    }

    #[rstest]
    #[case(OverlapperType::AiList)]
    #[case(OverlapperType::Bits)]
    fn test_empty_query(#[case] overlapper_type: OverlapperType) {
        let regions = vec![Region {
            chr: "chr1".to_string(),
            start: 100,
            end: 200,
            rest: None,
        }];
        let rs = RegionSet::from(regions);
        let gi = rs.into_genome_index(overlapper_type);

        let query = RegionSet::from(vec![]);
        let hits = gi.find_overlaps(&query);
        assert_eq!(hits.len(), 0);
    }

    #[rstest]
    #[case(OverlapperType::AiList)]
    #[case(OverlapperType::Bits)]
    fn test_query_nonexistent_chromosome(#[case] overlapper_type: OverlapperType) {
        let regions = vec![Region {
            chr: "chr1".to_string(),
            start: 100,
            end: 200,
            rest: None,
        }];
        let rs = RegionSet::from(regions);
        let gi = rs.into_genome_index(overlapper_type);

        let query = RegionSet::from(vec![Region {
            chr: "chr99".to_string(),
            start: 100,
            end: 200,
            rest: None,
        }]);

        let hits = gi.find_overlaps(&query);
        assert_eq!(hits.len(), 0);
    }

    #[rstest]
    #[case(OverlapperType::AiList)]
    #[case(OverlapperType::Bits)]
    fn test_with_metadata(#[case] overlapper_type: OverlapperType) {
        let regions = vec![
            Region {
                chr: "chr1".to_string(),
                start: 100,
                end: 200,
                rest: Some("gene_a".to_string()),
            },
            Region {
                chr: "chr1".to_string(),
                start: 300,
                end: 400,
                rest: Some("gene_b".to_string()),
            },
        ];
        let rs = RegionSet::from(regions);
        let gi = rs.into_genome_index(overlapper_type);

        let query = RegionSet::from(vec![Region {
            chr: "chr1".to_string(),
            start: 150,
            end: 250,
            rest: None,
        }]);

        let hits = gi.find_overlaps(&query);
        assert_eq!(hits.len(), 1);
        assert!(hits[0].1.val.is_some());
        assert_eq!(hits[0].1.val.as_ref().unwrap(), "gene_a");
    }

    #[rstest]
    #[case(OverlapperType::AiList)]
    #[case(OverlapperType::Bits)]
    fn test_overlapping_query_regions(#[case] overlapper_type: OverlapperType) {
        let regions = vec![
            Region {
                chr: "chr1".to_string(),
                start: 100,
                end: 200,
                rest: None,
            },
            Region {
                chr: "chr1".to_string(),
                start: 300,
                end: 400,
                rest: None,
            },
        ];
        let rs = RegionSet::from(regions);
        let gi = rs.into_genome_index(overlapper_type);

        // Two query regions that each hit different index regions
        let query = RegionSet::from(vec![
            Region {
                chr: "chr1".to_string(),
                start: 150,
                end: 250,
                rest: None,
            },
            Region {
                chr: "chr1".to_string(),
                start: 350,
                end: 450,
                rest: None,
            },
        ]);

        let hits = gi.find_overlaps(&query);
        assert_eq!(hits.len(), 2);
    }
}

use std::{collections::HashMap, fmt::Debug};

use gtars_core::models::{Interval, RegionSet, Region};
use thiserror::Error;
use num_traits::{PrimInt, Unsigned};

use crate::{Overlapper, OverlapperType, Bits, AiList};

#[derive(Debug, Error)]
pub enum GenomeIndexError {
    #[error("Error parsing region: {0}")]
    RegionParsingError(String),
    #[error("Error converting interval coordinates to u32: start={0}, end={1}")]
    CoordinateConversionError(String, String),
}

pub struct GenomeIndex<I, T> {
    index_maps: HashMap<String, Box<dyn Overlapper<I, T>>>,
    overlapper_type: OverlapperType,
}

impl<I, T> GenomeIndex<I, T>
where
    I: PrimInt + Unsigned + Send + Sync + Debug,
    T: Eq + Clone + Send + Sync + Debug,
{
    pub fn find_overlaps(&self, rs: &RegionSet) -> Result<Vec<Region>, GenomeIndexError> {
        let mut final_hits = Vec::new();
        for r in rs {
            let lapper = self.index_maps.get(&r.chr);
            match lapper {
                Some(lapper) => {
                    let start = I::from(r.start);
                    let end = I::from(r.end);
                    if let (Some(start), Some(end)) = (start, end) {
                        for iv in lapper.find_iter(start, end) {
                            let start_u32 = iv.start.to_u32().ok_or_else(|| {
                                GenomeIndexError::CoordinateConversionError(
                                    format!("{:?}", iv.start),
                                    format!("{:?}", iv.end)
                                )
                            })?;
                            let end_u32 = iv.end.to_u32().ok_or_else(|| {
                                GenomeIndexError::CoordinateConversionError(
                                    format!("{:?}", iv.start),
                                    format!("{:?}", iv.end)
                                )
                            })?;
                            let rest = format!("{:?}", iv.val);
                            final_hits.push(Region {
                                chr: r.chr.clone(),
                                start: start_u32,
                                end: end_u32,
                                rest: Some(rest)
                            });
                        }
                    } else {
                        return Err(GenomeIndexError::RegionParsingError(
                            format!("Could not parse region start and end: {r}")
                        ));
                    }
                    
                }
                None => continue
            }
        }
        Ok(final_hits)
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
    fn into_genome_index(self, overlapper_type: OverlapperType) -> GenomeIndex<u32, Option<String>> {
        // instantiate the tree and list of intervals
        let mut core: HashMap<String, Box<dyn Overlapper<u32, Option<String>>>> = HashMap::default();
        let mut intervals: HashMap<String, Vec<Interval<u32, Option<String>>>> = HashMap::default();

        // STEP 1: filter/organize/sort regions into vectors, one for each chrom
        for region in self.regions.into_iter() {

            // create interval
            let interval = Interval { start: region.start, end: region.end, val: region.rest };

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

    #[test]
    fn test_basic() {
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
        let gi = rs.into_genome_index(OverlapperType::AiList);

        
        let query = RegionSet::from(vec![
            Region {
                chr: "chr1".to_string(),
                start: 110,
                end: 210,
                rest: None,
            },
            Region {
                chr: "chr1".to_string(),
                start: 120,
                end: 220,
                rest: None,
            },
        ]);
        let hits = gi.find_overlaps(&query);
        println!("{:?}", hits.unwrap());
        // assert_eq!(hits.is_ok(), true);
        // assert_eq!(hits.unwrap().len(), 2);
        // let hits = hits.unwrap();
        // println!("{:?}", hits.clone());
    }
}

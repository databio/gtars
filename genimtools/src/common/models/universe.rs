use std::collections::HashMap;
use std::path::Path;

use crate::common::consts::{PAD_CHR, PAD_END, PAD_START, UNKNOWN_CHR, UNKNOWN_END, UNKNOWN_START};
use crate::common::models::region::Region;
use crate::common::utils::{extract_regions_from_bed_file, generate_region_to_id_map};

pub struct Universe {
    pub regions: Vec<Region>,
    pub region_to_id: HashMap<Region, u32>,
    length: u32,
}

impl Universe {
    pub fn convert_region_to_id(&self, region: &Region) -> u32 {
        let id = self.region_to_id.get(region);
        match id {
            Some(id) => id.to_owned(),
            None => self
                .region_to_id
                .get(&Region {
                    chr: UNKNOWN_CHR.to_string(),
                    start: UNKNOWN_START as u32,
                    end: UNKNOWN_END as u32,
                })
                .unwrap()
                .to_owned(),
        }
    }

    pub fn convert_chr_start_end_to_id(&self, chr: &str, start: u32, end: u32) -> u32 {
        let region = Region {
            chr: chr.to_string(),
            start,
            end,
        };
        self.convert_region_to_id(&region)
    }

    pub fn len(&self) -> u32 {
        self.length
    }

    pub fn is_empty(&self) -> bool {
        self.length == 0
    }

    pub fn padding_token(&self) -> Region {
        Region {
            chr: PAD_CHR.to_string(),
            start: PAD_START as u32,
            end: PAD_END as u32,
        }
    }

    pub fn unknown_token(&self) -> Region {
        Region {
            chr: UNKNOWN_CHR.to_string(),
            start: UNKNOWN_START as u32,
            end: UNKNOWN_END as u32,
        }
    }
}

impl From<Vec<Region>> for Universe {
    fn from(value: Vec<Region>) -> Self {
        // make a copy of the regions
        let mut regions = value.clone();

        // create the region to id map and add the Unk token if it doesn't exist
        let mut region_to_id = generate_region_to_id_map(&regions);
        let total_regions = region_to_id.len();

        // add Unk and Pad token if they doesn't exist
        // its possible the vocab file passed
        // in does have the Unk/Pad token, but
        // we don't know that here
        let unk = Region {
            chr: UNKNOWN_CHR.to_string(),
            start: UNKNOWN_START as u32,
            end: UNKNOWN_END as u32,
        };

        // do the same with the Pad token
        let pad = Region {
            chr: PAD_CHR.to_string(),
            start: PAD_START as u32,
            end: PAD_END as u32,
        };

        // notify if the Unk token is already in the vocab
        if region_to_id.contains_key(&unk) {
            // pass
        } else {
            regions.push(unk.to_owned());
        }

        // notify if the Pad token is already in the vocab
        if region_to_id.contains_key(&pad) {
            // pass
        } else {
            regions.push(pad.to_owned());
        }

        // add the Unk token to the region to id map
        region_to_id.entry(unk).or_insert(total_regions as u32);
        let total_regions = region_to_id.len();

        // add the Pad token to the region to id map
        region_to_id.entry(pad).or_insert(total_regions as u32);
        let total_regions = region_to_id.len();

        Universe {
            regions,
            region_to_id,
            length: total_regions as u32,
        }
    }
}

impl From<&Path> for Universe {
    fn from(value: &Path) -> Self {
        let regions = extract_regions_from_bed_file(value);

        let mut regions = match regions {
            Ok(r) => r,
            // should probably change this to something else,
            // but couldn't figure out how to return a `Result`
            // from a trait implementation
            Err(e) => panic!("{e}"),
        };

        let mut region_to_id = generate_region_to_id_map(&regions);
        let total_regions = region_to_id.len();

        // add Unk and Pad token if they doesn't exist
        // its possible the vocab file passed
        // in does have the Unk/Pad token, but
        // we don't know that here
        let unk = Region {
            chr: UNKNOWN_CHR.to_string(),
            start: UNKNOWN_START as u32,
            end: UNKNOWN_END as u32,
        };

        // do the same with the Pad token
        let pad = Region {
            chr: PAD_CHR.to_string(),
            start: PAD_START as u32,
            end: PAD_END as u32,
        };

        // notify if the Unk token is already in the vocab
        if region_to_id.contains_key(&unk) {
            // pass
        } else {
            regions.push(unk.to_owned());
        }

        // notify if the Pad token is already in the vocab
        if region_to_id.contains_key(&pad) {
            // pass
        } else {
            regions.push(pad.to_owned());
        }

        // add the Unk token to the region to id map
        region_to_id.entry(unk).or_insert(total_regions as u32);
        let total_regions = region_to_id.len();

        // add the Pad token to the region to id map
        region_to_id.entry(pad).or_insert(total_regions as u32);
        let total_regions = region_to_id.len();

        Universe {
            regions,
            region_to_id,
            length: total_regions as u32,
        }
    }
}

use std::str::FromStr;

use anyhow::Result;

use crate::common::models::Region;

#[allow(unused)]
pub struct Fragment {
    pub chr: String,
    pub start: u32,
    pub end: u32,
    pub barcode: String,
    pub read_support: u32,
}

impl FromStr for Fragment {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        let parts: Vec<&str> = s.split_whitespace().collect();
        // dont check file integrity right now
        // if parts.len() != 6 {
        //     anyhow::bail!(
        //         "Error parsing fragment file line: {}. Is your fragment file malformed? Found {} parts.",
        //         s,
        //         parts.len()
        //     )
        // }

        let start = parts[1].parse::<u32>()?;
        let end = parts[2].parse::<u32>()?;
        let read_support = parts[4].parse::<u32>()?;

        Ok(Fragment {
            chr: parts[0].to_string(),
            start,
            end,
            barcode: parts[3].to_string(),
            read_support,
        })
    }
}

impl From<Fragment> for Region {
    fn from(val: Fragment) -> Self {
        Region {
            chr: val.chr,
            start: val.start,
            end: val.end,
        }
    }
}

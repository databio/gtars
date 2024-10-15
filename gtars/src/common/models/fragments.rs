use std::str::FromStr;

use anyhow::Result;

pub struct Fragment {
    chr: String,
    start: u32,
    end: u32,
    barcode: String,
    read_support: u32,
}

impl FromStr for Fragment {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        let parts: Vec<&str> = s.split_whitespace().collect();
        if parts.len() != 5 {
            anyhow::bail!("Error parsing fragment file line: {}. Is your fragment file malformed?", s)
        }

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

use std::str::FromStr;

use anyhow::Error;

pub enum ScoringMode {
    Atac,
    Chip,
}

impl FromStr for ScoringMode {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "atac" => Ok(ScoringMode::Atac),
            "chip" => Ok(ScoringMode::Chip),
            _ => Err(Error::msg(format!("Invalid scoring mode: {}", s))),
        }
    }
}

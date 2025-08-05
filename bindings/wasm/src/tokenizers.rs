use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct BedEntry {
    pub chr: String,
    pub start: u32,
    pub end: u32,
}

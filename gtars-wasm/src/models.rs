use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub struct BedEntries(pub Vec<(String, u32, u32)>);

#[derive(Serialize, Deserialize)]
pub struct BedEntriesFull(pub Vec<(String, u32, u32, String)>);
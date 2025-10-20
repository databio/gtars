use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub struct BedEntries(pub Vec<(String, u32, u32)>);

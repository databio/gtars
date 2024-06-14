use std::collections::HashMap;

use rust_lapper::{Interval, Lapper};

use crate::common::models::Universe;

pub struct MetaTokenizer {
    pub universe: Universe,
    tree: HashMap<String, Lapper<u32, u32>>,
    secondary_trees: Option<Vec<HashMap<String, Lapper<u32, u32>>>>,
}

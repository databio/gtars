use std::collections::HashMap;

use rust_lapper::{Lapper, Interval};

use crate::common::models::Universe;

pub struct MetaTokenizer {
    pub universe: Universe,
    tree: HashMap<String, Lapper<u32,u32>>
}
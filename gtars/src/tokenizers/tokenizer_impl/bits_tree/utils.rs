use std::collections::HashMap;

use rust_lapper::{Interval, Lapper};

use crate::tokenizers::universe::Universe;

///
/// Simple wrapper function that will create a [Lapper] object (an interval tree)
/// from a [Universe] struct.
///
/// # Arguments:
/// - universe: the universe to create the interval tree for.
pub fn create_interval_tree_from_universe(
    universe: &Universe,
) -> HashMap<String, Lapper<u32, u32>> {
    // instantiate the tree and list of intervals
    let mut tree: HashMap<String, Lapper<u32, u32>> = HashMap::new();
    let mut intervals: HashMap<String, Vec<Interval<u32, u32>>> = HashMap::new();

    for region in universe.regions.iter() {
        // create interval
        let interval = Interval {
            start: region.start,
            stop: region.end,
            val: universe.convert_token_to_id(region).unwrap(),
        };

        // use chr to get the vector of intervals
        let chr_intervals = intervals.entry(region.chr.clone()).or_default();

        // push interval to vector
        chr_intervals.push(interval);
    }

    // build the tree
    for (chr, chr_intervals) in intervals.into_iter() {
        let lapper: Lapper<u32, u32> = Lapper::new(chr_intervals);
        tree.insert(chr.to_string(), lapper);
    }

    tree
}

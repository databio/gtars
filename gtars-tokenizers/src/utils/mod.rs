//!
//! This module contains utility functions for tokenizers. Basic things
//! like universe prapration and special token handling are done here.
//!
use std::path::Path;

use fxhash::FxHashMap as HashMap;

use special_tokens::SpecialTokens;

use gtars_core::models::Interval;
use gtars_overlaprs::{AIList, Bits, Overlapper};

use super::TokenizerError;
use super::config::TokenizerType;
use super::universe::Universe;

pub mod fragments;
pub mod special_tokens;

// Re-export fragment functions
pub use fragments::{count_fragments_by_barcode, tokenize_fragment_file};

///
/// Prepare the universe and special tokens. This function will build
/// the universe struct and prepare the special tokens if they are provided.
///
/// Doing these together is necessary, because the special tokens contribute
/// to the universe/vocab.
///
/// # Arguments:
/// - config: the tokenizer config
///
pub fn prepare_universe_and_special_tokens<P: AsRef<Path>>(
    universe_file: P,
    special_tokens: SpecialTokens,
) -> Result<(Universe, SpecialTokens), TokenizerError> {
    let mut universe = Universe::try_from(universe_file.as_ref())?;
    universe.add_special_tokens(&special_tokens);
    Ok((universe, special_tokens))
}

///
/// Simple wrapper function that will create a [Lapper] object (an interval tree)
/// from a [Universe] struct.
///
/// # Arguments:
/// - universe: the universe to create the interval tree for.
pub fn create_tokenize_core_from_universe(
    universe: &Universe,
    overlapper_type: TokenizerType,
) -> HashMap<String, Box<dyn Overlapper<u32, u32>>> {
    // instantiate the tree and list of intervals
    let mut core: HashMap<String, Box<dyn Overlapper<u32, u32>>> = HashMap::default();
    let mut intervals: HashMap<String, Vec<Interval<u32, u32>>> = HashMap::default();

    for region in universe.regions.iter() {
        // skip any special tokens that snuck into the regions
        if let Some(special_tokens) = &universe.special_tokens {
            if special_tokens.contains(&region.to_string()) {
                continue;
            }
        }

        let parts = region.split(":").collect::<Vec<&str>>();

        let chr = parts[0].to_string();
        let start_end = parts[1];

        let start_end_parts = start_end.split("-").collect::<Vec<&str>>();
        let start = start_end_parts[0];
        let end = start_end_parts[1];

        let start = start.parse::<u32>().unwrap();
        let end = end.parse::<u32>().unwrap();

        let val = universe.convert_token_to_id(region).unwrap();

        // create interval
        let interval = Interval { start, end, val };

        // use chr to get the vector of intervals
        let chr_intervals = intervals.entry(chr.clone()).or_default();

        // push interval to vector
        chr_intervals.push(interval);
    }

    // build the chromosome to tree mapping
    for (chr, chr_intervals) in intervals.into_iter() {
        let lapper: Box<dyn Overlapper<u32, u32>> = match overlapper_type {
            TokenizerType::Bits => Box::new(Bits::build(chr_intervals)),
            TokenizerType::AIList => Box::new(AIList::build(chr_intervals)),
        };
        core.insert(chr.to_string(), lapper);
    }

    core
}

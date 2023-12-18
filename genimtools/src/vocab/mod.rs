use std::collections::HashMap;
use std::fs;
use std::path::Path;

use indicatif::{ProgressBar, ProgressStyle};

use crate::common::consts::{UNKNOWN_CHR, UNKNOWN_END, UNKNOWN_START};
use crate::common::models::{Region, RegionSet};
use crate::tokenizers::{Tokenizer, TreeTokenizer};

pub mod cli;

pub mod consts {
    pub const FILE_EXTENSION: &str = ".bed";
    pub const DEFAULT_OUTPUT: &str = "output.bed";
    pub const DEFAULT_MIN_COUNT: u32 = 0;
    pub const VERSION: &str = env!("CARGO_PKG_VERSION");
    pub const PRUNE_CMD: &str = "prune";
}

///
/// Creates a count map for how many times each universe token appears in the training data.
///
/// # Arguments
/// - `data_path` - Path to the training data. This should be a folder of bed files.
/// - `universe_path` - Path to the universe file we want to prune.
/// - `output_path` - Path to the output file (the pruned universe).
pub fn create_count_map(
    data_path: &str,
    universe_path: &str,
) -> Result<HashMap<Region, u32>, Box<dyn std::error::Error>> {
    // set up the tokenizer
    let universe_path = Path::new(universe_path);
    let tokenizer = TreeTokenizer::from(universe_path);

    // set up the counter
    let mut counter: HashMap<Region, u32> = HashMap::new();

    // get a list of paths
    let data_path = Path::new(data_path);
    let paths = fs::read_dir(data_path)?;

    // iterate over the paths
    let num_paths = paths.count();
    let bar = ProgressBar::new(num_paths as u64);
    bar.set_style(
        ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
        )
        .unwrap()
        .progress_chars("##-"),
    );

    let paths = fs::read_dir(data_path)?;
    for path in paths {
        if let Some(extension) = path.as_ref().unwrap().path().extension() {
            if extension != consts::FILE_EXTENSION.trim_start_matches('.') {
                continue;
            }
        }

        // get regions from the file
        let path = path?.path();
        let path = Path::new(&path);

        let region_set = RegionSet::try_from(path)?;

        // skip if there are no regions
        if region_set.is_empty() {
            continue;
        }

        let tokens = tokenizer.tokenize_region_set(&region_set).unwrap();

        // count the tokens
        tokens.into_iter().for_each(|region| {
            let region = Region {
                chr: region.chr,
                start: region.start,
                end: region.end,
            };

            let cnt = counter.entry(region).or_insert_with(|| 0);
            *cnt += 1;
        });

        bar.inc(1);
    }

    // drop unknown token from hashmap, that is not needed
    counter.remove(&Region {
        chr: UNKNOWN_CHR.to_string(),
        start: UNKNOWN_START as u32,
        end: UNKNOWN_END as u32,
    });

    Ok(counter)
}

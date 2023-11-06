use std::path::Path;
use std::fs;
use std::collections::HashMap;

use indicatif::{ProgressBar, ProgressStyle};
use gtokenizers::tokenizers::TreeTokenizer;
use gtokenizers::models::region::Region;
use gtokenizers::models::region_set::RegionSet;
use gtokenizers::io::extract_regions_from_bed_file;
use gtokenizers::tokenizers::traits::{Tokenizer, UNKNOWN_CHR, UNKNOWN_START, UNKNOWN_END};

pub mod consts {
    pub const FILE_EXTENSION: &str = ".bed";
    pub const DEFAULT_OUTPUT: &str = "output.bed";
    pub const DEFAULT_MIN_COUNT: u32 = 0;
    pub const VERSION: &str = env!("CARGO_PKG_VERSION");
}

///
/// Creates a count map for how many times each universe token appears in the training data.
/// 
/// # Arguments
/// - `data_path` - Path to the training data. This should be a folder of bed files.
/// - `universe_path` - Path to the universe file we want to prune.
/// - `output_path` - Path to the output file (the pruned universe).
pub fn create_count_map(data_path: &str, universe_path: &str) -> Result<HashMap<Region, u32>, Box<dyn std::error::Error>> {

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
    bar.set_style(ProgressStyle::with_template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}")
        .unwrap()
        .progress_chars("##-"));

    let paths = fs::read_dir(data_path)?;
    for path in paths {

        if let Some(extension) = path.as_ref().unwrap().path().extension() {
            if extension != crate::consts::FILE_EXTENSION.trim_start_matches('.') {
                continue;
            }
        }

        // get regions from the file
        let path = path?.path();
        let path = Path::new(&path);
        let regions = extract_regions_from_bed_file(path)?;

        // skip if there are no regions
        if regions.is_empty() {
            continue;
        }
        
        let region_set = RegionSet::from(regions);
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
        chr: String::from(UNKNOWN_CHR),
        start: UNKNOWN_START as u32,
        end: UNKNOWN_END as u32,
    });

            
    Ok(counter)
}

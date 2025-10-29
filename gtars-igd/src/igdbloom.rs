use std::collections::HashMap;

use std::fs;
use std::io::{Error, Read, Write};
use std::path::Path;

use gtars_core::models::{Region, RegionSet};
use gtars_tokenizers::tokenizer::Tokenizer;

#[cfg(feature = "bloom")]
use bloomfilter::Bloom;

#[cfg(feature = "bloom")]
pub fn tokenize_then_create_bloom_for_each_file(universe_tokenizer: &Tokenizer, bed_file: &str, child_directory: &str, num_of_items: usize, false_positive_rate: f64){
    // we must first tokenize against a universe and then create a bloom filter for each chromosome
    // from that tokenization

    //TODO implement random generation of seed and create filters from this seed.
    let mut seed = [0u8; 32];

    // First load regions from bed file
    let bed = Path::new(&bed_file);
    let regions = RegionSet::try_from(bed).unwrap();

    let path = Path::new(bed_file);


    let filename = path.file_name()

        .and_then(|os_str| os_str.to_str()).unwrap();

    // Create bloom filter for each chromosome

    let bloom_filter_path = format!("{}/{}.bloom", child_directory, filename);

fn file_exists(path: &str) -> bool {
    Path::new(path).exists() && Path::new(path).is_file()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::path::PathBuf;

    #[rstest]
    fn test_bloom_filter() {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");

        let tempbedpath = PathBuf::from(path_to_crate)
            .parent()
            .unwrap()
            .join("tests/data/dir_of_files/dir_beds/dummy2.bed");
        let bed_path = tempbedpath.to_string_lossy();

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let child_directory = path.into_os_string().into_string().unwrap();
        let num_of_items = 1000;
        let false_positive_rate = 0.5;

        let tokenizer =
            Tokenizer::from_auto(bed_path.as_ref()).expect("Failed to create tokenizer from config.");

        tokenize_then_create_bloom_for_each_file(&tokenizer, &bed_path, &child_directory, num_of_items, false_positive_rate);
    }
}
}



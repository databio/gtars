use std::collections::HashMap;
use std::env::current_exe;
use std::fs;
use std::fs::File;
use std::io::{Error, Read, Write};
use std::path::Path;
use anyhow::Context;
// for now use simple implementation of bloom filter in rust
use bloomfilter;
use bloomfilter::Bloom;
use crate::common::models::{Region, RegionSet};
use crate::tokenizers::{Tokenizer, TreeTokenizer};

pub fn create_bloom_filter_main(){

    let universe_path  ="/home/drc/Downloads/bloom_testing/real_data/data/universe.merged.pruned.filtered100k.bed";

    // Just create the universe once
    let universe_path = Path::new(&universe_path);
    let universe_tree_tokenizer = TreeTokenizer::try_from(universe_path).unwrap();

    let num_of_items = 10000;
    let false_positive_rate = 0.001;
    let seed = 42;

    // make parent directory to hold sub directories
    let save_path ="/home/drc/Downloads/bloom_testing/test1/";
    let name = "test";
    let parent_directory = format!("{}{}/",save_path.clone(),name.clone());
    make_parent_directory(parent_directory.as_str()).unwrap();

    let bed_directory = "/home/drc/Downloads/bloom_testing/test1/all_bed_files/";
    let mut child_meta_data_map: HashMap<String, String> = HashMap::new(); // File name and directory to store files
    // For now we will store the individual bloom filters in child directories
    make_child_directories(parent_directory, bed_directory, &mut child_meta_data_map);

    //DEBUG CHECK TO ENSURE WE CHANGED HASHMAP
    for (key, value) in child_meta_data_map {  // &map to borrow, not take ownership
        println!("Key: {}, Value: {}", key, value); // key is child directory  value is the absolute path to the bed file

        tokenize_then_create_bloom(&universe_tree_tokenizer, value.as_str(), key.as_str(), num_of_items, false_positive_rate)




    }

    // TODO
    // TODO Based on collected BED files and assuming your universe was created using those same bed files
    // TODO tokenize the bed files and create individual bloom filters within each folder using the tokenized data

    //
    //
    // let example_region = "chr1|500|800";
    //
    // // Create Empty Bloom Filter
    // println!("Creating Bloom Filter");
    // //let igd_bloom_filter = bloomfilter::Bloom::new_for_fp_rate_with_seed(num_of_items, false_positive_rate, &[]);
    // let mut igd_bloom_filter = bloomfilter::Bloom::new_for_fp_rate(num_of_items, false_positive_rate).unwrap();
    // igd_bloom_filter.set(&example_region);   // insert 10 in the bloom filter
    // // igd_bloom_filter.check(&10); // return true
    // // igd_bloom_filter.check(&20); // return false
    //
    // let is_empty =  igd_bloom_filter.is_empty();
    // let num_hashes =  igd_bloom_filter.number_of_hash_functions();
    // let bloom_len =  igd_bloom_filter.len();
    //
    //
    // println!("Bloom empty: {}", is_empty);
    // println!("Num hashes: {}", num_hashes);
    // println!("Bloom length: {}", bloom_len);
    //
    // let final_save_path = format!("{}{}.bloom",save_path.clone(),name.clone());
    //
    // write_bloom_filter_to_disk(igd_bloom_filter, final_save_path.clone());
    // load_bloom_filter_from_disk(final_save_path.clone());

}

fn make_child_directories(parent_directory: String, bed_directory: &str, meta_data: &mut HashMap<String, String>) {

    //let mut bed_files = Vec::new();
    for entry in fs::read_dir(bed_directory).unwrap() {
        let entry = entry.unwrap();
        let path = entry.path();

        if path.is_file() {
            if let Some(extension) = path.extension() {
                if extension == "bed" {
                    let full_name = path.to_str().unwrap().to_string();
                    println!("Here is the full name: {}", full_name.clone());

                    let name_without_extension = path
                        .file_stem()
                        .unwrap()
                        .to_str()
                        .unwrap()
                        .to_string();
                    println!("Here is the name without extension: {}", name_without_extension.clone());


                    let single_parent_directory = format!("{}{}/",parent_directory,name_without_extension);
                    make_parent_directory(single_parent_directory.as_str());

                    meta_data.insert(single_parent_directory,full_name.clone());

                }
            }
        }
    }


}

pub fn tokenize_then_create_bloom(universe_tree_tokenizer: &TreeTokenizer, bed_file: &str, child_directory: &str, num_of_items: usize, false_positive_rate: f64){
    // we must first tokenize against a universe and then create a bloom filter for each chromosome
    // from that tokenization


    // First tokenize regions
    let bed = Path::new(&bed_file);
    let regions = RegionSet::try_from(bed)
        .with_context(|| "There was an error reading in the bedfile to be tokenized!").unwrap();
    let tokenized_regions = universe_tree_tokenizer.tokenize_region_set(&regions);

    let mut tokenized_regions_iter = tokenized_regions.into_iter();
    let first_region: Region  = tokenized_regions_iter.next().unwrap().into();
    let chr_copy = first_region.chr.clone();
    let bloom_filter_path = format!("{}{}.bloom",child_directory, chr_copy);
    println!("Initial Region: Searching for bloomfilter: {}", bloom_filter_path);

    let mut current_bloom_filter: Bloom<String> = bloomfilter::Bloom::new_for_fp_rate(num_of_items, false_positive_rate).unwrap();

    if file_exists(&*bloom_filter_path) {
        println!("File exists!");
        current_bloom_filter = load_bloom_filter_from_disk(bloom_filter_path);
    } else {
        println!("File does not exist! creating bloom filter in memory");
        current_bloom_filter = bloomfilter::Bloom::new_for_fp_rate(num_of_items, false_positive_rate).unwrap();
    }

    // add initial region
    let line = format!("{}|{}|{}", first_region.chr.clone(), first_region.start.clone(), first_region.end.clone());
    current_bloom_filter.set(&line);

    // set current to previous
    let mut previous_chrom: String = chr_copy.clone();


    // Create bloom filters
    for tokenized_region in tokenized_regions_iter {
        let region: Region = tokenized_region.into();
        let chr_copy = region.chr.clone();

        if chr_copy != previous_chrom{
            println!("Switching chromsomes...{}  {}", chr_copy, previous_chrom);

            // save any progress made to disk
            let previous_bloom_filter_path = format!("{}{}.bloom",child_directory, previous_chrom);
            write_bloom_filter_to_disk(current_bloom_filter, previous_bloom_filter_path);

            previous_chrom = chr_copy.clone();

            // Look for bloom filter at new location:
            let bloom_filter_path = format!("{}{}.bloom",child_directory, chr_copy);
            println!("Searching for bloomfilter: {}", bloom_filter_path);
            if file_exists(&*bloom_filter_path) {
                println!("File exists!");
                current_bloom_filter = load_bloom_filter_from_disk(bloom_filter_path);
            } else {
                println!("File does not exist! creating bloom filter in memory");
                current_bloom_filter = bloomfilter::Bloom::new_for_fp_rate(num_of_items, false_positive_rate).unwrap();
            }
        }

        let line = format!("{}|{}|{}", region.chr.clone(), region.start.clone(), region.end.clone());
        current_bloom_filter.set(&line);

       println!("{}",line);
    }




}

pub fn make_parent_directory(parent_directory: &str) -> Result<(), Error> {


    let parent_path = Path::new(&parent_directory);

    if !parent_path.exists() {

        match fs::create_dir_all(&parent_directory) {
            Ok(_) => {println!("Parent directory created successfully: {}", parent_directory);
                Ok(())
            }
            ,
            Err(e) => {
                eprintln!("Error creating parent directory: {}", e);
                return Err(e); // Example: Returning the error
            }
        }
    } else {
        println!("Parent directory already exists: {}", parent_directory);
        Ok(())
    }
}

pub fn write_bloom_filter_to_disk(igd_bloom_filter: Bloom<String>, save_path: String){

    println!("Here is final save path: {}", save_path);

    let slice = igd_bloom_filter.as_slice();

    let mut file = File::create(save_path).unwrap();
    file.write_all(slice);
    file.flush().unwrap();


}


pub fn search_bloom_filter(){

    let path_to_bloom = "";

    let universe_file = "";

    let query_bed_file = "";


}

pub fn load_bloom_filter_from_disk(load_path: String) -> Bloom<String>
{
    println!("Loading Bloom Filter");
    let mut file = File::open(load_path).unwrap();

    let mut buffer = Vec::new();

    file.read_to_end(&mut buffer).unwrap();
    if buffer.is_empty() {
       eprintln!("FILE IS EMPTY");
    }

    let loaded_igd_bloom_filter: Bloom<String> = Bloom::from_slice(&*buffer).unwrap();
    let is_empty = loaded_igd_bloom_filter.is_empty();
    let num_hashes = loaded_igd_bloom_filter.number_of_hash_functions();
    let bloom_len = loaded_igd_bloom_filter.len();

    println!("Bloom empty: {}", is_empty);
    println!("Num hashes: {}", num_hashes);
    println!("Bloom length: {}", bloom_len);

    loaded_igd_bloom_filter

}

fn file_exists(path: &str) -> bool {
    Path::new(path).exists() && Path::new(path).is_file()
}

#[cfg(test)]
mod tests{
    use super::*;
    use rstest::rstest;

    #[rstest]
    fn test_true_is_true(){

        create_bloom_filter_main();
        pretty_assertions::assert_eq!(true, true);

    }

}



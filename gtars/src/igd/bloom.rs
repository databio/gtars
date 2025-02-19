use std::collections::HashMap;
use std::env::current_exe;
use std::ffi::OsStr;
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


/// struct that contains the file path that was used to create the bloom filter as well as the parent
pub struct bloom_filter_local {
    pub parent_directory: String,
    pub bloom_file_path: String,
    pub bloom_filter: Bloom<String>
}

pub fn create_bloom_filter_main(action: &str){

    let universe_file  ="/home/drc/Downloads/bloom_testing/real_data/data/universe.merged.pruned.filtered100k.bed";

    //let query_bed = "/home/drc/Downloads/bloom_testing/test1/query1.bed";
    let query_bed = "/home/drc/Downloads/bloom_testing/test1/query2.bed";

    // Just create the universe once
    let universe_path = Path::new(&universe_file);
    let universe_tree_tokenizer = TreeTokenizer::try_from(universe_path).unwrap();

    let num_of_items = 10000;
    let false_positive_rate = 0.001;
    let seed = 42;

    // make parent directory to hold sub directories
    let save_path ="/home/drc/Downloads/bloom_testing/test1/";
    let name = "test2";
    let parent_directory = format!("{}{}/",save_path.clone(),name.clone());
    make_parent_directory(parent_directory.as_str()).unwrap();

    //let bed_directory = "/home/drc/Downloads/bloom_testing/test1/all_bed_files/";
    let bed_directory = "/home/drc/Downloads/bloom_testing/test1/all_bed_files_real/";

    match action{

        "create" => {create_bloom_filters(parent_directory.clone(),bed_directory,universe_tree_tokenizer, num_of_items,false_positive_rate);}
        "search" => {search_bloom_filter(parent_directory.as_str(), universe_file, query_bed);}


        _ => {println!("No action given, Select search or create")}
    }


}

fn create_bloom_filters(parent_directory: String, bed_directory: &str, universe_tree_tokenizer: TreeTokenizer, num_of_items: usize, false_positive_rate: f64) {
    println!("Entering create bloomfilters");
    let mut child_meta_data_map: HashMap<String, String> = HashMap::new(); // File name and directory to store files
    // For now we will store the individual bloom filters in child directories
    make_child_directories(parent_directory.clone(), bed_directory, &mut child_meta_data_map);

    //DEBUG CHECK TO ENSURE WE CHANGED HASHMAP
    for (key, value) in child_meta_data_map {  // &map to borrow, not take ownership
        println!("Key: {}, Value: {}", key, value); // key is child directory  value is the absolute path to the bed file

        tokenize_then_create_bloom(&universe_tree_tokenizer, value.as_str(), key.as_str(), num_of_items, false_positive_rate)


    }

}

fn make_child_directories(parent_directory: String, bed_directory: &str, meta_data: &mut HashMap<String, String>) {

    //let mut bed_files = Vec::new();
    for entry in fs::read_dir(bed_directory).unwrap() {
        let entry = entry.unwrap();
        let path = entry.path();

        if path.is_file() {
            println!("Found a file....");
            if let Some(extension) = path.extension() {
                if let Some(extension) = extension.to_str(){
                    println!("Found this extension: {}", extension);
                match extension{
                    "bed" | "gz" => {
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

                    _ => {}
                }}

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
    let mut chr_copy = first_region.chr.clone();
    let mut bloom_filter_path = format!("{}{}.bloom",child_directory, chr_copy);
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
        chr_copy = region.chr.clone();

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

    println!("What is this value at the end? {}", chr_copy);
    bloom_filter_path = format!("{}{}.bloom",child_directory, chr_copy);
    write_bloom_filter_to_disk(current_bloom_filter, bloom_filter_path);




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


pub fn search_bloom_filter(path_to_bloom_directory: &str, path_to_universe: &str, query_bed_file: &str){

    println!("SEARCH BLOOM FILTER");
    let bloom_files = find_bloom_files(path_to_bloom_directory).unwrap();

    // TODO ensure there ARE actually bloom files before bothering to proceed.

    // Create universe and tokenize the query using the universe
    // THIS SHOULD BE THE SAME UNIVERSE THAT THE BLOOMFILTER WAS CREATED FROM.
    let universe_path = Path::new(&path_to_universe);
    let universe_tree_tokenizer = TreeTokenizer::try_from(universe_path).unwrap();
    let bed = Path::new(&query_bed_file);
    let regions = RegionSet::try_from(bed)
        .with_context(|| "There was an error reading in the bedfile to be tokenized!").unwrap();
    let tokenized_regions = universe_tree_tokenizer.tokenize_region_set(&regions);
    let mut tokenized_regions_iter = tokenized_regions.into_iter();

    let mut hits: HashMap<String, HashMap<String, u32>> = Default::default(); // key is file path value is counts
    let mut bloom_path_maps : HashMap<Bloom<String>, String> = Default::default(); //TODO this might be a poor way to keep track of this relationship.

    let mut found = false;
    let first_region: Region  = tokenized_regions_iter.next().unwrap().into();
    let mut chr_copy = first_region.chr.clone();
    let mut start_copy = first_region.start.clone();
    let mut end_copy = first_region.end.clone();


    // skip regions if NO bloom filter exist
    while !found {
        if bloom_files.contains_key(&chr_copy) {
            found = true;
        } else {
            match tokenized_regions_iter.next() {
                Some(region_str) => {
                    let next_region: Region = region_str.into();
                    chr_copy = next_region.chr.clone();
                    start_copy = next_region.start.clone();
                    end_copy = next_region.end.clone();
                }
                None => {
                    println!("No more regions to process.");
                    break;
                }
            }
        }
    }

    // We have to do some work to check the very first region and ensure the proper filters are loaded
    // before engaging the main loop
    let mut previous_chrom = chr_copy.clone();
    println!("Found initial region with chromosome: {}", chr_copy);
    let mut paths_to_blooms =  bloom_files.get(&chr_copy).unwrap();
    let mut bloom_filter_structs: Vec<bloom_filter_local> = vec![];

    for path in paths_to_blooms.iter(){
        println!("Here are paths for loading bloom filters: {}", path);
        let path_copy = path.clone();
        let path = Path::new(path);
        let parent = path.parent().and_then(|p| p.to_str()).unwrap();
        let current_bloom_filter = load_bloom_filter_from_disk(path_copy.clone());

        let current_bloom_struct = bloom_filter_local{
            parent_directory: parent.to_string(),
            bloom_file_path: path_copy.to_string(),
            bloom_filter: current_bloom_filter,
        };

        bloom_filter_structs.push(current_bloom_struct);

    }

    // Check the very first region
    if bloom_filter_structs.len()>0{
        println!("Checking first value.....");
        // implies we found bloom filters related to the current chromosome
        let line = format!("{}|{}|{}", chr_copy.clone(), start_copy.clone(), end_copy.clone());

        for s in bloom_filter_structs.iter(){

            let result = s.bloom_filter.check(&line);

            println!("Found something: {}", result);
            if result {
                //*hits.entry(key).or_insert(0) += 1; // Increment by 1
                let outerkey = s.parent_directory.clone(); //will tell us the parent file this bloom was created form
                let innerkey = s.bloom_file_path.clone(); // will tell us the chromosome this came from
                *hits.entry(outerkey.clone()) // Clone key1 for entry()
                    .or_insert_with(HashMap::new) // If key1 doesn't exist, create it
                    .entry(innerkey) // Get the entry for key2 in the inner HashMap
                    .or_insert(0) += 1; // Increment the count (or initialize to 0 if it doesn't exist)
            }

        }

    }



    for tokenized_region in tokenized_regions_iter{

        // First check if there is even a bloom filter anywhere for this chromosome
        let region: Region = tokenized_region.into();

        chr_copy = region.chr.clone();
        start_copy = region.start.clone();
        end_copy = region.end.clone();

        if chr_copy != previous_chrom{
            // if switching chromosomes...load new bloom filters
            // best performance if query bed files is already sorted by chromosome.

            if !bloom_files.contains_key(&chr_copy){
                continue

            } else{
                previous_chrom = chr_copy.clone();
                println!("KEY EXISTS: {}", chr_copy);

                paths_to_blooms =  bloom_files.get(&chr_copy).unwrap();
                //bloom_filters.clear(); // clear to prepare pushing new bloom filters from disk into memory
                bloom_filter_structs.clear();

                for path in paths_to_blooms.iter(){
                    let path_copy = path.clone();
                    let path = Path::new(path);
                    let parent = path.parent().and_then(|p| p.to_str()).unwrap();
                    let current_bloom_filter = load_bloom_filter_from_disk(path_copy.clone());

                    let current_bloom_struct = bloom_filter_local{
                        parent_directory: parent.to_string(),
                        bloom_file_path: path_copy.to_string(),
                        bloom_filter: current_bloom_filter,
                    };

                    bloom_filter_structs.push(current_bloom_struct);


                }

                if bloom_filter_structs.len()>0{
                    // implies we found bloom filters related to the current chromosome
                    let line = format!("{}|{}|{}", chr_copy.clone(), start_copy.clone(), end_copy.clone());
                    for s in bloom_filter_structs.iter(){

                        let result = s.bloom_filter.check(&line);

                        if result {
                            let outerkey = s.parent_directory.clone(); //will tell us the parent file this bloom was created form
                            let innerkey = s.bloom_file_path.clone(); // will tell us the chromosome this came from
                            *hits.entry(outerkey.clone()) // Clone key1 for entry()
                                .or_insert_with(HashMap::new) // If key1 doesn't exist, create it
                                .entry(innerkey) // Get the entry for key2 in the inner HashMap
                                .or_insert(0) += 1; // Increment the count (or initialize to 0 if it doesn't exist)

                        }


                    }

                }

            }

        } else{

            let line = format!("{}|{}|{}", chr_copy.clone(), start_copy.clone(), end_copy.clone());
            for s in bloom_filter_structs.iter(){
                let result = s.bloom_filter.check(&line);
                if result {
                    let outerkey = s.parent_directory.clone(); //will tell us the parent file this bloom was created form
                    let innerkey = s.bloom_file_path.clone(); // will tell us the chromosome this came from
                    *hits.entry(outerkey.clone()) // Clone key1 for entry()
                        .or_insert_with(HashMap::new) // If key1 doesn't exist, create it
                        .entry(innerkey) // Get the entry for key2 in the inner HashMap
                        .or_insert(0) += 1; // Increment the count (or initialize to 0 if it doesn't exist)

                }
            }


        }

    }

    // Print the hashmap
    println!("HERE ARE THE FINAL RESULTS:");
    for (outer_key, inner_map) in &hits {
        println!("Outer key: {}", outer_key);
        for (inner_key, count) in inner_map {
            println!("  Inner key: {}, Count: {}", inner_key, count);
        }
    }

}
fn find_bloom_files(root_dir: &str) -> Result<HashMap<String, Vec<String>>, std::io::Error> {
    let mut all_files = HashMap::new();

    let root_path = Path::new(root_dir);
    if !root_path.is_dir() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::NotADirectory,
            "Root directory must be a directory",
        ));
    }

    fn traverse_directory(
        dir: &Path,
        all_files: &mut HashMap<String, Vec<String>>,
    ) -> Result<(), std::io::Error> {
        for entry in fs::read_dir(dir)? {
            let entry = entry?;
            let path = entry.path();

            if path.is_dir() {
                traverse_directory(&path, all_files)?;
            } else if path.is_file() {
                let file_name = entry.file_name().to_string_lossy();
                let base_name = match path.file_stem() {
                    Some(stem) => stem.to_string_lossy().to_string(),
                    None => continue,
                };
                let absolute_path = match path.canonicalize() {
                    Ok(abs_path) => abs_path.to_string_lossy().to_string(),
                    Err(e) => {
                        eprintln!("Error getting absolute path for {:?}: {}", path, e);
                        continue;
                    }
                };

                all_files.entry(base_name) // Get or create the entry for this base name
                    .or_insert_with(Vec::new)
                    .push(absolute_path);
            }
        }
        Ok(())
    }

    traverse_directory(root_path, &mut all_files)?;

    Ok(all_files)
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

        //create_bloom_filter_main("create");
        create_bloom_filter_main("search");
        pretty_assertions::assert_eq!(true, true);

    }

}



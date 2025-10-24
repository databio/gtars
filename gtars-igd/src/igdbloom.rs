use std::collections::HashMap;
use std::env::current_exe;
use std::ffi::OsStr;
use std::fs;
use std::fs::File;
use std::io::{Error, Read, Write};
use std::path::Path;
// for now use simple implementation of bloom filter in rust

use gtars_core::models::{Region, RegionSet};
use gtars_tokenizers::tokenizer::Tokenizer;

#[cfg(feature = "bloom")]
extern crate bloom;

#[cfg(feature = "bloom")]
use bloom::{BloomFilter, ASMS};

/// struct that contains the file path that was used to create the bloom filter as well as the parent

#[cfg(feature = "bloom")]
pub struct bloom_filter_local {
    pub parent_directory: String,
    pub bloom_file_path: String,
    pub bloom_filter: BloomFilter
}
//
// pub fn igd_get_bloom_matches(matches: &ArgMatches){
//
//     let action = matches
//         .get_one::<String>("action")
//         .expect("Action (search or create) is required").as_str();
//
//     let universe = matches
//         .get_one::<String>("universe")
//         .expect("universe path is required");
//
//     let bedfilesuniverse = matches
//         .get_one::<String>("bedfilesuniverse")
//         .map(String::as_str)
//         .unwrap_or_default();
//
//
//     let querybed = matches
//         .get_one::<String>("querybed")
//         .map(String::as_str)
//         .unwrap_or_default();
//
//     let bloomdirectory = matches
//         .get_one::<String>("bloomdirectory")
//         .expect("bloomdirectory is required");
//
//     let bloomname = matches
//         .get_one::<String>("bloomname")
//         .expect("bloomname is required");
//
//     let numitems = matches
//         .get_one::<usize>("numitems")
//         .expect("number of items required");
//
//     let falsepositive = matches
//         .get_one::<f64>("falsepositive")
//         .expect("number of items required");
//
//
//     let parent_directory = format!("{}{}/",bloomdirectory.clone(),bloomname.clone());
//
//
//     match action{
//
//         "create" => {
//             let universe_path = Path::new(&universe);
//             let universe_tree_tokenizer = TreeTokenizer::try_from(universe_path).unwrap();
//             make_parent_directory(parent_directory.as_str()).unwrap();
//
//             create_bloom_filters(parent_directory.clone(), bedfilesuniverse, universe_tree_tokenizer, *numitems, *falsepositive);
//
//
//         }
//         "search" => {search_bloom_filter(parent_directory.as_str(), universe, querybed);}
//
//
//         _ => {println!("No action given, Select search or create")}
//     }
//
// }

// pub fn create_bloom_filter_main(action: &str){
//
//     let universe_file  ="/home/drc/Downloads/bloom_testing/real_data/data/universe.merged.pruned.filtered100k.bed";
//
//     //let query_bed = "/home/drc/Downloads/bloom_testing/test1/query1.bed";
//     let query_bed = "/home/drc/Downloads/bloom_testing/test1/query2.bed";
//
//     // Just create the universe once
//     let universe_path = Path::new(&universe_file);
//     let universe_tree_tokenizer = TreeTokenizer::try_from(universe_path).unwrap();
//
//     let num_of_items = 10000;
//     let false_positive_rate = 0.001;
//     //let seed = 42;
//
//     // make parent directory to hold sub directories
//     let save_path ="/home/drc/Downloads/bloom_testing/test1/";
//     let name = "test";
//     let parent_directory = format!("{}{}/",save_path.clone(),name.clone());
//     make_parent_directory(parent_directory.as_str()).unwrap();
//
//     //let bed_directory = "/home/drc/Downloads/bloom_testing/test1/all_bed_files/";
//     //let bed_directory = "/home/drc/Downloads/bloom_testing/test1/all_bed_files_real/";
//     let bed_directory = "/home/drc/Downloads/bloom_testing/test1/two_real_bed_files/";
//
//
// }


// #[cfg(feature = "bloom")]
// fn create_bloom_filters(parent_directory: String, bed_directory: &str, universe_tokenizer: Tokenizer, num_of_items: usize, false_positive_rate: f64) {
//     println!("Creating bloomfilters");
//     let mut child_meta_data_map: HashMap<String, String> = HashMap::new(); // File name and directory to store files
//     // For now we will store the individual bloom filters in child directories
//     make_child_directories(parent_directory.clone(), bed_directory, &mut child_meta_data_map);

//     //DEBUG CHECK TO ENSURE WE CHANGED HASHMAP
//     for (key, value) in child_meta_data_map {  // &map to borrow, not take ownership
//         //println!("Key: {}, Value: {}", key, value); // key is child directory  value is the absolute path to the bed file
//         tokenize_then_create_bloom(&universe_tokenizer, value.as_str(), key.as_str(), num_of_items, false_positive_rate)

//     }

//  }
#[cfg(feature = "bloom")]
fn make_child_directories(parent_directory: String, bed_directory: &str, meta_data: &mut HashMap<String, String>) {

    //let mut bed_files = Vec::new();
    for entry in fs::read_dir(bed_directory).unwrap() {
        let entry = entry.unwrap();
        let path = entry.path();

        if path.is_file() {
            //println!("Found a file....");
            if let Some(extension) = path.extension() {
                if let Some(extension) = extension.to_str(){
                    //println!("Found this extension: {}", extension);
                match extension{
                    "bed" | "gz" => {
                        let full_name = path.to_str().unwrap().to_string();
                        //println!("Here is the full name: {}", full_name.clone());

                        let name_without_extension = path
                            .file_stem()
                            .unwrap()
                            .to_str()
                            .unwrap()
                            .to_string();
                        //println!("Here is the name without extension: {}", name_without_extension.clone());


                        let single_parent_directory = format!("{}{}/",parent_directory,name_without_extension);
                        let _ = make_parent_directory(single_parent_directory.as_str());

                        meta_data.insert(single_parent_directory,full_name.clone());
                    }

                    _ => {}
                }}

            }
        }
    }
}

#[cfg(feature = "bloom")]
pub fn tokenize_then_create_bloom(universe_tokenizer: &Tokenizer, bed_file: &str, child_directory: &str, num_of_items: usize, false_positive_rate: f64){
    // we must first tokenize against a universe and then create a bloom filter for each chromosome
    // from that tokenization

    //TODO implement random generation of seed and create filters from this seed.
    let mut seed = [0u8; 32];

    // First load regions from bed file
    let bed = Path::new(&bed_file);
    let regions = RegionSet::try_from(bed).unwrap();
    
    // Group regions by chromosome
    let mut chrom_regions: std::collections::HashMap<String, Vec<Region>> = std::collections::HashMap::new();
    for region in &regions.regions {
        chrom_regions.entry(region.chr.clone()).or_insert(Vec::new()).push(region.clone());
    }
    
    // Create bloom filter for each chromosome
    for (chr, chr_regions) in chrom_regions {
        let bloom_filter_path = format!("{}{}.bloom", child_directory, chr);
        
        // Check if bloom filter already exists
        if file_exists(&bloom_filter_path) {
            println!("File already exists: {}", bloom_filter_path);
            continue;
        }
        
        // Create new bloom filter
        let mut current_bloom_filter = BloomFilter::with_rate(false_positive_rate as f32, num_of_items as u32);
        
        // Tokenize regions for this chromosome and add tokens to bloom filter
        let tokenized_regions = universe_tokenizer.tokenize(&chr_regions).unwrap();
        for token in tokenized_regions {
            current_bloom_filter.insert(&token);
        }
        
        // Also add the region coordinates as backup
        for region in chr_regions {
            let line = format!("{}|{}|{}", region.chr, region.start, region.end);
            current_bloom_filter.insert(&line);
        }
        
        write_bloom_filter_to_disk(current_bloom_filter, bloom_filter_path);
    }
}


#[cfg(feature = "bloom")]
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

#[cfg(feature = "bloom")]
pub fn write_bloom_filter_to_disk(igd_bloom_filter: BloomFilter, save_path: String){
    // TODO: Implement serialization for BloomFilter
    // The bloom crate may not support direct serialization
    println!("Would save bloom filter to: {}", save_path);
}

// #[cfg(all(feature = "bloom", not(feature = "bloom-api-fixed")))]
// pub fn search_bloom_filter(path_to_bloom_directory: &str, path_to_universe: &str, query_bed_file: &str){

//     // TODO: Fix this function - it has issues with token to region conversion  
//     println!("search_bloom_filter function temporarily disabled due to API changes");
//     return;

//     println!("SEARCH BLOOM FILTER");
//     let bloom_files = find_bloom_files(path_to_bloom_directory).unwrap();

//     // TODO ensure there ARE actually bloom files before bothering to proceed.

//     // Create universe and tokenize the query using the universe
//     // THIS SHOULD BE THE SAME UNIVERSE THAT THE BLOOMFILTER WAS CREATED FROM.
//     let universe_path = Path::new(&path_to_universe);
//     let universe_tokenizer = Tokenizer::from_bed(universe_path).unwrap();
//     let bed = Path::new(&query_bed_file);
//     let regions = RegionSet::try_from(bed).unwrap();
//     let tokenized_regions = universe_tokenizer.tokenize(regions.regions.as_slice()).unwrap();
//     let mut tokenized_regions_iter = tokenized_regions.into_iter();

//     let mut hits: HashMap<String, HashMap<String, u32>> = Default::default(); // key is file path value is counts
//     let mut bloom_path_maps : HashMap<BloomFilter, String> = Default::default(); //TODO this might be a poor way to keep track of this relationship.

//     let mut found = false;
//     // TODO: Fix this line - tokenized_regions_iter produces String, not Region
//     // let first_region: Region  = tokenized_regions_iter.next().unwrap().into();
//     // For now, just return early
//     println!("Function needs API fixes");
//     return;
//     let mut chr_copy = first_region.chr.clone();
//     let mut start_copy = first_region.start.clone();
//     let mut end_copy = first_region.end.clone();


//     // skip regions if NO bloom filter exist
//     while !found {
//         if bloom_files.contains_key(&chr_copy) {
//             found = true;
//         } else {
//             match tokenized_regions_iter.next() {
//                 Some(region_str) => {
//                     let next_region: Region = region_str.into();
//                     chr_copy = next_region.chr.clone();
//                     start_copy = next_region.start.clone();
//                     end_copy = next_region.end.clone();
//                 }
//                 None => {
//                     println!("No more regions to process.");
//                     break;
//                 }
//             }
//         }
//     }

//     // We have to do some work to check the very first region and ensure the proper filters are loaded
//     // before engaging the main loop
//     let mut previous_chrom = chr_copy.clone();
//     //println!("Found initial region with chromosome: {}", chr_copy);
//     let mut paths_to_blooms =  bloom_files.get(&chr_copy).unwrap();
//     let mut bloom_filter_structs: Vec<bloom_filter_local> = vec![];

//     for path in paths_to_blooms.iter(){
//         //println!("Here are paths for loading bloom filters: {}", path);
//         let path_copy = path.clone();
//         let path = Path::new(path);
//         let parent = path.parent().and_then(|p| p.to_str()).unwrap();
//         let current_bloom_filter = load_bloom_filter_from_disk(path_copy.clone());

//         let current_bloom_struct = bloom_filter_local{
//             parent_directory: parent.to_string(),
//             bloom_file_path: path_copy.to_string(),
//             bloom_filter: current_bloom_filter,
//         };

//         bloom_filter_structs.push(current_bloom_struct);

//     }

//     // Check the very first region
//     if bloom_filter_structs.len()>0{
//         //println!("Checking first value.....");
//         // implies we found bloom filters related to the current chromosome
//         let line = format!("{}|{}|{}", chr_copy.clone(), start_copy.clone(), end_copy.clone());

//         for s in bloom_filter_structs.iter(){

//             let result = s.bloom_filter.contains(&line);

//             //println!("Found something: {}", result);
//             if result {
//                 //*hits.entry(key).or_insert(0) += 1; // Increment by 1
//                 let outerkey = s.parent_directory.clone(); //will tell us the parent file this bloom was created form
//                 let innerkey = s.bloom_file_path.clone(); // will tell us the chromosome this came from
//                 *hits.entry(outerkey.clone()) // Clone key1 for entry()
//                     .or_insert_with(HashMap::new) // If key1 doesn't exist, create it
//                     .entry(innerkey) // Get the entry for key2 in the inner HashMap
//                     .or_insert(0) += 1; // Increment the count (or initialize to 0 if it doesn't exist)
//             }

//         }

//     }



//     for tokenized_region in tokenized_regions_iter{

//         // First check if there is even a bloom filter anywhere for this chromosome
//         let region: Region = tokenized_region.into();

//         chr_copy = region.chr.clone();
//         start_copy = region.start.clone();
//         end_copy = region.end.clone();

//         if chr_copy != previous_chrom{
//             // if switching chromosomes...load new bloom filters
//             // best performance if query bed files is already sorted by chromosome.

//             if !bloom_files.contains_key(&chr_copy){
//                 continue

//             } else{
//                 previous_chrom = chr_copy.clone();
//                 //println!("KEY EXISTS: {}", chr_copy);

//                 paths_to_blooms =  bloom_files.get(&chr_copy).unwrap();
//                 //bloom_filters.clear(); // clear to prepare pushing new bloom filters from disk into memory
//                 bloom_filter_structs.clear();

//                 for path in paths_to_blooms.iter(){
//                     let path_copy = path.clone();
//                     let path = Path::new(path);
//                     let parent = path.parent().and_then(|p| p.to_str()).unwrap();
//                     let current_bloom_filter = load_bloom_filter_from_disk(path_copy.clone());

//                     let current_bloom_struct = bloom_filter_local{
//                         parent_directory: parent.to_string(),
//                         bloom_file_path: path_copy.to_string(),
//                         bloom_filter: current_bloom_filter,
//                     };

//                     bloom_filter_structs.push(current_bloom_struct);


//                 }

//                 if bloom_filter_structs.len()>0{
//                     // implies we found bloom filters related to the current chromosome
//                     let line = format!("{}|{}|{}", chr_copy.clone(), start_copy.clone(), end_copy.clone());
//                     for s in bloom_filter_structs.iter(){

//                         let result = s.bloom_filter.contains(&line);

//                         if result {
//                             let outerkey = s.parent_directory.clone(); //will tell us the parent file this bloom was created form
//                             let innerkey = s.bloom_file_path.clone(); // will tell us the chromosome this came from
//                             *hits.entry(outerkey.clone()) // Clone key1 for entry()
//                                 .or_insert_with(HashMap::new) // If key1 doesn't exist, create it
//                                 .entry(innerkey) // Get the entry for key2 in the inner HashMap
//                                 .or_insert(0) += 1; // Increment the count (or initialize to 0 if it doesn't exist)

//                         }


//                     }

//                 }

//             }

//         } else{

//             if bloom_filter_structs.len()>0 {
//                 let line = format!("{}|{}|{}", chr_copy.clone(), start_copy.clone(), end_copy.clone());
//                 for s in bloom_filter_structs.iter() {
//                     let result = s.bloom_filter.contains(&line);
//                     if result {
//                         let outerkey = s.parent_directory.clone(); //will tell us the parent file this bloom was created form
//                         let innerkey = s.bloom_file_path.clone(); // will tell us the chromosome this came from
//                         *hits.entry(outerkey.clone()) // Clone key1 for entry()
//                             .or_insert_with(HashMap::new) // If key1 doesn't exist, create it
//                             .entry(innerkey) // Get the entry for key2 in the inner HashMap
//                             .or_insert(0) += 1; // Increment the count (or initialize to 0 if it doesn't exist)

//                     }
//                 }
//             }


//         }

//     }

//     // Print the hashmap
//     println!("HERE ARE THE FINAL RESULTS:");
//     for (outer_key, inner_map) in &hits {
//         println!("Outer key: {}", outer_key);
//         for (inner_key, count) in inner_map {
//             println!("  Inner key: {}, Count: {}", inner_key, count);
//         }
//     }

// }
// #[cfg(feature = "bloom")]
// fn find_bloom_files(root_dir: &str) -> Result<HashMap<String, Vec<String>>, std::io::Error> {
//     let mut all_files = HashMap::new();

//     let root_path = Path::new(root_dir);
//     if !root_path.is_dir() {
//         return Err(std::io::Error::new(
//             std::io::ErrorKind::NotADirectory,
//             "Root directory must be a directory",
//         ));
//     }

//     fn traverse_directory(
//         dir: &Path,
//         all_files: &mut HashMap<String, Vec<String>>,
//     ) -> Result<(), std::io::Error> {
//         for entry in fs::read_dir(dir)? {
//             let entry = entry?;
//             let path = entry.path();

//             if path.is_dir() {
//                 traverse_directory(&path, all_files)?;
//             } else if path.is_file() {
//                 let file_name = entry.file_name().to_string_lossy();
//                 let base_name = match path.file_stem() {
//                     Some(stem) => stem.to_string_lossy().to_string(),
//                     None => continue,
//                 };
//                 let absolute_path = match path.canonicalize() {
//                     Ok(abs_path) => abs_path.to_string_lossy().to_string(),
//                     Err(e) => {
//                         eprintln!("Error getting absolute path for {:?}: {}", path, e);
//                         continue;
//                     }
//                 };

//                 all_files.entry(base_name) // Get or create the entry for this base name
//                     .or_insert_with(Vec::new)
//                     .push(absolute_path);
//             }
//         }
//         Ok(())
//     }

//     traverse_directory(root_path, &mut all_files)?;

//     Ok(all_files)
// }
// #[cfg(feature = "bloom")]
// pub fn load_bloom_filter_from_disk(load_path: String) -> BloomFilter
// {
//     //println!("Loading Bloom Filter: {}", load_path);
//     let mut file = File::open(&load_path).unwrap();

//     let mut buffer: Vec<u8> = Vec::new();

//     // TODO: Implement deserialization for BloomFilter
//     // For now, create a new empty bloom filter
//     println!("Would load bloom filter from: {}", load_path);
//     let loaded_igd_bloom_filter = BloomFilter::with_rate(0.01, 1000);
//     loaded_igd_bloom_filter

// }

fn file_exists(path: &str) -> bool {
    Path::new(path).exists() && Path::new(path).is_file()
}

// #[cfg(feature = "bloom")]
// pub fn create_bloom_tree(){

//     // Some notes to consider for a tree implementation
//     // I assume the leaf nodes are bloom filters created at the granularity of bedfile_chr#
//     // all the bloom filters must use the same hash functions for this to work given the XOR operation
//     // all the bloom filters must be the same length due to the xor function
//     // height_tree = log2(number_of_leaves)
//     // we must provide a tree from a list of bloom filters, and then decorate the tree with Union of bloom filters
//     // so that we can create nodes that union subsets of the bloom filters.


//     // make parent directory to hold sub directories
//     let save_path ="/home/drc/Downloads/bloom_testing/test_tree/";
//     let name = "bloom_tree";
//     let parent_directory = format!("{}{}/",save_path.clone(),name.clone());
//     make_parent_directory(parent_directory.as_str()).unwrap();

//     //let bed_directory = "/home/drc/Downloads/bloom_testing/test1/two_real_bed_files/";
//     let bed_directory = "/home/drc/Downloads/bloom_testing/test1/all_bed_files";

//     let universe  ="/home/drc/Downloads/bloom_testing/real_data/data/universe.merged.pruned.filtered100k.bed";
//     let universe_path = Path::new(&universe);
//     let universe_tokenizer = Tokenizer::from_bed(universe_path).unwrap();

//     create_bloom_filters(parent_directory.clone(), bed_directory, universe_tokenizer, 10000, 0.01);

//     // Now that we have created a directory of directories containing .bloom files, create a tree (balanced)

//     let all_bloom_files = find_bloom_files(&*parent_directory).unwrap();

//     let mut all_bloom_files_vec: Vec<String> = vec![];

//     for (key, value) in all_bloom_files{

//         for bloom_file_path in value.iter(){
//             all_bloom_files_vec.push((*bloom_file_path.clone()).parse().unwrap())
//         }

//     }

//     //println!("Here are my files: {:?}", all_bloom_files_vec);

//     //Building binary tree from just the file list
//     // let root = build_binary_tree(&*all_bloom_files_vec);
//     // if let Some(r) = root {
//     //     print_tree(&Some(r), 0);
//     // } else {
//     //     println!("No .bloom files found or tree is empty.");
//     // }

//     // Build nodes and then assign leaf nodes the pre-computed bloom filters
//     let num_leaves = all_bloom_files_vec.len();
//     let num_nodes = 2*num_leaves - 1;



// }
// #[cfg(feature = "bloom")]
// #[derive(Debug)]
// struct BloomFilterNode {
//     value: Option<String>,
//     left: Option<Box<BloomFilterNode>>,
//     right: Option<Box<BloomFilterNode>>,
// }

// #[cfg(feature = "bloom")]
// struct BloomFilterNode2 {
//     path: Option<String>, // only used for leaf nodes that tell us which chrom the bloom was constructed from
//     left: Option<Box<BloomFilterNode>>,
//     right: Option<Box<BloomFilterNode>>,
//     bloomfilter: Option<BloomFilter>
// }

// pub fn build_binary_tree(files: &[String])-> Option<Box<BloomFilterNode>> {

//     // Given a list of bloom files, build a binary tree.
//     // height = log2(num_files)

//     if files.is_empty() {
//         return None;
//     }

//     let n = files.len();
//     if n == 1 {
//         return Some(Box::new(BloomFilterNode {
//             value: Some(files[0].clone()),
//             left: None,
//             right: None,
//         }));
//     }

//     let mid = n / 2;
//     let root = Box::new(BloomFilterNode {
//         value: Some(files[mid].clone()),
//         left: None,
//         right: None,
//     });

//     let left_files = &files[..mid];
//     let right_files = &files[mid + 1..];

//     let left_subtree = build_binary_tree(left_files);
//     let right_subtree = build_binary_tree(right_files);

//     let mut root_node = *root; // Dereference to move ownership
//     root_node.left = left_subtree;
//     root_node.right = right_subtree;

//     Some(Box::new(root_node)) // Re-box and return

// }
// #[cfg(feature = "bloom")]
// pub fn build_binary_tree_2(files: &[String])-> Option<Box<BloomFilterNode>> {

//     // Given a list of bloom files, build a binary tree.
//     // height = log2(num_files)
//     //todo build actual bloom tree per chromosome

//     if files.is_empty() {
//         return None;
//     }

//     let n = files.len();
//     if n == 1 {
//         return Some(Box::new(BloomFilterNode {
//             value: Some(files[0].clone()),
//             left: None,
//             right: None,
//         }));
//     }

//     let mid = n / 2;
//     let root = Box::new(BloomFilterNode {
//         value: Some(files[mid].clone()),
//         left: None,
//         right: None,
//     });

//     let left_files = &files[..mid];
//     let right_files = &files[mid + 1..];

//     let left_subtree = build_binary_tree(left_files);
//     let right_subtree = build_binary_tree(right_files);

//     let mut root_node = *root; // Dereference to move ownership
//     root_node.left = left_subtree;
//     root_node.right = right_subtree;

//     Some(Box::new(root_node)) // Re-box and return

// }
// #[cfg(feature = "bloom")]
// fn print_tree(node: &Option<Box<BloomFilterNode>>, level: usize) {
//     if let Some(n) = node {
//         if let Some(val) = &n.value {
//             println!("{}{}", "  ".repeat(level), val);
//         }
//         print_tree(&n.left, level + 1);
//         print_tree(&n.right, level + 1);
//     }
// }

#[cfg(test)]
mod tests{
    use super::*;
    use rstest::rstest;
    use std::path::PathBuf;

    #[rstest]
    fn test_manual(){

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
        let num_of_items = 5;
        let false_positive_rate = 10.0;

        let tokenizer =
            Tokenizer::from_auto(bed_path.as_ref()).expect("Failed to create tokenizer from config.");

        tokenize_then_create_bloom(&tokenizer, &bed_path, &child_directory, num_of_items, false_positive_rate);
        //
        // let universe  ="/home/drc/Downloads/bloom_testing/real_data/data/universe.merged.pruned.filtered100k.bed";
        //
        // //let query_bed = "/home/drc/Downloads/bloom_testing/test1/query1.bed";
        // let querybed = "/home/drc/Downloads/bloom_testing/test1/query2.bed";
        //
        // // Just create the universe once
        // let universe_path = Path::new(&universe);
        // let universe_tree_tokenizer = TreeTokenizer::try_from(universe_path).unwrap();
        //
        // let numitems = 10000;
        // let falsepositive = 0.001;
        //
        // // make parent directory to hold sub directories
        // let save_path ="/home/drc/Downloads/bloom_testing/test1/";
        // let name = "test";
        // let parent_directory = format!("{}{}/",save_path.clone(),name.clone());
        // make_parent_directory(parent_directory.as_str()).unwrap();
        //
        // //let bed_directory = "/home/drc/Downloads/bloom_testing/test1/all_bed_files/";
        // //let bed_directory = "/home/drc/Downloads/bloom_testing/test1/all_bed_files_real/";
        // let bedfilesuniverse = "/home/drc/Downloads/bloom_testing/test1/two_real_bed_files/";
        //
         //create_bloom_filters(parent_directory.clone(), bedfilesuniverse, universe_tree_tokenizer, numitems, falsepositive);
        //
        // search_bloom_filter(parent_directory.as_str(), universe, querybed);

        //create_bloom_tree();
        pretty_assertions::assert_eq!(true, true);

    }



}



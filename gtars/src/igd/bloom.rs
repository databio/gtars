use std::fs::File;
use std::io::{Read, Write};
// for now use simple implementation of bloom filter in rust
use bloomfilter;
use bloomfilter::Bloom;

pub fn create_bloom_filter_main(){

    let num_of_items = 10000;
    let false_positive_rate = 0.001;
    let seed = 42;

    let save_path ="/home/drc/Downloads/bloom_testing/test1/";
    let name = "test";

    let universe_file = "";

    let query_bed_file = "";

    let example_region = "chr1|500|800";

    // Create Empty Bloom Filter
    println!("Creating Bloom Filter");
    //let igd_bloom_filter = bloomfilter::Bloom::new_for_fp_rate_with_seed(num_of_items, false_positive_rate, &[]);
    let mut igd_bloom_filter = bloomfilter::Bloom::new_for_fp_rate(num_of_items, false_positive_rate).unwrap();
    igd_bloom_filter.set(&example_region);   // insert 10 in the bloom filter
    // igd_bloom_filter.check(&10); // return true
    // igd_bloom_filter.check(&20); // return false

    let is_empty =  igd_bloom_filter.is_empty();
    let num_hashes =  igd_bloom_filter.number_of_hash_functions();
    let bloom_len =  igd_bloom_filter.len();


    println!("Bloom empty: {}", is_empty);
    println!("Num hashes: {}", num_hashes);
    println!("Bloom length: {}", bloom_len);

    let final_save_path = format!("{}{}.bloom",save_path.clone(),name.clone());

    write_bloom_filter_to_disk(igd_bloom_filter, final_save_path.clone());
    load_bloom_filter_from_disk(final_save_path.clone());

}

pub fn write_bloom_filter_to_disk(igd_bloom_filter: Bloom<&str>, save_path: String){

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

pub fn load_bloom_filter_from_disk(load_path: String)
{
    println!("Loading Bloom Filter");
    let mut file = File::open(load_path).unwrap();

    let mut buffer = Vec::new();

    file.read_to_end(&mut buffer).unwrap();
    if buffer.is_empty() {
       eprintln!("FILE IS EMPTY");
    }

    let loaded_igd_bloom_filter: Bloom<&str> = Bloom::from_slice(&*buffer).unwrap();

    let is_empty = loaded_igd_bloom_filter.is_empty();
    let num_hashes = loaded_igd_bloom_filter.number_of_hash_functions();
    let bloom_len = loaded_igd_bloom_filter.len();


    println!("Bloom empty: {}", is_empty);
    println!("Num hashes: {}", num_hashes);
    println!("Bloom length: {}", bloom_len);


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



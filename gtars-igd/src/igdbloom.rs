use std::collections::HashMap;
use std::fs;
use std::io::{Error, Read, Write};
use std::path::Path;
use gtars_core::models::{Region, RegionSet};
use gtars_tokenizers::tokenizer::Tokenizer;


#[cfg(feature = "bloom")]
use bloomfilter::Bloom;

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

    // Check if bloom filter already exists
    if file_exists(&bloom_filter_path) {
        println!("File already exists: {}", bloom_filter_path);

    } else{

        // Create new bloom filter
        let mut current_bloom_filter = Bloom::new_for_fp_rate( num_of_items as usize, false_positive_rate as f64).unwrap();

        // Tokenize regions for this chromosome and add tokens to bloom filter
        let tokenized_regions = universe_tokenizer.tokenize(&regions.regions).unwrap();
        for token in tokenized_regions {
            current_bloom_filter.set(&token);
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
pub fn write_bloom_filter_to_disk(igd_bloom_filter: Bloom<String>, save_path: String) -> Result<(), std::io::Error>{

    let bytes = igd_bloom_filter.to_bytes();

    match std::fs::write(&save_path, bytes) {
        Ok(_) => {
            println!("Successfully saved bloom filter to: {}", &save_path);
            Ok(())
        }
        Err(e) => {
            eprintln!("Failed to save bloom filter: {}", e);
            Err(e)
        }
    }


}

#[cfg(feature = "bloom")]
pub fn load_bloom_filter_from_disk(load_path: &str) -> Result<Bloom<String>, Box<dyn std::error::Error>> {
    let bytes = std::fs::read(load_path)?;

    let filter = Bloom::from_bytes(bytes).map_err(|e| format!("Bloom filter error: {}", e))?;

    println!("Successfully loaded bloom filter from: {}", load_path);
    Ok(filter)
}


fn file_exists(path: &str) -> bool {
    Path::new(path).exists() && Path::new(path).is_file()
}



#[cfg(test)]
mod tests{
    use super::*;
    use rstest::rstest;
    use std::path::PathBuf;

    #[rstest]
    fn test_bloom_filter(){

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

        // Can we load the bloom filter?

        let path = Path::new(bed_path.as_ref());
        let filename = path.file_name()
            .and_then(|os_str| os_str.to_str()).unwrap();
        let bloom_filter_path = format!("{}/{}.bloom", child_directory, filename);

        let loaded_filter = load_bloom_filter_from_disk(bloom_filter_path.as_str()).unwrap();

        let result = loaded_filter.check(&"chr1:22-30".to_string());
        pretty_assertions::assert_eq!(true, result);

        let result = loaded_filter.check(&"chr1:23-31".to_string());
        pretty_assertions::assert_eq!(false, result);

    }



}



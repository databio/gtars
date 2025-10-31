use gtars_core::models::RegionSet;
use gtars_tokenizers::tokenizer::Tokenizer;
use std::collections::HashMap;
use std::fs;
use std::io::{Error, Read, Write};
use std::path::Path;

#[cfg(feature = "bloom")]
use bloomfilter::Bloom;

#[cfg(feature = "bloom")]
pub fn tokenize_then_create_bloom_for_each_file(
    universe_tokenizer: &Tokenizer,
    bed_file: &str,
    child_directory: &str,
    num_of_items: usize,
    false_positive_rate: f64,
) {
    //TODO this function does a couple of things, potentially better to disentangle them in the future

    //TODO implement random generation of seed and create filters from this seed.
    //let mut seed = [0u8; 32];

    // First load regions from bed file
    let bed = Path::new(&bed_file);
    let regions = RegionSet::try_from(bed).unwrap();

    let path = Path::new(bed_file);

    let filename = path.file_name().and_then(|os_str| os_str.to_str()).unwrap();

    // Create bloom filter for the current BED file

    let bloom_filter_path = format!("{}/{}.bloom", child_directory, filename);

    if file_exists(&bloom_filter_path) {
        println!("File already exists: {}", bloom_filter_path);
    } else {
        let mut current_bloom_filter =
            Bloom::new_for_fp_rate(num_of_items as usize, false_positive_rate as f64).unwrap();

        // Tokenize regions for this chromosome and add these regions to bloom filter
        let tokenized_regions = universe_tokenizer.tokenize(&regions.regions).unwrap();
        for token in tokenized_regions {
            current_bloom_filter.set(&token);
        }

        let _ = write_bloom_filter_to_disk(current_bloom_filter, bloom_filter_path);
    }
}

#[cfg(feature = "bloom")]
pub fn make_parent_directory(parent_directory: &str) -> Result<(), Error> {
    let parent_path = Path::new(&parent_directory);

    if !parent_path.exists() {
        match fs::create_dir_all(&parent_directory) {
            Ok(_) => {
                println!(
                    "Parent directory created successfully: {}",
                    parent_directory
                );
                Ok(())
            }
            Err(e) => {
                eprintln!("Error creating parent directory: {}", e);
                Err(e)
            }
        }
    } else {
        println!("Parent directory already exists: {}", parent_directory);
        Ok(())
    }
}

#[cfg(feature = "bloom")]
pub fn write_bloom_filter_to_disk(
    igd_bloom_filter: Bloom<String>,
    save_path: String,
) -> Result<(), std::io::Error> {
    let bytes = igd_bloom_filter.to_bytes();

    match fs::write(&save_path, bytes) {
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
pub fn load_bloom_filter_from_disk(
    load_path: &str,
) -> Result<Bloom<String>, Box<dyn std::error::Error>> {
    let bytes = fs::read(load_path)?;

    let filter = Bloom::from_bytes(bytes).map_err(|e| format!("Bloom filter error: {}", e))?;

    println!("Successfully loaded bloom filter from: {}", load_path);
    Ok(filter)
}

fn file_exists(path: &str) -> bool {
    Path::new(path).exists() && Path::new(path).is_file()
}

#[cfg(feature = "bloom")]
pub fn process_bed_directory(
    universe_tokenizer: &Tokenizer,
    input_directory: &str,
    output_directory: &str,
    num_of_items: usize,
    false_positive_rate: f64,
) -> Result<usize, Box<dyn std::error::Error>> {
    // Create output directory if it doesn't exist
    make_parent_directory(output_directory)?;

    let input_path = Path::new(input_directory);

    if !input_path.exists() || !input_path.is_dir() {
        return Err(format!(
            "Input directory does not exist or is not a directory: {}",
            input_directory
        )
        .into());
    }

    let mut processed_count = 0;

    for entry in fs::read_dir(input_path)? {
        let entry = entry?;
        let path = entry.path();

        if path.is_file() {
            if let Some(extension) = path.extension() {
                if extension == "bed" {
                    let bed_file = path.to_str().ok_or("Invalid path")?;
                    println!("Processing BED file: {}", bed_file);

                    tokenize_then_create_bloom_for_each_file(
                        universe_tokenizer,
                        bed_file,
                        output_directory,
                        num_of_items,
                        false_positive_rate,
                    );

                    processed_count += 1;
                }
            }
        }
    }

    println!(
        "Processed {} BED files from directory: {}",
        processed_count, input_directory
    );
    Ok(processed_count)
}

#[cfg(feature = "bloom")]
pub fn load_bloom_directory(
    bloom_directory: &str,
) -> Result<HashMap<String, Bloom<String>>, Box<dyn std::error::Error>> {
    let bloom_path = Path::new(bloom_directory);

    if !bloom_path.exists() || !bloom_path.is_dir() {
        return Err(format!(
            "Bloom directory does not exist or is not a directory: {}",
            bloom_directory
        )
        .into());
    }

    let mut bloom_filters: HashMap<String, Bloom<String>> = HashMap::new();

    for entry in fs::read_dir(bloom_path)? {
        let entry = entry?;
        let path = entry.path();

        if path.is_file() {
            if let Some(extension) = path.extension() {
                if extension == "bloom" {
                    let bloom_file = path.to_str().ok_or("Invalid path")?;

                    let filename = path
                        .file_stem()
                        .and_then(|os_str| os_str.to_str())
                        .ok_or("Invalid filename")?
                        .to_string();

                    println!("Loading bloom filter: {}", bloom_file);

                    match load_bloom_filter_from_disk(bloom_file) {
                        Ok(filter) => {
                            bloom_filters.insert(filename, filter);
                        }
                        Err(e) => {
                            eprintln!("Failed to load bloom filter from {}: {}", bloom_file, e);
                            // Continue loading other files even if one fails
                        }
                    }
                }
            }
        }
    }

    println!(
        "Loaded {} bloom filters from directory: {}",
        bloom_filters.len(),
        bloom_directory
    );
    Ok(bloom_filters)
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

        let child_directory = path.into_os_string().into_string().unwrap();
        let num_of_items = 1000;
        let false_positive_rate = 0.5;

        let tokenizer = Tokenizer::from_auto(bed_path.as_ref())
            .expect("Failed to create tokenizer from config.");

        // Can we create the bloom filter and save to disk?
        tokenize_then_create_bloom_for_each_file(
            &tokenizer,
            &bed_path,
            &child_directory,
            num_of_items,
            false_positive_rate,
        );

        // Can we load the saved bloom filter and query?
        let path = Path::new(bed_path.as_ref());
        let filename = path.file_name().and_then(|os_str| os_str.to_str()).unwrap();
        let bloom_filter_path = format!("{}/{}.bloom", child_directory, filename);

        let loaded_filter = load_bloom_filter_from_disk(bloom_filter_path.as_str()).unwrap();

        let result = loaded_filter.check(&"chr1:22-30".to_string());
        pretty_assertions::assert_eq!(true, result);

        let result = loaded_filter.check(&"chr1:23-31".to_string());
        pretty_assertions::assert_eq!(false, result);
    }

    #[rstest]
    fn test_process_bed_directory() {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");

        let input_dir = PathBuf::from(path_to_crate)
            .parent()
            .unwrap()
            .join("tests/data/dir_of_files/dir_beds");

        let tempdir = tempfile::tempdir().unwrap();
        let output_dir = tempdir.path();

        let num_of_items = 1000;
        let false_positive_rate = 0.5;

        let sample_bed = input_dir.join("dummy2.bed");
        let tokenizer = Tokenizer::from_auto(sample_bed.to_str().unwrap())
            .expect("Failed to create tokenizer from config.");

        let processed_count = process_bed_directory(
            &tokenizer,
            input_dir.to_str().unwrap(),
            output_dir.to_str().unwrap(),
            num_of_items,
            false_positive_rate,
        )
        .unwrap();

        pretty_assertions::assert_ne!(processed_count, 0, "Should process at least one BED file");

        // Verify bloom files were created
        let bloom_files: Vec<_> = fs::read_dir(output_dir)
            .unwrap()
            .filter_map(|e| e.ok())
            .filter(|e| e.path().extension().and_then(|s| s.to_str()) == Some("bloom"))
            .collect();

        pretty_assertions::assert_eq!(bloom_files.len(), processed_count);

        // Now load them back into memory
        let loaded_filters = load_bloom_directory(output_dir.to_str().unwrap()).unwrap();

        pretty_assertions::assert_eq!(loaded_filters.len(), processed_count);

        // Test that we can query the loaded filters
        for (filename, filter) in loaded_filters.iter() {
            println!("Loaded filter for: {}", filename);
            let _ = filter.check(&"chr1:22-30".to_string());
        }
    }
}

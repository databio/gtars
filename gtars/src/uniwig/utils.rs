use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;
use std::str::FromStr;
use crate::uniwig::reading::{create_chrom_vec_default_score, create_chrom_vec_scores};
use crate::uniwig::{Chromosome, FileType};
use serde_json::Value;
use std::collections::HashMap;
use std::env;
use std::fs::{self, File};
use std::path::Path;
use std::io::Write;
use std::io::Read;
use byteorder::{LittleEndian, ReadBytesExt};
use ndarray::{Array, Ix1};
use ndarray_npy::read_npy;


pub struct FileInfo {
    pub file_type: FileType,
    pub is_gzipped: bool,
}

pub fn get_file_info(path: &PathBuf) -> FileInfo {
    let mut file_type = FileType::UNKNOWN;
    let mut is_gzipped = false;

    if let Some(os_str_filename) = path.file_name() {
        if let Some(filename) = os_str_filename.to_str() {
            // Check for .gz first
            if filename.ends_with(".gz") {
                is_gzipped = true;
                if let Some(base_filename) = filename.strip_suffix(".gz") {
                    // Try to get the extension before .gz
                    if let Some(ext) = PathBuf::from(base_filename).extension().and_then(|e| e.to_str()) {
                        file_type = FileType::from_str(ext).unwrap_or(FileType::UNKNOWN);
                    } else {
                        // If there's no extension before .gz (e.g., "my_data.gz"),
                        // you might want to handle this specifically or leave as UNKNOWN.
                        // For now, we'll try to parse the whole base_filename as a type
                        file_type = FileType::from_str(base_filename).unwrap_or(FileType::UNKNOWN);
                    }
                }
            } else {
                // Not gzipped, just get the direct extension
                if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
                    file_type = FileType::from_str(ext).unwrap_or(FileType::UNKNOWN);
                }
            }
        }
    }

    FileInfo {
        file_type,
        is_gzipped,
    }
}


/// Attempt to compress counts before writing to bedGraph
pub fn compress_counts(
    count_results: &mut (Vec<u32>, Vec<i32>),
    start_position: i32,
) -> (Vec<u32>, Vec<u32>, Vec<u32>) {
    let mut final_starts: Vec<u32> = Vec::new();
    let mut final_ends: Vec<u32> = Vec::new();
    let mut final_counts: Vec<u32> = Vec::new();

    // .0 are the counts, .1 are the positions to track
    let mut previous_count = count_results.0[0];

    let previous_start = start_position as u32;
    let mut current_start = previous_start;

    let mut current_end = start_position as u32;

    for (u, _i) in count_results.0.iter().zip(count_results.1.iter()) {
        let current_count = *u;
        current_end += 1;

        if current_count != previous_count {
            final_starts.push(current_start);
            final_ends.push(current_end);
            final_counts.push(previous_count);
            current_start = current_end;
            previous_count = current_count;
        } else {
            previous_count = current_count;
        }
    }

    // Must add these lines else we will not get the closing interval (since previous count will be = current count at the close).
    final_starts.push(current_start);
    final_ends.push(current_end);
    final_counts.push(previous_count);

    // println!("Final Starts:{:?}", final_starts);
    // println!("Final Ends:{:?}", final_ends);
    // println!("Final Counts:{:?}", final_counts);

    (final_starts, final_ends, final_counts)
}

/// Determine if there is a size associated with a Chromosome
/// Only return chromosomes that have an associated size.
pub fn get_final_chromosomes(
    ft: &Result<FileType, String>,
    filepath: &str,
    chrom_sizes: &std::collections::HashMap<String, u32>,
    score: bool,
) -> Vec<Chromosome> {


    let mut chromosomes: Vec<Chromosome> = Vec::new();
    
    let path = PathBuf::from(filepath);

    if path.is_dir() {
        
            let mut combined_chromosome_map: HashMap<String, Chromosome> = HashMap::new();
        
            for entry_result in fs::read_dir(path).unwrap() {
                let entry = entry_result.unwrap();
                let single_file_path = entry.path();

                
        
                if single_file_path.is_file() {
                    let file_info = get_file_info(&single_file_path);
                        let single_file_path = single_file_path.to_str().unwrap();
                    match file_info.file_type {
                        
                            FileType::BED | FileType::NARROWPEAK => {
                                println!("Processing file: {}", single_file_path);
                                //let chromosomes_from_file = read_bed_vec(single_file_path);
                                let chromosomes_from_file = if score {
                                    create_chrom_vec_scores(single_file_path)
                                } else {
                                    create_chrom_vec_default_score(single_file_path) // Use bed reader if score flag is false
                                };
                                // Merge chromosomes from this file into the combined map
                                for chrom_data in chromosomes_from_file {
                                    let entry = combined_chromosome_map
                                        .entry(chrom_data.chrom.clone())
                                        .or_insert_with(|| Chromosome {
                                            chrom: chrom_data.chrom,
                                            starts: Vec::new(),
                                            ends: Vec::new(),
                                        });
                                    entry.starts.extend(chrom_data.starts);
                                    entry.ends.extend(chrom_data.ends);
                                }
                            }
                            FileType::BAM => {
                                println!("WARNING: Skipping BAM file ({}). Not supported at this time for direct parsing.", single_file_path);
                            }
                            FileType::UNKNOWN => {
                                println!("WARNING: Skipping file with unknown extension: {}", single_file_path);
                            }
                        }
                }
            }
        
            // Convert the combined map back to a Vec<Chromosome> and ensure final sorting
            //let mut final_chromosomes: Vec<Chromosome> = combined_chromosome_map.into_values().collect();
        //chromosomes = combined_chromosome_map.into_values().collect();
        let mut final_chromosomes: Vec<Chromosome> = combined_chromosome_map.into_values().collect();

        for chromosome in &mut final_chromosomes {
            // Sort the starts and ends for each chromosome after all data for it has been aggregated
            chromosome.starts.sort_unstable_by_key(|&(pos, _)| pos);
            chromosome.ends.sort_unstable_by_key(|&(pos, _)| pos);
        }
        // And finally, sort the chromosomes themselves by name for consistent output
        final_chromosomes.sort_unstable_by(|a, b| a.chrom.cmp(&b.chrom));
        
        chromosomes = final_chromosomes;
    }
    
    else{

        chromosomes = match ft {
            Ok(FileType::BED) => create_chrom_vec_default_score(filepath),
            Ok(FileType::NARROWPEAK) => {
                if score {
                    println!("FileType is NarrowPeak and Score = True...Counting based on Score");
                    create_chrom_vec_scores(filepath) // if score flag enabled, this will attempt to read narrowpeak scores
                } else {
                    create_chrom_vec_default_score(filepath)
                }
            }
            _ => create_chrom_vec_default_score(filepath),
        };
        
        
    }
    

    let num_chromosomes = chromosomes.len();

    println!("PreProcessing each chromosome...");
    let mut final_chromosomes: Vec<Chromosome> = Vec::with_capacity(num_chromosomes);
    for chromosome in chromosomes.iter() {
        if chromosome.starts.len() != chromosome.ends.len() {
            break;
        }

        // Check if there is an available chrom size, if not exclude it from our final list
        let _current_chrom_size = match chrom_sizes.get(&chromosome.chrom) {
            Some(size) => *size as i32, // Dereference to get the i32 value
            None => {
                continue; // Or handle the error differently
            }
        };

        final_chromosomes.push(chromosome.clone())
    }

    println!(
        "Initial chroms: {}  vs Final chroms: {}",
        chromosomes.len(),
        final_chromosomes.len()
    );
    if chromosomes.len() != final_chromosomes.len() {
        println!("Some chromosomes were not found in chrom.sizes file and will be skipped...")
    }

    final_chromosomes
}


/// Custom comparator for version sorting
pub fn version_sort(a: &String, b: &String) -> std::cmp::Ordering {
    use std::cmp::Ordering;

    let mut split_a = a.split(|c: char| !c.is_numeric()).filter_map(|s| s.parse::<usize>().ok());
    let mut split_b = b.split(|c: char| !c.is_numeric()).filter_map(|s| s.parse::<usize>().ok());

    loop {
        match (split_a.next(), split_b.next()) {
            (Some(x), Some(y)) => match x.cmp(&y) {
                Ordering::Equal => continue,
                ord => return ord,
            },
            (Some(_), None) => return Ordering::Greater,
            (None, Some(_)) => return Ordering::Less,
            (None, None) => return a.cmp(b), // Fallback to lexicographical if needed
        }
    }
}

pub fn read_u32_npy(npy_file_path: &Path) -> Result<Vec<u32>, Box<dyn std::error::Error>> {
    // Open the file
    let mut file = File::open(npy_file_path)?;

    // Read the entire file into a buffer
    let mut buffer = vec![];
    file.read_to_end(&mut buffer)?;

    // Skip the header
    let header_end = buffer
        .iter()
        .position(|&b| b == b'\n') // Find the end of the header
        .ok_or("Invalid NPY file: missing header newline")?
        + 1;

    let mut cursor = &buffer[header_end..]; // Skip to the data section
    let mut values = vec![];

    // Read remaining bytes as `u32` in little-endian order
    while let Ok(value) = cursor.read_u32::<LittleEndian>() {
        values.push(value);
    }

    Ok(values)
}


// two inputs are equivalent to --fileheader when running uniwig
pub fn npy_to_wig(npy_header: &Path, wig_header: &Path) -> Result<(), Box<dyn std::error::Error>> {
    //TODO add test
    // Create output/wiggle folder
    std::fs::create_dir_all(wig_header)?;
    // Read the JSON file
    let input_file_path = npy_header.join("npy_meta.json");
    let json_data = fs::read_to_string(&input_file_path)?;

    // Deserialize JSON into a HashMap (unordered)
    let dictionary: HashMap<String, HashMap<String, i32>> = serde_json::from_str(&json_data)?;

    // Sort outer keys using version sorting
    let mut sorted_outer_keys: Vec<String> = dictionary.keys().cloned().collect();
    sorted_outer_keys.sort_by(version_sort);

    // Define the list of inner keys to include
    let inner_keys_filter = vec!["start", "core", "end"];
    let step_key = "stepsize";

    // Iterate through the list of inner keys
    for target_inner_key in &inner_keys_filter {
        println!("Preparing {} wiggle file", target_inner_key);
        // Construct the output file name
        let output_file_path = wig_header.join(format!("{}_{}.wig", wig_header.display(), target_inner_key));
        let mut output_file = File::create(&output_file_path)?;

        // Check this inner key across all sorted outer dictionaries
        for outer_key in &sorted_outer_keys {
            let inner_dict = dictionary.get(outer_key).unwrap();
            let mut value = *inner_dict.get(*target_inner_key).unwrap();
            if *target_inner_key == "start" || *target_inner_key == "core" {
                value += 1;
            }
            let step_value = inner_dict.get(step_key).unwrap();
            writeln!(
                output_file,
                "fixedStep chrom={} start={} step={}",
                outer_key, value, step_value
            )?;

            let npy_file_path = npy_header.join(format!("{}_{}.npy", outer_key, target_inner_key));
            let array = read_u32_npy(&npy_file_path)?;

            // Write the array values row by row
            for value in array.iter() {
                writeln!(output_file, "{}", value)?;
            }
        }
    }

    Ok(())
}



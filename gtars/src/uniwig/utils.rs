use crate::uniwig::reading::{read_bed_vec, read_narrow_peak_vec};
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
        current_end = current_end + 1;

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
    let chromosomes: Vec<Chromosome> = match ft {
        Ok(FileType::BED) => read_bed_vec(filepath),
        Ok(FileType::NARROWPEAK) => {
            if score {
                println!("FileType is NarrowPeak and Score = True...Counting based on Score");
                read_narrow_peak_vec(filepath) // if score flag enabled, this will attempt to read narrowpeak scores
            } else {
                read_bed_vec(filepath)
            }
        }
        _ => read_bed_vec(filepath),
    };

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


// how npy_to_wig was run in local test:
// fn main() {
//     let args: Vec<String> = env::args().collect();
//     if args.len() < 3 {
//         eprintln!("Usage: cargo run <npy_file_header> <wiggle_file_header>");
//         std::process::exit(1);
//     }

//     let npy_header = Path::new(&args[1]);
//     let wig_header = Path::new(&args[2]);

//     if let Err(e) = npy_to_wig(npy_header, wig_header) {
//         eprintln!("Error: {}", e);
//         std::process::exit(1);
//     }
// }

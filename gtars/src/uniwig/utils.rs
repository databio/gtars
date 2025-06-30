use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;
use std::str::FromStr;
use crate::uniwig::reading::{read_bed_vec, read_narrow_peak_vec};
use crate::uniwig::{Chromosome, FileType};

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
                    if let Some(ext_str) = single_file_path.extension().and_then(|ext| ext.to_str()) {
                        let ft_string = FileType::from_str(ext_str);
                        let single_file_path = single_file_path.to_str().unwrap();
                        match ft_string {
                            Ok(FileType::BED) => {
                                println!("Processing BED file: {}", single_file_path);
                                let chromosomes_from_file = read_bed_vec(single_file_path);
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
                            Ok(FileType::NARROWPEAK) => {
                                println!("Processing NARROWPEAK file: {}", single_file_path);
                                let chromosomes_from_file = if score {
                                    read_narrow_peak_vec(single_file_path)
                                } else {
                                    read_bed_vec(single_file_path) // Use bed reader if score flag is false
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
                            Ok(FileType::BAM) => {
                                println!("WARNING: Skipping BAM file ({}). Not supported at this time for direct parsing.", single_file_path);
                            }
                            Err(_) => {
                                println!("WARNING: Skipping file with unknown extension: {}", single_file_path);
                            }
                        }
                    } else {
                        println!("WARNING: Skipping file with no extension: {}", single_file_path.display());
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

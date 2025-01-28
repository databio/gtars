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

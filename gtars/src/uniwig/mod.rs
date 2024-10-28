use clap::ArgMatches;

use indicatif::ProgressBar;

use rayon::prelude::*;
use std::error::Error;

use std::io::{BufWriter, Write};

use crate::uniwig::counting::{core_counts, start_end_counts};
use crate::uniwig::reading::{
    get_seq_reads_bam, read_bam_header, read_bed_vec, read_chromosome_sizes, read_narrow_peak_vec,
};
use crate::uniwig::utils::{compress_counts, get_final_chromosomes};
use crate::uniwig::writing::{
    write_bw_files, write_combined_files, write_to_bed_graph_file, write_to_npy_file,
    write_to_wig_file,
};
use std::str::FromStr;
// use noodles::sam as sam;
//use bstr::BString;

pub mod cli;
pub mod counting;
pub mod reading;
mod utils;
pub mod writing;

pub mod consts {
    pub const UNIWIG_CMD: &str = "uniwig";
}

#[derive(Debug)]
enum FileType {
    BED,
    BAM,
    NARROWPEAK,
}

impl FromStr for FileType {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "bed" => Ok(FileType::BED),
            "bam" => Ok(FileType::BAM),
            "narrowpeak" => Ok(FileType::NARROWPEAK),
            _ => Err(format!("Invalid file type: {}", s)),
        }
    }
}

// Chromosome representation for Bed File Inputs
#[derive(Debug)]
pub struct Chromosome {
    pub chrom: String,
    pub starts: Vec<(i32, i32)>,
    pub ends: Vec<(i32, i32)>,
}
impl Clone for Chromosome {
    fn clone(&self) -> Self {
        Self {
            chrom: self.chrom.clone(),
            starts: self.starts.clone(),
            ends: self.ends.clone(),
        }
    }
}

/// Matches items from CLAP args before running uniwig_main
pub fn run_uniwig(matches: &ArgMatches) {
    //println!("I am running. Here are the arguments: {:?}", matches);

    let filepath = matches
        .get_one::<String>("file")
        .expect("file path is required");

    let filetype = matches
        .get_one::<String>("filetype")
        .expect("file type is required");

    let chromsizerefpath = matches
        .get_one::<String>("chromref")
        .cloned()
        .unwrap_or_else(|| filepath.clone());

    let bwfileheader = matches
        .get_one::<String>("fileheader")
        .expect("fileheader is required");

    let smoothsize = matches
        .get_one::<i32>("smoothsize")
        .expect("smoothsize required");

    let output_type = matches
        .get_one::<String>("outputtype")
        .expect("output type is required");

    let num_threads = matches
        .get_one::<i32>("threads")
        .expect("requires integer value");

    let score = matches.get_one::<bool>("score").unwrap_or_else(|| &false);

    let stepsize = matches
        .get_one::<i32>("stepsize")
        .expect("requires integer value");

    let zoom = matches
        .get_one::<i32>("zoom")
        .expect("requires integer value");

    uniwig_main(
        *smoothsize,
        filepath,
        chromsizerefpath.as_str(),
        bwfileheader,
        output_type,
        filetype,
        *num_threads,
        *score,
        *stepsize,
        *zoom,
    )
    .expect("Uniwig failed.");
}

/// Ensures that the start position for every wiggle file is at a minimum equal to `1`
fn clamped_start_position(start: i32, smoothsize: i32) -> i32 {
    std::cmp::max(1, start - smoothsize)
}

/// Main function
pub fn uniwig_main(
    smoothsize: i32,
    filepath: &str,
    chromsizerefpath: &str,
    bwfileheader: &str,
    output_type: &str,
    filetype: &str,
    num_threads: i32,
    score: bool,
    stepsize: i32,
    zoom: i32,
) -> Result<(), Box<dyn Error>> {
    // Must create a Rayon thread pool in which to run our iterators
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads as usize)
        .build()
        .unwrap();

    // Determine File Type
    let ft = FileType::from_str(filetype.to_lowercase().as_str());
    // Set up output file names

    let og_output_type = output_type; // need this later for conversion
    let mut output_type = output_type;
    if output_type == "bedgraph" || output_type == "bw" || output_type == "bigwig" {
        output_type = "bedGraph" // we must create bedgraphs first before creating bigwig files
    }

    let mut meta_data_file_names: [String; 3] = [
        "placeholder1".to_owned(),
        "placeholder2".to_owned(),
        "placeholder3".to_owned(),
    ];

    meta_data_file_names[0] = format!("{}{}.{}", bwfileheader, "start", "meta");
    meta_data_file_names[1] = format!("{}{}.{}", bwfileheader, "end", "meta");
    meta_data_file_names[2] = format!("{}{}.{}", bwfileheader, "core", "meta");

    let chrom_sizes = match read_chromosome_sizes(chromsizerefpath) {
        // original program gets chromosome size from a .sizes file, e.g. chr1 248956422
        // the original program simply pushes 0's until the end of the chromosome length and writes these to file.
        // can we instead just use the last endsite for each chromosome to save space in th wiggle file?
        Ok(chrom_sizes) => chrom_sizes,
        Err(err) => {
            println!("Error reading chromosome sizes: {}", err);
            return Err(Box::from("An error occurred")); // Exit the main function on error
        }
    };



    match ft {
        //BED AND NARROWPEAK WORKFLOW
        Ok(FileType::BED) | Ok(FileType::NARROWPEAK) => {

            let mut final_chromosomes = get_final_chromosomes(&ft, filepath, &chrom_sizes,score);

            let bar = ProgressBar::new(final_chromosomes.len() as u64);

            // Pool installs iterator
            pool.install(|| {
                final_chromosomes
                    .par_iter_mut()
                    .for_each(|chromosome: &mut Chromosome| {
                        bar.inc(1);
                        match ft {
                            Ok(FileType::BAM) => {
                                get_seq_reads_bam(chromosome, filepath);
                            }
                            _ => {}
                        };

                        let primary_start = chromosome.starts[0].clone();
                        let primary_end = chromosome.ends[0].clone();

                        let current_chrom_size = *chrom_sizes.get(&chromosome.chrom).unwrap() as i32;
                        let chrom_name = chromosome.chrom.clone();

                        // Iterate 3 times to output the three different files.
                        for j in 0..3 {
                            // Original code uses:
                            // bwOpen, then bwCreateChromList, then bwWriteHdr

                            let mut _success_count = 0;
                            let mut _failure_count = 0;

                            if smoothsize != 0 {
                                match j {
                                    0 => {
                                        let mut count_result = match ft {
                                            Ok(FileType::BED) => start_end_counts(
                                                &chromosome.starts,
                                                current_chrom_size,
                                                smoothsize,
                                                stepsize,
                                            ),
                                            _ => start_end_counts(
                                                &chromosome.starts,
                                                current_chrom_size,
                                                smoothsize,
                                                stepsize,
                                            ),
                                        };

                                        match output_type {
                                            "file" => {
                                                //print!("Writing to CLI");
                                                let handle = &std::io::stdout();
                                                let mut buf = BufWriter::new(handle);
                                                for count in &count_result.0 {
                                                    writeln!(buf, "{}", count)
                                                        .expect("failed to write line");
                                                }
                                                buf.flush().unwrap();
                                            }
                                            "wig" => {
                                                //println!("Writing to wig file!");
                                                let file_name = format!(
                                                    "{}{}_{}.{}",
                                                    bwfileheader, chrom_name, "start", output_type
                                                );
                                                write_to_wig_file(
                                                    &count_result.0,
                                                    file_name.clone(),
                                                    chrom_name.clone(),
                                                    clamped_start_position(primary_start.0, smoothsize),
                                                    stepsize,
                                                );
                                            }
                                            "bedGraph" => {
                                                let file_name = format!(
                                                    "{}{}_{}.{}",
                                                    bwfileheader, chrom_name, "start", output_type
                                                );
                                                let count_info: (Vec<u32>, Vec<u32>, Vec<u32>) =
                                                    compress_counts(
                                                        &mut count_result,
                                                        clamped_start_position(primary_start.0, smoothsize),
                                                    );
                                                write_to_bed_graph_file(
                                                    &count_info,
                                                    file_name.clone(),
                                                    chrom_name.clone(),
                                                    stepsize,
                                                );
                                            }
                                            "csv" => {
                                                panic!("Write to CSV. Not Implemented");
                                            }
                                            "npy" => {
                                                let file_name = format!(
                                                    "{}{}_{}.{}",
                                                    bwfileheader, chrom_name, "start", output_type
                                                );
                                                write_to_npy_file(
                                                    &count_result.0,
                                                    file_name.clone(),
                                                    chrom_name.clone(),
                                                    clamped_start_position(primary_start.0, smoothsize),
                                                    stepsize,
                                                    meta_data_file_names[0].clone(),
                                                );
                                            }
                                            _ => {
                                                println!("Defaulting to npy file...");
                                                let file_name = format!(
                                                    "{}{}_{}.{}",
                                                    bwfileheader, chrom_name, "start", output_type
                                                );
                                                write_to_npy_file(
                                                    &count_result.0,
                                                    file_name.clone(),
                                                    chrom_name.clone(),
                                                    clamped_start_position(primary_start.0, smoothsize),
                                                    stepsize,
                                                    meta_data_file_names[0].clone(),
                                                );
                                            }
                                        }
                                    }
                                    1 => {
                                        let mut count_result = match ft {
                                            Ok(FileType::BED) => start_end_counts(
                                                &chromosome.ends,
                                                current_chrom_size,
                                                smoothsize,
                                                stepsize,
                                            ),
                                            _ => start_end_counts(
                                                &chromosome.ends,
                                                current_chrom_size,
                                                smoothsize,
                                                stepsize,
                                            ),
                                        };

                                        match output_type {
                                            "file" => {
                                                let handle = &std::io::stdout();
                                                let mut buf = BufWriter::new(handle);
                                                for count in &count_result.0 {
                                                    writeln!(buf, "{}", count)
                                                        .expect("failed to write line");
                                                }
                                                buf.flush().unwrap();
                                            }
                                            "bedGraph" => {
                                                let file_name = format!(
                                                    "{}{}_{}.{}",
                                                    bwfileheader, chrom_name, "end", output_type
                                                );

                                                let count_info: (Vec<u32>, Vec<u32>, Vec<u32>) =
                                                    compress_counts(
                                                        &mut count_result,
                                                        clamped_start_position(primary_end.0, smoothsize),
                                                    );
                                                write_to_bed_graph_file(
                                                    &count_info,
                                                    file_name.clone(),
                                                    chrom_name.clone(),
                                                    stepsize,
                                                );
                                            }
                                            "wig" => {
                                                let file_name = format!(
                                                    "{}{}_{}.{}",
                                                    bwfileheader, chrom_name, "end", output_type
                                                );
                                                write_to_wig_file(
                                                    &count_result.0,
                                                    file_name.clone(),
                                                    chrom_name.clone(),
                                                    clamped_start_position(primary_end.0, smoothsize),
                                                    stepsize,
                                                );
                                            }
                                            "csv" => {
                                                panic!("Write to CSV. Not Implemented");
                                            }
                                            "npy" => {
                                                let file_name = format!(
                                                    "{}{}_{}.{}",
                                                    bwfileheader, chrom_name, "end", output_type
                                                );
                                                write_to_npy_file(
                                                    &count_result.0,
                                                    file_name.clone(),
                                                    chrom_name.clone(),
                                                    clamped_start_position(primary_start.0, smoothsize),
                                                    stepsize,
                                                    meta_data_file_names[1].clone(),
                                                );
                                            }
                                            _ => {
                                                println!("Defaulting to npy file...");
                                                let file_name = format!(
                                                    "{}{}_{}.{}",
                                                    bwfileheader, chrom_name, "end", output_type
                                                );
                                                write_to_npy_file(
                                                    &count_result.0,
                                                    file_name.clone(),
                                                    chrom_name.clone(),
                                                    clamped_start_position(primary_start.0, smoothsize),
                                                    stepsize,
                                                    meta_data_file_names[1].clone(),
                                                );
                                            }
                                        }
                                    }
                                    2 => {
                                        let mut core_results = match ft {
                                            Ok(FileType::BED) => core_counts(
                                                &chromosome.starts,
                                                &chromosome.ends,
                                                current_chrom_size,
                                                stepsize,
                                            ),
                                            _ => core_counts(
                                                &chromosome.starts,
                                                &chromosome.ends,
                                                current_chrom_size,
                                                stepsize,
                                            ),
                                        };

                                        match output_type {
                                            "file" => {
                                                let handle = &std::io::stdout();
                                                let mut buf = BufWriter::new(handle);
                                                for count in &core_results.0 {
                                                    writeln!(buf, "{}", count)
                                                        .expect("failed to write line");
                                                }
                                                buf.flush().unwrap();
                                            }
                                            "bedGraph" => {
                                                let file_name = format!(
                                                    "{}{}_{}.{}",
                                                    bwfileheader, chrom_name, "core", output_type
                                                );

                                                let count_info: (Vec<u32>, Vec<u32>, Vec<u32>) =
                                                    compress_counts(&mut core_results, primary_start.0);
                                                write_to_bed_graph_file(
                                                    &count_info,
                                                    file_name.clone(),
                                                    chrom_name.clone(),
                                                    stepsize,
                                                );
                                            }
                                            "wig" => {
                                                let file_name = format!(
                                                    "{}{}_{}.{}",
                                                    bwfileheader, chrom_name, "core", output_type
                                                );
                                                write_to_wig_file(
                                                    &core_results.0,
                                                    file_name.clone(),
                                                    chrom_name.clone(),
                                                    primary_start.0,
                                                    stepsize,
                                                );
                                            }
                                            "csv" => {
                                                panic!("Write to CSV. Not Implemented");
                                            }
                                            "npy" => {
                                                let file_name = format!(
                                                    "{}{}_{}.{}",
                                                    bwfileheader, chrom_name, "core", output_type
                                                );
                                                write_to_npy_file(
                                                    &core_results.0,
                                                    file_name.clone(),
                                                    chrom_name.clone(),
                                                    primary_start.0,
                                                    stepsize,
                                                    meta_data_file_names[2].clone(),
                                                );
                                            }
                                            _ => {
                                                println!("Defaulting to npy file...");
                                                let file_name = format!(
                                                    "{}{}_{}.{}",
                                                    bwfileheader, chrom_name, "core", output_type
                                                );
                                                write_to_npy_file(
                                                    &core_results.0,
                                                    file_name.clone(),
                                                    chrom_name.clone(),
                                                    primary_start.0,
                                                    stepsize,
                                                    meta_data_file_names[2].clone(),
                                                );
                                            }
                                        }
                                    }
                                    _ => panic!("Unexpected value: {}", j), // Handle unexpected values
                                }
                            }
                        }
                    })
            });

            bar.finish();

            let vec_strings = vec!["start", "core", "end"];
            //let vec_strings = vec!["start"];

            let bar = ProgressBar::new(vec_strings.len() as u64);
            match output_type {
                "wig" | "bedGraph" => {
                    println!("Combining {} Files", output_type);

                    for location in vec_strings.iter() {
                        bar.inc(1);
                        write_combined_files(*location, output_type, bwfileheader, &final_chromosomes);
                    }
                }
                _ => {}
            }
            bar.finish();

            match og_output_type {
                "bw" | "bigWig" => {
                    println!("Writing bigWig files");
                    write_bw_files(bwfileheader, chromsizerefpath, num_threads, zoom);
                }

                _ => {}
            }


        }
        //BAM REQUIRES DIFFERENT WORKFLOW
        Ok(FileType::BAM) => {

            println!("Do nothing for bam right now");


                },

        _ => {
            panic!("Unknown File Type provided");
        },
    };


    println!("FINISHED");

    Ok(())
}

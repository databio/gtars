use clap::ArgMatches;
use std::collections::HashMap;

use indicatif::ProgressBar;

use rayon::prelude::*;
use std::error::Error;
use std::fs::{remove_file, File};
use std::io::{BufRead, BufReader, BufWriter, Write};

use crate::uniwig::counting::{
    bam_to_bed_no_counts, core_counts, start_end_counts, variable_core_counts_bam_to_bw,
    variable_shifted_bam_to_bw, variable_start_end_counts_bam_to_bw, BAMRecordError,
};
use crate::uniwig::reading::read_chromosome_sizes;
use crate::uniwig::utils::{compress_counts, get_final_chromosomes};
use crate::uniwig::writing::{
    write_bw_files, write_combined_files, write_to_bed_graph_file, write_to_npy_file,
    write_to_wig_file,
};
use bigtools::beddata::BedParserStreamingIterator;
use bigtools::utils::cli::bedgraphtobigwig::BedGraphToBigWigArgs;
use bigtools::utils::cli::bigwigmerge::{get_merged_vals, ChromGroupReadImpl};
use bigtools::utils::cli::BBIWriteArgs;
use bigtools::utils::reopen::ReopenableFile;
use bigtools::{BigWigRead, BigWigWrite, InputSortType};
use noodles::bam;
use noodles::bam::io::reader::Query;
use noodles::bgzf::Reader;
use os_pipe::PipeWriter;
use rayon::ThreadPool;
use std::path::PathBuf;
use std::str::FromStr;
use std::sync::{Arc, Mutex};
use std::thread;
use tokio::runtime;

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

// #[derive(Debug)]
// enum OutSelection {
//     STARTS,
//     ENDS,
//     CORE,
// }

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

    //let default_vec = &vec!["start", "end", "core"];
    let count_types = matches
        .get_one::<String>("counttype")
        .expect("output type is required");

    // let mut vec_count_type: Vec<&str> = Vec::new();
    let vec_count_type = match count_types.as_str() {
        "all" => {
            vec!["start", "end", "core"]
        }
        "start" => {
            vec!["start"]
        }
        "end" => {
            vec!["end"]
        }
        "core" => {
            vec!["core"]
        }
        "shift" => {
            vec!["shift"]
        }

        _ => {
            vec!["start", "end", "core"]
        }
    };

    //println!("FOUND count_type {:?}", vec_count_type);

    let num_threads = matches
        .get_one::<i32>("threads")
        .expect("requires integer value");

    let bam_scale = matches
        .get_one::<f32>("bamscale")
        .expect("requires int value");

    let score = matches.get_one::<bool>("score").unwrap_or_else(|| &false);
    let bam_shift = matches
        .get_one::<bool>("no-bamshift")
        .unwrap_or_else(|| &true);

    let debug = matches.get_one::<bool>("debug").unwrap_or_else(|| &false);

    let stepsize = matches
        .get_one::<i32>("stepsize")
        .expect("requires integer value");

    let zoom = matches
        .get_one::<i32>("zoom")
        .expect("requires integer value");

    uniwig_main(
        vec_count_type,
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
        *debug,
        *bam_shift,
        *bam_scale,
    )
    .expect("Uniwig failed.");
}

/// Ensures that the start position is at a minimum equal to `1`
fn clamped_start_position(start: i32, smoothsize: i32, wig_shift: i32) -> i32 {
    std::cmp::max(1, start - smoothsize + wig_shift)
}
/// Ensure that the start position is at a minimum equal to `0`
fn clamped_start_position_zero_pos(start: i32, smoothsize: i32) -> i32 {
    std::cmp::max(0, start - smoothsize)
}

/// Main function
pub fn uniwig_main(
    vec_count_type: Vec<&str>,
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
    debug: bool,
    bam_shift: bool,
    bam_scale: f32,
) -> Result<(), Box<dyn Error>> {
    // Must create a Rayon thread pool in which to run our iterators
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads as usize)
        .build()
        .unwrap();

    // Determine Input File Type
    let input_filetype = FileType::from_str(filetype.to_lowercase().as_str());
    // Set up output file names

    let mut meta_data_file_names: [String; 3] = [
        "placeholder1".to_owned(),
        "placeholder2".to_owned(),
        "placeholder3".to_owned(),
    ];

    meta_data_file_names[0] = format!("{}{}.{}", bwfileheader, "start", "meta");
    meta_data_file_names[1] = format!("{}{}.{}", bwfileheader, "end", "meta");
    meta_data_file_names[2] = format!("{}{}.{}", bwfileheader, "core", "meta");

    let mut npy_meta_data_map: HashMap<String, HashMap<String, i32>> = HashMap::new();

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

    match input_filetype {
        //BED AND NARROWPEAK WORKFLOW
        Ok(FileType::BED) | Ok(FileType::NARROWPEAK) => {
            // Pare down chromosomes if necessary
            let mut final_chromosomes =
                get_final_chromosomes(&input_filetype, filepath, &chrom_sizes, score);

            // Some housekeeping depending on output type
            let og_output_type = output_type; // need this later for conversion
            let mut output_type = output_type;
            if output_type == "bedgraph" || output_type == "bw" || output_type == "bigwig" {
                output_type = "bedGraph" // we must create bedgraphs first before creating bigwig files
            }

            let bar = ProgressBar::new(final_chromosomes.len() as u64);

            // Pool installs iterator via rayon crate
            pool.install(|| {
                final_chromosomes
                    .par_iter_mut()
                    .for_each(|chromosome: &mut Chromosome| {
                        bar.inc(1);

                        let primary_start = chromosome.starts[0].clone();
                        let primary_end = chromosome.ends[0].clone();

                        let current_chrom_size =
                            *chrom_sizes.get(&chromosome.chrom).unwrap() as i32;
                        let chrom_name = chromosome.chrom.clone();

                        // Iterate 3 times to output the three different files.
                        for j in 0..3 {
                            // todo change these to be ooptional based on vec_count_type
                            // Original code uses:
                            // bwOpen, then bwCreateChromList, then bwWriteHdr

                            let mut _success_count = 0;
                            let mut _failure_count = 0;

                            if smoothsize != 0 {
                                match j {
                                    0 => {
                                        let mut count_result = start_end_counts(
                                            &chromosome.starts,
                                            current_chrom_size,
                                            smoothsize,
                                            stepsize,
                                        );

                                        match output_type {
                                            "file" => {
                                                panic!("Writing to file currently not supported");
                                            }
                                            "csv" => {
                                                panic!("Write to CSV. Not Implemented");
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
                                                    clamped_start_position(
                                                        primary_start.0,
                                                        smoothsize,
                                                        1, //must shift wiggle starts and core by 1 since it is 1 based
                                                    ),
                                                    stepsize,
                                                    current_chrom_size,
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
                                                        clamped_start_position_zero_pos(
                                                            primary_start.0,
                                                            smoothsize,
                                                        ),
                                                    );
                                                write_to_bed_graph_file(
                                                    &count_info,
                                                    file_name.clone(),
                                                    chrom_name.clone(),
                                                    stepsize,
                                                );
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
                                                    clamped_start_position_zero_pos(
                                                        primary_start.0,
                                                        smoothsize,
                                                    ),
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
                                                    clamped_start_position_zero_pos(
                                                        primary_start.0,
                                                        smoothsize,
                                                    ),
                                                    stepsize,
                                                    meta_data_file_names[0].clone(),
                                                );
                                            }
                                        }
                                    }
                                    1 => {
                                        let mut count_result = start_end_counts(
                                            &chromosome.ends,
                                            current_chrom_size,
                                            smoothsize,
                                            stepsize,
                                        );
                                        match output_type {
                                            "file" => {
                                                panic!("Writing to file not currently supported.")
                                            }
                                            "csv" => {
                                                panic!("Write to CSV. Not Implemented");
                                            }
                                            "bedGraph" => {
                                                let file_name = format!(
                                                    "{}{}_{}.{}",
                                                    bwfileheader, chrom_name, "end", output_type
                                                );

                                                let count_info: (Vec<u32>, Vec<u32>, Vec<u32>) =
                                                    compress_counts(
                                                        &mut count_result,
                                                        clamped_start_position(
                                                            primary_end.0,
                                                            smoothsize,
                                                            0,
                                                        ),
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
                                                    clamped_start_position(
                                                        primary_end.0,
                                                        smoothsize,
                                                        0, // ends already 1 based, do not shift further
                                                    ),
                                                    stepsize,
                                                    current_chrom_size,
                                                );
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
                                                    clamped_start_position(
                                                        primary_end.0,
                                                        smoothsize,
                                                        0,
                                                    ),
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
                                                    clamped_start_position(
                                                        primary_end.0,
                                                        smoothsize,
                                                        0,
                                                    ),
                                                    stepsize,
                                                    meta_data_file_names[1].clone(),
                                                );
                                            }
                                        }
                                    }
                                    2 => {
                                        let mut core_results = core_counts(
                                            &chromosome.starts,
                                            &chromosome.ends,
                                            current_chrom_size,
                                            stepsize,
                                        );
                                        match output_type {
                                            "file" => {
                                                panic!("Writing to file not supported.")
                                            }
                                            "csv" => {
                                                panic!("Write to CSV. Not Implemented");
                                            }
                                            "bedGraph" => {
                                                let file_name = format!(
                                                    "{}{}_{}.{}",
                                                    bwfileheader, chrom_name, "core", output_type
                                                );

                                                let count_info: (Vec<u32>, Vec<u32>, Vec<u32>) =
                                                    compress_counts(
                                                        &mut core_results,
                                                        clamped_start_position_zero_pos(
                                                            primary_start.0,
                                                            0,
                                                        ),
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
                                                    bwfileheader, chrom_name, "core", output_type
                                                );
                                                write_to_wig_file(
                                                    &core_results.0,
                                                    file_name.clone(),
                                                    chrom_name.clone(),
                                                    clamped_start_position(primary_start.0, 0, 1), //starts are 1 based must be shifted by 1
                                                    stepsize,
                                                    current_chrom_size,
                                                );
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
                                                    clamped_start_position_zero_pos(
                                                        primary_start.0,
                                                        0,
                                                    ),
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
                                                    clamped_start_position_zero_pos(
                                                        primary_start.0,
                                                        0,
                                                    ),
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

            //let vec_strings = vec!["start", "core", "end"];
            //let vec_strings = vec!["start"];

            let bar = ProgressBar::new(vec_count_type.len() as u64);
            match output_type {
                "wig" | "bedGraph" => {
                    println!("Combining {} Files", output_type);

                    for location in vec_count_type.iter() {
                        bar.inc(1);
                        write_combined_files(
                            *location,
                            output_type,
                            bwfileheader,
                            &final_chromosomes,
                        );
                    }
                }
                "npy" => {
                    // populate hashmap for the npy meta data
                    for chromosome in final_chromosomes.iter() {
                        let chr_name = chromosome.chrom.clone();
                        let current_chrom_size =
                            *chrom_sizes.get(&chromosome.chrom).unwrap() as i32;
                        npy_meta_data_map.insert(
                            chr_name,
                            HashMap::from([
                                ("stepsize".to_string(), stepsize),
                                ("reported_chrom_size".to_string(), current_chrom_size),
                            ]),
                        );
                    }

                    for location in vec_count_type.iter() {
                        let temp_meta_file_name =
                            format!("{}{}.{}", bwfileheader, *location, "meta");

                        if let Ok(file) = File::open(&temp_meta_file_name) {
                            let reader = BufReader::new(file);

                            for line in reader.lines() {
                                let line = line.unwrap();
                                let parts: Vec<&str> = line.split_whitespace().collect();
                                if parts.len() >= 3 {
                                    let chrom = parts[1].split('=').nth(1).expect(
                                        "Processing npy metadata file: Missing chromosome in line",
                                    );
                                    let start_str = parts[2].split('=')
                                        .nth(1)
                                        .expect("Processing npy metadata file: Missing start position in line");
                                    let starting_position: i32 = start_str.parse().expect(
                                        "Processing npy metadata file: Invalid start position",
                                    );

                                    if let Some(current_chr_data) = npy_meta_data_map.get_mut(chrom)
                                    {
                                        current_chr_data.insert(
                                            (*location.to_string()).parse().unwrap(),
                                            starting_position,
                                        );
                                    }
                                }
                            }
                            // Remove the file after it is used.
                            let path = std::path::Path::new(&temp_meta_file_name);
                            let _ = remove_file(path).unwrap();
                        }
                    }
                    //write combined metadata as json
                    let json_string = serde_json::to_string_pretty(&npy_meta_data_map).unwrap();
                    let combined_npy_meta_file_path =
                        format!("{}{}.{}", bwfileheader, "npy_meta", "json");
                    let mut file = File::create(combined_npy_meta_file_path).unwrap();
                    file.write_all(json_string.as_bytes()).unwrap();
                }
                _ => {}
            }
            bar.finish();

            match og_output_type {
                "bw" | "bigWig" => {
                    println!("Writing bigWig files");
                    if zoom != 1 {
                        println!(
                            "Only zoom level 1 is supported at this time, zoom level supplied {}",
                            zoom
                        );
                    }
                    let zoom = 1; //overwrite zoom
                    write_bw_files(bwfileheader, chromsizerefpath, num_threads, zoom);
                }

                _ => {}
            }
        }
        //BAM REQUIRES DIFFERENT WORKFLOW
        Ok(FileType::BAM) => {
            if chromsizerefpath == filepath {
                panic!("Must provide a valid chrom.sizes file for processing bam files. Provided file: {}", chromsizerefpath);
            }

            let _ = process_bam(
                vec_count_type,
                filepath,
                bwfileheader,
                chrom_sizes,
                chromsizerefpath,
                num_threads,
                zoom,
                pool,
                smoothsize,
                stepsize,
                output_type,
                debug,
                bam_shift,
                bam_scale,
            );
        }

        _ => {
            panic!("Unknown File Type provided");
        }
    };

    println!("FINISHED");

    Ok(())
}

/// This is for bam workflows where bam is the input file.
/// Currently, supports bam -> bigwig (start, end, core) and bam -> bed (shifted core values only).
/// You must provide a .bai file alongside the bam file! Create one: `samtools index your_file.bam`
fn process_bam(
    mut vec_count_type: Vec<&str>,
    filepath: &str,
    bwfileheader: &str,
    chrom_sizes: HashMap<String, u32>,
    chrom_sizes_ref_path: &str,
    num_threads: i32,
    zoom: i32,
    pool: ThreadPool,
    smoothsize: i32,
    stepsize: i32,
    output_type: &str,
    debug: bool,
    bam_shift: bool,
    bam_scale: f32,
) -> Result<(), Box<dyn Error>> {
    println!("Begin bam processing workflow...");
    let fp_string = filepath.to_string();
    let chrom_sizes_ref_path_string = chrom_sizes_ref_path.to_string();

    let list_of_valid_chromosomes: Vec<String> = chrom_sizes.keys().cloned().collect(); //taken from chrom.sizes as source of truth
    let mut final_chromosomes: Vec<String> = Vec::with_capacity(list_of_valid_chromosomes.len());

    // pre-process chromosomes that are actually in the bam file BEFORE spawning threads.
    for chromosome in list_of_valid_chromosomes.iter() {
        let region = chromosome.parse().unwrap();
        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(filepath)
            .unwrap();
        let header = reader.read_header().unwrap();
        match reader.query(&header, &region).map(Box::new) {
            Err(..) => {
                if debug {
                    eprintln!("Region not found, skipping region {}", region); //TODO only print if a debug mode is set?
                }

                continue;
            }

            Ok(mut records) => {
                // TODO does this pre-processing make downstream error handling redundant? No, because the functions are public.
                let first_record_option = records.next();

                match first_record_option {
                    Some(Ok(..)) => final_chromosomes.push(chromosome.clone()), // Extract the record
                    Some(Err(err)) => {
                        // Handle the error no first record
                        if debug {
                            eprintln!(
                                "Error reading the first record for chrom: {} {:?} Skipping...",
                                chromosome, err
                            );
                        }
                    }
                    None => {
                        // Handle no records
                        if debug {
                            eprintln!("No records exist for chrom: {} Skipping...", chromosome);
                        }
                    }
                };
            }
        }
    }

    //let out_selection_vec: Vec<&str>;

    if !bam_shift {
        //do nothing, just keep user output selection for starts, ends, core
    } else {
        if vec_count_type.len() > 1 {
            println!("bam_shift defaults to true for bam processing, but more than one count_type was selected. Defaulting to shift workflow which will produce a single file count file.");
        }
        vec_count_type = vec!["shift"];
    }

    match output_type {
        // Must merge all individual CHRs bw files...
        "bw" => {
            // TODO Add progress bars...
            pool.install(|| {
                final_chromosomes
                    .par_iter()
                    .for_each(|chromosome_string: &String| {
                        let out_selection_vec = vec_count_type.clone();

                        //let out_selection_vec = vec![OutSelection::STARTS];

                        for selection in out_selection_vec.iter() {
                            match selection {
                                &"start" => {
                                    process_bw_in_threads(
                                        &chrom_sizes,
                                        chromosome_string,
                                        smoothsize,
                                        stepsize,
                                        num_threads,
                                        zoom,
                                        bwfileheader,
                                        &fp_string,
                                        &chrom_sizes_ref_path_string,
                                        "start",
                                        bam_shift,
                                        bam_scale,
                                    );
                                }
                                &"end" => {
                                    process_bw_in_threads(
                                        &chrom_sizes,
                                        chromosome_string,
                                        smoothsize,
                                        stepsize,
                                        num_threads,
                                        zoom,
                                        bwfileheader,
                                        &fp_string,
                                        &chrom_sizes_ref_path_string,
                                        "end",
                                        bam_shift,
                                        bam_scale,
                                    );
                                }
                                &"core" => {
                                    process_bw_in_threads(
                                        &chrom_sizes,
                                        chromosome_string,
                                        smoothsize,
                                        stepsize,
                                        num_threads,
                                        zoom,
                                        bwfileheader,
                                        &fp_string,
                                        &chrom_sizes_ref_path_string,
                                        "core",
                                        bam_shift,
                                        bam_scale,
                                    );
                                }
                                &"shift" => {
                                    process_bw_in_threads(
                                        &chrom_sizes,
                                        chromosome_string,
                                        smoothsize,
                                        stepsize,
                                        num_threads,
                                        zoom,
                                        bwfileheader,
                                        &fp_string,
                                        &chrom_sizes_ref_path_string,
                                        "shift",
                                        bam_shift,
                                        bam_scale,
                                    );
                                }
                                _ => {
                                    println!("Must specify start, end, or core.")
                                }
                            }
                        }
                    })
            });

            println!("Merging all bigwig files...");
            //let out_selection_vec = vec!["start", "end", "core"];
            //let out_selection_vec = vec!["start"];

            for selection in vec_count_type.iter() {
                let combined_bw_file_name =
                    format!("{}_{}.{}", bwfileheader, selection, output_type);

                let final_file_path = combined_bw_file_name.clone();

                let mut inputs: Vec<String> = Vec::new();

                for chrom in final_chromosomes.iter() {
                    let file_name =
                        format!("{}_{}_{}.{}", bwfileheader, chrom, selection, output_type);
                    let result = File::open(&file_name);
                    match result {
                        Ok(_) => {
                            // File exists, add it to the input list
                            inputs.push(file_name);
                        }
                        Err(error) => {
                            // Just pass for now, this could happen if there are chroms in the bam header but no .bw files were created for those chroms
                            eprintln!("Error opening file: {}", error);
                        }
                    }
                    //inputs.push(file_name);
                }

                let mut bigwigs: Vec<BigWigRead<ReopenableFile>> = vec![];

                let inputs_clone = inputs.clone();

                for input in inputs {
                    match BigWigRead::open_file(&input) {
                        Ok(bw) => bigwigs.push(bw),
                        Err(e) => {
                            eprintln!(
                                "Error when opening bigwig {}. Skipping due to error: {:?}",
                                input, e
                            );
                        }
                    }
                }

                let threshold = 0.0; // default
                let adjust = Some(0.0); // default
                let clip = Some(100000000.0); // arbitrary but large because we don't want to clip
                let (iter, chrom_map) = get_merged_vals(bigwigs, 10, threshold, adjust, clip)?;

                let outb = BigWigWrite::create_file(combined_bw_file_name, chrom_map)?;
                let runtime = if num_threads == 1 {
                    runtime::Builder::new_current_thread().build().unwrap()
                } else {
                    runtime::Builder::new_multi_thread()
                        .worker_threads(num_threads as usize)
                        .build()
                        .unwrap()
                };
                let all_values = ChromGroupReadImpl {
                    iter: Box::new(iter),
                };

                //println!("WRITING COMBINED BW FILE: {}", combined_bw_file_name.clone());
                // outb.write(all_values, runtime)?;

                match outb.write(all_values, runtime) {
                    Ok(_) => {
                        eprintln!("Successfully wrote file: {}", final_file_path);
                    }
                    Err(err) => {
                        eprintln!("Error writing to BigWig file: {}", err);
                        // Delete the partially written file
                        std::fs::remove_file(final_file_path).unwrap_or_else(|e| {
                            eprintln!("Error deleting file: {}", e);
                        });
                    }
                }

                // CLean up after writing merged bigwig
                for input in inputs_clone.iter() {
                    std::fs::remove_file(input).unwrap_or_else(|e| {
                        eprintln!("Error deleting file: {}", e);
                    });
                }
            }
        }

        "bed" => {
            pool.install(|| {
                final_chromosomes
                    .par_iter()
                    .for_each(|chromosome_string: &String| {
                        let out_selection_vec = vec_count_type.clone();
                        //let out_selection_vec = vec![OutSelection::STARTS];

                        for selection in out_selection_vec.iter() {
                            match selection {
                                &"start" => {
                                    println!(
                                        "Only shift output is implemented for bam to BED file. (bamshift must be set to true)"
                                    );
                                }
                                &"end" => {
                                    println!(
                                        "Only shift output is implemented for bam to BED file. (bamshift must be set to true)"
                                    );
                                }
                                &"core" => {
                                    println!(
                                        "Only shift output is implemented for bam to BED file. (bamshift must be set to true)"
                                    );
                                }
                                &"shift" => {
                                    process_bed_in_threads(
                                        chromosome_string,
                                        smoothsize,
                                        bwfileheader,
                                        &fp_string,
                                        "shift",
                                    );
                                }
                                _ => {
                                    println!("Must specify start, end, or core")
                                }
                            }
                        }
                    })
            });

            // Combine bed files
            let out_selection_vec = vec_count_type.clone();
            for location in out_selection_vec.iter() {
                // this is a work around since we need to make a String to Chrom
                // so that we can re-use write_combined_files
                // use vec of Strings to make vec of empty chrom structs
                let mut chromosome_vec: Vec<Chromosome> = Vec::new();
                for chrom_string in final_chromosomes.iter() {
                    let chrom_name = chrom_string.clone();

                    let chromosome = Chromosome {
                        chrom: chrom_name,
                        starts: vec![],
                        ends: vec![],
                    };
                    chromosome_vec.push(chromosome);
                }

                write_combined_files(*location, output_type, bwfileheader, &chromosome_vec);
            }
        }

        _ => {

            // todo combine files for non bw outputs
        }
    }

    Ok(())
}

/// This option is for outputting BAM counts to any other file type that is not BW
/// Currently this will use FIXED step counting while outputting to bw uses variable step counting
// fn output_bam_counts_non_bw(    chrom_sizes: &HashMap<String, u32>,
//                                 chromosome_string: &String,
//                                 smoothsize: i32,
//                                 stepsize: i32,
//                                 num_threads: i32,
//                                 zoom: i32,
//                                 bwfileheader: &str,
//                                 fp_String: &String,
//                                 chrom_sizes_ref_path_String: &String,
//                                 sel: &str,) {
//
//     let region = chromosome_string.parse().unwrap();
//     let mut reader = bam::io::indexed_reader::Builder::default()
//         .build_from_path(fp_String)
//         .unwrap();
//     let header = reader.read_header().unwrap();
//
//     let mut records = reader.query(&header, &region).map(Box::new).unwrap();
//
//
//     match sel {
//         "start" | "end" => {
//             println!("fixed_core_counts for bam to other file file type (not bw or BED) currently not implemented.");
//             // fixed_start_end_counts_bam(
//             //     &mut records,
//             //     current_chrom_size,
//             //     smoothsize,
//             //     stepsize,
//             //     output_type,
//             //     chromosome_string,
//             //     bwfileheader,
//             //     "end",
//             //     false,
//             // );
//         }
//
//         "core" => {
//             println!("fixed_core_counts for bam to other file file type (not bw) currently not implemented.");
//         }
//
//         _ => {eprintln!("improper selection: {}", sel)}
//     }
//
//
//
// }

/// Creates a Producer/Consumer workflow for reading bam sequences and outputting to Bed files across threads.
fn process_bed_in_threads(
    chromosome_string: &String,
    smoothsize: i32,
    bwfileheader: &str,
    fp_string: &String,
    sel: &str,
) {
    let (reader, writer) = os_pipe::pipe().unwrap();
    let write_fd = Arc::new(Mutex::new(writer));
    let read_fd = Arc::new(Mutex::new(reader));

    let smoothsize_cloned = smoothsize.clone();

    let chromosome_string_cloned = chromosome_string.clone();

    let file_name = format!("{}{}_{}", bwfileheader, chromosome_string, sel);

    let fpclone = fp_string.clone(); // we must clone this string here, not before, else we get lifetime issues.

    let producer_handle = thread::spawn(move || {
        let region = chromosome_string_cloned.parse().unwrap();
        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(fpclone)
            .unwrap();
        let header = reader.read_header().unwrap();

        let mut records = reader.query(&header, &region).map(Box::new).unwrap();

        match bam_to_bed_no_counts(
            &mut records,
            smoothsize_cloned,
            &chromosome_string_cloned,
            write_fd,
        ) {
            Ok(_) => {
                eprintln!("Processing successful for {}", chromosome_string_cloned);
            }
            Err(err) => {
                eprintln!("Error processing records: {:?}", err);
            }
        }
    });

    let consumer_handle = thread::spawn(move || {
        let mut file_lock = read_fd.lock().unwrap(); // Acquire lock for writing
        let reader = std::io::BufReader::new(&mut *file_lock);

        let file_path = PathBuf::from(file_name);
        let new_file_path = file_path.with_extension("bed");

        let new_file_path = new_file_path.to_str().unwrap();

        // Create a new file
        let mut writer = std::fs::File::create(new_file_path).unwrap();

        // Read data from the reader and write it to the file
        for line in reader.lines() {
            let line = line.unwrap();
            writeln!(&mut writer, "{}", line).unwrap();
        }
    });

    producer_handle.join().unwrap();
    consumer_handle.join().unwrap();
}

/// Creates a Producer/Consumer workflow for reading bam sequences and outputting to bigwig files across threads.
fn process_bw_in_threads(
    chrom_sizes: &HashMap<String, u32>,
    chromosome_string: &String,
    smoothsize: i32,
    stepsize: i32,
    num_threads: i32,
    zoom: i32,
    bwfileheader: &str,
    fp_string: &String,
    chrom_sizes_ref_path_string: &String,
    sel: &str,
    bam_shift: bool,
    bam_scale: f32,
) {
    let (reader, writer) = os_pipe::pipe().unwrap();
    let write_fd = Arc::new(Mutex::new(writer));
    let read_fd = Arc::new(Mutex::new(reader));

    let current_chrom_size = *chrom_sizes.get(&chromosome_string.clone()).unwrap() as i32;

    let current_chrom_size_cloned = current_chrom_size.clone();
    let smoothsize_cloned = smoothsize.clone();
    let stepsize_cloned = stepsize.clone();
    let chromosome_string_cloned = chromosome_string.clone();
    let sel_clone = String::from(sel); // for some reason, even cloning a &str will lead to errors below when sel is moved to a new thread.

    let file_name = format!("{}_{}_{}", bwfileheader, chromosome_string, sel);

    let fpclone = fp_string.clone(); // we must clone this string here, not before, else we get lifetime issues.
    let chr_sz_ref_clone = chrom_sizes_ref_path_string.clone();

    let producer_handle = thread::spawn(move || {
        let region = chromosome_string_cloned.parse().unwrap();
        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(fpclone)
            .unwrap();
        let header = reader.read_header().unwrap();

        let records = reader.query(&header, &region).map(Box::new).unwrap();

        match determine_counting_func(
            records,
            current_chrom_size_cloned,
            smoothsize_cloned,
            stepsize_cloned,
            &chromosome_string_cloned,
            sel_clone.as_str(),
            write_fd,
            bam_shift,
            bam_scale,
        ) {
            Ok(_) => {
                //eprintln!("Processing successful for {}", chromosome_string_cloned);
            }
            Err(err) => {
                eprintln!("Error processing records: {:?}", err);
            }
        }
    });

    let consumer_handle = thread::spawn(move || {
        let mut file_lock = read_fd.lock().unwrap(); // Acquire lock for writing
        let mut reader = std::io::BufReader::new(&mut *file_lock);

        let file_path = PathBuf::from(file_name);
        let new_file_path = file_path.with_extension("bw");

        let new_file_path = new_file_path.to_str().unwrap();

        let mut outb = create_bw_writer(&*chr_sz_ref_clone, new_file_path, num_threads, zoom);

        let runtime = if num_threads == 1 {
            outb.options.channel_size = 0;
            runtime::Builder::new_current_thread().build().unwrap()
        } else {
            runtime::Builder::new_multi_thread()
                .worker_threads(num_threads as usize)
                .build()
                .unwrap()
        };
        let allow_out_of_order_chroms = !matches!(outb.options.input_sort_type, InputSortType::ALL);

        let vals =
            BedParserStreamingIterator::from_bedgraph_file(&mut reader, allow_out_of_order_chroms);
        match outb.write(vals, runtime) {
            Ok(_) => {
                //eprintln!("Successfully wrote file: {}", new_file_path);
            }
            Err(err) => {
                eprintln!("Error writing to BigWig file: {}", err);
                // Delete the partially written file
                std::fs::remove_file(new_file_path).unwrap_or_else(|e| {
                    eprintln!("Error deleting file: {}", e);
                });
            }
        }
    });

    producer_handle.join().unwrap();
    consumer_handle.join().unwrap();
}

/// This function determines if the starts/end counting function should be selected or the core counting function
/// Currently only variable step is supported, however, fixed_step has been written and can be added or replaced below if the user wishes.
/// Replacing the variable funcs with fixed step funcs will result in performance loss and greater processing times.
fn determine_counting_func(
    mut records: Box<Query<Reader<File>>>,
    current_chrom_size_cloned: i32,
    smoothsize_cloned: i32,
    stepsize_cloned: i32,
    chromosome_string_cloned: &String,
    sel_clone: &str,
    write_fd: Arc<Mutex<PipeWriter>>,
    bam_shift: bool,
    bam_scale: f32,
) -> Result<(), BAMRecordError> {
    //let bam_shift: bool = true; // This is to ensure a shifted position workflow is used when doing bams

    let count_result: Result<(), BAMRecordError> = match bam_shift {
        true => {
            match variable_shifted_bam_to_bw(
                &mut records,
                current_chrom_size_cloned,
                smoothsize_cloned,
                stepsize_cloned,
                &chromosome_string_cloned,
                sel_clone,
                write_fd,
                bam_scale,
            ) {
                Ok(_) => Ok(()),
                Err(err) => {
                    //eprintln!("Error processing records for {} {:?}", sel_clone,err);
                    Err(err)
                }
            }
        }
        false => {
            match sel_clone {
                "start" | "end" => {
                    match variable_start_end_counts_bam_to_bw(
                        &mut records,
                        current_chrom_size_cloned,
                        smoothsize_cloned,
                        stepsize_cloned,
                        &chromosome_string_cloned,
                        sel_clone,
                        write_fd,
                    ) {
                        Ok(_) => Ok(()),
                        Err(err) => {
                            //eprintln!("Error processing records for {} {:?}", sel_clone,err);
                            Err(err)
                        }
                    }
                }

                "core" => {
                    match variable_core_counts_bam_to_bw(
                        &mut records,
                        current_chrom_size_cloned,
                        stepsize_cloned,
                        &chromosome_string_cloned,
                        write_fd,
                    ) {
                        Ok(_) => {
                            //eprintln!("Processing successful for {}", chromosome_string_cloned);
                            Ok(())
                        }
                        Err(err) => {
                            //eprintln!("Error processing records for {}: {:?}", sel_clone,err);
                            Err(err)
                        }
                    }
                }

                &_ => {
                    eprintln!(
                        "Error processing records, improper selection: {}",
                        sel_clone
                    );
                    Err(BAMRecordError::IncorrectSel)
                }
            }
        }
    };

    count_result
}

/// Creates the bigwig writer struct for use with the BigTools crate
pub fn create_bw_writer(
    chrom_sizes_ref_path: &str,
    new_file_path: &str,
    num_threads: i32,
    zoom: i32,
) -> BigWigWrite<File> {
    //TODO do we need to force zooms? Related to https://github.com/jackh726/bigtools/issues/63
    let bedgraphargstruct = BedGraphToBigWigArgs {
        bedgraph: String::from("-"),
        chromsizes: chrom_sizes_ref_path.to_string(),
        output: new_file_path.to_string(),
        parallel: "auto".to_string(),
        single_pass: false,
        write_args: BBIWriteArgs {
            nthreads: num_threads as usize,
            nzooms: zoom as u32, // this does NOT force zooms
            zooms: None,         // this will force zooms
            uncompressed: false,
            sorted: "start".to_string(),
            block_size: 256,      //default
            items_per_slot: 1024, //default
            inmemory: false,
        },
    };
    let chrom_map: HashMap<String, u32> =
        BufReader::new(File::open(bedgraphargstruct.chromsizes).unwrap())
            .lines()
            .filter(|l| match l {
                Ok(s) => !s.is_empty(),
                _ => true,
            })
            .map(|l| {
                let words = l.expect("Split error");
                let mut split = words.split_whitespace();
                (
                    split.next().expect("Missing chrom").to_owned(),
                    split.next().expect("Missing size").parse::<u32>().unwrap(),
                )
            })
            .collect();

    let mut outb: BigWigWrite<File> =
        BigWigWrite::create_file(bedgraphargstruct.output, chrom_map).unwrap();
    outb.options.max_zooms = bedgraphargstruct.write_args.nzooms;
    outb.options.manual_zoom_sizes = bedgraphargstruct.write_args.zooms;
    outb.options.compress = !bedgraphargstruct.write_args.uncompressed;
    outb.options.input_sort_type = InputSortType::START;
    outb.options.block_size = bedgraphargstruct.write_args.block_size;
    outb.options.inmemory = bedgraphargstruct.write_args.inmemory;

    outb
}

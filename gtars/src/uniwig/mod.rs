use clap::ArgMatches;
use std::collections::HashMap;

use indicatif::ProgressBar;

use rayon::prelude::*;
use std::error::Error;
use std::fs::{create_dir_all, File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Write};

use crate::uniwig::counting::{
    core_counts, fixed_start_end_counts_bam, fixed_start_end_counts_bam_to_bw, start_end_counts,
};
use crate::uniwig::reading::{
    get_seq_reads_bam, read_bam_header, read_bed_vec, read_chromosome_sizes, read_narrow_peak_vec,
};
use crate::uniwig::utils::{compress_counts, get_final_chromosomes};
use crate::uniwig::writing::{
    write_bw_files, write_combined_files, write_to_bed_graph_file, write_to_npy_file,
    write_to_wig_file,
};
use bigtools::utils::cli::bedgraphtobigwig::{bedgraphtobigwig, BedGraphToBigWigArgs};
use bigtools::utils::cli::BBIWriteArgs;
use noodles::bam;
use noodles::sam::alignment::Record;
use rayon::ThreadPool;
use std::ops::Deref;
use std::path::PathBuf;
use std::str::FromStr;
use bigtools::beddata::BedParserStreamingIterator;
use bigtools::{BigWigWrite, InputSortType};
use tokio::runtime;
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

#[derive(Debug)]
enum OutSelection {
    STARTS,
    ENDS,
    CORE,
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
    let fixed = true;

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
            let og_output_type = output_type; // need this later for conversion
            let mut output_type = output_type;

            if output_type == "bedgraph" || output_type == "bw" || output_type == "bigwig" {
                output_type = "bedGraph" // we must create bedgraphs first before creating bigwig files
            }

            let mut final_chromosomes = get_final_chromosomes(&ft, filepath, &chrom_sizes, score);

            let bar = ProgressBar::new(final_chromosomes.len() as u64);

            // Pool installs iterator
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
                                                    clamped_start_position(
                                                        primary_start.0,
                                                        smoothsize,
                                                    ),
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
                                                        clamped_start_position(
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
                                                    clamped_start_position(
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
                                                    clamped_start_position(
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
                                                        clamped_start_position(
                                                            primary_end.0,
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
                                                    ),
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
                                                    clamped_start_position(
                                                        primary_start.0,
                                                        smoothsize,
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
                                                        primary_start.0,
                                                        smoothsize,
                                                    ),
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
                                                    compress_counts(
                                                        &mut core_results,
                                                        primary_start.0,
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
                        write_combined_files(
                            *location,
                            output_type,
                            bwfileheader,
                            &final_chromosomes,
                        );
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
            if chromsizerefpath == filepath {
                panic!("Must provide a valid chrom.sizes file for processing bam files. Provided file: {}", chromsizerefpath);
            }

            // Read sequences in chunks, do counts, send to bigTools via streamer.
            // Check that bam is sorted? Can noodles do that ahead of time? Error if not sorted.
            // Check for associated .bai file, if it does not exist create it
            //print!("Writing to CLI");
            // let handle = &std::io::stdout();
            // let mut buf = BufWriter::new(handle);
            // for count in &count_result.0 {
            //     writeln!(buf, "{}", count)
            //         .expect("failed to write line");
            // }
            // buf.flush().unwrap();
            let _ = process_bam(
                filepath,
                bwfileheader,
                chrom_sizes,
                chromsizerefpath,
                num_threads,
                zoom,
                pool,
                smoothsize,
                stepsize,
                fixed,
                output_type,
            );
            // match og_output_type {
            //     "bw" | "bigWig" => {
            //         println!("Writing bigWig files");
            //
            //         process_bam(filepath, bwfileheader,chrom_sizes, num_threads, zoom, pool, smoothsize, stepsize, fixed)
            //     }
            //     &_ => Ok({})
            // }
        }

        _ => {
            panic!("Unknown File Type provided");
        }
    };

    println!("FINISHED");

    Ok(())
}

fn process_bam(
    filepath: &str,
    bwfileheader: &str,
    chrom_sizes: HashMap<String, u32>,
    chrom_sizes_ref_path: &str,
    num_threads: i32,
    zoom: i32,
    pool: ThreadPool,
    smoothsize: i32,
    stepsize: i32,
    fixed: bool,
    output_type: &str,
) -> Result<(), Box<dyn Error>> {
    println!("Begin Process bam");

    let mut reader = bam::io::indexed_reader::Builder::default().build_from_path(filepath)?;
    let header = reader.read_header()?;

    let list_of_valid_chromosomes: Vec<String> = chrom_sizes.keys().cloned().collect(); //taken from chrom.sizes as source of truth

    pool.install(|| {
        list_of_valid_chromosomes
            .par_iter()
            .for_each(|chromosome_string: &String| {
                let region = chromosome_string.parse().unwrap(); // can this be coordinate?
                let current_chrom_size =
                    *chrom_sizes.get(&chromosome_string.clone()).unwrap() as i32;

                // let out_selection_vec =
                //     vec![OutSelection::STARTS, OutSelection::ENDS, OutSelection::CORE];
                let out_selection_vec = vec![OutSelection::STARTS];

                for selection in out_selection_vec.iter() {
                    match selection {
                        OutSelection::STARTS => {
                            let mut reader = bam::io::indexed_reader::Builder::default()
                                .build_from_path(filepath)
                                .unwrap();
                            let header = reader.read_header().unwrap();

                            match reader.query(&header, &region).map(Box::new) {
                                Err(_) => {} //Do nothing. //println!("Region not found in bam file, skipping region {}", region),

                                Ok(mut records) => {
                                    // let first = records.next().unwrap();
                                    // let first_start= first.unwrap().alignment_start().unwrap().unwrap().get();
                                    // You could get the first value and shift setting up the file headers BEFORE the counting

                                    match output_type {
                                        "bw" => {
                                            let file_name = format!(
                                                "{}_{}_{}",
                                                chromosome_string, bwfileheader, "start"
                                            );
                                            let file_path = PathBuf::from(file_name);
                                            let new_file_path = file_path.with_extension("bw");
                                            let new_file_path = new_file_path.to_str().unwrap();

                                            //TODO remove local path
                                            let new_file_path = "/home/drc/Downloads/refactor_test_gtars/example.bw";

                                            let bedgraphargstruct = BedGraphToBigWigArgs {
                                                bedgraph: String::from("-"),
                                                chromsizes: chrom_sizes_ref_path.to_string(),
                                                output: new_file_path.to_string(),
                                                parallel: "auto".to_string(),
                                                single_pass: false,
                                                write_args: BBIWriteArgs {
                                                    nthreads: num_threads as usize,
                                                    nzooms: zoom as u32,
                                                    zooms: None,
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

                                            let mut outb = BigWigWrite::create_file(bedgraphargstruct.output, chrom_map).unwrap();
                                            outb.options.max_zooms = bedgraphargstruct.write_args.nzooms;
                                            outb.options.compress = !bedgraphargstruct.write_args.uncompressed;
                                            outb.options.input_sort_type = InputSortType::START;
                                            outb.options.block_size = bedgraphargstruct.write_args.block_size;
                                            outb.options.inmemory = bedgraphargstruct.write_args.inmemory;
                                            let runtime = if bedgraphargstruct.write_args.nthreads == 1 {
                                                outb.options.channel_size = 0;
                                                runtime::Builder::new_current_thread().build().unwrap()
                                            } else {
                                                runtime::Builder::new_multi_thread()
                                                    .worker_threads(bedgraphargstruct.write_args.nthreads)
                                                    .build()
                                                    .unwrap()
                                            };
                                            let allow_out_of_order_chroms = !matches!(outb.options.input_sort_type, InputSortType::ALL);

                                            let bedgraph_line = fixed_start_end_counts_bam_to_bw(
                                                &mut records,
                                                current_chrom_size,
                                                smoothsize,
                                                stepsize,
                                                chromosome_string,
                                                bwfileheader,
                                                "start",
                                                true,
                                            );
                                            println!("after_fixed_start");
                                            let vals = BedParserStreamingIterator::from_bedgraph_file(bedgraph_line, allow_out_of_order_chroms);
                                            outb.write(vals, runtime).unwrap();

                                        }
                                        _ => {
                                            fixed_start_end_counts_bam(
                                                &mut records,
                                                current_chrom_size,
                                                smoothsize,
                                                stepsize,
                                                output_type,
                                                chromosome_string,
                                                bwfileheader,
                                                "start",
                                                false,
                                            );
                                        }
                                    }
                                }
                            }
                        }
                        OutSelection::ENDS => {
                            let mut reader = bam::io::indexed_reader::Builder::default()
                                .build_from_path(filepath)
                                .unwrap();
                            let header = reader.read_header().unwrap();
                            match reader.query(&header, &region).map(Box::new) {
                                Err(_) => {} //Do nothing. //println!("Region not found in bam file, skipping region {}", region),

                                Ok(mut records) => {
                                    // fixed_start_end_counts_bam(
                                    //     &mut records,
                                    //     current_chrom_size,
                                    //     smoothsize,
                                    //     stepsize,
                                    //     output_type,
                                    //     chromosome_string,
                                    //     bwfileheader,
                                    //     "end",
                                    //     false,
                                    // );
                                }
                            }
                        }
                        OutSelection::CORE => {
                            let mut reader = bam::io::indexed_reader::Builder::default()
                                .build_from_path(filepath)
                                .unwrap();
                            let header = reader.read_header().unwrap();
                            match reader.query(&header, &region).map(Box::new) {
                                Err(_) => {} //Do nothing. //println!("Region not found in bam file, skipping region {}", region),

                                Ok(mut records) => {}
                            }
                        }
                    }
                }

                // let mut reader = bam::io::indexed_reader::Builder::default().build_from_path(filepath).unwrap();
                // let header = reader.read_header().unwrap();
                //
                // let region = chromosome_string.parse().unwrap(); // can this be coordinate?
                // let current_chrom_size = *chrom_sizes.get(&chromosome_string.clone()).unwrap() as i32;

                // match reader.query(&header, &region).map(Box::new){
                //     Err(_) =>{},//Do nothing. //println!("Region not found in bam file, skipping region {}", region),
                //
                //     Ok(mut records) => {
                //
                //             for selection in out_selection_vec.iter() {
                //
                //                 match selection {
                //
                //                     OutSelection::STARTS =>{
                //
                //                         match fixed {
                //
                //                             true => {
                //                                 println!("Counting starts");
                //                                 //todo matching output type here might be redundandt if we need to do it anyway later for file writing...
                //                                 // match output_type {
                //                                 //
                //                                 //     "wig" => {
                //                                 //         //DETERMINE HEADER
                //                                 //         // can't do this
                //                                 //         //let iter = records.copied().peekable();
                //                                 //
                //                                 //     }
                //                                 //
                //                                 //     _ =>{println!("Unknown output type");
                //                                 //
                //                                 //     }
                //                                 //
                //                                 //
                //                                 // }
                //                                 fixed_start_end_counts_bam(&mut records,current_chrom_size,smoothsize,stepsize, output_type, chromosome_string, bwfileheader, "start");
                //
                //                                 //fixed_start_end_counts_bam(&mut records,current_chrom_size,smoothsize,stepsize);
                //
                //                             }
                //                             _ => {println!("Variable step not implemented")}
                //
                //
                //                         }
                //
                //
                //
                //
                //
                //                     }
                //
                //                     OutSelection::ENDS =>{
                //                         //TODO
                //                         match fixed {
                //
                //                             true => {
                //                                 println!("Counting ends");
                //                                 fixed_start_end_counts_bam(&mut records,current_chrom_size,smoothsize,stepsize, output_type, chromosome_string, bwfileheader, "end");
                //                                 //println!("Variable step not implemented")
                //
                //                             }
                //                             _ => {println!("Variable step not implemented")}
                //
                //
                //                         }
                //
                //                     }
                //
                //                     OutSelection::CORE =>{
                //                         //TODO
                //                         match fixed {
                //
                //                             true => {
                //                                 //fixed_start_end_counts_bam(&mut records,current_chrom_size,smoothsize,stepsize);
                //                                 println!("CORE NOT IMPLEMENTED")
                //
                //                             }
                //                             _ => {println!("Variable step not implemented")}
                //
                //
                //                         }
                //
                //                     }
                //                     _ => panic!("Unexpected value: {:?}", selection), // Handle unexpected values
                //
                //
                //                 }
                //
                //             }
                //
                //     },
                //
                // }
            })
    });

    Ok(())
}

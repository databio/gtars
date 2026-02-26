pub mod counting;
pub mod reading;
pub mod utils;
pub mod writing;

use indicatif::ProgressBar;
use rayon::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::fs::{File, remove_file};
use std::io::{BufRead, BufReader, Write};

use self::counting::{
    BAMRecordError, bam_to_bed_no_counts, core_counts, start_end_counts,
    variable_core_counts_bam_to_bw, variable_shifted_bam_to_bw,
    variable_start_end_counts_bam_to_bw,
};
use self::reading::read_chromosome_sizes;
use self::utils::{
    Chromosome, clamped_start_position, clamped_start_position_zero_pos, compress_counts,
    get_final_chromosomes,
};
use self::writing::{
    write_bw_files, write_combined_files, write_to_bed_graph_file, write_to_npy_file,
    write_to_wig_file, write_to_wig_file_variable,
};
use bigtools::beddata::BedParserStreamingIterator;
use bigtools::utils::cli::BBIWriteArgs;
use bigtools::utils::cli::bedgraphtobigwig::BedGraphToBigWigArgs;
use bigtools::utils::cli::bigwigmerge::{ChromGroupReadImpl, get_merged_vals};
use bigtools::utils::reopen::ReopenableFile;
use bigtools::{BigWigRead, BigWigWrite, InputSortType};
use gtars_core::utils::FileType;
use noodles::bam;
use noodles::bam::io::reader::Query;
use noodles::bgzf::Reader;
use noodles::sam::alignment::Record as SamRecord;
use os_pipe::PipeWriter;
use rayon::ThreadPool;
use std::path::PathBuf;
use std::str::FromStr;
use std::sync::{Arc, Mutex};
use std::thread;
use tokio::runtime;

/// Main function
#[allow(clippy::too_many_arguments)]
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
    wigstep: &str,
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

                        let primary_start = chromosome.starts[0];
                        let primary_end = chromosome.ends[0];

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
                                                if wigstep == "variable" {
                                                    write_to_wig_file_variable(
                                                        &count_result.0,
                                                        file_name.clone(),
                                                        chrom_name.clone(),
                                                        clamped_start_position(
                                                            primary_start.0,
                                                            smoothsize,
                                                            0, // no shift needed - coordinates already 1-based from BED conversion
                                                        ),
                                                        stepsize,
                                                        current_chrom_size,
                                                    );
                                                } else {
                                                    write_to_wig_file(
                                                        &count_result.0,
                                                        file_name.clone(),
                                                        chrom_name.clone(),
                                                        clamped_start_position(
                                                            primary_start.0,
                                                            smoothsize,
                                                            0, // no shift needed - coordinates already 1-based from BED conversion
                                                        ),
                                                        stepsize,
                                                        current_chrom_size,
                                                    );
                                                }
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
                                                    count_result.0,
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
                                                    count_result.0,
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
                                                if wigstep == "variable" {
                                                    write_to_wig_file_variable(
                                                        &count_result.0,
                                                        file_name.clone(),
                                                        chrom_name.clone(),
                                                        clamped_start_position(
                                                            primary_end.0,
                                                            smoothsize,
                                                            0,
                                                        ),
                                                        stepsize,
                                                        current_chrom_size,
                                                    );
                                                } else {
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
                                            }

                                            "npy" => {
                                                let file_name = format!(
                                                    "{}{}_{}.{}",
                                                    bwfileheader, chrom_name, "end", output_type
                                                );
                                                write_to_npy_file(
                                                    count_result.0,
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
                                                    count_result.0,
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
                                                if wigstep == "variable" {
                                                    write_to_wig_file_variable(
                                                        &core_results.0,
                                                        file_name.clone(),
                                                        chrom_name.clone(),
                                                        clamped_start_position(primary_start.0, 0, 0), // no shift - already 1-based
                                                        stepsize,
                                                        current_chrom_size,
                                                    );
                                                } else {
                                                    write_to_wig_file(
                                                        &core_results.0,
                                                        file_name.clone(),
                                                        chrom_name.clone(),
                                                        clamped_start_position(primary_start.0, 0, 0), // no shift - already 1-based
                                                        stepsize,
                                                        current_chrom_size,
                                                    );
                                                }
                                            }
                                            "npy" => {
                                                let file_name = format!(
                                                    "{}{}_{}.{}",
                                                    bwfileheader, chrom_name, "core", output_type
                                                );
                                                write_to_npy_file(
                                                    core_results.0,
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
                                                    core_results.0,
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
                            location,
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
                            remove_file(path).unwrap();
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
#[allow(clippy::too_many_arguments)]
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
        //TODO if no .bai file exists, the below line will fail and won't properly tell you WHY it failed, issue #57
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
            println!(
                "bam_shift defaults to true for bam processing, but more than one count_type was selected. Defaulting to shift workflow which will produce a single file count file."
            );
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
                            match *selection {
                                "start" => {
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
                                "end" => {
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
                                "core" => {
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
                                "shift" => {
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
                            match *selection {
                                "start" => {
                                    println!(
                                        "Only shift output is implemented for bam to BED file. (bamshift must be set to true)"
                                    );
                                }
                                "end" => {
                                    println!(
                                        "Only shift output is implemented for bam to BED file. (bamshift must be set to true)"
                                    );
                                }
                                "core" => {
                                    println!(
                                        "Only shift output is implemented for bam to BED file. (bamshift must be set to true)"
                                    );
                                }
                                "shift" => {
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

                write_combined_files(location, output_type, bwfileheader, &chromosome_vec);
            }
        }

        "wig" | "bedgraph" => {
            // Process BAM reads for wig/bedGraph output.
            // Collects BAM alignment positions into memory and uses the same
            // start_end_counts algorithm as the BED path, so output is equivalent
            // for reads at the same genomic positions.
            let norm_output_type = if output_type == "bedgraph" { "bedGraph" } else { output_type };

            for chromosome_string in final_chromosomes.iter() {
                let current_chrom_size =
                    *chrom_sizes.get(chromosome_string).unwrap() as i32;
                let region = chromosome_string.parse().unwrap();

                for selection in vec_count_type.iter() {
                    let mut reader = bam::io::indexed_reader::Builder::default()
                        .build_from_path(filepath)
                        .unwrap();
                    let header = reader.read_header().unwrap();
                    let mut records = reader.query(&header, &region).map(Box::new).unwrap();

                    // Collect positions from BAM records (score=1 per read)
                    let mut positions: Vec<(i32, i32)> = Vec::new();
                    for record in records.by_ref() {
                        let record = record.unwrap();
                        let pos: i32 = match *selection {
                            "start" => record.alignment_start().unwrap().unwrap().get() as i32,
                            "end" => SamRecord::alignment_end(&record).unwrap().unwrap().get() as i32,
                            "core" => {
                                eprintln!("Core counts for BAM non-BW output not yet implemented. Skipping.");
                                break;
                            }
                            _ => break,
                        };
                        positions.push((pos, 1));
                    }

                    if positions.is_empty() || *selection == "core" {
                        continue;
                    }

                    // Use same counting algorithm as BED path
                    let mut count_result =
                        start_end_counts(&positions, current_chrom_size, smoothsize, stepsize);
                    let primary_start = positions[0];

                    match norm_output_type {
                        "wig" => {
                            let file_name = format!(
                                "{}{}_{}.{}",
                                bwfileheader, chromosome_string, selection, norm_output_type
                            );
                            write_to_wig_file(
                                &count_result.0,
                                file_name,
                                chromosome_string.clone(),
                                clamped_start_position(primary_start.0, smoothsize, 0),
                                stepsize,
                                current_chrom_size,
                            );
                        }
                        "bedGraph" => {
                            let file_name = format!(
                                "{}{}_{}.{}",
                                bwfileheader, chromosome_string, selection, norm_output_type
                            );
                            let count_info = compress_counts(
                                &mut count_result,
                                clamped_start_position_zero_pos(primary_start.0, smoothsize),
                            );
                            write_to_bed_graph_file(
                                &count_info,
                                file_name,
                                chromosome_string.clone(),
                                current_chrom_size,
                            );
                        }
                        _ => {}
                    }
                }
            }

            // Build chromosome vec for write_combined_files
            let chromosome_vec: Vec<Chromosome> = final_chromosomes
                .iter()
                .map(|chrom_name| Chromosome {
                    chrom: chrom_name.clone(),
                    starts: vec![],
                    ends: vec![],
                })
                .collect();

            let norm_output_type = if output_type == "bedgraph" { "bedGraph" } else { output_type };
            for location in vec_count_type.iter() {
                if *location != "core" {
                    write_combined_files(
                        location,
                        norm_output_type,
                        bwfileheader,
                        &chromosome_vec,
                    );
                }
            }
        }

        _ => {

            // todo combine files for non bw outputs
        }
    }

    Ok(())
}

#[allow(clippy::empty_line_after_doc_comments)]
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
#[allow(clippy::ptr_arg)]
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

    let smoothsize_cloned = smoothsize;

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

        let file_path_with_ext = format!("{}.bed", file_name);
        let file_path = PathBuf::from(file_path_with_ext);

        let new_file_path = file_path.to_str().unwrap();

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
#[allow(clippy::too_many_arguments)]
#[allow(clippy::ptr_arg)]
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

    let current_chrom_size_cloned = current_chrom_size;
    let smoothsize_cloned = smoothsize;
    let stepsize_cloned = stepsize;
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

        let mut outb = create_bw_writer(&chr_sz_ref_clone, new_file_path, num_threads, zoom);

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
#[allow(clippy::too_many_arguments)]
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
                chromosome_string_cloned,
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
                        chromosome_string_cloned,
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
                        chromosome_string_cloned,
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

#[cfg(test)]
mod tests {

    use rstest::{fixture, rstest};
    use std::fs;
    use std::fs::File;
    use std::fs::read_dir;
    use std::io::{BufRead, BufReader, Read};
    use std::path::{Path, PathBuf};

    use super::{Chromosome, uniwig_main};
    use gtars_core::utils::parse_bedlike_file;

    use super::counting::{core_counts, start_end_counts};
    use super::reading::{
        create_chrom_vec_default_score, create_chrom_vec_scores, read_bam_header,
        read_chromosome_sizes,
    };

    use super::utils::npy_to_wig;

    //use super::utils::npy_to_wig;
    use super::writing::write_bw_files;

    // use gtars::bbcache::client::BBClient;

    //FIXTURES
    #[fixture]
    fn path_to_data() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/data")
    }

    #[fixture]
    fn path_to_bed_file() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/data/tokenizers/peaks.bed")
    }

    #[fixture]
    fn path_to_sorted_small_bed_file() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/data/test_sorted_small.bed")
    }

    #[fixture]
    fn path_to_small_bam_file() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/data/test_chr22_small.bam")
        //"/home/drc/Downloads/bam files for rust test/test1_sort_dedup.bam"
    }

    #[fixture]
    fn path_to_dummy_bam_file() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/data/dummy.bam")
    }

    #[fixture]
    fn path_to_chrom_sizes_file() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/hg38.chrom.sizes")
    }

    #[fixture]
    fn path_to_bed_file_gzipped() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/data/tokenizers/peaks.bed.gz")
    }

    #[fixture]
    fn path_to_dummy_bed_file() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/data/dummy.bed")
    }

    #[fixture]
    fn path_to_dummy_chromsizes() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/data/dummy.chrom.sizes")
    }

    #[fixture]
    fn path_to_dummy_narrowpeak() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/data/dummy.narrowPeak")
    }

    #[fixture]
    fn path_to_start_wig_output() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/data/out/_start.wig")
    }

    #[fixture]
    fn path_to_core_wig_output() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/data/out/_core.wig")
    }

    #[fixture]
    fn path_to_start_bedgraph_output() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/data/out/_start.bedGraph")
    }

    #[fixture]
    fn path_to_core_bedgraph_output() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/data/out/_core.bedGraph")
    }

    #[fixture]
    fn path_to_bed_gz_from_bb() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/data/6b2e163a1d4319d99bd465c6c78a9741.bed.gz")
    }

    #[fixture]
    fn bbid() -> &'static str {
        "6b2e163a1d4319d99bd465c6c78a9741"
    }

    #[fixture]
    fn bsid() -> &'static str {
        "gse127562"
    }

    #[fixture]
    fn path_to_bedset() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/data/bedset")
    }

    // UNIWIG TESTS
    #[rstest]
    fn test_uniwig_parsed_bed_file(path_to_bed_file: PathBuf) {
        let path = &path_to_bed_file;
        let file = File::open(path).unwrap();

        let mut reader = BufReader::new(file);
        let first_line = reader.by_ref().lines().next().unwrap().expect("expect");
        println!("{:?}", first_line);

        let result = parse_bedlike_file(&first_line);

        if let Some((ctg, st, en)) = result {
            println!("ctg: {}", ctg);
            println!("st: {}", st);
            println!("en: {}", en);
            assert_eq!(st, 7915738);
        } else {
            panic!("Failed to parse BED record");
        }
    }

    #[rstest]
    fn test_create_chrom_vec_default_score(
        path_to_bed_file: PathBuf,
        path_to_bed_file_gzipped: PathBuf,
    ) {
        let result1 = create_chrom_vec_default_score(&path_to_bed_file.to_string_lossy());
        assert_eq!(result1.len(), 20);

        let result2 = create_chrom_vec_default_score(&path_to_bed_file_gzipped.to_string_lossy());
        assert_eq!(result2.len(), 20);
    }

    #[rstest]
    fn test_create_chrom_vec_scores() {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let path_to_narrow_peak = PathBuf::from(path_to_crate)
            .parent()
            .unwrap()
            .join("tests/data/dummy.narrowPeak");
        let result1 = create_chrom_vec_scores(&path_to_narrow_peak.to_string_lossy());
        assert_eq!(result1.len(), 1);

        let path_to_narrow_peak_gzipped = PathBuf::from(path_to_crate)
            .parent()
            .unwrap()
            .join("tests/data/dummy.narrowPeak.gz");

        let result2 = create_chrom_vec_scores(&path_to_narrow_peak_gzipped.to_string_lossy());
        assert_eq!(result2.len(), 1);
    }

    #[rstest]
    fn test_read_scored_core_counts() {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let path_to_narrow_peak = PathBuf::from(path_to_crate)
            .parent()
            .unwrap()
            .join("tests/data/dummy.narrowPeak");
        let chrom_sizes = read_chromosome_sizes(&path_to_narrow_peak.to_string_lossy()).unwrap();
        let narrow_peak_vec: Vec<Chromosome> =
            create_chrom_vec_scores(&path_to_narrow_peak.to_string_lossy());
        let stepsize = 1;

        for chromosome in narrow_peak_vec.iter() {
            let current_chrom_size = *chrom_sizes.get(&chromosome.chrom).unwrap() as i32;
            let _result = core_counts(
                &chromosome.starts,
                &chromosome.ends,
                current_chrom_size,
                stepsize,
            );
        }
    }

    #[rstest]
    fn test_read_scored_starts_counts() {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let path_to_narrow_peak = PathBuf::from(path_to_crate)
            .parent()
            .unwrap()
            .join("tests/data/dummy.narrowPeak");
        let chrom_sizes = read_chromosome_sizes(&path_to_narrow_peak.to_string_lossy()).unwrap();
        let narrow_peak_vec: Vec<Chromosome> =
            create_chrom_vec_scores(&path_to_narrow_peak.to_string_lossy());
        let stepsize = 1;
        let smooth_size = 1;

        for chromosome in narrow_peak_vec.iter() {
            let current_chrom_size = *chrom_sizes.get(&chromosome.chrom).unwrap() as i32;
            let _result = start_end_counts(
                &chromosome.starts,
                current_chrom_size,
                smooth_size,
                stepsize,
            );
        }
    }

    #[rstest]
    fn test_read_bed_vec_length(path_to_sorted_small_bed_file: PathBuf) {
        let chromosomes: Vec<Chromosome> =
            create_chrom_vec_default_score(&path_to_sorted_small_bed_file.to_string_lossy());
        let num_chromosomes = chromosomes.len();

        assert_eq!(num_chromosomes, 5);
    }

    #[rstest]
    fn test_read_bam_header(path_to_small_bam_file: PathBuf) {
        let chromosomes: Vec<Chromosome> =
            read_bam_header(&path_to_small_bam_file.to_string_lossy());
        let num_chromosomes = chromosomes.len();
        println!("Number of chroms: {}", num_chromosomes);
        assert_eq!(num_chromosomes, 1);
    }

    #[rstest]
    fn test_process_bam(
        path_to_small_bam_file: PathBuf,
    ) -> Result<(), Box<dyn std::error::Error + 'static>> {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let chromsizerefpath: String = format!("{}{}", path_to_crate, "/../tests/hg38.chrom.sizes");
        let chromsizerefpath = chromsizerefpath.as_str();
        let combinedbedpath = &path_to_small_bam_file.to_string_lossy();

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let bwfileheader_path = path.into_os_string().into_string().unwrap();
        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 1;
        let output_type = "bw";
        let filetype = "bam";
        let num_threads = 2;
        let score = false;
        let stepsize = 1;
        let zoom = 0;

        let vec_count_type = vec!["start", "end", "core"];

        uniwig_main(
            vec_count_type,
            smoothsize,
            combinedbedpath,
            chromsizerefpath,
            bwfileheader,
            output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
            "fixed",
        )
        .expect("Uniwig main failed!");

        Ok(())
    }

    #[rstest]
    fn test_process_bam_to_bed(
        path_to_small_bam_file: PathBuf,
    ) -> Result<(), Box<dyn std::error::Error + 'static>> {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let chromsizerefpath: String = format!("{}{}", path_to_crate, "/../tests/hg38.chrom.sizes");
        let chromsizerefpath = chromsizerefpath.as_str();
        let combinedbedpath = &path_to_small_bam_file.to_string_lossy();

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let bwfileheader_path = path.into_os_string().into_string().unwrap();
        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 1;
        let output_type = "bed";
        let filetype = "bam";
        let num_threads = 2;
        let score = false;
        let stepsize = 1;
        let zoom = 0;
        let vec_count_type = vec!["start", "end", "core"];

        uniwig_main(
            vec_count_type,
            smoothsize,
            combinedbedpath,
            chromsizerefpath,
            bwfileheader,
            output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
            "fixed",
        )
        .expect("Uniwig main failed!");

        Ok(())
    }

    #[rstest]
    fn test_run_uniwig_main_wig_type() -> Result<(), Box<dyn std::error::Error + 'static>> {
        // This test uses the bed file to determine chromsizes for speed
        let path_to_crate = env!("CARGO_MANIFEST_DIR");

        let tempbedpath = format!("{}{}", path_to_crate, "/../tests/data/test5.bed");
        println!("{}", tempbedpath);
        let combinedbedpath = tempbedpath.as_str();

        let chromsizerefpath = combinedbedpath;

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let bwfileheader_path = path.into_os_string().into_string().unwrap();
        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 5;
        let output_type = "wig";
        let filetype = "bed";
        let num_threads = 6;
        let score = false;
        let stepsize = 1;
        let zoom = 0;
        let vec_count_type = vec!["start", "end", "core"];

        uniwig_main(
            vec_count_type,
            smoothsize,
            combinedbedpath,
            chromsizerefpath,
            bwfileheader,
            output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
            "fixed",
        )
        .expect("Uniwig main failed!");

        Ok(())
    }

    #[rstest]
    fn test_run_uniwig_main_npy_type() -> Result<(), Box<dyn std::error::Error + 'static>> {
        // This test uses the bed file to determine chromsizes for speed
        let path_to_crate = env!("CARGO_MANIFEST_DIR");

        let tempbedpath = format!("{}{}", path_to_crate, "/../tests/data/test5.bed");
        let combinedbedpath = tempbedpath.as_str();

        let chromsizerefpath = combinedbedpath;

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let bwfileheader_path = path.into_os_string().into_string().unwrap();
        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 2;
        let output_type = "npy";
        let filetype = "bed";
        let num_threads = 6;
        let score = false;
        let stepsize = 1;
        let zoom = 0;
        let vec_count_type = vec!["start", "end", "core"];

        uniwig_main(
            vec_count_type,
            smoothsize,
            combinedbedpath,
            chromsizerefpath,
            bwfileheader,
            output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
            "fixed",
        )
        .expect("Uniwig main failed!");
        Ok(())
    }

    #[rstest]
    fn test_run_uniwig_main_directory_type() -> Result<(), Box<dyn std::error::Error + 'static>> {
        // This test uses the bed file to determine chromsizes for speed
        let path_to_crate = env!("CARGO_MANIFEST_DIR");

        let tempbedpath = PathBuf::from(path_to_crate)
            .parent()
            .unwrap()
            .join("tests/data/dir_of_files/dir_beds/");
        let combinedbedpath = tempbedpath.to_string_lossy();

        //let chromsizerefpath = combinedbedpath;

        let chromsizerefpath = PathBuf::from(path_to_crate)
            .parent()
            .unwrap()
            .join("tests/data/dir_of_files/dummy.chrom.sizes");
        let chromsizerefpath = chromsizerefpath.to_string_lossy();
        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let bwfileheader_path = path.into_os_string().into_string().unwrap();
        let bwfileheader = bwfileheader_path.as_str();

        //let bwfileheader = "/home/drc/Downloads/gtars_uniwig_30june2025/output/";

        let smoothsize: i32 = 2;
        let output_type = "wig";
        let filetype = "bed";
        let num_threads = 6;
        let score = false;
        let stepsize = 1;
        let zoom = 0;
        let vec_count_type = vec!["start", "end", "core"];

        uniwig_main(
            vec_count_type,
            smoothsize,
            &combinedbedpath,
            &chromsizerefpath,
            bwfileheader,
            output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
            "fixed",
        )
        .expect("Uniwig main failed!");
        Ok(())
    }

    #[rstest]
    fn test_run_uniwig_main_directory_narrowpeaks_type()
    -> Result<(), Box<dyn std::error::Error + 'static>> {
        // This test uses the bed file to determine chromsizes for speed
        let path_to_crate = env!("CARGO_MANIFEST_DIR");

        let tempbedpath = format!(
            "{}{}",
            path_to_crate, "/../tests/data/dir_of_files/dir_narrowpeaks/"
        );
        let combinedbedpath = tempbedpath.as_str();

        //let chromsizerefpath = combinedbedpath;

        let chromsizerefpath = format!(
            "{}{}",
            path_to_crate, "/../tests/data/dir_of_files/dummy.chrom.sizes"
        );
        let chromsizerefpath = chromsizerefpath.as_str();
        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let bwfileheader_path = path.into_os_string().into_string().unwrap();
        let bwfileheader = bwfileheader_path.as_str();

        //let bwfileheader = "/home/drc/Downloads/gtars_uniwig_30june2025/output/";

        let smoothsize: i32 = 2;
        let output_type = "wig";
        let filetype = "narrowpeak";
        let num_threads = 6;
        let score = true;
        let stepsize = 1;
        let zoom = 0;
        let vec_count_type = vec!["start", "end", "core"];

        uniwig_main(
            vec_count_type,
            smoothsize,
            combinedbedpath,
            chromsizerefpath,
            bwfileheader,
            output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
            "fixed",
        )
        .expect("Uniwig main failed!");
        Ok(())
    }

    #[rstest]
    fn test_reading_chrom_sizes() {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");

        // Read from sizes file
        let chromsizerefpath: String = format!("{}{}", path_to_crate, "/../tests/hg38.chrom.sizes");
        let chrom_sizes = read_chromosome_sizes(chromsizerefpath.as_str()).unwrap();
        let chrom_name = String::from("chr13");
        let current_chrom_size = chrom_sizes[&chrom_name.clone()] as i32;
        assert_eq!(current_chrom_size, 114364328);

        // Read from BED file
        let tempbedpath = format!("{}{}", path_to_crate, "/../tests/data/test5.bed");
        let combinedbedpath = tempbedpath.as_str();
        let chrom_sizes = read_chromosome_sizes(combinedbedpath).unwrap();
        let chrom_name = String::from("chr1");
        let current_chrom_size = chrom_sizes[&chrom_name.clone()] as i32;
        assert_eq!(current_chrom_size, 32);
    }

    #[rstest]
    fn test_uniwig_mismatched_chrom_sizes(_path_to_bed_file: PathBuf) {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");

        // Read from sizes file
        let chromsizerefpath: String = format!("{}{}", path_to_crate, "/../tests/hg38.chrom.sizes");

        // Read from BED file that contains chromosomes not in size file
        let tempbedpath = format!(
            "{}{}",
            path_to_crate, "/../tests/data/test_unknown_chrom.bed"
        );
        let combinedbedpath = tempbedpath.as_str();
        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let bwfileheader_path = path.into_os_string().into_string().unwrap();
        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 5;
        let output_type = "npy";
        let filetype = "bed";
        let num_threads: i32 = 6;
        let score = false;
        let stepsize = 1;
        let zoom = 0;
        let vec_count_type = vec!["start", "end", "core"];

        let result = uniwig_main(
            vec_count_type,
            smoothsize,
            combinedbedpath,
            &chromsizerefpath,
            bwfileheader,
            output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
            "fixed",
        );

        assert!(result.is_ok());
    }

    #[rstest]
    fn test_uniwig_write_bw(_path_to_bed_file: PathBuf) {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");

        let original_bedgraph_path = format!("{}/../tests/data/test1.bedGraph", path_to_crate);
        let chrom_sizes_path = format!("{}/../tests/data/dummy.chrom.sizes", path_to_crate);

        let temp_dir = tempfile::tempdir().unwrap();
        let temp_bedgraph_path = temp_dir.path().join("test1.bedGraph");

        fs::copy(&original_bedgraph_path, &temp_bedgraph_path)
            .expect("Failed to copy .bedGraph file to temporary directory");

        let num_threads = 2;
        let zoom = 0;

        write_bw_files(
            temp_bedgraph_path.to_str().expect("Invalid temp path"), // Use the path in the temp directory
            chrom_sizes_path.as_str(),
            num_threads,
            zoom,
        );
    }

    #[rstest]
    fn test_uniwig_wiggle_output(
        _path_to_dummy_bed_file: PathBuf,
        _path_to_dummy_chromsizes: PathBuf,
        _path_to_start_wig_output: PathBuf,
        _path_to_core_wig_output: PathBuf,
    ) {
        let chromsizerefpath = &_path_to_dummy_chromsizes.to_string_lossy();
        let combinedbedpath = &_path_to_dummy_bed_file.to_string_lossy();
        let test_output_path = _path_to_start_wig_output.as_path();
        let core_test_output_path = _path_to_core_wig_output;

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let mut bwfileheader_path = path.into_os_string().into_string().unwrap();
        bwfileheader_path.push_str("/final/");

        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 1;
        let output_type = "wig";
        let filetype = "bed";
        let num_threads: i32 = 2;
        let score = false;
        let stepsize = 1;
        let zoom = 0;
        let vec_count_type = vec!["start", "end", "core"];

        let result = uniwig_main(
            vec_count_type,
            smoothsize,
            combinedbedpath,
            chromsizerefpath,
            bwfileheader,
            output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
            "fixed",
        );

        assert!(result.is_ok());

        // Test _start.wig output
        let path = PathBuf::from(&tempdir.path());
        let mut final_start_file_path = path.into_os_string().into_string().unwrap();
        final_start_file_path.push_str("/final/_start.wig");
        let final_start_file_path = final_start_file_path.as_str();

        let file1 = File::open(final_start_file_path).unwrap();
        let file2 = File::open(test_output_path).unwrap();

        let reader1 = BufReader::new(file1);
        let reader2 = BufReader::new(file2);

        let mut lines1 = reader1.lines();
        let mut lines2 = reader2.lines();

        loop {
            let line1 = lines1.next().transpose().unwrap();
            let line2 = lines2.next().transpose().unwrap();

            match (line1, line2) {
                (Some(line1), Some(line2)) => {
                    assert_eq!(line1, line2);
                }
                (None, None) => {
                    break; // Both files reached the end
                }
                _ => {
                    panic!("FILES ARE NOT EQUAL!!!")
                }
            }
        }

        // Test _core.wig output
        let path = PathBuf::from(&tempdir.path());
        let mut final_core_file_path = path.into_os_string().into_string().unwrap();
        final_core_file_path.push_str("/final/_core.wig");
        let final_core_file_path = final_core_file_path.as_str();

        let file1 = File::open(final_core_file_path).unwrap();
        let file2 = File::open(core_test_output_path).unwrap();

        let reader1 = BufReader::new(file1);
        let reader2 = BufReader::new(file2);

        let mut lines1 = reader1.lines();
        let mut lines2 = reader2.lines();

        loop {
            let line1 = lines1.next().transpose().unwrap();
            let line2 = lines2.next().transpose().unwrap();

            match (line1, line2) {
                (Some(line1), Some(line2)) => {
                    assert_eq!(line1, line2);
                }
                (None, None) => {
                    break; // Both files reached the end
                }
                _ => {
                    panic!("FILES ARE NOT EQUAL!!!")
                }
            }
        }
    }

    #[rstest]
    fn test_uniwig_bedgraph_output(
        _path_to_dummy_bed_file: PathBuf,
        _path_to_dummy_chromsizes: PathBuf,
        _path_to_start_bedgraph_output: PathBuf,
        _path_to_core_bedgraph_output: PathBuf,
    ) {
        let chromsizerefpath = &_path_to_dummy_chromsizes.to_string_lossy();
        let combinedbedpath = &_path_to_dummy_bed_file.to_string_lossy();
        let test_output_path = _path_to_start_bedgraph_output.as_path();
        let core_test_output_path = _path_to_core_bedgraph_output.as_path();

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let mut bwfileheader_path = path.into_os_string().into_string().unwrap();
        bwfileheader_path.push_str("/final/");

        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 1;
        let output_type = "bedgraph";
        let filetype = "bed";
        let num_threads: i32 = 2;
        let score = false;
        let stepsize = 1;
        let zoom = 0;
        let vec_count_type = vec!["start", "end", "core"];

        let result = uniwig_main(
            vec_count_type,
            smoothsize,
            combinedbedpath,
            chromsizerefpath,
            bwfileheader,
            output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
            "fixed",
        );

        assert!(result.is_ok());

        // Test _start.wig output
        let path = PathBuf::from(&tempdir.path());
        let mut final_start_file_path = path.into_os_string().into_string().unwrap();
        final_start_file_path.push_str("/final/_start.bedGraph");
        let final_start_file_path = final_start_file_path.as_str();

        let file1 = File::open(final_start_file_path).unwrap();
        let file2 = File::open(test_output_path).unwrap();

        let reader1 = BufReader::new(file1);
        let reader2 = BufReader::new(file2);

        let mut lines1 = reader1.lines();
        let mut lines2 = reader2.lines();

        loop {
            let line1 = lines1.next().transpose().unwrap();
            let line2 = lines2.next().transpose().unwrap();

            match (line1, line2) {
                (Some(line1), Some(line2)) => {
                    assert_eq!(line1, line2);
                }
                (None, None) => {
                    break; // Both files reached the end
                }
                _ => {
                    panic!("FILES ARE NOT EQUAL!!!")
                }
            }
        }

        // Test _core.wig output
        let path = PathBuf::from(&tempdir.path());
        let mut final_core_file_path = path.into_os_string().into_string().unwrap();
        final_core_file_path.push_str("/final/_core.bedGraph");
        let final_core_file_path = final_core_file_path.as_str();

        let file1 = File::open(final_core_file_path).unwrap();
        let file2 = File::open(core_test_output_path).unwrap();

        let reader1 = BufReader::new(file1);
        let reader2 = BufReader::new(file2);

        let mut lines1 = reader1.lines();
        let mut lines2 = reader2.lines();

        loop {
            let line1 = lines1.next().transpose().unwrap();
            let line2 = lines2.next().transpose().unwrap();

            match (line1, line2) {
                (Some(line1), Some(line2)) => {
                    assert_eq!(line1, line2);
                }
                (None, None) => {
                    break; // Both files reached the end
                }
                _ => {
                    panic!("FILES ARE NOT EQUAL!!!")
                }
            }
        }
    }

    #[rstest]
    fn test_uniwig_bam_wig_output(
        path_to_dummy_bam_file: PathBuf,
        path_to_dummy_chromsizes: PathBuf,
        path_to_start_wig_output: PathBuf,
    ) {
        let chromsizerefpath = &path_to_dummy_chromsizes.to_string_lossy();
        let combinedbedpath = &path_to_dummy_bam_file.to_string_lossy();
        let test_output_path = path_to_start_wig_output.as_path();

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        let mut bwfileheader_path = path.into_os_string().into_string().unwrap();
        bwfileheader_path.push_str("/final/");
        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 1;
        let output_type = "wig";
        let filetype = "bam";
        let num_threads: i32 = 2;
        let score = false;
        let stepsize = 1;
        let zoom = 0;
        let vec_count_type = vec!["start"];

        let result = uniwig_main(
            vec_count_type,
            smoothsize,
            combinedbedpath,
            chromsizerefpath,
            bwfileheader,
            output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            false, // bam_shift=false to produce start/end/core outputs
            1.0,
            "fixed",
        );
        assert!(result.is_ok());

        // Test _start.wig output
        let path = PathBuf::from(&tempdir.path());
        let mut final_start_file_path = path.into_os_string().into_string().unwrap();
        final_start_file_path.push_str("/final/_start.wig");
        let final_start_file_path = final_start_file_path.as_str();

        let file1 = File::open(final_start_file_path).unwrap();
        let file2 = File::open(test_output_path).unwrap();

        let reader1 = BufReader::new(file1);
        let reader2 = BufReader::new(file2);

        let mut lines1 = reader1.lines();
        let mut lines2 = reader2.lines();

        loop {
            let line1 = lines1.next().transpose().unwrap();
            let line2 = lines2.next().transpose().unwrap();

            match (line1, line2) {
                (Some(line1), Some(line2)) => {
                    assert_eq!(line1, line2);
                }
                (None, None) => {
                    break; // Both files reached the end
                }
                _ => {
                    panic!("FILES ARE NOT EQUAL!!!")
                }
            }
        }
    }

    #[rstest]
    fn test_uniwig_bam_bedgraph_output(
        path_to_dummy_bam_file: PathBuf,
        path_to_dummy_chromsizes: PathBuf,
        path_to_start_bedgraph_output: PathBuf,
    ) {
        let chromsizerefpath = &path_to_dummy_chromsizes.to_string_lossy();
        let combinedbedpath = &path_to_dummy_bam_file.to_string_lossy();
        let test_output_path = path_to_start_bedgraph_output.as_path();

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        let mut bwfileheader_path = path.into_os_string().into_string().unwrap();
        bwfileheader_path.push_str("/final/");
        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 1;
        let output_type = "bedgraph";
        let filetype = "bam";
        let num_threads: i32 = 2;
        let score = false;
        let stepsize = 1;
        let zoom = 0;
        let vec_count_type = vec!["start"];

        let result = uniwig_main(
            vec_count_type,
            smoothsize,
            combinedbedpath,
            chromsizerefpath,
            bwfileheader,
            output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            false, // bam_shift=false to produce start/end/core outputs
            1.0,
            "fixed",
        );
        assert!(result.is_ok());

        // Test _start.bedGraph output
        let path = PathBuf::from(&tempdir.path());
        let mut final_start_file_path = path.into_os_string().into_string().unwrap();
        final_start_file_path.push_str("/final/_start.bedGraph");
        let final_start_file_path = final_start_file_path.as_str();

        let file1 = File::open(final_start_file_path).unwrap();
        let file2 = File::open(test_output_path).unwrap();

        let reader1 = BufReader::new(file1);
        let reader2 = BufReader::new(file2);

        let mut lines1 = reader1.lines();
        let mut lines2 = reader2.lines();

        loop {
            let line1 = lines1.next().transpose().unwrap();
            let line2 = lines2.next().transpose().unwrap();

            match (line1, line2) {
                (Some(line1), Some(line2)) => {
                    assert_eq!(line1, line2);
                }
                (None, None) => {
                    break; // Both files reached the end
                }
                _ => {
                    panic!("FILES ARE NOT EQUAL!!!")
                }
            }
        }
    }

    #[rstest]
    fn test_process_narrowpeak(
        path_to_dummy_narrowpeak: PathBuf,
    ) -> Result<(), Box<dyn std::error::Error + 'static>> {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let chromsizerefpath = PathBuf::from(path_to_crate)
            .parent()
            .unwrap()
            .join("tests/data/dummy.chrom.sizes");
        let chromsizerefpath = chromsizerefpath.to_string_lossy();
        let combinedbedpath = &path_to_dummy_narrowpeak.to_string_lossy();

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let mut bwfileheader_path = path.into_os_string().into_string().unwrap();
        bwfileheader_path.push_str("/final/");
        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 1;
        let output_type = "bw";
        let filetype = "narrowpeak";
        let num_threads = 2;
        let score = true;
        let stepsize = 1;
        let zoom = 2;
        let vec_count_type = vec!["start", "end", "core"];

        uniwig_main(
            vec_count_type,
            smoothsize,
            combinedbedpath,
            &chromsizerefpath,
            bwfileheader,
            output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
            "fixed",
        )
        .expect("Uniwig main failed!");

        Ok(())
    }

    #[rstest]
    fn test_process_bed_to_bw(
        _path_to_dummy_bed_file: PathBuf,
    ) -> Result<(), Box<dyn std::error::Error + 'static>> {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let chromsizerefpath = PathBuf::from(path_to_crate)
            .parent()
            .unwrap()
            .join("tests/data/dummy.chrom.sizes");
        let chromsizerefpath = chromsizerefpath.to_string_lossy();
        let combinedbedpath = &_path_to_dummy_bed_file.to_string_lossy();

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        let mut bwfileheader_path = path.into_os_string().into_string().unwrap();
        bwfileheader_path.push_str("/final/");
        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 1;
        let output_type = "bw";
        let filetype = "bed";
        let num_threads = 2;
        let score = true;
        let stepsize = 1;
        let zoom = 1;
        let vec_count_type = vec!["start", "end", "core"];

        uniwig_main(
            vec_count_type,
            smoothsize,
            combinedbedpath,
            &chromsizerefpath,
            bwfileheader,
            output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
            "fixed",
        )
        .expect("Uniwig main failed!");

        Ok(())
    }

    #[rstest]
    fn test_npy_to_wig(
        _path_to_dummy_bed_file: PathBuf,
        _path_to_dummy_chromsizes: PathBuf,
    ) -> Result<(), Box<dyn std::error::Error + 'static>> {
        let chromsizerefpath = &_path_to_dummy_chromsizes.to_string_lossy();
        let combinedbedpath = &_path_to_dummy_bed_file.to_string_lossy();
        let tempdir = tempfile::tempdir()?; // use `?` for idiomatic error handling
        let path = PathBuf::from(tempdir.path());

        let smoothsize = 1;
        let wig_output_type = "wig";
        let npy_output_type = "npy";
        let filetype = "bed";
        let num_threads = 6;
        let score = false;
        let stepsize = 1;
        let zoom = 0;
        let vec_count_type = vec!["start", "end", "core"];

        // Generate npy output
        let npyfileheader_path = format!("{}/npyfinal/", path.display());
        let npyfileheader = npyfileheader_path.as_str();

        let _ = uniwig_main(
            vec_count_type.clone(),
            smoothsize,
            combinedbedpath,
            chromsizerefpath,
            npyfileheader,
            npy_output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
            "fixed",
        );

        // Generate wig output
        let wigfileheader_path = format!("{}/wigfinal/", path.display());
        let wigfileheader = wigfileheader_path.as_str();

        let _ = uniwig_main(
            vec_count_type.clone(),
            smoothsize,
            combinedbedpath,
            chromsizerefpath,
            wigfileheader,
            wig_output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
            "fixed",
        );

        // Run npy_to_wig
        let genwigfileheader_path = format!("{}/genwigfinal/", path.display());
        let genwigfileheader = genwigfileheader_path.as_str();

        let npy_header_path = Path::new(npyfileheader);
        let gen_wig_header_path = Path::new(genwigfileheader);
        let _ = npy_to_wig(npy_header_path, gen_wig_header_path);

        // Compare output directories
        let ref_wig_header_path = Path::new(wigfileheader);

        let mut files1: Vec<_> = read_dir(ref_wig_header_path)?
            .map(|entry| entry.unwrap().file_name().into_string().unwrap())
            .collect();
        let mut files2: Vec<_> = read_dir(gen_wig_header_path)?
            .map(|entry| entry.unwrap().file_name().into_string().unwrap())
            .collect();

        files1.sort();
        files2.sort();

        assert_eq!(files1, files2, "Directory file names differ");

        for file_name in files1 {
            let path1 = gen_wig_header_path.join(&file_name);
            let path2 = ref_wig_header_path.join(&file_name);

            let mut f1 = File::open(&path1)?;
            let mut f2 = File::open(&path2)?;

            let mut buf1 = Vec::new();
            let mut buf2 = Vec::new();

            f1.read_to_end(&mut buf1)?;
            f2.read_to_end(&mut buf2)?;

            assert_eq!(
                buf1,
                buf2,
                "File contents differ between:\n  {}\nand\n  {}",
                path1.display(),
                path2.display()
            );
        }

        Ok(())
    }

    /// Test smoothing near chromosome start (clamping behavior).
    ///
    /// A single read at position 3 with smoothsize=5 should produce coverage
    /// at positions 1-8 (8 positions). The smoothing window extends from
    /// position 3-5=-2 to 3+5=8, clamped to chromosome bounds: 1-8.
    ///
    /// Input: single read at BED 0-based position 2-3 (1-based position 3)
    /// Smoothsize: 5
    /// Expected window: 3±5 = -2 to 8, clamped to 1-8
    #[rstest]
    fn test_smoothing_clamp_at_chromosome_start() {
        use std::io::Write;

        let tempdir = tempfile::tempdir().unwrap();
        let temp_path = tempdir.path();

        // Create single-read BED file: chr1 2 3 (0-based, = position 3 in 1-based)
        let bed_path = temp_path.join("single.bed");
        let mut bed_file = File::create(&bed_path).unwrap();
        writeln!(bed_file, "chr1\t2\t3").unwrap();

        // Create chrom.sizes file
        let chrom_path = temp_path.join("chrom.sizes");
        let mut chrom_file = File::create(&chrom_path).unwrap();
        writeln!(chrom_file, "chr1\t20").unwrap();

        // Run uniwig
        let output_path = temp_path.join("output");
        std::fs::create_dir_all(&output_path).unwrap();
        let output_header = format!("{}/", output_path.display());

        uniwig_main(
            vec!["start"],
            5, // smoothsize
            bed_path.to_str().unwrap(),
            chrom_path.to_str().unwrap(),
            &output_header,
            "wig",
            "bed",
            1,
            false,
            1,
            0,
            false,
            true,
            1.0,
            "fixed",
        )
        .expect("uniwig_main failed");

        // Read output and count positions with value 1
        let wig_path = output_path.join("_start.wig");
        let content = std::fs::read_to_string(&wig_path).unwrap();
        let ones_count = content.lines().filter(|line| *line == "1").count();

        // Smoothing window 3±5 clamped to chromosome bounds = positions 1-8
        assert_eq!(
            ones_count, 8,
            "Expected 8 positions with value 1 (window clamped to 1-8), got {}",
            ones_count
        );
    }
}

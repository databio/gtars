use clap::builder::OsStr;
use clap::ArgMatches;
use flate2::read::GzDecoder;
use indicatif::ProgressBar;
use ndarray::Array;
use ndarray_npy::write_npy;
use noodles::bam;
use rayon::prelude::*;
use std::error::Error;
use std::fs::{create_dir_all, remove_file, File, OpenOptions};
use std::io;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::ops::Deref;
use std::path::Path;
use std::str::FromStr;
// use noodles::sam as sam;
//use bstr::BString;

pub mod cli;

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
pub struct Chromosome {
    chrom: String,
    starts: Vec<i32>,
    ends: Vec<i32>,
}
impl Clone for Chromosome {
    fn clone(&self) -> Self {
        Self {
            chrom: self.chrom.clone(),   // Clone the string
            starts: self.starts.clone(), // Clone the vector
            ends: self.ends.clone(),     // Clone the vector
        }
    }
}

// Chromosome representation for NarrowPeak Inputs
pub struct NarrowPeakChromosome {
    pub chrom: String,
    pub starts: Vec<(i32, i32)>, // first value of tuple is coordinate, 2nd is the narrowpeak score
    pub ends: Vec<(i32, i32)>,   // first value of tuple is coordinate, 2nd is the narrowpeak score
}
impl Clone for NarrowPeakChromosome {
    fn clone(&self) -> Self {
        Self {
            chrom: self.chrom.clone(),   // Clone the string
            starts: self.starts.clone(), // Clone the vector
            ends: self.ends.clone(),     // Clone the vector
        }
    }
}

/// Reads combined bed file from a given path.
/// Returns Vec of Chromosome struct
pub fn read_bed_vec(combinedbedpath: &str) -> Vec<Chromosome> {
    let path = Path::new(combinedbedpath);

    let file = File::open(path).unwrap();

    let is_gzipped = path.extension().unwrap_or(&OsStr::from("bed")) == "gz";

    // We must encapsulate in a box and use a dynamic Read trait so that either case could continue.
    let reader: Box<dyn Read> = match is_gzipped {
        true => Box::new(GzDecoder::new(file)),
        false => Box::new(file),
    };

    let reader = BufReader::new(reader);

    let mut chromosome = Chromosome {
        chrom: "".to_string(),
        starts: vec![],
        ends: vec![],
    };

    let mut chromosome_vec: Vec<Chromosome> = Vec::new();

    let mut chrom = String::new();

    for line in reader.lines() {
        //println!("Here is line{:?}", line);

        // Must use a 2nd let statement to appease the borrow-checker
        let line_string = line.unwrap();
        let s = line_string.as_str();

        let (parsed_chr, parsed_start, parsed_end) = parse_bed_file(s).unwrap();

        if chrom.is_empty() {
            // Initial chromosome
            chromosome.chrom = String::from(parsed_chr.trim());
            chrom = String::from(parsed_chr.trim());
            chromosome.starts.push(parsed_start);
            chromosome.ends.push(parsed_end);
            continue;
        }

        if String::from(parsed_chr.trim()) != chrom {
            // If the parsed chrom is not the same as the current, sort, and then push to vector
            // then reset chromosome struct using the newest parsed_chr
            chromosome.starts.sort_unstable();
            chromosome.ends.sort_unstable();

            chromosome_vec.push(chromosome.clone());

            chromosome.chrom = String::from(parsed_chr.trim());
            chrom = String::from(parsed_chr.trim());

            chromosome.starts = vec![];
            chromosome.ends = vec![]
        }

        chromosome.starts.push(parsed_start);
        chromosome.ends.push(parsed_end);
    }

    // Is this final sort and push actually necessary?
    chromosome.starts.sort_unstable();
    chromosome.ends.sort_unstable();
    chromosome_vec.push(chromosome.clone());

    println!("Reading Bed file complete.");

    chromosome_vec
}

pub fn read_narrow_peak_vec(combinedbedpath: &str) -> Vec<NarrowPeakChromosome> {
    let path = Path::new(combinedbedpath);

    let file = File::open(path).unwrap();

    let is_gzipped = path.extension().unwrap_or(&OsStr::from("narrowpeak")) == "gz";

    // We must encapsulate in a box and use a dynamic Read trait so that either case could continue.
    let reader: Box<dyn Read> = match is_gzipped {
        true => Box::new(GzDecoder::new(file)),
        false => Box::new(file),
    };

    let reader = BufReader::new(reader);

    let mut npchromosome = NarrowPeakChromosome {
        chrom: "".to_string(),
        starts: vec![],
        ends: vec![],
    };

    let mut chromosome_vec: Vec<NarrowPeakChromosome> = Vec::new();

    let mut chrom = String::new();

    for line in reader.lines() {
        //println!("Here is line{:?}", line);

        // Must use a 2nd let statement to appease the borrow-checker
        let line_string = line.unwrap();
        let s = line_string.as_str();

        let (parsed_chr, parsed_start, parsed_end, parsed_score) =
            parse_narrow_peak_file(s).unwrap();

        if chrom.is_empty() {
            // Initial chromosome
            npchromosome.chrom = String::from(parsed_chr.trim());
            chrom = String::from(parsed_chr.trim());
            npchromosome.starts.push((parsed_start, parsed_score));
            npchromosome.ends.push((parsed_end, parsed_score));
            continue;
        }

        if String::from(parsed_chr.trim()) != chrom {
            // If the parsed chrom is not the same as the current, sort, and then push to vector
            // then reset chromosome struct using the newest parsed_chr
            //npchromosome.starts.sort_unstable();
            //npchromosome.ends.sort_unstable();
            npchromosome.starts.sort_unstable_by(|a, b| a.0.cmp(&b.0));
            npchromosome.ends.sort_unstable_by(|a, b| a.0.cmp(&b.0));

            chromosome_vec.push(npchromosome.clone());

            npchromosome.chrom = String::from(parsed_chr.trim());
            chrom = String::from(parsed_chr.trim());

            npchromosome.starts = vec![];
            npchromosome.ends = vec![]
        }

        npchromosome.starts.push((parsed_start, parsed_score));
        npchromosome.ends.push((parsed_end, parsed_score));
    }

    // Is this final sort and push actually necessary?
    // npchromosome.starts.sort_unstable();
    // npchromosome.ends.sort_unstable();
    npchromosome.starts.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    npchromosome.ends.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    chromosome_vec.push(npchromosome.clone());

    println!("Reading narrowPeak file complete.");

    chromosome_vec
}
pub fn parse_narrow_peak_file(line: &str) -> Option<(String, i32, i32, i32)> {
    let mut fields = line.split('\t');
    // Get the first field which should be chromosome.
    let ctg = fields.next()?;
    // Parse 2nd and 3rd string as integers or return -1 if failure
    let st = fields
        .next()
        .and_then(|s| s.parse::<i32>().ok())
        .unwrap_or(-1);
    let en = fields
        .next()
        .and_then(|s| s.parse::<i32>().ok())
        .unwrap_or(-1);

    let _ = fields.next();

    let narrow_peak_score = fields
        .next()
        .and_then(|s| s.parse::<i32>().ok())
        .unwrap_or(-1);

    // Original code had a remainder of the line, r, but it does not appear to have been used
    // in any way

    Some((ctg.parse().unwrap(), st, en, narrow_peak_score))
}
/// Parses each line of given bed file into a contig (chromosome), starts and ends
pub fn parse_bed_file(line: &str) -> Option<(String, i32, i32)> {
    let mut fields = line.split('\t');
    // Get the first field which should be chromosome.
    let ctg = fields.next()?;
    // Parse 2nd and 3rd string as integers or return -1 if failure
    let st = fields
        .next()
        .and_then(|s| s.parse::<i32>().ok())
        .unwrap_or(-1);
    let en = fields
        .next()
        .and_then(|s| s.parse::<i32>().ok())
        .unwrap_or(-1);

    // Original code had a remainder of the line, r, but it does not appear to have been used
    // in any way

    Some((ctg.parse().unwrap(), st, en))
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

    uniwig_main(
        *smoothsize,
        filepath,
        chromsizerefpath.as_str(),
        bwfileheader,
        output_type,
        filetype,
        *num_threads,
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
) -> Result<(), Box<dyn Error>> {
    // Must create a Rayon thread pool in which to run our iterators
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads as usize)
        .build()
        .unwrap();

    // Determine File Type
    let ft = FileType::from_str(filetype.to_lowercase().as_str());

    let stepsize = 1;
    // Set up output file names

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

    let chromosomes: Vec<Chromosome> = match ft {
        Ok(FileType::BED) => read_bed_vec(filepath),
        Ok(FileType::BAM) => read_bam_header(filepath),
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
    let bar = ProgressBar::new(final_chromosomes.len() as u64);

    // Pool installs iterator
    pool.install(|| {
        final_chromosomes
            .par_iter()
            .for_each(|chromosome: &Chromosome| {
                // Need these for setting wiggle header
                bar.inc(1);
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
                                let count_result = match ft {
                                    Ok(FileType::BED) => smooth_fixed_start_end_wiggle(
                                        &chromosome.starts,
                                        current_chrom_size,
                                        smoothsize,
                                        stepsize,
                                    ),
                                    Ok(FileType::BAM) => smooth_fixed_start_end_wiggle_bam(
                                        &chromosome.starts,
                                        current_chrom_size,
                                        smoothsize,
                                        stepsize,
                                    ),
                                    _ => smooth_fixed_start_end_wiggle(
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
                                            clamped_start_position(primary_start, smoothsize),
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
                                            clamped_start_position(primary_start, smoothsize),
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
                                            clamped_start_position(primary_start, smoothsize),
                                            stepsize,
                                            meta_data_file_names[0].clone(),
                                        );
                                    }
                                }
                            }
                            1 => {
                                let count_result = match ft {
                                    Ok(FileType::BED) => smooth_fixed_start_end_wiggle(
                                        &chromosome.ends,
                                        current_chrom_size,
                                        smoothsize,
                                        stepsize,
                                    ),
                                    Ok(FileType::BAM) => smooth_fixed_start_end_wiggle_bam(
                                        &chromosome.ends,
                                        current_chrom_size,
                                        smoothsize,
                                        stepsize,
                                    ),
                                    _ => smooth_fixed_start_end_wiggle(
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
                                    "wig" => {
                                        let file_name = format!(
                                            "{}{}_{}.{}",
                                            bwfileheader, chrom_name, "end", output_type
                                        );
                                        write_to_wig_file(
                                            &count_result.0,
                                            file_name.clone(),
                                            chrom_name.clone(),
                                            clamped_start_position(primary_end, smoothsize),
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
                                            clamped_start_position(primary_start, smoothsize),
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
                                            clamped_start_position(primary_start, smoothsize),
                                            stepsize,
                                            meta_data_file_names[1].clone(),
                                        );
                                    }
                                }
                            }
                            2 => {
                                let core_results = match ft {
                                    Ok(FileType::BED) => fixed_core_wiggle(
                                        &chromosome.starts,
                                        &chromosome.ends,
                                        current_chrom_size,
                                        stepsize,
                                    ),
                                    Ok(FileType::BAM) => fixed_core_wiggle_bam(
                                        &chromosome.starts,
                                        &chromosome.ends,
                                        current_chrom_size,
                                        stepsize,
                                    ),
                                    _ => fixed_core_wiggle(
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
                                    "wig" => {
                                        let file_name = format!(
                                            "{}{}_{}.{}",
                                            bwfileheader, chrom_name, "core", output_type
                                        );
                                        write_to_wig_file(
                                            &core_results.0,
                                            file_name.clone(),
                                            chrom_name.clone(),
                                            primary_start,
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
                                            primary_start,
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
                                            primary_start,
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

    let bar = ProgressBar::new(vec_strings.len() as u64);
    match output_type {
        "wig" => {
            println!("Combining Wig Files");

            for location in vec_strings.iter() {
                bar.inc(1);
                write_combined_wig_files(*location, output_type, bwfileheader, &final_chromosomes);
            }
        }
        _ => {}
    }
    bar.finish();
    println!("FINISHED");

    Ok(())
}

fn fixed_core_wiggle_bam(
    _p0: &Vec<i32>,
    _p1: &Vec<i32>,
    _p2: i32,
    _p3: i32,
) -> (Vec<u32>, Vec<i32>) {
    println!("smooth_fixed_start_end_wiggle_bam");

    let v_coordinate_positions: Vec<i32> = Vec::new(); // these are the final coordinates after any adjustments
    let v_coord_counts: Vec<u32> = Vec::new(); // u8 stores 0:255 This may be insufficient. u16 max is 65535

    (v_coord_counts, v_coordinate_positions)
}

fn smooth_fixed_start_end_wiggle_bam(
    _p0: &Vec<i32>,
    _p1: i32,
    _p2: i32,
    _p3: i32,
) -> (Vec<u32>, Vec<i32>) {
    println!("smooth_fixed_start_end_wiggle_bam");

    let v_coordinate_positions: Vec<i32> = Vec::new(); // these are the final coordinates after any adjustments
    let v_coord_counts: Vec<u32> = Vec::new(); // u8 stores 0:255 This may be insufficient. u16 max is 65535

    (v_coord_counts, v_coordinate_positions)
}

pub fn read_bam_header(filepath: &str) -> Vec<Chromosome> {
    // BAM and SAM format specification https://samtools.github.io/hts-specs/SAMv1.pdf
    println!("READ BAM HEADER PLACE HOLDER");

    let mut reader = bam::io::reader::Builder.build_from_path(filepath).unwrap();
    let header = reader.read_header();

    let references = header.unwrap();
    let references = references.reference_sequences();

    let mut chromosome = Chromosome {
        chrom: "".to_string(),
        starts: vec![],
        ends: vec![],
    };
    let mut chromosome_vec: Vec<Chromosome> = Vec::new();

    for ref_key in references {
        let chrom_name_vec = ref_key.0.deref().clone();
        let chrom_name = String::from_utf8((*chrom_name_vec).to_owned()).unwrap();

        //For later
        // use bstr::BString;
        //
        // let s = BString::from("Hello, world!");
        chromosome.chrom = chrom_name;
        chromosome.starts.push(0); //default values for now, less important for bam
        chromosome.ends.push(0); //default values for now, less important for bam
        chromosome_vec.push(chromosome.clone());
    }

    chromosome_vec
}

fn write_to_npy_file(
    counts: &Vec<u32>,
    filename: String,
    chromname: String,
    start_position: i32,
    stepsize: i32,
    metafilename: String,
) {
    // For future reference `&Vec<u32>` is a SLICE and thus we must use the `to_vec` function below when creating an array
    // https://users.rust-lang.org/t/why-does-std-to-vec-exist/45893/9

    // Write the NumPy Files
    let arr = Array::from_vec(counts.to_vec());
    let _ = write_npy(filename, &arr);

    // Write to the metadata file.
    // Note: there should be a single metadata file for starts, ends and core

    let path = std::path::Path::new(&metafilename).parent().unwrap();
    let _ = create_dir_all(path);

    let mut file = OpenOptions::new()
        .create(true) // Create the file if it doesn't exist
        .append(true) // Append data to the existing file if it does exist
        .open(metafilename)
        .unwrap();

    // The original wiggle file header. This can be anything we wish it to be. Currently space delimited.
    let mut wig_header = "fixedStep chrom=".to_string()
        + chromname.as_str()
        + " start="
        + start_position.to_string().as_str()
        + " step="
        + stepsize.to_string().as_str();
    wig_header.push_str("\n");
    file.write_all(wig_header.as_ref()).unwrap();
}

fn write_combined_wig_files(
    location: &str,
    output_type: &str,
    bwfileheader: &str,
    chromosomes: &Vec<Chromosome>,
) {
    let combined_wig_file_name = format!("{}_{}.{}", bwfileheader, location, output_type);
    let path = std::path::Path::new(&combined_wig_file_name)
        .parent()
        .unwrap();
    let _ = create_dir_all(path);

    let mut combined_file = OpenOptions::new()
        .create(true) // Create the file if it doesn't exist
        .append(true) // Append data to the existing file if it does exist
        .open(combined_wig_file_name)
        .unwrap();

    let mut inputs: Vec<String> = Vec::new();

    for chrom in chromosomes.iter() {
        let file_name = format!(
            "{}{}_{}.{}",
            bwfileheader, chrom.chrom, location, output_type
        );
        inputs.push(file_name);
    }

    for input_file in inputs {
        // copy single file to the combined file
        let mut input = File::open(&input_file).unwrap();
        io::copy(&mut input, &mut combined_file).expect("cannot copy file!!");

        // Remove the file after it is combined.
        let path = std::path::Path::new(&input_file);
        let _ = remove_file(path).unwrap();
    }
}

#[allow(unused_variables)]
fn write_to_wig_file(
    counts: &Vec<u32>,
    filename: String,
    chromname: String,
    start_position: i32,
    stepsize: i32,
) {
    let path = std::path::Path::new(&filename).parent().unwrap();
    let _ = create_dir_all(path);

    let mut file = OpenOptions::new()
        .create(true) // Create the file if it doesn't exist
        .append(true) // Append data to the existing file if it does exist
        .open(filename)
        .unwrap();

    let wig_header = "fixedStep chrom=".to_string()
        + chromname.as_str()
        + " start="
        + start_position.to_string().as_str()
        + " step="
        + stepsize.to_string().as_str();
    file.write_all(wig_header.as_ref()).unwrap();
    file.write_all(b"\n").unwrap();

    let mut buf = BufWriter::new(file);

    for count in counts.iter() {
        writeln!(&mut buf, "{}", count).unwrap();
    }
    buf.flush().unwrap();
}

/// Reads chromosome size file from path and returns chromosome sizes hash map
pub fn read_chromosome_sizes(
    chrom_size_path: &str,
) -> Result<std::collections::HashMap<String, u32>, Box<dyn Error>> {
    let chrom_size_file = File::open(Path::new(chrom_size_path))?;

    // Get FIle extension
    let path = Path::new(chrom_size_path);
    let extension = path.extension().and_then(|ext| ext.to_str());

    let mut chrom_sizes = std::collections::HashMap::new();
    let reader = BufReader::new(chrom_size_file);

    match extension {
        //TODO what if the user provides a zipped bed file or a zipped narrowPeak and not a .sizes file? This will probably fail.
        Some("bed") => {
            // Read BED file
            //println!("Processing BED file: {}", chrom_size_path);
            for line in reader.lines() {
                let line = line?; // Propagate the potential error
                let mut iter = line.split('\t');
                let chrom_name = iter.next().unwrap().to_owned();
                let _ = iter.next().unwrap();
                let size_str = iter.next().unwrap();
                let size = size_str.parse::<u32>()?;

                chrom_sizes.insert(chrom_name, size);
            }
        }
        Some("narrowPeak") => {
            // TODO refactor the above case and this case to simply call a function
            // Read narrowPeak
            for line in reader.lines() {
                let line = line?; // Propagate the potential error
                let mut iter = line.split('\t');
                let chrom_name = iter.next().unwrap().to_owned();
                let _ = iter.next().unwrap();
                let size_str = iter.next().unwrap();
                let size = size_str.parse::<u32>()?;

                chrom_sizes.insert(chrom_name, size);
            }
        }
        Some("sizes") => {
            // Read sizes file
            // Note this may lead to slower performance as uniwig will pad the remaining chromosome with zeros
            // this is a remainder from legacy uniwig for creating wiggle files and bigwigs
            // It could potentially be removed in future versions if deemed unnecessary.
            for line in reader.lines() {
                let line = line?; // Propagate the potential error
                let mut iter = line.split_whitespace();
                let chrom_name = iter.next().unwrap().to_owned();
                let size_str = iter.next().unwrap();
                let size = size_str.parse::<u32>()?;

                chrom_sizes.insert(chrom_name, size);
            }
        }
        _ => {
            panic!("Unsupported file type: {}", chrom_size_path);
        }
    }

    Ok(chrom_sizes)
}

/// This function is a more direct port of smoothFixedStartEndBW from uniwig written in CPP.
/// It allows the user to accumulate reads of either starts or ends.
/// Counts occur between a start coordinate (cutSite) and an end site (endSite) where the endsite is determined based on
/// the level of smoothing.
/// counts are reported over a stepsize (with a default of stepsize = 1).
/// Unlike the original function, it does not write to disk in chunks. it simply returns a vector of accumulated reads.
#[allow(unused_variables)]
pub fn smooth_fixed_start_end_wiggle(
    starts_vector: &Vec<i32>,
    chrom_size: i32,
    smoothsize: i32,
    stepsize: i32,
) -> (Vec<u32>, Vec<i32>) {
    //println!("BEGIN smooth_Fixed_Start_End_Wiggle");

    let vin_iter = starts_vector.iter();

    let mut v_coordinate_positions: Vec<i32> = Vec::new(); // these are the final coordinates after any adjustments
    let mut v_coord_counts: Vec<u32> = Vec::new(); // u8 stores 0:255 This may be insufficient. u16 max is 65535

    let mut coordinate_position = 1;

    let mut count: u32 = 0;

    let mut coordinate_value: i32;
    let mut prev_coordinate_value = 0;

    let mut adjusted_start_site: i32;
    let mut current_end_site: i32;

    let mut collected_end_sites: Vec<i32> = Vec::new();

    //println!("DEBUG: START SITE BEFORE ADJUSTMENT -> {}",starts_vector[0].clone());

    adjusted_start_site = starts_vector[0].clone(); // get first coordinate position
    adjusted_start_site = adjusted_start_site - smoothsize; // adjust based on smoothing
                                                            //println!("DEBUG: START SITE AFTER ADJUSTMENT -> {}",adjusted_start_site.clone());
                                                            //Check endsite generation
    current_end_site = adjusted_start_site + 1 + smoothsize * 2;

    //println!("DEBUG: INITIAL ENDSITE -> {}", current_end_site.clone());

    if adjusted_start_site < 1 {
        adjusted_start_site = 1;
    }

    //println!("DEBUG: SKIPPING UNTIL COORDINATE_POSITION < ADJUSTEDSTARTSITE -> {}  {}", coordinate_position.clone(), adjusted_start_site.clone());
    while coordinate_position < adjusted_start_site {
        // Just skip until we reach the initial adjusted start position
        // Note that this function will not return 0s at locations before the initial start site
        coordinate_position = coordinate_position + stepsize;
    }

    //println!("DEBUG: SKIPPING UNTIL COORDINATE_POSITION < ADJUSTEDSTARTSITE -> {}  {}", coordinate_position.clone(), adjusted_start_site.clone());

    //prev_coordinate_value = adjusted_start_site;

    for coord in vin_iter.skip(0) {
        //println!("DEBUG: BEGIN COORDINATE ITERATION");
        coordinate_value = *coord;
        //println!("DEBUG: COORDINATE VALUE {}", coordinate_value.clone());
        adjusted_start_site = coordinate_value - smoothsize;
        count += 1;

        if adjusted_start_site < 1 {
            adjusted_start_site = 1;
        }

        //current_end_site = adjusted_start_site + 1 + smoothsize*2; //

        collected_end_sites.push(adjusted_start_site + 1 + smoothsize * 2);

        //println!("DEBUG: Coordinate Value: {}, Adjusted Start Site: {}, New Endsite: {} ", coordinate_value.clone(), adjusted_start_site.clone(), adjusted_start_site + 1 + smoothsize*2);

        if adjusted_start_site == prev_coordinate_value {
            continue;
        }

        while coordinate_position < adjusted_start_site {
            while current_end_site == coordinate_position {
                count = count - 1;

                if collected_end_sites.last() == None {
                    current_end_site = 0; // From original code. Double check this is the proper way.
                } else {
                    current_end_site = collected_end_sites.remove(0)
                }
            }

            if coordinate_position % stepsize == 0 {
                // Step size defaults to 1, so report every value
                v_coord_counts.push(count);
                v_coordinate_positions.push(coordinate_position);
                //println!("DEBUG: Reporting count: {} at position: {} for adjusted start site: {}",count, coordinate_position, adjusted_start_site);
            }

            //println!("DEBUG: Incrementing coordinate_position: {}  -> {}", coordinate_position,  coordinate_position +1);
            coordinate_position = coordinate_position + 1;
        }

        prev_coordinate_value = adjusted_start_site;
    }

    count = count + 1; // We must add 1 extra value here so that our calculation during the tail as we close out the end sites does not go negative.
                       // this is because the code above subtracts twice during the INITIAL end site closure. So we are missing one count and need to make it up else we go negative.
                       //

    while coordinate_position < chrom_size {
        // Apply an bound to push the final coordinates otherwise it will become truncated.

        while current_end_site == coordinate_position {
            count = count - 1;

            if collected_end_sites.last() == None {
                current_end_site = 0; // From original code. Double check this is the proper way.
            } else {
                current_end_site = collected_end_sites.remove(0)
            }
        }

        if coordinate_position % stepsize == 0 {
            // Step size defaults to 1, so report every value
            v_coord_counts.push(count);
            v_coordinate_positions.push(coordinate_position); // This is ONLY the starts
                                                              //println!("DEBUG: Reporting count: {} at start position: {} and end position: {}", count, coordinate_position, current_end_site);
        }

        //println!("DEBUG: Incrementing coordinate_position: {}  -> {}", coordinate_position,  coordinate_position +1);
        coordinate_position = coordinate_position + 1;
    }

    //println!("DEBUG: FINAL LENGTHS... Counts: {:?}  Positions: {:?}", v_coord_counts, v_coordinate_positions);
    (v_coord_counts, v_coordinate_positions)
}

/// This function is a more direct port of fixedCoreBW from uniwig written in CPP
/// It allows the user to accumulate reads across paired starts and ends.
/// Counts occur between a start coordinate (cutSite) and an end site (endSite) where the endsite is determined based on
/// the paired ends.
/// Counts are reported over a stepsize (with a default of stepsize = 1)
/// Unlike the original function, it does not write to disk in chunks. it simply returns a vector of accumulated reads.
#[allow(unused_variables)]
pub fn fixed_core_wiggle(
    starts_vector: &Vec<i32>,
    ends_vector: &Vec<i32>,
    chrom_size: i32,
    stepsize: i32,
) -> (Vec<u32>, Vec<i32>) {
    //println!("BEGIN Fixed_Core_Wiggle");

    //println!("STARTS VECTOR LENGTH: {}  END VECTORS LENGTH: {}", starts_vector.len().clone(), ends_vector.len().clone());

    let mut v_coordinate_positions: Vec<i32> = Vec::new(); // these are the final coordinates after any adjustments
    let mut v_coord_counts: Vec<u32> = Vec::new(); // u8 stores 0:255 This may be insufficient. u16 max is 65535

    let mut coordinate_position = 1;

    let mut count = 0;

    let mut coordinate_value: i32;
    let mut prev_coordinate_value = 0;

    let mut current_start_site: i32;
    let mut current_end_site: i32;

    let mut collected_end_sites: Vec<i32> = Vec::new();

    current_start_site = starts_vector[0].clone(); // get first coordinate position
    current_end_site = ends_vector[0];

    //Check endsite generation
    //current_end_site = adjusted_start_site + 1 + smoothsize*2;

    if current_start_site < 1 {
        current_start_site = 1;
    }

    while coordinate_position < current_start_site {
        // Just skip until we reach the initial adjusted start position
        // Note that this function will not return 0s at locations before the initial start site
        coordinate_position = coordinate_position + stepsize;
    }

    //prev_coordinate_value = current_start_site;

    for (index, coord) in starts_vector.iter().enumerate().skip(0) {
        coordinate_value = *coord;

        current_start_site = coordinate_value;

        count += 1;

        if current_start_site < 1 {
            current_start_site = 1;
        }

        let current_index = index;

        //current_end_site = ends_vector[current_index];

        collected_end_sites.push(ends_vector[current_index]);

        if current_start_site == prev_coordinate_value {
            continue;
        }

        while coordinate_position < current_start_site {
            while current_end_site == coordinate_position {
                count = count - 1;

                if collected_end_sites.last() == None {
                    current_end_site = 0; // From original code. Double check this is the proper way.
                } else {
                    current_end_site = collected_end_sites.remove(0)
                }
            }

            if coordinate_position % stepsize == 0 {
                // Step size defaults to 1, so report every value
                v_coord_counts.push(count);
                v_coordinate_positions.push(coordinate_position); // This is ONLY the starts
                                                                  //println!("DEBUG: Reporting count: {} at start position: {} and end position: ",count, coordinate_position);
            }

            //println!("DEBUG: Incrementing coordinate_position: {}  -> {}", coordinate_position,  coordinate_position +1);
            coordinate_position = coordinate_position + 1;
        }

        prev_coordinate_value = current_start_site;
    }

    count = count + 1; // We must add 1 extra value here so that our calculation during the tail as we close out the end sites does not go negative.
                       // this is because the code above subtracts twice during the INITIAL end site closure. So we are missing one count and need to make it up else we go negative.
                       //

    while coordinate_position < chrom_size {
        while current_end_site == coordinate_position {
            count = count - 1;

            if collected_end_sites.last() == None {
                current_end_site = 0; // From original code. Double check this is the proper way.
            } else {
                current_end_site = collected_end_sites.remove(0)
            }
        }

        if coordinate_position % stepsize == 0 {
            // Step size defaults to 1, so report every value
            v_coord_counts.push(count);
            v_coordinate_positions.push(coordinate_position); // This is ONLY the starts
                                                              //println!("DEBUG: Reporting count: {} at start position: {} and end position: ",count, coordinate_position);
        }

        //println!("DEBUG: Incrementing coordinate_position: {}  -> {}", coordinate_position,  coordinate_position +1);
        coordinate_position = coordinate_position + 1;
    }

    //println!("DEBUG: FINAL LENGTHS... Counts: {}  Positions: {}", v_coord_counts.len(), v_coordinate_positions.len());
    (v_coord_counts, v_coordinate_positions)
}


#[allow(unused_variables)]
pub fn smooth_fixed_start_end_narrow_peak(
    starts_vector: &Vec<(i32,i32)>,
    chrom_size: i32,
    smoothsize: i32,
    stepsize: i32,
) -> (Vec<u32>, Vec<i32>) {

    let mut v_coordinate_positions: Vec<i32> = Vec::new(); // these are the final coordinates after any adjustments
    let mut v_coord_counts: Vec<u32> = Vec::new(); // u8 stores 0:255 This may be insufficient. u16 max is 65535

    let mut coordinate_position = 1;

    let mut count = 0;

    let mut coordinate_value: (i32,i32);
    let mut prev_coordinate_value = 0;

    let mut adjusted_start_site: (i32,i32);
    let mut current_end_site: (i32,i32);

    let mut collected_end_sites: Vec<(i32,i32)> = Vec::new();

    adjusted_start_site = starts_vector[0].clone(); // get first coordinate position

    adjusted_start_site.0 = adjusted_start_site.0 - smoothsize; // adjust based on smoothing

    current_end_site = adjusted_start_site;
    current_end_site.0 = adjusted_start_site.0 + 1 + smoothsize * 2;

    if adjusted_start_site.0 < 1 {
        adjusted_start_site.0 = 1;
    }

    while coordinate_position < adjusted_start_site.0 {
        // Just skip until we reach the initial adjusted start position
        // Note that this function will not return 0s at locations before the initial start site
        coordinate_position = coordinate_position + stepsize;
    }
    // prev_coordinate_value = adjusted_start_site.0;

    for (index, coord) in starts_vector.iter().enumerate().skip(0) {
        coordinate_value = *coord;

        adjusted_start_site = coordinate_value;
        adjusted_start_site.0 = coordinate_value.0 - smoothsize;

        let current_score = adjusted_start_site.1;

        count += current_score;

        if adjusted_start_site.0 < 1 {
            adjusted_start_site.0 = 1;
        }

        let current_index = index;

        if current_index != 0{ // this is already added at the beginning of the functions
            current_end_site = adjusted_start_site;
            current_end_site.0 = adjusted_start_site.0 + 1 + smoothsize*2;
            collected_end_sites.push(current_end_site);
        }

        if adjusted_start_site.0 == prev_coordinate_value {
            continue;
        }

        while coordinate_position < adjusted_start_site.0 {

            while current_end_site.0 == coordinate_position {
                count = count - current_score;

                if collected_end_sites.last() == None {
                    current_end_site.0 = 0; // From original code. Double check this is the proper way.
                } else {
                    current_end_site = collected_end_sites.remove(0);
                }
            }

            if coordinate_position % stepsize == 0 {
                // Step size defaults to 1, so report every value
                v_coord_counts.push(count as u32);
                v_coordinate_positions.push(coordinate_position); // This is ONLY the starts

            }

            coordinate_position = coordinate_position + 1;
        }

        prev_coordinate_value = adjusted_start_site.0;
    }


    while coordinate_position < chrom_size {
        while current_end_site.0 == coordinate_position {
            let current_score = adjusted_start_site.1;

            count = count - current_score;

            if collected_end_sites.last() == None {
                current_end_site.0 = 0; // From original code. Double check this is the proper way.
            } else {
                current_end_site = collected_end_sites.remove(0)
            }
        }

        if coordinate_position % stepsize == 0 {
            // Step size defaults to 1, so report every value
            v_coord_counts.push(count as u32);
            v_coordinate_positions.push(coordinate_position); // This is ONLY the starts
        }

        coordinate_position = coordinate_position + 1;
    }

    (v_coord_counts, v_coordinate_positions)

}

//Counts based on NarrowPeak Scores
pub fn fixed_core_narrow_peak(
    starts_vector: &Vec<(i32,i32)>,
    ends_vector: &Vec<(i32,i32)>,
    chrom_size: i32,
    stepsize: i32,
) -> (Vec<u32>, Vec<i32>) {

    let mut v_coordinate_positions: Vec<i32> = Vec::new(); // these are the final coordinates after any adjustments
    let mut v_coord_counts: Vec<u32> = Vec::new(); // u8 stores 0:255 This may be insufficient. u16 max is 65535

    let mut coordinate_position = 1;

    let mut count = 0;

    let mut coordinate_value: (i32,i32);
    let mut prev_coordinate_value = 0;

    let mut current_start_site: (i32,i32);
    let mut current_end_site: (i32,i32);

    let mut collected_end_sites: Vec<(i32,i32)> = Vec::new();

    current_start_site = starts_vector[0].clone(); // get first coordinate position
    current_end_site = ends_vector[0].clone();


    if current_start_site.0 < 1 {
        current_start_site.0 = 1;
    }

    while coordinate_position < current_start_site.0 {
        // Just skip until we reach the initial adjusted start position
        // Note that this function will not return 0s at locations before the initial start site
        coordinate_position = coordinate_position + stepsize;
    }


    for (index, coord) in starts_vector.iter().enumerate().skip(0) {

        coordinate_value = *coord;

        current_start_site = coordinate_value;

        let current_score = current_start_site.1;

        count += current_score;

        if current_start_site.0 < 1 {
            current_start_site.0 = 1;
        }

        let current_index = index;

        if current_index != 0{ // this is already added at the beginning of the functions
            collected_end_sites.push(ends_vector[current_index]);
        }


        if current_start_site.0 == prev_coordinate_value {
            continue;
        }

        while coordinate_position < current_start_site.0 {
            while current_end_site.0 == coordinate_position {
                count = count - current_score;

                if collected_end_sites.last() == None {
                    current_end_site.0 = 0; // From original code. Double check this is the proper way.
                } else {
                    current_end_site = collected_end_sites.remove(0);
                }
            }

            if coordinate_position % stepsize == 0 {
                // Step size defaults to 1, so report every value
                v_coord_counts.push(count as u32);
                v_coordinate_positions.push(coordinate_position); // This is ONLY the starts
            }

            coordinate_position = coordinate_position + 1;
        }

        prev_coordinate_value = current_start_site.0;
    }


    while coordinate_position < chrom_size {
        while current_end_site.0 == coordinate_position {
            let current_score = current_start_site.1;

            count = count - current_score;

            if collected_end_sites.last() == None {
                current_end_site.0 = 0; // From original code. Double check this is the proper way.
            } else {
                current_end_site = collected_end_sites.remove(0)
            }
        }

        if coordinate_position % stepsize == 0 {
            // Step size defaults to 1, so report every value
            v_coord_counts.push(count as u32);
            v_coordinate_positions.push(coordinate_position); // This is ONLY the starts
        }
        coordinate_position = coordinate_position + 1;
    }

    (v_coord_counts, v_coordinate_positions)
}

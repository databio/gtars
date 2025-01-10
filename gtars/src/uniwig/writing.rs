use crate::uniwig::Chromosome;
use bigtools::utils::cli::bedgraphtobigwig::{bedgraphtobigwig, BedGraphToBigWigArgs};
use bigtools::utils::cli::BBIWriteArgs;
use indicatif::ProgressBar;
use ndarray::Array;
use ndarray_npy::write_npy;
use std::fs::{create_dir_all, remove_file, File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::{fs, io};

/// Write output to npy files
pub fn write_to_npy_file(
    counts: &[u32],
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
    wig_header.push('\n');
    file.write_all(wig_header.as_ref()).unwrap();
}

/// Write either combined bedGraph, wiggle files, and bed files
/// Requires a list of Chromosomes
pub fn write_combined_files(
    location: &str,
    output_type: &str,
    bwfileheader: &str,
    chromosomes: &[Chromosome], // TODO make this a vec of Strings instead? Since we only care about the names.
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

/// Write output  to a wiggle file
pub fn write_to_wig_file(
    counts: &[u32],
    filename: String,
    chromname: String,
    start_position: i32,
    stepsize: i32,
    chrom_size: i32,
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

    for count in counts.iter().take(chrom_size as usize) {
        // must set upper bound for wiggles based on reported chromsize, this is for downstream tool interoperability
        writeln!(&mut buf, "{}", count).unwrap();
    }
    buf.flush().unwrap();
}

/// Write output to bedgraph file
pub fn write_to_bed_graph_file(
    count_info: &(Vec<u32>, Vec<u32>, Vec<u32>),
    filename: String,
    chromname: String,
    _stepsize: i32,
) {
    let path = std::path::Path::new(&filename).parent().unwrap();
    let _ = create_dir_all(path);

    if count_info.0.len() != count_info.1.len() || count_info.0.len() != count_info.2.len() {
        panic!("count info vectors are not equal!")
    }

    let n_index = count_info.0.len();

    let file = OpenOptions::new()
        .create(true) // Create the file if it doesn't exist
        .append(true) // Append data to the existing file if it does exist
        .open(filename)
        .unwrap();

    let mut buf = BufWriter::new(file);

    for i in 0..n_index {
        writeln!(
            &mut buf,
            "{}\t{}\t{}\t{}",
            chromname, count_info.0[i], count_info.1[i], count_info.2[i]
        )
        .unwrap();
    }
    buf.flush().unwrap();
}

/// Converts uniwig generated bedGraphs to bigWig files
pub fn write_bw_files(location: &str, chrom_sizes: &str, num_threads: i32, zoom_level: i32) {
    //Collect all bedGraph files in the given location/directory
    let mut bed_graph_files = Vec::new();

    let mut location_path = location;

    if !location_path.ends_with("/") {
        let mut temp_path = Path::new(location_path);
        let parent_location_path = temp_path.parent().unwrap();
        location_path = parent_location_path.to_str().unwrap();
    }

    for entry in fs::read_dir(location_path).unwrap() {
        let entry = entry.unwrap();
        let path = entry.path();

        if path.is_file() {
            let extension = path.extension().unwrap();
            let extension = extension.to_str().unwrap().to_lowercase();
            let extension = extension.as_str();

            match extension {
                "bedgraph" => {
                    bed_graph_files.push(path.to_str().unwrap().to_string());
                }
                _ => {
                    continue;
                }
            }
        }
    }

    let bar = ProgressBar::new(bed_graph_files.len() as u64);
    for file in bed_graph_files.iter() {
        bar.inc(1);
        let file_path = PathBuf::from(file);
        let new_file_path = file_path.with_extension("bw");
        let new_file_path = new_file_path.to_str().unwrap();

        let current_arg_struct = BedGraphToBigWigArgs {
            bedgraph: file.to_string(),
            chromsizes: chrom_sizes.to_string(),
            output: new_file_path.to_string(),
            parallel: "auto".to_string(),
            single_pass: false,
            write_args: BBIWriteArgs {
                nthreads: num_threads as usize,
                nzooms: zoom_level as u32,
                zooms: None,
                uncompressed: false,
                sorted: "start".to_string(),
                block_size: 256,      //default
                items_per_slot: 1024, //default
                inmemory: false,
            },
        };

        let _ = bedgraphtobigwig(current_arg_struct);
    }
    bar.finish();
}

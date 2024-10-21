use crate::uniwig::Chromosome;
use ndarray::Array;
use ndarray_npy::write_npy;
use std::fs::{create_dir_all, remove_file, File, OpenOptions};
use std::io;
use std::io::{BufWriter, Write};

pub fn write_to_npy_file(
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
    wig_header.push('\n');
    file.write_all(wig_header.as_ref()).unwrap();
}

pub fn write_combined_wig_files(
    location: &str,
    output_type: &str,
    bwfileheader: &str,
    chromosomes: &[Chromosome],
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
pub fn write_to_wig_file(
    counts: &[u32],
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

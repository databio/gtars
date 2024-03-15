use clap::ArgMatches;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::fs::{File};
use std::error::Error;


pub mod cli;

pub mod consts {
    pub const UNIWIG_CMD: &str = "uniwig";

}

pub struct Chromosome {
    chrom: String,
    starts: Vec<u32>,
    ends: Vec<u32>,
}

pub fn show_chromosomes_map(){

}

pub fn show_chromosomes_vec(){


}

pub fn read_bed_map(){

}

pub fn read_bed_vec(){

}

pub fn parse_bed_file(line: &str) -> Option<(String, i32, i32)> {
    // TODO Eventually refactor all bed file parsing to a single shared function

    let mut fields = line.split('\t');
    // Get the first field which should be chromosome.
    let ctg = fields.next()?;
    // Parse 2nd and 3rd string as integers or return -1 if failure
    let st = fields.next().and_then(|s| s.parse::<i32>().ok()).unwrap_or(-1);
    let en = fields.next().and_then(|s| s.parse::<i32>().ok()).unwrap_or(-1);

    // Original code had a remainder of the line, r, but it does not appear to have been used
    // in any way

    Some((ctg.parse().unwrap(), st, en))

}


pub fn run_uniwig(matches: &ArgMatches) {
    println!("Im running. Here are the arguments: {:?}", matches);


    // Placeholder Arguments

    let sorted: bool = true;
    let smoothsize: i32 = 5;
    let writesize: i32 = 1;
    let combinedbedpath: &str = "/home/drc/GITHUB/genimtools/genimtools/tests/data/peaks.bed";
    let chromsizerefpath: String = "/home/drc/GITHUB/genimtools/genimtools/tests/hg38.chrom.sizes".to_string();
    let bwfileheader: &str = "test";


    uniwig_main(sorted, smoothsize, writesize, combinedbedpath,chromsizerefpath,bwfileheader)


}

pub fn uniwig_main(sorted: bool, smoothsize:i32, writesize:i32, combinedbedpath: &str,chromsizerefpath:String,bwfileheader: &str){
    // Main Function

    println!("Hello from Uniwig main");

    match read_chromosome_sizes(combinedbedpath) {
        Ok(chrom_sizes) => {
            println!("Chromosome sizes:");
            for (chrom, size) in chrom_sizes.iter() {
                println!("{}: {}", chrom, size);
            }
        }
        Err(err) => println!("Error reading chromosome sizes: {}", err),
    }


}

fn read_chromosome_sizes(chrom_size_path: &str) -> Result<std::collections::HashMap<String, i32>, Box<dyn Error>> {
    let chrom_size_file = File::open(Path::new(chrom_size_path))?;
    let mut chrom_sizes = std::collections::HashMap::new();
    let reader = BufReader::new(chrom_size_file);

    for line in reader.lines() {
        let line = line?; // Propagate the potential error
        let mut iter = line.split('\t');
        let chrom_name = iter.next().unwrap().to_owned();
        let size_str = iter.next().unwrap();
        let size = size_str.parse::<i32>()?;

        chrom_sizes.insert(chrom_name, size);
    }

    Ok(chrom_sizes)
}

use clap::ArgMatches;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::fs::{File};
use std::error::Error;
use clap::builder::OsStr;
use flate2::read::GzDecoder;


pub mod cli;

pub mod consts {
    pub const UNIWIG_CMD: &str = "uniwig";

}

pub struct Chromosome {
    chrom: String,
    starts: Vec<i32>,
    ends: Vec<i32>,
}
impl Clone for Chromosome {
    fn clone(&self) -> Self {
        Self {
            chrom: self.chrom.clone(), // Clone the string
            starts: self.starts.clone(), // Clone the vector
            ends: self.ends.clone(),   // Clone the vector
        }
    }
}


pub fn read_bed_map(combinedbedpath: &str){


}

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

    let mut chromosome = Chromosome{
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

        if chrom.is_empty(){
            // Initial chromosome
            chromosome.chrom = parsed_chr.clone();
            chrom = parsed_chr.clone();
            chromosome.starts.push(parsed_start);
            chromosome.ends.push(parsed_end);
        }


        if parsed_chr != chrom{

            // If the parsed chrom is not the same as the current, sort, and then push to vector
            // then reset chromosome struct using the newest parsed_chr
            chromosome.starts.sort_unstable();
            chromosome.ends.sort_unstable();

            chromosome_vec.push(chromosome.clone());

            chromosome.chrom =parsed_chr;

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

    return chromosome_vec

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
    println!("I am running. Here are the arguments: {:?}", matches);


    // Placeholder Arguments

    let sorted: bool = true;
    let smoothsize: i32 = 5;
    let writesize: i32 = 1;
    let combinedbedpath: &str = "/home/drc/GITHUB/genimtools/genimtools/tests/data/peaks.bed";
    let chromsizerefpath: String = "/home/drc/GITHUB/genimtools/genimtools/tests/hg38.chrom.sizes".to_string();
    let bwfileheader: &str = "test";


    uniwig_main(sorted, smoothsize, writesize, combinedbedpath,chromsizerefpath,bwfileheader)


}

pub fn uniwig_main(sorted: bool, _smoothsize:i32, _writesize:i32, combinedbedpath: &str, _chromsizerefpath:String, bwfileheader: &str){
    // Main Function

    println!("Hello from Uniwig main");

    // Set up output file names

    let mut file_names: [String; 3] = ["placeholder1".to_owned(), "placeholder2".to_owned(), "placeholder3".to_owned()];

    file_names[0] = format!("{}_{}", bwfileheader, "start.bw");
    file_names[1] = format!("{}_{}", bwfileheader, "end.bw");
    file_names[2] = format!("{}_{}", bwfileheader, "core.bw");

    let _chrom_sizes = match read_chromosome_sizes(combinedbedpath) {
        Ok(chrom_sizes) => chrom_sizes,
        Err(err) => {
            println!("Error reading chromosome sizes: {}", err);
            return; // Exit the main function on error
        }
    };


    if sorted {

        println!("Sorted is true");

        let mut _chromosomes: Vec<Chromosome> = read_bed_vec(combinedbedpath);


    } else{
        println!("read_bed_map goes here if sorted is untrue");
        // std::map<std::string, chromosome> chromosomes;
        read_bed_map(combinedbedpath);


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

use std::collections::HashMap;
use clap::ArgMatches;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::fs::{File};
use std::error::Error;
use bigtools::BBIFile::BigWig;
use clap::builder::OsStr;
use flate2::read::GzDecoder;
use bigtools::BigWigWrite;
use bigtools::bedchromdata::BedParserStreamingIterator;
use bigtools::bed::bedparser::BedParser;


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
    let bwfileheader: &str = "/home/drc/Downloads/test";


    uniwig_main(sorted, smoothsize, writesize, combinedbedpath,chromsizerefpath,bwfileheader)


}

pub fn uniwig_main(sorted: bool, smoothsize:i32, _writesize:i32, combinedbedpath: &str, _chromsizerefpath:String, bwfileheader: &str){
    // Main Function

    println!("Hello from Uniwig main");

    // Set up output file names

    let mut file_names: [String; 3] = ["placeholder1".to_owned(), "placeholder2".to_owned(), "placeholder3".to_owned()];

    file_names[0] = format!("{}_{}", bwfileheader, "start.bw");
    file_names[1] = format!("{}_{}", bwfileheader, "end.bw");
    file_names[2] = format!("{}_{}", bwfileheader, "core.bw");

    let chrom_sizes = match read_chromosome_sizes(combinedbedpath) {
        Ok(chrom_sizes) => chrom_sizes,
        Err(err) => {
            println!("Error reading chromosome sizes: {}", err);
            return; // Exit the main function on error
        }
    };


    if sorted {

        println!("Sorted is true");

        let mut chromosomes: Vec<Chromosome> = read_bed_vec(combinedbedpath);

        let num_chromosomes = chromosomes.len();
        // Preallocate memory based on number of chromsomes from previous step
        let mut chroms: Vec<String> = Vec::with_capacity(num_chromosomes);
        let mut chr_lens: Vec<i32> = Vec::with_capacity(num_chromosomes);

        for chromosome in chromosomes.iter(){
            let chrom_name = chromosome.chrom.clone();
            chroms.push(chrom_name);
            chr_lens.push(chrom_sizes[&chromosome.chrom] as i32); // retrieve size from hashmap
        }

        // Original Steps
        // Create bigwig file
        // Create header from chroms and chr lens
        // write to bigwig file with smoothing IF smoothsize is set
        // original code skips this if smoothsize is not set
        // Close bigwig file

        // Using BigTools
        let runtime = tokio::runtime::Builder::new_multi_thread()
            .worker_threads(6)
            .build()
            .expect("Unable to create runtime.");

        // let vals_iter = BedParser::from_bed_file(combinedbedpath);
        // let vals = BedParserStreamingIterator::new(vals_iter, false);

        // Iterate 3 times to output the three different files.
        for j in 0..3 {
            // Original code uses:
            // bwOpen, then bwCreateChromList, then bwWriteHdr

            let mut success_count  = 0;
            let mut failure_count = 0;

            println!("Processing each chromosome...");

            let mut out = BigWigWrite::create_file(file_names[j].clone());

            if smoothsize!=0 {
                match j {
                    0 => {
                        println!("Write Starts Here");
                    },
                    1 => {
                        println!("Write Ends Here");
                    },
                    2 => {
                        println!("Write Core Here");
                    },
                    _ => println!("Unexpected value: {}", j), // Handle unexpected values
                }
            }



        }



        // Using BigTools Bed Parsing as Alternative

        //let path = Path::new(combinedbedpath);
        // let path = PathBuf::from(combinedbedpath);
        //
        // let file = File::open(path).unwrap();
        //
        // // let is_gzipped = path.extension().unwrap_or(&OsStr::from("bed")) == "gz";
        // //
        // // // We must encapsulate in a box and use a dynamic Read trait so that either case could continue.
        // // let reader: Box<dyn Read> = match is_gzipped {
        // //     true => Box::new(GzDecoder::new(file)),
        // //     false => Box::new(file),
        // // };
        //
        // //let reader = BufReader::new(file);
        //
        // let vals_iter = BedParser::from_bed_file(file);
        //
        // let vals = BedParserStreamingIterator::new(vals_iter, true);
        //
        // println!("DONE");
        // let mut out = BigWigWrite::create_file(file_names[0].clone());
        //
        // out.options.block_size = 5;

        // WHen opening bed file using the bed parser, the func returns Ok((chrom, BedEntry { start, end, rest })
        // from the testing case, the bigtools crate opens from a bedgraph which returns Some(Ok((chrom, Value { start, end, value })))
        // Value is required (not BedEntry) when writing to a BigWig file (it throws a compiler error).
        // out.write(chrom_sizes, vals, runtime).unwrap();
        // let mut chrom_map = HashMap::new();
        // chrom_map.insert("chr17".to_string(), 83257441);

        //out.write(chrom_map, vals, runtime).unwrap();





    } else{
        println!("read_bed_map goes here if sorted is untrue");
        // std::map<std::string, chromosome> chromosomes;
        read_bed_map(combinedbedpath);


    }



}

fn read_chromosome_sizes(chrom_size_path: &str) -> Result<std::collections::HashMap<String, u32>, Box<dyn Error>> {
    let chrom_size_file = File::open(Path::new(chrom_size_path))?;
    let mut chrom_sizes = std::collections::HashMap::new();
    let reader = BufReader::new(chrom_size_file);

    for line in reader.lines() {
        let line = line?; // Propagate the potential error
        let mut iter = line.split('\t');
        let chrom_name = iter.next().unwrap().to_owned();
        let size_str = iter.next().unwrap();
        let size = size_str.parse::<u32>()?;

        chrom_sizes.insert(chrom_name, size);
    }

    Ok(chrom_sizes)
}

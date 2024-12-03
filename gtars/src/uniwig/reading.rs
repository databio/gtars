use crate::uniwig::Chromosome;
use clap::builder::OsStr;
use flate2::read::GzDecoder;
use noodles::bam;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::ops::Deref;
use std::path::Path;

//const UNMAPPED: &str = "*";

/// Reads combined bed file from a given path.
/// Returns Vec of Chromosome struct
pub fn read_bed_vec(combinedbedpath: &str) -> Vec<Chromosome> {
    let default_score = 1; // this will later be used for the count, which, by default, was originally = 1
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
            chromosome.starts.push((parsed_start, default_score));
            chromosome.ends.push((parsed_end, default_score));
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

        chromosome.starts.push((parsed_start, default_score));
        chromosome.ends.push((parsed_end, default_score));
    }

    // Is this final sort and push actually necessary?
    chromosome.starts.sort_unstable();
    chromosome.ends.sort_unstable();
    chromosome_vec.push(chromosome.clone());

    println!("Reading Bed file complete.");

    // println!("Here are chrom starts");
    //
    // for start in chromosome.starts.iter(){
    //     println!("{}",start.0);
    // }
    //
    // println!("Here are chrom ends");
    //
    // for end in chromosome.ends.iter(){
    //     println!("{}",end.0);
    // }

    chromosome_vec
}

/// Reads narrowPeak files and returns a `Vec<Chromosome>`
/// Pushes narrowPeak scores along with chrom coordinates.
/// Differs from read_bed_vec in that read_bed_vec pushes a default score (1).
pub fn read_narrow_peak_vec(combinedbedpath: &str) -> Vec<Chromosome> {
    // For narrowpeak there is no default score, we attempt to parse it from the file
    //
    let path = Path::new(combinedbedpath);

    let file = File::open(path).unwrap();

    let is_gzipped = path.extension().unwrap_or(&OsStr::from("narrowpeak")) == "gz";

    // We must encapsulate in a box and use a dynamic Read trait so that either case could continue.
    let reader: Box<dyn Read> = match is_gzipped {
        true => Box::new(GzDecoder::new(file)),
        false => Box::new(file),
    };

    let reader = BufReader::new(reader);

    let mut npchromosome = Chromosome {
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

    // println!("Here are chrom starts");
    //
    // for start in npchromosome.starts.iter(){
    //     println!("{}",start.0);
    // }
    //
    // println!("Here are chrom ends");
    //
    // for end in npchromosome.ends.iter(){
    //     println!("{}",end.0);
    // }

    chromosome_vec
}

/// Parses narrowPeak file grabbing chrom (ctg), start, end, and a numerical score in column 5
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
/// This ignores any other columns beyond start and ends.
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
        Some("bed") | Some("narrowPeak") => {
            // Read BED file
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

/// A wrapper around Noodles package to retrieve information from the bam header.
/// Returns a `Vec<Chromosome>`
pub fn read_bam_header(filepath: &str) -> Vec<Chromosome> {
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
        //let length: NonZeroUsize = ref_key.1.clone().length();
        // let length = ref_key.1.other_fields().clone().into_values();
        // let length = length.len();
        // for (key, value) in ref_key.1.other_fields().iter(){
        //     println!("new iteration");
        //     println!("here is key = {:?}", key);
        //     println!("here is value = {:?}", value);
        //     println!("Done");
        //
        // }
        // println!("here is length = {:?}", length);
        chromosome.chrom = chrom_name;
        chromosome.starts.push((0, 0)); //default values for now, less important for bam
        chromosome.ends.push((0, 0)); //default values for now, less important for bam
        chromosome_vec.push(chromosome.clone());
    }
    //
    // for c in &chromosome_vec{
    //     println!("chromsome= {:?}", c);
    // }
    //TODO this could just as easily be a Vec<String>?
    // In fact I think we later convert to Vec<String> after assessing the final chromosomes.
    chromosome_vec
}

// pub fn get_seq_reads_bam(chromosome: &mut Chromosome, filepath: &str) {
//     // read bam seq info into the current Chromosome
//
//     // TODO this function requires there to be an associated .bai file in the same directory as the .bam file
//     // And the error message if it does not exist is not very helpful.
//     let src = String::from(filepath);
//     let raw_region = String::from(chromosome.chrom.clone());
//     //let raw_region = String::from("chr1");
//
//     let mut reader = bam::io::indexed_reader::Builder::default()
//         .build_from_path(src)
//         .unwrap();
//     let header = reader.read_header().unwrap();
//
//     let records: Box<dyn Iterator<Item = io::Result<bam::Record>>> = if raw_region == UNMAPPED {
//         reader.query_unmapped().map(Box::new).unwrap()
//     } else {
//         let region = raw_region.parse().unwrap();
//         reader.query(&header, &region).map(Box::new).unwrap()
//     };
//
//     // remove the placeholder (0,0 ))
//     chromosome.starts.remove(0);
//     chromosome.ends.remove(0);
//     let default_score = 1;
//
//     for result in records {
//         let record = result.unwrap();
//         //TODO Determine position shift via what flags are set
//         let start_position = record.alignment_start().unwrap().unwrap();
//         let start = start_position.get();
//         let end_position = record.alignment_end().unwrap().unwrap();
//         let end = end_position.get();
//         chromosome.starts.push((start as i32, default_score));
//         chromosome.ends.push((end as i32, default_score));
//     }
//
//     chromosome.starts.sort_unstable_by(|a, b| a.0.cmp(&b.0));
//     chromosome.ends.sort_unstable_by(|a, b| a.0.cmp(&b.0));
//
//     println!(
//         "Finished reading seq for chrom: {}",
//         chromosome.chrom.clone()
//     );
// }

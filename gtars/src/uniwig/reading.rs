use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::ops::Deref;
use std::path::Path;
use clap::builder::OsStr;
use flate2::read::GzDecoder;
use noodles::bam;
use crate::uniwig::Chromosome;

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
        starts_with_scores: vec![],
        ends_with_scores: vec![],
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

pub fn read_narrow_peak_vec(combinedbedpath: &str) -> Vec<Chromosome> {
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
        starts_with_scores: vec![],
        ends_with_scores: vec![],
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
            npchromosome
                .starts_with_scores
                .push((parsed_start, parsed_score));
            npchromosome
                .ends_with_scores
                .push((parsed_end, parsed_score));
            continue;
        }

        if String::from(parsed_chr.trim()) != chrom {
            // If the parsed chrom is not the same as the current, sort, and then push to vector
            // then reset chromosome struct using the newest parsed_chr
            //npchromosome.starts.sort_unstable();
            //npchromosome.ends.sort_unstable();
            npchromosome
                .starts_with_scores
                .sort_unstable_by(|a, b| a.0.cmp(&b.0));
            npchromosome
                .ends_with_scores
                .sort_unstable_by(|a, b| a.0.cmp(&b.0));

            chromosome_vec.push(npchromosome.clone());

            npchromosome.chrom = String::from(parsed_chr.trim());
            chrom = String::from(parsed_chr.trim());

            npchromosome.starts_with_scores = vec![];
            npchromosome.ends_with_scores = vec![]
        }

        npchromosome
            .starts_with_scores
            .push((parsed_start, parsed_score));
        npchromosome
            .ends_with_scores
            .push((parsed_end, parsed_score));
    }

    // Is this final sort and push actually necessary?
    // npchromosome.starts.sort_unstable();
    // npchromosome.ends.sort_unstable();
    npchromosome
        .starts_with_scores
        .sort_unstable_by(|a, b| a.0.cmp(&b.0));
    npchromosome
        .ends_with_scores
        .sort_unstable_by(|a, b| a.0.cmp(&b.0));
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
        starts_with_scores: vec![],
        ends_with_scores: vec![],
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
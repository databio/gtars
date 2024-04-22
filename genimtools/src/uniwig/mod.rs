use clap::ArgMatches;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;
use std::fs::{File, OpenOptions};
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
            chromosome.chrom = String::from(parsed_chr.trim());
            chrom = String::from(parsed_chr.trim());
            chromosome.starts.push(parsed_start);
            chromosome.ends.push(parsed_end);
            continue;
        }


        if String::from(parsed_chr.trim()) != chrom{

            // If the parsed chrom is not the same as the current, sort, and then push to vector
            // then reset chromosome struct using the newest parsed_chr
            chromosome.starts.sort_unstable();
            chromosome.ends.sort_unstable();

            chromosome_vec.push(chromosome.clone());

            chromosome.chrom =String::from(parsed_chr.trim());
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

    //chromosome_vec.sort_by_key(|c| c.chrom.clone());

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

    let combinedbedpath = matches
        .get_one::<String>("bed")
        .expect("combined bed path is required");

    let chromsizerefpath = matches
        .get_one::<String>("chromref")
        .expect("chromref path path is required");

    let bwfileheader = matches
        .get_one::<String>("fileheader")
        .expect("fileheader is required");

    let smoothsize = matches
        .get_one::<i32>("smoothsize")
        .expect("smoothsize required");

    let output_type = matches
        .get_one::<String>("outputtype")
        .expect("output type is required");


    uniwig_main(*smoothsize, combinedbedpath, chromsizerefpath, bwfileheader, output_type)


}

pub fn uniwig_main(smoothsize:i32, combinedbedpath: &str, _chromsizerefpath: &String, bwfileheader: &str, output_type: &str){
    // Main Function

    //println!("Hello from Uniwig main");

    let stepsize = 1;

    // Set up output file names

    let mut file_names: [String; 3] = ["placeholder1".to_owned(), "placeholder2".to_owned(), "placeholder3".to_owned()];

    // TODO determine potential file types
    file_names[0] = format!("{}_{}.{}", bwfileheader, "start", output_type);
    file_names[1] = format!("{}_{}.{}", bwfileheader, "end", output_type);
    file_names[2] = format!("{}_{}.{}", bwfileheader, "core", output_type);



    let chrom_sizes = match read_chromosome_sizes(combinedbedpath) {
        // original program gets chromosome size from a .sizes file, e.g. chr1 248956422
        // the original program simply pushes 0's until the end of the chromosome length and writes these to file.
        // can we instead just use the last endsite for each chromosome to save space in th wiggle file?

        Ok(chrom_sizes) => chrom_sizes,
        Err(err) => {
            println!("Error reading chromosome sizes: {}", err);
            return; // Exit the main function on error
        }
    };





    let mut chromosomes: Vec<Chromosome> = read_bed_vec(combinedbedpath);

    let num_chromosomes = chromosomes.len();

    println!(" DEBUG Number of Chromosomes{:?}", num_chromosomes);

    // Preallocate memory based on number of chromsomes from previous step
    let mut chroms: Vec<String> = Vec::with_capacity(num_chromosomes);
    let mut chr_lens: Vec<i32> = Vec::with_capacity(num_chromosomes);

    println!("Processing each chromosome...");
    for chromosome in chromosomes.iter() {


        if chromosome.starts.len() != chromosome.ends.len(){
            println!("Chromosome starts and ends are not equal!");
            break
        }

        let chrom_name = chromosome.chrom.clone();
        //println!("DEBUG: CHROM NAME -> {}",chromosome.chrom.clone());
        chroms.push(chrom_name.clone());

        //chr_lens.push(chrom_sizes[&chromosome.chrom] as i32); // retrieve size from hashmap
        let current_chrom_size =chrom_sizes[&chromosome.chrom] as i32;


        // Iterate 3 times to output the three different files.
        for j in 0..3 {
            // Original code uses:
            // bwOpen, then bwCreateChromList, then bwWriteHdr

            let mut _success_count = 0;
            let mut _failure_count = 0;



            if smoothsize != 0 {
                match j {
                    0 => {
                        println!("Write Starts Here");
                        //println!("DEBUG: HERE is Initial VEC FOR STARTS:{:?}", chromosome.starts.clone());
                        //let count_result = count_coordinate_reads(&chromosome.starts);
                        //println!("DEBUG: HERE is COUNT VEC FOR STARTS:{:?}", result);

                        let count_result = smooth_Fixed_Start_End_Wiggle(&chromosome.starts,current_chrom_size,smoothsize, stepsize);

                        match output_type {
                            "wig" => {

                                println!("Writing to wig file!");
                                write_to_wig_file(&count_result.1, &count_result.0, file_names[0].clone(), chrom_name.clone());


                            },
                            "csv" => {println!("Write to CSV. Not Implemented");},
                            _ => {println!("Default to wig file.")},
                        }
                    },
                    1 => {
                        println!("Write Ends Here");
                        //let count_result = count_coordinate_reads(&chromosome.ends);
                        let count_result = smooth_Fixed_Start_End_Wiggle(&chromosome.ends,current_chrom_size,smoothsize, stepsize);
                        //println!("DEBUG: HERE is COUNT VEC FOR STARTS:{:?}", result);

                        match output_type {
                            "wig" => {

                                println!("Writing to wig file!");
                                write_to_wig_file(&count_result.1, &count_result.0, file_names[1].clone(), chrom_name.clone());

                            },
                            "csv" => {println!("Write to CSV. Not Implemented");},
                            _ => {println!("Default to wig file.")},
                        }
                    },
                    2 => {

                            println!("Write Core Here");

                            let core_results = Fixed_Core_Wiggle(&chromosome.starts,&chromosome.ends,current_chrom_size, stepsize);

                            match output_type {
                                "wig" => {

                                    println!("Writing to CORE RESULTS wig file!");
                                    //write_to_wig_file(&chromosome.starts, &count_result, file_names[0].clone(), chrom_name.clone());
                                    write_to_wig_file(&core_results.1, &core_results.0, file_names[2].clone(), chrom_name.clone());


                                },
                                "csv" => {println!("Write to CSV. Not Implemented");},
                                _ => {println!("Default to wig file.")},
                            }

                    },
                    _ => println!("Unexpected value: {}", j), // Handle unexpected values
                }
            }
        }
    }







}

fn write_to_wig_file(coordinates: &Vec<i32>, counts: &Vec<u32>, filename: String, chromname: String) {

    let mut file = OpenOptions::new()
        .create(true)  // Create the file if it doesn't exist
        .append(true)  // Append data to the existing file if it does exist
        .open(filename).unwrap();

    //println!("DEBUG: fixedStep chrom={}",chromname.clone());
    let wig_header = "fixedStep chrom=".to_string() + chromname.as_str() + " start=1 step=1";
    file.write_all(wig_header.as_ref()).unwrap();
    file.write_all(b"\n").unwrap();

    let mut position = 0;

    for count in counts.iter(){
        //TODO THis is inefficient to iterate over ALL counts when the above coordinate vecs could act as an index
        if *count == 0 {
            position += 1;
            continue
        } else{

            //println!("DEBUG COORDINATE = {} COUNTS= {}",position, count);
            //let wig_line = position.to_string() + " " + count.to_string().as_str();
            let wig_line = count.to_string();
            file.write_all(wig_line.as_ref()).unwrap();
            file.write_all(b"\n").unwrap();
            position+=1;
        }

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
        let _ = iter.next().unwrap();
        let size_str = iter.next().unwrap(); // we really want the 3rd column which is the end column.
        let size = size_str.parse::<u32>()?;

        chrom_sizes.insert(chrom_name, size);
    }

    Ok(chrom_sizes)
}

pub fn count_coordinate_reads(input_vector: &Vec<i32>) -> Vec<u8> {
    // Take a pre-sorted vector of potentially repeated positions and count the repeats for each position
    // else place a 0 at the position if no counts exist.

    // based on smoothFixedStarEndBW from orig uniwig but does not use a stepsize nor perform any smoothing

    //println!("DEBUG: Executing count_coordinate_reads");

    let vin_iter = input_vector.iter();
    let mut v_coord_counts: Vec<u8> = Vec::new(); // u8 stores 0:255 This may be insufficient. u16 max is 65535

    let mut coordinate_position = 1;
    let mut count = 0;

    let mut coordinate_value = 0;
    let mut prev_coordinate_value = 0;

    for coord in vin_iter{

        coordinate_value = *coord;

        if coordinate_value == prev_coordinate_value
        {
            count +=1;
            continue;

        }
        while prev_coordinate_value > coordinate_position {
            // add zeros in-between reads and increment until we "catch up" to the next coordinate position in the vector
            v_coord_counts.push(0);
            coordinate_position +=1;
        }

        v_coord_counts.push(count);
        prev_coordinate_value = coordinate_value;
        count = 1;
        coordinate_position +=1;
    }

    // Must finish out final value
    while coordinate_value > coordinate_position{
        v_coord_counts.push(0);
        coordinate_position += 1;
    }

    v_coord_counts.push(count);

    return v_coord_counts
}

pub fn count_coordinate_reads_start_end(starts_vector: &Vec<i32>, ends_vector: &Vec<i32>) -> Vec<u8> {
    // Take a pre-sorted vector of potentially repeated positions and count the repeats for each position
    // within a window based on the end point
    // else place a 0 at the position if no counts exist.

    // based on fixedCoreBW from orig uniwig but does not use a stepsize

    // TODO in progress

    //println!("DEBUG: Executing count_coordinate_reads");

    let vin_iter = starts_vector.iter();
    let mut v_coord_counts: Vec<u8> = Vec::new(); // u8 stores 0:255 This may be insufficient. u16 max is 65535

    let mut coordinate_position = 1;
    let mut count = 0;

    let mut coordinate_value = 0;
    let mut prev_coordinate_value = 0;

    for coord in vin_iter{

        coordinate_value = *coord;

        if coordinate_value == prev_coordinate_value
        {
            count +=1;
            continue;

        }
        while prev_coordinate_value > coordinate_position {
            // add zeros in-between reads and increment until we "catch up" to the next coordinate position in the vector
            v_coord_counts.push(0);
            coordinate_position +=1;
        }

        v_coord_counts.push(count);
        prev_coordinate_value = coordinate_value;
        count = 1;
        coordinate_position +=1;
    }

    // Must finish out final value
    while coordinate_value > coordinate_position{
        v_coord_counts.push(0);
        coordinate_position += 1;
    }

    v_coord_counts.push(count);

    return v_coord_counts
}

pub fn smooth_Fixed_Start_End_Wiggle(starts_vector: &Vec<i32>, chrom_size: i32, smoothsize: i32, stepsize:i32) -> (Vec<u32>, Vec<i32>) {
    // This function is a more direct port of smoothFixedStartEndBW from uniwig written in CPP
    // It allows the user to accumulate reads of either starts or ends
    // Counts occur between a start coordinate (cutSite) and an end site (endSite) where the endsite is determined based on
    // the level of smoothing.
    // counts are reported over a stepsize (with a default of stepsize = 1)
    // Unlike the original function, it does not write to disk in chunks. it simply returns a vector of accumulated reads.



    //println!("BEGIN smooth_Fixed_Start_End_Wiggle");

    let vin_iter = starts_vector.iter();

    let mut v_coordinate_positions: Vec<i32> = Vec::new(); // these are the final coordinates after any adjustments
    let mut v_coord_counts: Vec<u32> = Vec::new(); // u8 stores 0:255 This may be insufficient. u16 max is 65535

    let mut coordinate_position = 1;

    let mut count:u32 = 0;

    let mut coordinate_value = 0;
    let mut prev_coordinate_value = 0;


    let mut adjusted_start_site =0;
    let mut current_end_site = 0;

    let mut collected_end_sites: Vec<i32> = Vec::new();

    adjusted_start_site = starts_vector[0].clone(); // get first coordinate position
    adjusted_start_site = adjusted_start_site - smoothsize; // adjust based on smoothing

    //Check endsite generation
    current_end_site = adjusted_start_site + 1 + smoothsize*2;

    if adjusted_start_site < 1{
        adjusted_start_site = 1;
    }

    while coordinate_position < adjusted_start_site{
        // Just skip until we reach the initial adjusted start position
        // Note that this function will not return 0s at locations before the initial start site
        coordinate_position = coordinate_position + stepsize;
    }

    //prev_coordinate_value = adjusted_start_site;

    for coord in vin_iter {
        coordinate_value = *coord;
        adjusted_start_site = coordinate_value - smoothsize;
        count += 1;

        if adjusted_start_site < 1{
            adjusted_start_site = 1;
        }

        //current_end_site = adjusted_start_site + 1 + smoothsize*2; //

        collected_end_sites.push(adjusted_start_site + 1 + smoothsize*2);

        if adjusted_start_site == prev_coordinate_value
        {
            count +=1;
            continue;

        }

        while coordinate_position < adjusted_start_site{

            while current_end_site==coordinate_position{

                count = count - 1;

                if collected_end_sites.last() == None {
                    current_end_site = 0; // From original code. Double check this is the proper way.
                } else {
                    current_end_site = collected_end_sites.remove(0)
                }

            }

            if coordinate_position%stepsize == 0{
                // Step size defaults to 1, so report every value
                v_coord_counts.push(count);
                v_coordinate_positions.push(coordinate_position);
                //println!("DEBUG: Reporting count: {} at position: {}",count, coordinate_position);

            }

            //println!("DEBUG: Incrementing coordinate_position: {}  -> {}", coordinate_position,  coordinate_position +1);
            coordinate_position = coordinate_position + 1;


        }

        prev_coordinate_value = adjusted_start_site;

    }

    count = count + 1; // We must add 1 extra value here so that our calculation during the tail as we close out the end sites does not go negative.
    // this is because the code above subtracts twice during the INITIAL end site closure. So we are missing one count and need to make it up else we go negative.
    //

    while coordinate_position <= chrom_size{

        while current_end_site==coordinate_position{

            count = count - 1;

            if collected_end_sites.last() == None {
                current_end_site = 0; // From original code. Double check this is the proper way.
            } else {
                current_end_site = collected_end_sites.remove(0)
            }

        }

        if coordinate_position % stepsize == 0{
            // Step size defaults to 1, so report every value
            v_coord_counts.push(count);
            v_coordinate_positions.push(coordinate_position); // This is ONLY the starts
            //println!("DEBUG: Reporting count: {} at start position: {} and end position: ",count, coordinate_position);

        }

        //println!("DEBUG: Incrementing coordinate_position: {}  -> {}", coordinate_position,  coordinate_position +1);
        coordinate_position = coordinate_position + 1;


    }


    //println!("DEBUG: FINAL LENGTHS... Counts: {}  Positions: {}", v_coord_counts.len(), v_coordinate_positions.len());
    return (v_coord_counts, v_coordinate_positions)
}

pub fn Fixed_Core_Wiggle(starts_vector: &Vec<i32>, ends_vector: &Vec<i32>, chrom_size: i32, stepsize:i32) -> (Vec<u32>, Vec<i32>) {
    // This function is a more direct port of fixedCoreBW from uniwig written in CPP
    // It allows the user to accumulate reads of across paired starts and ends.
    // Counts occur between a start coordinate (cutSite) and an end site (endSite) where the endsite is determined based on
    // the paired ends.
    // Counts are reported over a stepsize (with a default of stepsize = 1)
    // Unlike the original function, it does not write to disk in chunks. it simply returns a vector of accumulated reads.

    //println!("BEGIN Fixed_Core_Wiggle");

    //println!("STARTS VECTOR LENGTH: {}  END VECTORS LENGTH: {}", starts_vector.len().clone(), ends_vector.len().clone());

    // TODO STARTS AND ENDS MUST BE EQUAL

    let mut v_coordinate_positions: Vec<i32> = Vec::new(); // these are the final coordinates after any adjustments
    let mut v_coord_counts: Vec<u32> = Vec::new(); // u8 stores 0:255 This may be insufficient. u16 max is 65535

    let mut coordinate_position = 1;

    let mut count = 0;

    let mut coordinate_value = 0;
    let mut prev_coordinate_value = 0;


    let mut current_start_site =0;
    let mut current_end_site = 0;

    let mut collected_end_sites: Vec<i32> = Vec::new();

    current_start_site = starts_vector[0].clone(); // get first coordinate position
    current_end_site = ends_vector[0];

    //Check endsite generation
    //current_end_site = adjusted_start_site + 1 + smoothsize*2;

    if current_start_site < 1{
        current_start_site = 1;
    }

    while coordinate_position < current_start_site{
        // Just skip until we reach the initial adjusted start position
        // Note that this function will not return 0s at locations before the initial start site
        coordinate_position = coordinate_position + stepsize;
    }

    //prev_coordinate_value = current_start_site;

    for (index, coord) in starts_vector.iter().enumerate() {
        coordinate_value = *coord;

        current_start_site = coordinate_value;

        count += 1;

        if current_start_site < 1{
            current_start_site = 1;
        }

        let current_index =  index;

        //current_end_site = ends_vector[current_index];

        collected_end_sites.push(ends_vector[current_index]);

        if current_start_site == prev_coordinate_value
        {
            count +=1;
            continue;

        }

        while coordinate_position < current_start_site{

            while current_end_site==coordinate_position{

                count = count - 1;

                if collected_end_sites.last() == None {
                    current_end_site = 0; // From original code. Double check this is the proper way.
                } else {
                    current_end_site = collected_end_sites.remove(0)
                }

            }

            if coordinate_position % stepsize == 0{
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

    while coordinate_position <= chrom_size{

        while current_end_site==coordinate_position{

            count = count - 1;

            if collected_end_sites.last() == None {
                current_end_site = 0; // From original code. Double check this is the proper way.
            } else {
                current_end_site = collected_end_sites.remove(0)
            }

        }

        if coordinate_position % stepsize == 0{
            // Step size defaults to 1, so report every value
            v_coord_counts.push(count);
            v_coordinate_positions.push(coordinate_position); // This is ONLY the starts
            //println!("DEBUG: Reporting count: {} at start position: {} and end position: ",count, coordinate_position);

        }

        //println!("DEBUG: Incrementing coordinate_position: {}  -> {}", coordinate_position,  coordinate_position +1);
        coordinate_position = coordinate_position + 1;


    }



    //println!("DEBUG: FINAL LENGTHS... Counts: {}  Positions: {}", v_coord_counts.len(), v_coordinate_positions.len());
    return (v_coord_counts, v_coordinate_positions)
}
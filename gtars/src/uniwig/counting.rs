use bigtools::beddata::BedParserStreamingIterator;
use bigtools::utils::cli::bedgraphtobigwig::BedGraphToBigWigArgs;
use bigtools::{BigWigWrite, InputSortType};
use noodles::bam;
use noodles::bam::io::reader::Query;
use noodles::bam::io::Reader;
use noodles::bgzf;
use noodles::sam::alignment::Record;
use std::collections::HashMap;
use std::fs::{create_dir_all, File, OpenOptions};
use std::io;
use std::io::{stdout, BufRead, BufReader, BufWriter, Cursor, Error, Write};
use std::sync::Arc;
use std::os::unix::io::{AsRawFd, FromRawFd};
use tokio::runtime;

#[derive(Debug)]
pub enum BAMRecordError {
    IoError(std::io::Error),
    NoFirstRecord,
    IncorrectSel,
}

impl From<std::io::Error> for BAMRecordError {
    fn from(err: std::io::Error) -> Self {
        BAMRecordError::IoError(err)
    }
}


/// This function is a more direct port of smoothFixedStartEndBW from uniwig written in CPP.
/// It allows the user to accumulate reads of either starts or ends.
/// Counts occur between a start coordinate (cutSite) and an end site (endSite) where the endsite is determined based on
/// the level of smoothing.
/// counts are reported over a stepsize (with a default of stepsize = 1).
/// Unlike the original function, it does not write to disk in chunks. it simply returns a vector of accumulated reads.
#[allow(unused_variables)]
pub fn start_end_counts(
    starts_vector: &[(i32, i32)],
    chrom_size: i32,
    smoothsize: i32,
    stepsize: i32,
) -> (Vec<u32>, Vec<i32>) {
    //let vin_iter = starts_vector.iter();

    let mut v_coordinate_positions: Vec<i32> = Vec::new(); // these are the final coordinates after any adjustments
    let mut v_coord_counts: Vec<u32> = Vec::new(); // u8 stores 0:255 This may be insufficient. u16 max is 65535

    let mut coordinate_position = 1;

    let mut count: i32 = 0;

    let mut coordinate_value: (i32, i32);
    let mut prev_coordinate_value = 0;

    let mut adjusted_start_site: (i32, i32);
    let mut current_end_site: (i32, i32);

    let mut collected_end_sites: Vec<(i32, i32)> = Vec::new();

    adjusted_start_site = starts_vector[0]; // get first coordinate position

    adjusted_start_site.0 = adjusted_start_site.0 - smoothsize;

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

    for (index, coord) in starts_vector.iter().enumerate() {
        coordinate_value = *coord;

        adjusted_start_site = coordinate_value;
        adjusted_start_site.0 = coordinate_value.0 - smoothsize;

        let current_score = adjusted_start_site.1;

        count += current_score;

        if adjusted_start_site.0 < 1 {
            adjusted_start_site.0 = 1;
        }

        let current_index = index;

        let mut new_end_site = adjusted_start_site;
        new_end_site.0 = adjusted_start_site.0 + 1 + smoothsize * 2;
        collected_end_sites.push(new_end_site);

        if adjusted_start_site.0 == prev_coordinate_value {
            continue;
        }

        while coordinate_position < adjusted_start_site.0 {
            while current_end_site.0 == coordinate_position {
                count = count - current_score;

                if count < 0 {
                    count = 0;
                }

                if collected_end_sites.last() == None {
                    current_end_site.0 = 0;
                } else {
                    current_end_site = collected_end_sites.remove(0)
                }
            }

            if coordinate_position % stepsize == 0 {
                // Step size defaults to 1, so report every value
                v_coord_counts.push(count as u32);
                v_coordinate_positions.push(coordinate_position);
            }

            coordinate_position = coordinate_position + 1;
        }

        prev_coordinate_value = adjusted_start_site.0;
    }

    count = count + 1; // We must add 1 extra value here so that our calculation during the tail as we close out the end sites does not go negative.
                       // this is because the code above subtracts twice during the INITIAL end site closure. So we are missing one count and need to make it up else we go negative.

    while coordinate_position < chrom_size {
        // Apply a bound to push the final coordinates otherwise it will become truncated.

        while current_end_site.0 == coordinate_position {
            let current_score = adjusted_start_site.1;
            count = count - current_score;
            if count < 0 {
                count = 0;
            }

            if collected_end_sites.last() == None {
                current_end_site.0 = 0;
            } else {
                current_end_site = collected_end_sites.remove(0)
            }
        }

        if coordinate_position % stepsize == 0 {
            // Step size defaults to 1, so report every value
            v_coord_counts.push(count as u32);
            v_coordinate_positions.push(coordinate_position);
        }

        coordinate_position = coordinate_position + 1;
    }

    (v_coord_counts, v_coordinate_positions)
}

/// This function is a more direct port of fixedCoreBW from uniwig written in CPP
/// It allows the user to accumulate reads across paired starts and ends.
/// Counts occur between a start coordinate (cutSite) and an end site (endSite) where the endsite is determined based on
/// the paired ends.
/// Counts are reported over a stepsize (with a default of stepsize = 1)
/// Unlike the original function, it does not write to disk in chunks. it simply returns a vector of accumulated reads.
#[allow(unused_variables)]
pub fn core_counts(
    starts_vector: &[(i32, i32)],
    ends_vector: &[(i32, i32)],
    chrom_size: i32,
    stepsize: i32,
) -> (Vec<u32>, Vec<i32>) {
    let mut v_coordinate_positions: Vec<i32> = Vec::new(); // these are the final coordinates after any adjustments
    let mut v_coord_counts: Vec<u32> = Vec::new(); // u8 stores 0:255 This may be insufficient. u16 max is 65535

    let mut coordinate_position = 1;

    let mut count = 0;

    let mut coordinate_value: (i32, i32);
    let mut prev_coordinate_value = 0;

    let mut current_start_site: (i32, i32);
    let mut current_end_site: (i32, i32);

    let mut collected_end_sites: Vec<(i32, i32)> = Vec::new();

    current_start_site = starts_vector[0]; // get first coordinate position
    current_end_site = ends_vector[0];

    if current_start_site.0 < 1 {
        current_start_site.0 = 1;
    }

    while coordinate_position < current_start_site.0 {
        // Just skip until we reach the initial adjusted start position
        // Note that this function will not return 0s at locations before the initial start site
        coordinate_position = coordinate_position + stepsize;
    }

    for (index, coord) in starts_vector.iter().enumerate() {
        coordinate_value = *coord;

        current_start_site = coordinate_value;

        let current_score = current_start_site.1;
        count += current_score;

        if current_start_site.0 < 1 {
            current_start_site.0 = 1;
        }

        let current_index = index;

        collected_end_sites.push(ends_vector[current_index]);

        if current_start_site.0 == prev_coordinate_value {
            continue;
        }

        while coordinate_position < current_start_site.0 {
            while current_end_site.0 == coordinate_position {
                count = count - current_score;
                if count < 0 {
                    count = 0;
                }

                if collected_end_sites.last() == None {
                    current_end_site.0 = 0;
                } else {
                    current_end_site = collected_end_sites.remove(0)
                }
            }

            if coordinate_position % stepsize == 0 {
                // Step size defaults to 1, so report every value
                v_coord_counts.push(count as u32);
                v_coordinate_positions.push(coordinate_position);
            }

            coordinate_position = coordinate_position + 1;
        }

        prev_coordinate_value = current_start_site.0;
    }

    count = count + 1; // We must add 1 extra value here so that our calculation during the tail as we close out the end sites does not go negative.

    while coordinate_position < chrom_size {
        while current_end_site.0 == coordinate_position {
            let current_score = current_start_site.1;
            count = count - current_score;
            if count < 0 {
                count = 0;
            }
            if collected_end_sites.last().is_none() {
                current_end_site.0 = 0;
            } else {
                current_end_site = collected_end_sites.remove(0)
            }
        }

        if coordinate_position % stepsize == 0 {
            // Step size defaults to 1, so report every value
            v_coord_counts.push(count as u32);
            v_coordinate_positions.push(coordinate_position);
        }

        coordinate_position = coordinate_position + 1;
    }

    (v_coord_counts, v_coordinate_positions)
}

///Instead of counting based on in-memory chromosomes, this method takes a buffered reader and iterates
/// Primarily for use to count sequence reads in bam files.
pub fn fixed_start_end_counts_bam(
    records: &mut Box<Query<noodles::bgzf::reader::Reader<std::fs::File>>>,
    chrom_size: i32,
    smoothsize: i32,
    stepsize: i32,
    output_type: &str,
    chromosome_name: &String,
    bwfileheader: &str,
    out_sel: &str,
    std_out_sel: bool,
) -> (Vec<u32>, Vec<i32>) {
    //let vin_iter = starts_vector.iter();

    let mut v_coordinate_positions: Vec<i32> = Vec::new(); // these are the final coordinates after any adjustments
    let mut v_coord_counts: Vec<u32> = Vec::new(); // u8 stores 0:255 This may be insufficient. u16 max is 65535

    let mut coordinate_position = 1;

    let mut count: i32 = 0;

    let mut coordinate_value: i32;
    let mut prev_coordinate_value = 0;

    let mut adjusted_start_site: i32;
    let mut current_end_site: i32;

    let mut collected_end_sites: Vec<i32> = Vec::new();

    let first_record = records.next().unwrap().unwrap();

    let mut adjusted_start_site: i32 = match out_sel {
        "start" => first_record.alignment_start().unwrap().unwrap().get() as i32,
        "end" => first_record.alignment_end().unwrap().unwrap().get() as i32,
        _ => {
            panic!("unknown output selection must be either 'start', 'end', 'core'")
        }
    };

    //adjusted_start_site = first_record.alignment_start().unwrap().unwrap().get() as i32; // get first coordinate position

    adjusted_start_site = adjusted_start_site - smoothsize;

    //SETUP OUTPUT FILE HERE BECAUSE WE NEED TO KNOW INITIAL VALUES
    let file = set_up_file_output(
        output_type,
        adjusted_start_site,
        chromosome_name,
        bwfileheader,
        stepsize,
        out_sel,
        std_out_sel,
    );
    let file = file.unwrap();
    let mut buf = BufWriter::new(file);

    current_end_site = adjusted_start_site;
    current_end_site = adjusted_start_site + 1 + smoothsize * 2;

    if adjusted_start_site < 1 {
        adjusted_start_site = 1;
    }

    while coordinate_position < adjusted_start_site {
        // Just skip until we reach the initial adjusted start position
        // Note that this function will not return 0s at locations before the initial start site
        coordinate_position = coordinate_position + stepsize;
    }

    for coord in records {
        let mut coordinate_value: i32 = match out_sel {
            "start" => coord.unwrap().alignment_start().unwrap().unwrap().get() as i32,
            "end" => coord.unwrap().alignment_end().unwrap().unwrap().get() as i32,
            _ => {
                panic!("unknown output selection must be either 'start', 'end', 'core'")
            }
        };

        // coordinate_value = coord.unwrap().alignment_start().unwrap().unwrap().get() as i32;

        adjusted_start_site = coordinate_value;
        adjusted_start_site = coordinate_value - smoothsize;

        let current_score = adjusted_start_site;

        count += current_score;

        if adjusted_start_site < 1 {
            adjusted_start_site = 1;
        }

        //let current_index = index;

        let mut new_end_site = adjusted_start_site;
        new_end_site = adjusted_start_site + 1 + smoothsize * 2;
        collected_end_sites.push(new_end_site);

        if adjusted_start_site == prev_coordinate_value {
            continue;
        }

        while coordinate_position < adjusted_start_site {
            while current_end_site == coordinate_position {
                count = count - current_score;

                if count < 0 {
                    count = 0;
                }

                if collected_end_sites.last() == None {
                    current_end_site = 0;
                } else {
                    current_end_site = collected_end_sites.remove(0)
                }
            }

            if coordinate_position % stepsize == 0 {
                // Step size defaults to 1, so report every value
                //v_coord_counts.push(count as u32);

                match output_type {
                    "wig" => {
                        writeln!(&mut buf, "{}", count).unwrap();
                    }
                    "bedgraph" => {
                        writeln!(
                            &mut buf,
                            "{}\t{}\t{}\t{}",
                            chromosome_name, adjusted_start_site, current_end_site, count
                        )
                        .unwrap();
                    }
                    _ => {}
                }

                v_coordinate_positions.push(coordinate_position);
            }

            coordinate_position = coordinate_position + 1;
        }

        prev_coordinate_value = adjusted_start_site;
    }

    count = count + 1; // We must add 1 extra value here so that our calculation during the tail as we close out the end sites does not go negative.
                       // this is because the code above subtracts twice during the INITIAL end site closure. So we are missing one count and need to make it up else we go negative.

    while coordinate_position < chrom_size {
        // Apply a bound to push the final coordinates otherwise it will become truncated.

        while current_end_site == coordinate_position {
            let current_score = adjusted_start_site;
            count = count - current_score;
            if count < 0 {
                count = 0;
            }

            if collected_end_sites.last() == None {
                current_end_site = 0;
            } else {
                current_end_site = collected_end_sites.remove(0)
            }
        }

        if coordinate_position % stepsize == 0 {
            // Step size defaults to 1, so report every value
            //v_coord_counts.push(count as u32);
            match output_type {
                "wig" => {
                    writeln!(&mut buf, "{}", count).unwrap();
                }

                "bedgraph" => {
                    writeln!(
                        &mut buf,
                        "{}\t{}\t{}\t{}",
                        chromosome_name, adjusted_start_site, current_end_site, count
                    )
                    .unwrap();
                }

                _ => {}
            }
            v_coordinate_positions.push(coordinate_position);
        }

        coordinate_position = coordinate_position + 1;
    }

    buf.flush().unwrap();
    //println!("FInished with fixed_start_end_counts_bam");
    (v_coord_counts, v_coordinate_positions)
}


pub fn fixed_core_counts_bam_to_bw(
    records: &mut Box<Query<noodles::bgzf::reader::Reader<std::fs::File>>>,
    chrom_size: i32,
    stepsize: i32,
    chromosome_name: &String,
) -> Result<Cursor<String>, BAMRecordError> {

    let mut bedgraphlines = String::new();
    let mut coordinate_position = 1;
    let mut count: i32 = 0;
    let mut prev_coordinate_value = 0;
    let mut current_end_site: i32;
    let mut collected_end_sites: Vec<i32> = Vec::new();

    let first_record_option = records.next();

    let first_record = match first_record_option {
        Some(Ok(record)) => record,  // Extract the record
        Some(Err(err)) => {
            // Handle the error
            eprintln!("Error reading the first record for chrom: {} {:?} Skipping...", chromosome_name,err);
            return Err(BAMRecordError::NoFirstRecord);  // Example error handling
        }
        None => {
            // Handle no records
            eprintln!("Error reading the first record for chrom: {} Skipping...", chromosome_name);
            return Err(BAMRecordError::NoFirstRecord);
        }
    };

    let mut current_start_site = first_record.alignment_start().unwrap().unwrap().get() as i32;
    let mut current_end_site = first_record.alignment_end().unwrap().unwrap().get() as i32;

    if current_start_site < 1 {
        current_start_site = 1;
    }

    while coordinate_position < current_start_site {
        // Just skip until we reach the initial adjusted start position
        // Note that this function will not return 0s at locations before the initial start site
        coordinate_position = coordinate_position + stepsize;
    }

    for coord in records {

        let unwrapped_coord = coord.unwrap().clone();
        let mut current_start_site = unwrapped_coord.alignment_start().unwrap().unwrap().get() as i32;
        let new_end_site = unwrapped_coord.alignment_end().unwrap().unwrap().get() as i32;

        count += 1;

        if current_start_site < 1 {
            current_start_site = 1;
        }

        collected_end_sites.push(new_end_site);

        if current_start_site == prev_coordinate_value {
            continue;
        }

        while coordinate_position < current_start_site {
            while current_end_site == coordinate_position {
                count = count - 1;
                if count < 0 {
                    count = 0;
                }

                if collected_end_sites.last() == None {
                    current_end_site = 0;
                } else {
                    current_end_site = collected_end_sites.remove(0)
                }
            }

            if coordinate_position % stepsize == 0 {
                let single_line = format!("{}\t{}\t{}\t{}\n",
                                          chromosome_name, coordinate_position, coordinate_position+1, count);
                bedgraphlines.push_str(&*single_line);
            }

            coordinate_position = coordinate_position + 1;
        }

        prev_coordinate_value = current_start_site;
    }
    count = count + 1; // We must add 1 extra value here so that our calculation during the tail as we close out the end sites does not go negative.
    // this is because the code above subtracts twice during the INITIAL end site closure. So we are missing one count and need to make it up else we go negative.

    while coordinate_position < chrom_size {
        // Apply a bound to push the final coordinates otherwise it will become truncated.

        while current_end_site == coordinate_position {
            count = count - 1;
            if count < 0 {
                count = 0;
            }

            if collected_end_sites.last() == None {
                current_end_site = 0;
            } else {
                current_end_site = collected_end_sites.remove(0)
            }
        }

        if coordinate_position % stepsize == 0 {
            // Step size defaults to 1, so report every value
            let single_line = format!("{}\t{}\t{}\t{}\n",
                                      chromosome_name, coordinate_position, coordinate_position+1, count);
            bedgraphlines.push_str(&*single_line);
        }

        coordinate_position = coordinate_position + 1;
    }

    let cursor = Cursor::new(bedgraphlines);

    Ok(cursor)

}


///Instead of counting based on in-memory chromosomes, this method takes a buffered reader and iterates
/// Primarily for use to count sequence reads in bam files.
pub fn fixed_start_end_counts_bam_to_bw(
    records: &mut Box<Query<noodles::bgzf::reader::Reader<std::fs::File>>>,
    chrom_size: i32,
    smoothsize: i32,
    stepsize: i32,
    chromosome_name: &String,
    out_sel: &str,
    write_fd: Arc<dyn AsRawFd + Send + Sync>,
) -> Result<(), BAMRecordError> {
    let mut writer = std::io::BufWriter::new(unsafe { std::fs::File::from_raw_fd(write_fd.as_raw_fd()) });
    //let vin_iter = starts_vector.iter();

    //let mut vec_lines: Vec<String> = Vec::new();
    //let mut bedgraphlines = String::new();

    let mut v_coordinate_positions: Vec<i32> = Vec::new(); // these are the final coordinates after any adjustments
    let mut v_coord_counts: Vec<u32> = Vec::new(); // u8 stores 0:255 This may be insufficient. u16 max is 65535

    let mut coordinate_position = 1;

    let mut count: i32 = 0;

    let mut coordinate_value: i32;
    let mut prev_coordinate_value = 0;

    let mut adjusted_start_site: i32;
    let mut current_end_site: i32;

    let mut collected_end_sites: Vec<i32> = Vec::new();

    let first_record_option = records.next();

    let first_record = match first_record_option {
        Some(Ok(record)) => record,  // Extract the record
        Some(Err(err)) => {
            // Handle the error
            eprintln!("Error reading the first record for chrom: {} {:?} Skipping...", chromosome_name,err);
            writer.write_all(b"").unwrap();
            writer.flush().unwrap();
            return Err(BAMRecordError::NoFirstRecord);  // Example error handling
        }
        None => {
            // Handle no records
            eprintln!("Error reading the first record for chrom: {} Skipping...", chromosome_name);
            writer.write_all(b"").unwrap();
            writer.flush().unwrap();
            return Err(BAMRecordError::NoFirstRecord);
        }
    };


    let mut adjusted_start_site: i32 = match out_sel {
        "start" => first_record.alignment_start().unwrap().unwrap().get() as i32,
        "end" => first_record.alignment_end().unwrap().unwrap().get() as i32,
        _ => {
            writer.write_all(b"").unwrap();
            writer.flush().unwrap();
            return Err(BAMRecordError::IncorrectSel);  // Example error handling
            //panic!("unknown output selection must be either 'start', 'end', 'core'")
        }
    };


    adjusted_start_site = adjusted_start_site - smoothsize;


    current_end_site = adjusted_start_site;
    current_end_site = adjusted_start_site + 1 + smoothsize * 2;

    if adjusted_start_site < 1 {
        adjusted_start_site = 1;
    }

    while coordinate_position < adjusted_start_site {
        // Just skip until we reach the initial adjusted start position
        // Note that this function will not return 0s at locations before the initial start site
        coordinate_position = coordinate_position + stepsize;
    }

    for coord in records {
        let mut coordinate_value: i32 = match out_sel {
            "start" => coord.unwrap().alignment_start().unwrap().unwrap().get() as i32,
            "end" => coord.unwrap().alignment_end().unwrap().unwrap().get() as i32,
            _ => {
                writer.write_all(b"").unwrap();
                writer.flush().unwrap();
                return Err(BAMRecordError::IncorrectSel);
                //panic!("unknown output selection must be either 'start', 'end', 'core'")
            }
        };

        // coordinate_value = coord.unwrap().alignment_start().unwrap().unwrap().get() as i32;

        adjusted_start_site = coordinate_value;
        adjusted_start_site = coordinate_value - smoothsize;

        //let current_score = adjusted_start_site;

        count += 1;

        if adjusted_start_site < 1 {
            adjusted_start_site = 1;
        }

        //let current_index = index;

        let mut new_end_site = adjusted_start_site;
        new_end_site = adjusted_start_site + 1 + smoothsize * 2;
        collected_end_sites.push(new_end_site);

        if adjusted_start_site == prev_coordinate_value {
            continue;
        }

        while coordinate_position < adjusted_start_site {
            while current_end_site == coordinate_position {
                count = count - 1;

                if count < 0 {
                    count = 0;
                }

                if collected_end_sites.last() == None {
                    current_end_site = 0;
                } else {
                    current_end_site = collected_end_sites.remove(0)
                }
            }

            if coordinate_position % stepsize == 0 {
                let single_line = format!("{}\t{}\t{}\t{}\n",
                                          chromosome_name, coordinate_position, coordinate_position+1, count);
                writer.write_all(single_line.as_bytes())?;
                writer.flush()?;
                //eprintln!("{}",single_line);

            }

            coordinate_position = coordinate_position + 1;
        }

        prev_coordinate_value = adjusted_start_site;
    }

    count = count + 1; // We must add 1 extra value here so that our calculation during the tail as we close out the end sites does not go negative.
                       // this is because the code above subtracts twice during the INITIAL end site closure. So we are missing one count and need to make it up else we go negative.

    while coordinate_position < chrom_size {
        // Apply a bound to push the final coordinates otherwise it will become truncated.

        while current_end_site == coordinate_position {
            let current_score = adjusted_start_site;
            count = count - 1;
            if count < 0 {
                count = 0;
            }

            if collected_end_sites.last() == None {
                current_end_site = 0;
            } else {
                current_end_site = collected_end_sites.remove(0)
            }
        }

        if coordinate_position % stepsize == 0 {
            // Step size defaults to 1, so report every value
            let single_line = format!("{}\t{}\t{}\t{}\n",
                                      chromosome_name, coordinate_position, coordinate_position+1, count);
            writer.write_all(single_line.as_bytes())?;
            writer.flush()?;
            //eprintln!("{}",single_line);
        }

        coordinate_position = coordinate_position + 1;
    }

    Ok(())
}

fn set_up_file_output(
    output_type: &str,
    adjusted_start_site: i32,
    chromosome_name: &String,
    bwfileheader: &str,
    stepsize: i32,
    out_sel: &str,
    std_out_sel: bool,
) -> Result<Box<dyn Write>, io::Error> {
    if !std_out_sel {
        // SET UP FILE BASED ON NAME
        let filename = format!(
            "{}{}_{}.{}",
            bwfileheader, chromosome_name, out_sel, output_type
        );
        let path = std::path::Path::new(&filename).parent().unwrap();
        let _ = create_dir_all(path);
        //
        let mut file = OpenOptions::new()
            .create(true) // Create the file if it doesn't exist
            .append(true) // Append data to the existing file if it does exist
            .open(filename)
            .unwrap();

        match output_type {
            "wig" => {
                let wig_header = "fixedStep chrom=".to_string()
                    + chromosome_name.as_str()
                    + " "
                    + out_sel
                    + "="
                    + adjusted_start_site.to_string().as_str()
                    + " step="
                    + stepsize.to_string().as_str();
                file.write_all(wig_header.as_ref()).unwrap();
                file.write_all(b"\n").unwrap();
            }
            "bedgraph" => { // do nothing, no header needed
            }
            _ => {
                panic!("output type not recognized during file set up for writing!")
            }
        }

        Ok(Box::new(file))
    } else {
        Ok(Box::new(io::stdout()))
        // write to std_out, this will be useful for sending input to bigtools to create bw files
    }
}

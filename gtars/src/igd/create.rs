use crate::common::consts::{BED_FILE_EXTENSION, GZ_FILE_EXTENSION};
use crate::common::utils::get_dynamic_reader;
use anyhow::Result;
use byteorder::{LittleEndian, ReadBytesExt};
use clap::ArgMatches;
use std::collections::HashMap;
use std::fs;
use std::fs::{create_dir_all, File, OpenOptions};
use std::io::{BufRead, BufReader, Error, Read, Write};
use std::path::{Path, PathBuf};

pub const maxCount: i64 = 268435456; //16* = 4GB memory  // original code had this as i32

// Assuming a maximum length of 40 characters, original program constraint
pub const MAX_CHROM_NAME_LEN: usize = 40;

#[derive(Default, Clone)]
pub struct gdata_t {
    pub idx: i32,   //genomic object--data set index
    pub start: i32, //region start
    pub end: i32,   //region end
    pub value: i32,
}
impl gdata_t {
    /// Constructs new instance of a gdata_t
    pub fn new() -> Self {
        Self::default()
    }
}
#[derive(Default, Clone, Copy)]
pub struct gdata0_t {
    pub idx: i32,   //genomic object--data set index
    pub start: i32, //region start
    pub end: i32,   //region end
}
impl gdata0_t {
    /// Constructs new instance of a gdata0_t
    pub fn new() -> Self {
        Self::default()
    }
}

#[derive(Default, Clone)]
pub struct tile_t {
    pub ncnts: i32,          // batch counts
    pub nCnts: i32,          // total (batch) counts
    pub mcnts: i32,          //  max counts
    pub gList: Vec<gdata_t>, //genomic data
}
#[derive(Default)]
pub struct ctg_t {
    pub name: String,       //name of the contig
    pub mTiles: i32,        //determined by the interval start and end
    pub gTile: Vec<tile_t>, //tile data
}
impl ctg_t {
    /// Constructs new instance of a ctg
    pub fn new() -> Self {
        Self::default()
    }
}

#[derive(Default)]
pub struct igd_t {
    // this struct is used for SAVING to disk
    pub nbp: i32,        //data type: 0, 1, 2 etc; size differs
    pub gType: i32,      //data type: 0, 1, 2 etc; size differs
    pub nctg: i32,       //data type: 0, 1, 2 etc; size differs
    pub mctg: i32,       //data type: 0, 1, 2 etc; size differs
    pub total: i64,      // total region in each ctg
    pub ctg: Vec<ctg_t>, // this is the list of contigs (of size n-ctg)  // this might need to be a reference
    pub total_regions: i32,
    pub total_average: f32,
    pub average_length: f32,
}

impl igd_t {
    /// Constructs new instance of IGD
    pub fn new() -> Self {
        Self::default()
    }
}

impl tile_t {
    /// Constructs new instance of tile
    pub fn new() -> Self {
        Self::default()
    }
}

pub fn igd_get_create_matches(matches: &ArgMatches) {
    //println!("HELLO FROM IGD CREATE SUBMODULE!");

    let output_path = matches
        .get_one::<String>("output")
        .expect("Output path is required");

    let filelist = matches
        .get_one::<String>("filelist")
        .expect("File list path is required");

    let db_output_name = matches
        .get_one::<String>("dbname")
        .expect("File list path is required");

    let _igd = create_igd_f(output_path, filelist, db_output_name);
}

/// Creates IGD database from a directory of bed files.
pub fn create_igd_f(output_path: &String, filelist: &String, db_output_name: &String) -> igd_t {
    //println!("{}",db_output_name);
    //Initialize IGD into Memory
    let mut igd = igd_t::new();

    // create hash table
    let mut hash_table: HashMap<String, i32> = HashMap::new();

    igd.gType = 1;
    igd.nbp = 16384; // from og code tile_size = 16384;  -> this is the bin size (2^14) from the original paper
    igd.nctg = 0;
    igd.mctg = 32;
    igd.total = 0;

    //Check that file path exists and get number of files
    let mut all_bed_files: Vec<PathBuf> = Vec::new();

    let mut ix = 0;
    let (mut start, mut end) = (0, 0);
    let mut va: i32 = 0;

    // create Path obj from filepath
    let input_filepaths = if filelist.ends_with(".txt") {
        // if txt input, read paths from file
        let mut paths = Vec::new();
        if let Ok(file) = File::open(filelist) {
            let reader = BufReader::new(file);
            for line in reader.lines() {
                if let Ok(path) = line {
                    paths.push(PathBuf::from(path.trim()));
                }
            }
        }
        paths
    } else if filelist == "-" || filelist == "stdin" {
        // if you pass "-" assume you want to read files list from stdin
        let stdin = std::io::stdin();
        let locked = stdin.lock();
        let reader = BufReader::new(locked);

        let mut paths: Vec<PathBuf> = Vec::new();

        for line in reader.lines() {
            match line {
                Ok(line) => {
                    let path = PathBuf::from(line);
                    paths.push(path);
                }
                Err(e) => {
                    eprintln!("Error reading line: {}", e);
                }
            }
        }
        paths
    } else {
        // if dir input, get directory entries directly
        let entries = fs::read_dir(filelist).unwrap();
        let mut paths = Vec::new();

        for entry in entries {
            let p = entry.as_ref().unwrap().path();
            paths.push(p)
        }
        paths
    };

    //--------------------
    // Check each file and only keep the validated BED files
    //
    // -------------------
    for path in input_filepaths {
        // For now only take .bed files
        if let Some(extension) = path.extension() {
            if extension != BED_FILE_EXTENSION.trim_start_matches('.')
                && extension != GZ_FILE_EXTENSION.trim_start_matches('.')
            {
                continue;
            }
        } else {
            continue;
        } // This will skip files that do not have an extension

        let metadata = fs::metadata(&path).unwrap();
        let file_type = metadata.file_type();

        if file_type.is_file() {
            // open bed file
            // TODO original code uses gzopen (I assume for .gz files?)
            // let file = File::open(entry.path()).unwrap();
            // let mut reader = BufReader::new(file);

            let mut reader = get_dynamic_reader(&path).unwrap();

            // Read the very first line and see if it meets our criteria
            // MUST USE by_ref() otherwise borrow checker won't let code compile
            // ALSO bec careful to call by_ref() BEFORE .lines()
            //
            let first_line = reader.by_ref().lines().next().unwrap().expect("expect");

            //TODO Need to do error handling to ensure we gracefully continue if there is no data in the file.
            //let mut lines = reader.lines();

            // TODO Better name for og function?
            // TODO parse_bed -> parse_bed_file_line
            let ctg = parse_bed(&first_line, &mut start, &mut end, &mut va);
            // if it parses, add it to collected lines, increment ix
            match ctg {
                Some(_ctg) => {
                    //println!("ctg successfully parsed {}", ctg);
                    all_bed_files.push(path);
                    ix += 1;
                }
                None => continue,
            }
        }
    }

    //println!("ALL PARSED Lines from BED FILES:\n{:?}", all_bed_files);

    let n_files = ix; //all_bed_files.len();
    let _nf10 = n_files / 10;

    println!("Number of Bed Files found:\n{}", n_files);

    //Check that there is more than 0 files?

    //Prep memory allocation in a Rust-like manner
    // TODO original code checks that the bed file can be parsed BEFORE memory allocation
    // TODO but then re-parses the bed file again later.
    // og C code:
    //    int32_t *nr = calloc(n_files, sizeof(int32_t));
    //     double *avg = calloc(n_files, sizeof(double));
    let mut avg: Vec<i32> = Vec::with_capacity(n_files); //Can we use arrays? Is this an array? no, can we put an array on files.
    avg.resize(n_files, 0);

    let mut nr: Vec<i32> = Vec::with_capacity(n_files);
    nr.resize(n_files, 0);

    //--------------------
    // READ VALIDATED FILES
    // Note: this seems wasteful to load the file *again* using BufReader
    // Is there a better way than below?
    // -------------------
    // Initialize required variables
    let mut i0 = 0;

    let mut ig: usize;
    let mut m: i32;
    let (mut va, nf10) = (0, n_files / 10);

    while i0 < n_files {
        //from og code: 2.1 Start from (i0, L0): read till (i1, L1)
        ig = i0;
        m = 0;
        //from og code: 2.2 Read ~4GB data from files
        // og code skips first line (since its already in the vec but we need to reread the file.
        while m == 0 && ig < n_files {
            //og comment: m>0 defines breaks when reading maxCount

            let file_path_buf = &all_bed_files[ig]; // could not move all_bed_files, so using reference to the PathBuf
            let fp = file_path_buf.clone();
            // let file = File::open(fp).unwrap();
            // let mut reader = BufReader::new(file);

            let reader = get_dynamic_reader(&fp).unwrap();

            // let mut buffer = String::new();

            for line in reader.lines() {
                let line = line.expect("Error reading line"); // Handle errors
                if m != 0 {
                    break;
                }
                // TODO original code: if(nCols>4) va = atol(splits[4]);
                // assumes that 5th value it numeric from original .gz file. Is this valid?
                // va = score  ----> https://genome.ucsc.edu/FAQ/FAQformat.html#format1

                // for line in reader.lines() {
                //     let line = line.expect("Error reading line"); // Handle errors

                let ctg = parse_bed(&line, &mut start, &mut end, &mut va);

                match ctg {
                    Some(ctg) => {
                        // check that st>=0 and end <321000000   NOTE: these values taken from og code.
                        if start >= 0 && end < 321000000 {
                            igd_add(&mut igd, &mut hash_table, ctg, start, end, va, ig);
                            nr[ig] += 1;
                            avg[ig] += end - start;
                            //println!("DEBUG: after igd add");
                        }
                    }
                    None => continue,
                }

                if igd.total > maxCount {
                    m = 1;
                }
                //endpoint
            }

            if m == 0 {
                ig += 1;
            }

            if nf10 > 1 && ig % nf10 == 0 {
                println!(".") // Show progress for every 10 files
            }
        }

        //og: 2.3 save/append temp tiles to disc, add cnts to Cnts
        //
        igd_saveT(&mut igd, output_path);

        i0 = ig;
    }

    let tsv_save_path = format!("{}{}{}", output_path, db_output_name, ".tsv");
    let tsv_parent_path = tsv_save_path.clone();
    let path = std::path::Path::new(&tsv_parent_path).parent().unwrap();
    let result = create_file_with_parents(path);

    match result {
        Ok(_file) => (),
        Err(err) => println!("Error creating file: {}", err),
    }
    let mut file = OpenOptions::new()
        .create(true) // Create the file if it doesn't exist
        .append(true) // Append data to the existing file if it does exist
        .open(tsv_save_path)
        .unwrap();

    let initial_line = format!("Index\tFile\tNumber of Regions\t Avg size\n");
    let mut buffer = Vec::new();
    buffer.write_all((&initial_line).as_ref()).unwrap();

    let mut total_regions = 0;
    let mut total_avg_size = 0.0;

    for i in 0..n_files {
        let file_path = &all_bed_files[i].to_str().unwrap();

        let file_path = Path::new(file_path);
        //let filename = file_path.rsplitn(1, '/').next().unwrap_or(file_path);

        let filename = file_path.file_name().unwrap();
        let filename = filename.to_str().unwrap();

        total_regions += nr[i];

        total_avg_size += avg[i] as f32;

        // Write file summary
        //writeln!(fpi, "{} \t {} \t {} \t {}", i, filename, nr[i], avg[i] / nr[i]).expect("Couldn't write to file");
        let current_line = format!(
            "{} \t {} \t {} \t {} \n",
            i,
            filename,
            nr[i],
            avg[i] / nr[i]
        );
        buffer.write_all((&current_line).as_ref()).unwrap();
    }

    file.write_all(&buffer).unwrap();

    // Sort tile data and save into single files per ctg
    igd_save_db(&mut igd, output_path, db_output_name);

    igd.total_regions = total_regions;
    igd.total_average = total_avg_size;
    igd.average_length = total_avg_size / total_regions as f32;

    let save_path = format!("{}{}{}", output_path, db_output_name, ".igd");
    println!("IGD saved to: {}", save_path);
    println!(
        "Total Intervals: {}, l_avg: {}",
        igd.total_regions, igd.average_length
    );
    println!("nctg:{}  nbp:{}", igd.nctg, igd.nbp);

    igd // return for testing purposes
}

/// Saves the primary .igd database file by reading the temp_tiles, sorting them, and then writing the sorted tiles to disk.
pub fn igd_save_db(igd: &mut igd_t, output_path: &String, db_output_name: &String) {
    //println!("HELLO from igd_save_db");
    // this is the igd_save func from the original c code

    // sprintf(idFile, "%s%s%s_%i", oPath, "data0/", ctg->name, j);
    let save_path = format!("{}{}{}", output_path, db_output_name, ".igd");
    let parent_path = save_path.clone();

    let path = std::path::Path::new(&parent_path).parent().unwrap();
    let result = create_file_with_parents(path);

    match result {
        Ok(_file) => (),
        Err(err) => println!("Error creating file: {}", err),
    }

    let mut main_db_file = OpenOptions::new()
        .create(true) // Create the file if it doesn't exist
        .append(true) // Append data to the existing file if it does exist
        .open(save_path)
        .unwrap();

    let mut buffer = Vec::new();

    buffer.write_all(&igd.nbp.to_le_bytes()).unwrap();
    buffer.write_all(&igd.gType.to_le_bytes()).unwrap();
    buffer.write_all(&igd.nctg.to_le_bytes()).unwrap();

    for i in 0..igd.nctg {
        let idx = i.clone() as usize;
        let current_ctg = &igd.ctg[idx];

        buffer.write_all(&current_ctg.mTiles.to_le_bytes()).unwrap();

        //println!("writing current_ctg.mTile to databse: {} ", current_ctg.mTiles);
    }

    for i in 0..igd.nctg {
        let idx = i.clone() as usize;
        let current_ctg = &igd.ctg[idx];

        //let j = igd.nctg;

        let n = current_ctg.mTiles;

        //println!("iterating current_ctg.mTile to database: {} ", current_ctg.mTiles);

        for j in 0..n {
            let jdx = j.clone() as usize;

            if current_ctg.gTile[jdx].nCnts != 0 {
                //println!(" nCnts >0:  {} > 0, contig number: {}, mTile number: {}", current_ctg.gTile[jdx].nCnts, i ,j);
            }
            buffer
                .write_all(&current_ctg.gTile[jdx].nCnts.to_le_bytes())
                .unwrap();
            //}
        }
    }

    for i in 0..igd.nctg {
        let idx = i.clone() as usize;
        let current_ctg = &igd.ctg[idx];

        let mut name_bytes = current_ctg.name.as_bytes().to_vec();

        //40 bytes might actually be overkill?
        name_bytes.resize(MAX_CHROM_NAME_LEN, 0);

        //let len = std::cmp::min(name_bytes.len(), MAX_CHROM_NAME_LEN);
        //buffer.write_all(&name_bytes[..len]).unwrap();

        buffer.write_all(&name_bytes).unwrap();

        //println!("writing chromosome name, {}", current_ctg.name);
        //buffer.write_all((&current_ctg.name).as_ref()).unwrap();
    }

    main_db_file.write_all(&buffer).unwrap();

    //2. Sort and save tiles data

    let _k: i32;

    for i in 0..igd.nctg {
        let idx = i.clone() as usize;

        let current_ctg = &mut igd.ctg[idx];
        let n = current_ctg.mTiles;
        //println!("\ndebug mTiles for current contig: {}", current_ctg.mTiles);
        for j in 0..n {
            let jdx = j.clone() as usize;

            //current tile
            let q = &mut current_ctg.gTile[jdx];

            let nrec = q.nCnts;

            if nrec > 0 {
                // println!("nrec greater than 0: {}   Here is j index: {}", nrec, j);
                let save_path = format!(
                    "{}{}{}_{}{}",
                    output_path, "data0/", current_ctg.name, j, ".igd"
                );
                //println!("DEBUG retrieved saveT path:{}", save_path);
                let parent_path = save_path.clone();

                let path = std::path::Path::new(&parent_path);

                let mut temp_tile_file = match OpenOptions::new()
                    .create(true)
                    .append(true)
                    .read(true)
                    .open(path)
                {
                    Ok(temp_tile_file) => temp_tile_file,
                    Err(err) => {
                        println!("Error opening file: {}", err);
                        return;
                    }
                };

                // TODO we should delete the temp files after processing...
                //println!(" Reading from tempfile {:?}", temp_tile_file);

                // Read from Temp File

                let mut gdata: Vec<gdata_t> = Vec::new();
                //
                loop {
                    //TODO check that 16 is the right value when reading back the gdata_t structs
                    let mut buf = [0u8; 16];

                    let n = temp_tile_file.read(&mut buf).unwrap();

                    if n == 0 {
                        //println!("Breaking loop while reading tempfile");
                        break;
                    } else if n != 16 {
                        //panic!("Cannot read temp file.");
                        return;
                    }

                    let mut rdr = &buf[..] as &[u8];
                    let idx = rdr.read_i32::<LittleEndian>().unwrap();
                    let start = rdr.read_i32::<LittleEndian>().unwrap();
                    let end = rdr.read_i32::<LittleEndian>().unwrap();
                    let value = rdr.read_i32::<LittleEndian>().unwrap();

                    //println!("Looping through g_datat in temp files\n");
                    //println!("idx: {}  start: {} end: {}\n", idx,start,end);

                    gdata.push(gdata_t {
                        idx: idx,
                        start,
                        end,
                        value,
                    });
                }

                // Sort Data
                gdata.sort_by_key(|d| d.start); // Sort by start value

                // Write to database after sorting
                let mut temp_buffer = Vec::new();

                for data in gdata {
                    temp_buffer.write_all(&data.idx.to_le_bytes()).unwrap();
                    temp_buffer.write_all(&data.start.to_le_bytes()).unwrap();
                    temp_buffer.write_all(&data.end.to_le_bytes()).unwrap();
                    temp_buffer.write_all(&data.value.to_le_bytes()).unwrap();
                }

                let _ = main_db_file.write_all(&temp_buffer);
            }

            //q.nCnts = 0;
        }
    }

    //file.write_all(&buffer).unwrap();
}

/// Saves temporary tiles to disc to later be sorted before collating into main .igd file
pub fn igd_saveT(igd: &mut igd_t, output_file_path: &String) {
    //println!("HELLO from igd_saveT");

    // From OG COde:
    // TEMPORARILY save/append tiles to disc, add cnts to Cnts; reset tile.gList

    let mut nt = 0;
    //println!("Number of contigs to be saved {}", igd.nctg);
    for i in 0..igd.nctg {
        let idx = i.clone() as usize;
        let idx_2 = idx;
        let current_ctg = &mut igd.ctg[idx_2];
        nt = nt + current_ctg.mTiles;

        //println!("Number of mTiles to be saved {}", current_ctg.mTiles);
        for j in 0..current_ctg.mTiles {
            let jdx = j.clone() as usize;
            let jdx_2 = jdx;

            let current_tile = &mut current_ctg.gTile[jdx_2];

            if current_tile.ncnts > 0 {
                // Construct specific temp file on disk using this information

                // OG code
                // sprintf(idFile, "%s%s%s_%i", oPath, "data0/", ctg->name, j);
                let save_path = format!(
                    "{}{}{}_{}{}",
                    output_file_path, "data0/", current_ctg.name, j, ".igd"
                );

                let parent_path = save_path.clone();
                //println!("Saving saveT path, because current_tile.ncnts > 0:{} {}", current_tile.ncnts,save_path);

                //println!("{}", save_path);

                //todo this needs to create the path if it does not already exist!!!

                let path = std::path::Path::new(&parent_path).parent().unwrap();
                let result = create_file_with_parents(path);

                match result {
                    Ok(_file) => (),
                    Err(err) => println!("Error creating file: {}", err),
                }

                let mut file = OpenOptions::new()
                    .create(true) // Create the file if it doesn't exist
                    .append(true) // Append data to the existing file if it does exist
                    .open(save_path)
                    .unwrap();

                // Because gList is a Vector of structs, we must take each field
                // and convert it to byte representation before writing to a file...
                let mut buffer = Vec::new();
                for data in &current_tile.gList[..current_tile.ncnts as usize] {
                    buffer.write_all(&data.idx.to_le_bytes()).unwrap();
                    buffer.write_all(&data.start.to_le_bytes()).unwrap();
                    buffer.write_all(&data.end.to_le_bytes()).unwrap();
                    buffer.write_all(&data.value.to_le_bytes()).unwrap();
                }
                file.write_all(&buffer).unwrap();

                current_tile.nCnts = current_tile.nCnts + current_tile.ncnts;

                if current_tile.ncnts > 8 {
                    current_tile.mcnts = 8;
                } else {
                    current_tile.mcnts = 2;
                }
                current_tile.ncnts = 0;
            }
        }
    }
    println!(
        "Temporary Tiles:\n nCtgs (igd.nctg): {}, nRegions (igd.total): {}, nTiles (nt): {}",
        igd.nctg, igd.total, nt
    );
    igd.total = 0; // batch total
}

/// Creates file and any parent directories if they do not already exist.
fn create_file_with_parents(path: &Path) -> Result<File, Error> {
    // Create all parent directories if they don't exist (ignore errors)
    let _ = create_dir_all(path); // Discard the result (success or error)

    // Open the file for creation or append, ignoring errors if it exists
    let file = OpenOptions::new().create(true).append(true).open(path);

    match file {
        Ok(file) => {
            println!("File created or opened successfully!");
            Ok(file)
        }
        Err(_) => Ok(File::open(path).unwrap_or_else(|_| File::create(path).unwrap())), // Handle existing file or create new one
    }
}

/// Adds genomic interval to the igd struct
pub fn igd_add(
    igd: &mut igd_t,
    hash_table: &mut HashMap<String, i32>,
    chrm: String,
    start: i32,
    end: i32,
    v: i32,
    idx: usize,
) {
    //Add an interval
    // og code: layers: igd->ctg->gTile->gdata(list)
    //println!("HELLO from igd_add");
    // println!(
    //     "Entering IGD ADD Chrm {}, start {}, end {}, v {}, idx {}",
    //     chrm, start, end, v, idx
    // );
    if start >= end {
        // println!(
        //     "Start: {0} greater than End: {1}, returning from igd_add",
        //     start, end
        // );
        return;
    }
    let _absent: i32;
    let _i: i32;

    // Cloning chrm String because the hash table will own the key after insertion
    let key = chrm.clone();

    let n1 = start / igd.nbp;
    let n2 = (end - 1) / igd.nbp;

    // // create hash table
    // let mut hash_table: HashMap<String, i32> = HashMap::new();

    let key_check = hash_table.contains_key(&key);

    if key_check == false {
        // println!(
        //     "Key does not exist in hash map, creating for {}",
        //     key.clone()
        // );

        // Insert key and value (igd.nctg)
        hash_table.insert(key.clone(), igd.nctg);
        igd.nctg += 1;

        // initialize ctg
        let mut p = ctg_t::new();
        p.name = chrm;
        p.mTiles = 1 + n2;
        //p.gTile original code mallocs mTiles*sizeof title_t
        // however in Rust, structs have 0 size: https://doc.rust-lang.org/nomicon/exotic-sizes.html#zero-sized-types-zsts
        //p.gTile = Vec::with_capacity((p.mTiles as usize)*size_of(tile_t()));
        p.gTile = Vec::with_capacity(p.mTiles as usize);

        for _i in 0..p.mTiles {
            //println!("iterating of p.Mtiles");
            let mut new_tile: tile_t = tile_t::new();

            new_tile.ncnts = 0; //each batch
            new_tile.nCnts = 0; //total
            new_tile.mcnts = 2;

            for _j in 0..new_tile.mcnts {
                new_tile.gList.push(gdata_t::new());
            }

            p.gTile.push(new_tile);
        }

        igd.ctg.push(p);
    }
    // else {
    //     println!(
    //         "Key exists in hash map, skipping creation, key: {}",
    //         key.clone()
    //     )
    // }

    // Retrieve values from Hash Map

    let keycloned = key.clone();

    let index = hash_table.get(&keycloned).unwrap();
    let cloned_index = index.clone();

    let p = &mut igd.ctg[cloned_index as usize];

    if n2 + 1 >= p.mTiles {
        //println!("TRUE:{} vs {}", (n2 + 1), p.mTiles.clone());
        let tt = p.mTiles;

        p.mTiles = n2 + 1;
        p.gTile
            .resize(p.mTiles as usize, crate::igd::create::tile_t::default());
        // original code: p->gTile = realloc(p->gTile, p->mTiles*sizeof(tile_t));
        // Supposedly we may not need to do this ...  p.gTile = Vec::resize()   ???

        for i in tt..p.mTiles {
            let idx = i.clone() as usize;
            let idx_2 = idx as usize;

            let existing_tile: &mut tile_t = &mut p.gTile[idx_2];

            existing_tile.ncnts = 0;
            existing_tile.nCnts = 0;
            existing_tile.mcnts = 2;
            for _j in 0..existing_tile.mcnts {
                existing_tile.gList.push(gdata_t::new());
            }
        }
    }

    for i in n1..=n2 {
        //println!("iterating n1..n2");
        //println!("Adding data elements, iteration: {}", i);
        //this is inclusive of n1 and n2
        // Get index as usize
        let idx_1 = i.clone() as usize;
        let idx_2 = idx_1 as usize;
        // get the tile for the contig
        let existing_tile: &mut tile_t = &mut p.gTile[idx_2];

        if existing_tile.ncnts == existing_tile.mcnts {
            // println!(
            //     "DEBUG Existing tile: ncnts == mcnts {} vs {}",
            //     existing_tile.ncnts, existing_tile.mcnts
            // );

            // Expand number of elements by doubling, og used a realloc macro to achieve this...
            existing_tile
                .gList
                .resize((existing_tile.mcnts * 2) as usize, Default::default());
            existing_tile.mcnts *= 2;
        }

        let tile_idx = existing_tile.ncnts.clone() as usize;
        let gdata = &mut existing_tile.gList[tile_idx];
        existing_tile.ncnts = existing_tile.ncnts + 1;

        gdata.start = start;
        gdata.end = end;
        gdata.value = v;
        //println!("Adding to igd, start {}, idx {}", start,idx);
        gdata.idx = idx as i32;

        igd.total += 1;
    }

    //println!("DEBUG: Here is igd.total:  {}", igd.total);

    //println!("Finished from igd_add");
    return;
}

#[derive(PartialEq)] // So that we can do comparisons with equality operator
pub enum ParseBedResult {
    Str(String),
    Int(i32),
}

/// Reads bed file, returning contig and modifying borrowed start and end coordinate
pub fn parse_bed(line: &String, start: &mut i32, end: &mut i32, score: &mut i32) -> Option<String> {
    //println!("HERE IS THE LINE TO PARSE: {}", line);
    let mut fields = line.split('\t');
    // Get the first field which should be chromosome.
    let ctg = fields.next()?; // Why is ctg used as variable name in og code?
                              //println!("GOT CHR: {}", ctg);
                              // Parse 2nd and 3rd string as integers or return -1 if failure
    let st = fields
        .next()
        .and_then(|s| s.parse::<i32>().ok())
        .unwrap_or(-1);
    //println!("GOT st: {}", st);
    let en = fields
        .next()
        .and_then(|s| s.parse::<i32>().ok())
        .unwrap_or(-1);
    //println!("GOT en: {}", en);

    let _ = fields.next(); // throwaway the 4th bed column

    // 5th column is score, if it exists, else -1
    let va = fields
        .next()
        .and_then(|s| s.parse::<i32>().ok())
        .unwrap_or(-1);

    if !ctg.starts_with("chr") || ctg.len() >= 40 || en <= 0 {
        println!("RETURNING NONE, {}", ctg);
        return None;
    }

    *start = st;
    *end = en;
    *score = va;

    //println!("SUCCESSFULLY FINISHING PARSE");
    Some(ctg.parse().unwrap())
}

use crate::common::consts::{BED_FILE_EXTENSION, IGD_FILE_EXTENSION};
use crate::common::utils::get_dynamic_reader;
use crate::igd::create::{gdata0_t, gdata_t, parse_bed};

use byteorder::{LittleEndian, ReadBytesExt};
use clap::ArgMatches;

use std::collections::HashMap;

use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, Error, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

#[derive(Default)]
pub struct igd_t_from_disk {
    // int32_t nFiles;
    // info_t *finfo;
    // char fname[64];
    // int32_t nbp, gType, nCtg;			//data type: 0, 1, 2 etc; size differs
    // char **cName;						//name of ctgs
    // int32_t *nTile;						//num of tiles in each ctg
    // int32_t **nCnt;						//num of counts in each tile
    // int64_t **tIdx;
    pub nFiles: i32,
    pub file_info: Vec<info_t>,
    pub filename: String,
    pub nbp: i32,   //data type: 0, 1, 2 etc; size differs
    pub gType: i32, //data type: 0, 1, 2 etc; size differs
    pub nCtg: i32,  //data type: 0, 1, 2 etc; size differs
    // Original code uses pointer to pointers
    pub cName: Vec<String>,
    pub nTile: Vec<i32>,

    //pub nCnt: i32,
    pub nCnt: Vec<Vec<i32>>,
    //pub tIdx: i32,
    pub tIdx: Vec<Vec<i64>>,
}

impl igd_t_from_disk {
    /// Constructs new instance of igd_t_from_disk
    pub fn new() -> Self {
        Self::default()
    }
}

#[derive(Default)]
pub struct info_t {
    pub fileName: String, //dataset file
    pub nr: i32,          // number of regions in dataset
    pub md: f64,          // average width of the regions
}

// typedef struct{
//     char* fileName;				//dataset file
//     int32_t nr;					//number regions/dataset
//     double md;    				//average width of the regions
// } info_t;

/// Searches IGD database
pub fn igd_get_search_matches(matches: &ArgMatches) {
    let database_path = matches
        .get_one::<String>("database")
        .expect("Database path is required");

    let query = matches
        .get_one::<String>("query")
        .expect("Query bed file path is required");

    igd_search(database_path, query).expect("Error:");
}

#[allow(unused_variables)]
pub fn igd_search(database_path: &String, query_file_path: &String) -> Result<Vec<String>, String> {
    // First check that BOTH the igd database and the query are the proper file types
    // else raise error

    let mode = 1;
    let mut hash_table: HashMap<String, i32> = HashMap::new();

    let mut final_string_vec = Vec::new();

    match check_file_extension(database_path, IGD_FILE_EXTENSION) {
        Ok(_) => (),
        Err(e) => return Err(e),
    }

    match check_file_extension(query_file_path, BED_FILE_EXTENSION) {
        Ok(_) => (),
        Err(e) => {
            return Err(e);
        }
    }

    //println!("\n {} \n {}", database_path, query_file_path);

    //Get file info from the associated TSV

    //DEBUG
    //read_and_print_numbers(database_path.as_str());

    // Create IGD Struct from database
    let mut IGD: igd_t_from_disk =
        get_igd_info(database_path, &mut hash_table).expect("Could not open IGD");

    // If query "-q" used set to mode 1

    let tsv_path = get_tsv_path(database_path).unwrap();

    get_file_info_tsv(tsv_path, &mut IGD).unwrap(); //sets igd.finfo

    let nfiles = IGD.nFiles;
    //let mut hits: Vec<i64> = Vec::with_capacity(nfiles as usize);
    let mut hits: Vec<i64> = vec![0; nfiles as usize];

    //Open IGD database

    match mode {
        1 => {
            // Querying a bedfile

            if IGD.gType == 0 {
                //getOverlaps0(query_file_path, hits);
                println!("gType = 0");
            } else {
                getOverlaps(
                    &mut IGD,
                    database_path,
                    query_file_path,
                    &mut hits,
                    &mut hash_table,
                );
            }

            println!("index\t number of regions\t number of hits\t File_name");
            let format_string = format!("index\t number of regions\t number of hits\t File_name");
            final_string_vec.push(format_string);

            let mut total: i64 = 0;
            for (i, hit) in hits.iter().enumerate() {
                if *hit > 0 {
                    println!(
                        "{}\t{}\t{}\t{}",
                        i, IGD.file_info[i].nr, hit, IGD.file_info[i].fileName
                    );
                    let format_string = format!(
                        "{}\t{}\t{}\t{}",
                        i, IGD.file_info[i].nr, hit, IGD.file_info[i].fileName
                    );
                    final_string_vec.push(format_string);
                }
                total += hit;
            }

            println!("Total: {}", total);
        }

        _ => {
            println!("Invalid mode selected, exiting");
            return Err(String::from("Invalid mode selected, exiting"));
        }
    }

    println!("FINISHED");

    Ok(final_string_vec)
}
#[allow(unused_variables)]
pub fn getOverlaps(
    IGD: &igd_t_from_disk,
    database_path: &String,
    query_file: &String,
    hits: &mut Vec<i64>,
    hash_table: &mut HashMap<String, i32>,
) -> i32 {
    let mut start = 0;
    let mut end = 0;
    let mut va = 0;
    let mut ols = 0;

    let mut preChr = -6;
    let mut preIdx = -8;

    // Get Reader for QUERY FILE dynamically
    let path = Path::new(query_file);
    let reader = get_dynamic_reader(path).unwrap();

    // Also get Reader for database file (.igd)
    let parent_path = database_path.clone();

    let dbpath = std::path::Path::new(&parent_path);

    let db_file = OpenOptions::new()
        .create(true)
        .append(true)
        .read(true)
        .open(dbpath)
        .unwrap();

    let mut db_reader = BufReader::new(db_file);

    for line in reader.lines() {
        let line = line.unwrap();
        let ctg = parse_bed(&line, &mut start, &mut end, &mut va);
        // if it parses, add it to collected lines, increment ix
        match ctg {
            Some(ctg) => {
                //println!("ctg successfully parsed {}", ctg);

                let nl = get_overlaps(
                    &IGD,
                    ctg,
                    start,
                    end,
                    hits,
                    hash_table,
                    &mut preChr,
                    &mut preIdx,
                    path,
                    &mut db_reader,
                );

                ols += nl;
            }
            None => continue,
        }
    }

    return ols;
}

// trait ReaderSeeker: Read + Seek {
//
// }

#[allow(unused_variables)]
fn get_overlaps(
    IGD: &igd_t_from_disk,
    ctg: String,
    query_start: i32,
    query_end: i32,
    hits: &mut Vec<i64>,
    hash_table: &mut HashMap<String, i32>,
    preChr: &mut i32,
    preIdx: &mut i32,
    query_path: &Path,
    db_reader: &mut BufReader<File>,
) -> i32 {
    //println!("get overlaps main func");

    let ichr = get_id(ctg, hash_table);
    //println!("ichr from get_overlaps {}", ichr);

    if ichr < 0 {
        return 0;
    }

    // Define Boundary
    let n1 = query_start / IGD.nbp;
    let mut n2 = (query_end - 1) / IGD.nbp;
    let i: i32;
    let j: i32;
    let ni: i32;

    //int32_t tE, tS, tL, tR, tM, tmpi, tmpi1, mlen, mTile = IGD->nTile[ichr]-1;
    //int32_t nols = 0;

    let tE: i32;
    let mut tS: i32;
    let mut tL: i32;
    let mut tR: i32;
    let mut tM: i32;
    let mut tmpi: i32;
    let mut tmpi1: i32;
    let mlen: i32;

    let nols = 0; //number of overlaps

    let mTile = IGD.nTile[ichr as usize] - 1;

    if n1 > mTile {
        return 0;
    }

    // Min between n2 and mTile
    n2 = if n2 < mTile { n2 } else { mTile };

    tmpi = IGD.nCnt[ichr as usize][n1 as usize];
    tmpi1 = tmpi - 1;

    // println!(
    //     "prechr and preidx at the  begining of get_overlaps {}  {} \n",
    //     preChr, preIdx
    // );

    if tmpi > 0 {
        // println!(
        //     "n1 != *preIdx || ichr!= *preChr {} vs {}  {} vs {} \n",
        //     n1, preIdx, ichr, preChr
        // );

        //println!("Seek start here: {}",IGD.tIdx[ichr as usize][n1 as usize]);
        //let ichr = 1;
        db_reader
            .seek(SeekFrom::Start(IGD.tIdx[ichr as usize][n1 as usize] as u64))
            .unwrap();

        let mut gData: Vec<gdata_t> = Vec::new();
        for j in 0..tmpi {
            gData.push(gdata_t::default())
        }
        //let mut gData: Vec<gdata_t> = Vec::with_capacity(tmpi as usize);

        for i in 0..tmpi {
            let mut buf = [0u8; 16];

            let n = db_reader.read(&mut buf).unwrap();

            if n == 0 {
                //println!("Breaking loop while reading tempfile");
                break;
            } else if n != 16 {
                //panic!("Cannot read temp file.");
                break;
            }

            let mut rdr = &buf[..] as &[u8];
            let idx = rdr.read_i32::<LittleEndian>().unwrap();
            let start = rdr.read_i32::<LittleEndian>().unwrap();
            let end = rdr.read_i32::<LittleEndian>().unwrap();
            let value = rdr.read_i32::<LittleEndian>().unwrap();

            //println!("for tmpi>0 where tmpi = {}", tmpi);
            //println!("Looping through g_datat in temp files\n");
            //println!("idx: {}  start: {} end: {}\n", idx,start,end);

            gData[i as usize] = gdata_t {
                idx: idx,
                start,
                end,
                value,
            };

            *preIdx = n1;
            *preChr = ichr;
        }

        // check this code block. original code has outside this first check but that would potentially cause access to wrong
        // object in memory if it was not de-allocated?

        if query_end > gData[0].start {
            // sorted by start
            //println!("n1 != *preIdx || ichr != *preChr query_end > gData[0].start:  {} > {}", query_end,gData[0].start);
            // find the 1st rs<qe
            tL = 0;
            tR = tmpi1;

            while tL < tR - 1 {
                tM = (tL + tR) / 2; //result: tR=tL+1, tL.s<qe
                                    //println!("What is tM? {}", tM);
                if gData[tM as usize].start < query_end {
                    tL = tM; //right side
                } else {
                    tR = tM; //left side
                }
            }
            if gData[tR as usize].start < query_end {
                tL = tR;
            }
            //--------------------------
            for i in (0..=tL).rev() {
                //println!("Countdownfrom TL");
                // count down from tL (inclusive to tL)
                //println!("iterate over i: {} from tL {}", i, tL);
                //println!("gdata[i].end {} vs query start {}",gData[i as usize].end,query_start);
                if gData[i as usize].end > query_start {
                    //println!("ADDING TO HITS");
                    //println!(" > gData[i].end > query_start  {} > {}", gData[i as usize].end, query_start);
                    hits[gData[i as usize].idx as usize] = hits[gData[i as usize].idx as usize] + 1;
                }
            }
        }

        if n2 > n1 {
            //println!("n2>n1  {} vs {} ", n2, n1);

            let mut bd = IGD.nbp * (n1 + 1); // only keep the first
            for j in (n1 + 1)..=n2 {
                //n2 inclusive
                tmpi = IGD.nCnt[ichr as usize][j as usize];
                tmpi1 = tmpi - 1;
                if tmpi > 0 {
                    let mut gData: Vec<gdata_t> = Vec::with_capacity(tmpi as usize);

                    if j != *preIdx || ichr != *preChr {
                        // println!(
                        //     "j != *preIdx || ichr!= *preChr {} vs {}  {} vs {} \n",
                        //     j, preIdx, ichr, preChr
                        // );

                        db_reader
                            .seek(SeekFrom::Start(IGD.tIdx[ichr as usize][j as usize] as u64))
                            .unwrap();

                        for i in 0..tmpi {
                            let mut buf = [0u8; 16];

                            let n = db_reader.read(&mut buf).unwrap();

                            if n == 0 {
                                //println!("Breaking loop while reading tempfile");
                                break;
                            } else if n != 16 {
                                //panic!("Cannot read temp file.");
                                break;
                            }

                            let mut rdr = &buf[..] as &[u8];
                            let idx = rdr.read_i32::<LittleEndian>().unwrap();
                            let start = rdr.read_i32::<LittleEndian>().unwrap();
                            let end = rdr.read_i32::<LittleEndian>().unwrap();
                            let value = rdr.read_i32::<LittleEndian>().unwrap();

                            //println!("Looping through g_datat in temp files\n");
                            // println!("idx: {}  start: {} end: {}\n", idx,start,end);

                            gData.push(gdata_t {
                                idx: idx,
                                start,
                                end,
                                value,
                            });

                            *preIdx = j;
                            *preChr = ichr;
                        }
                    }

                    if query_end > gData[0].start {
                        //println!("n2>n1 query_end > gData[0].start:  {} > {}", query_end,gData[0].start);
                        tS = 0;

                        while tS < tmpi && gData[tS as usize].start < bd {
                            //query start < bd
                            tS = tS + 1;
                        }

                        tL = 0;
                        tR = tmpi1;

                        while tL < tR - 1 {
                            //result: tR=tL+1, tL.s<qe
                            tM = (tL + tR) / 2;

                            if gData[tM as usize].start < query_end {
                                // right side
                                tL = tM;
                            } else {
                                tR = tM;
                            }
                        }
                        if gData[tR as usize].start < query_end {
                            tL = tR;
                        }
                        //--------------------------
                        for i in (tS..=tL).rev() {
                            //println!("* gdata[i].end {} vs query start {}",gData[i as usize].end,query_start);
                            if gData[i as usize].end > query_start {
                                //println!("* gData[i].end > query_start  {} > {}", gData[i as usize].end, query_start);
                                hits[gData[i as usize].idx as usize] =
                                    hits[gData[i as usize].idx as usize] + 1;
                            }
                        }
                    }
                }
                bd = bd + IGD.nbp;
            }
        }
    }
    //println!("here are the hits {:?}", hits);
    return nols; //TODO this is from the original code but its not actually being used for anything. hits vec IS the main thing.
}

#[allow(unused_variables)]
fn get_id(ctg: String, hash_table: &mut HashMap<String, i32>) -> i32 {
    let key_check = hash_table.contains_key(&ctg);

    if key_check == false {
        -1
    } else {
        let value = hash_table.get(&ctg).unwrap();
        value.clone()
    }
}

// #[allow(unused_variables)]
// fn getOverlaps0(p0: &String, p1: Vec<i64>) {
//     println!("getoverlaps0");
// }

/// Given an igd path, simple give the .tsv path that is parallel to the  .igd path
pub fn get_tsv_path(igd_path: &str) -> Option<PathBuf> {
    let igd_path = Path::new(igd_path);
    let stem = igd_path.file_stem()?;
    let mut tsv_path = igd_path.with_file_name(stem);
    tsv_path.set_extension("tsv");
    Some(tsv_path)
}
// fn read_and_print_numbers(filename: &str) -> std::io::Result<()> {
//     // Just a debug function to determine what was actually written to a file.
//     let file = File::open(filename)?;
//     let mut reader = BufReader::new(file);
//
//     let mut buffer = [0u8; 4];
//
//     loop {
//         match reader.read_exact(&mut buffer) {
//             Ok(_) => {
//                 let number = u32::from_le_bytes(buffer);
//                 println!("{}", number);
//             }
//             Err(ref e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
//             Err(e) => return Err(e),
//         }
//     }
//
//     Ok(())
// }
#[allow(unused_variables)]
pub fn get_igd_info(
    database_path: &String,
    hash_table: &mut HashMap<String, i32>,
) -> Result<igd_t_from_disk, Error> {
    //println!("hello from get_igd_info");

    let mut igd = igd_t_from_disk::new();

    // Open file

    let parent_path = database_path.clone();

    let dbpath = std::path::Path::new(&parent_path);

    let temp_tile_file = match OpenOptions::new()
        .create(true)
        .append(true)
        .read(true)
        .open(dbpath)
    {
        Ok(temp_tile_file) => temp_tile_file,
        Err(err) => {
            println!("Error opening file: {}", err);
            return Err(err);
        }
    };

    let mut reader = BufReader::new(temp_tile_file);

    // TODO is this the correct buffer size given the way it was written to disk?
    //let mut buffer = [0u8; std::mem::size_of::<i32>()];
    let mut buffer = [0u8; 4];

    reader.read_exact(&mut buffer)?;
    let nbp = i32::from_le_bytes(buffer);
    reader.read_exact(&mut buffer)?;
    let gType = i32::from_le_bytes(buffer);
    reader.read_exact(&mut buffer)?;
    let nCtg = i32::from_le_bytes(buffer);

    println!("Found:\n nbp:{} gtype: {} nCtg: {}", nbp, gType, nCtg);

    igd.nbp = nbp;
    igd.gType = gType;
    igd.nCtg = nCtg;

    //println!("Found:\n nbp:{} gtype: {} nCtg: {}", nbp, gType, nCtg);
    let gdsize = if gType == 0 {
        std::mem::size_of::<gdata0_t>()
    } else {
        std::mem::size_of::<gdata_t>()
    };

    let tileS = igd.nCtg;
    let m = igd.nCtg; //the idx of a tile in the chrom

    let mut n_Tile: Vec<i32> = Vec::with_capacity(m as usize);
    for _ in 0..m {
        reader.read_exact(&mut buffer)?;
        let tile_value = i32::from_le_bytes(buffer);
        //n_Tile.push(reader.read_i32::<LittleEndian>()?);
        n_Tile.push(tile_value);
    }

    igd.nTile = n_Tile.clone();

    //println!("here is m: {}",m);
    // This calculation is from og code.
    // TODO The above buffer size might throw it off and should be double checked
    let mut chr_loc = (12 + 44 * m) as i64; // originally this is the header size in bytes

    //println!("Initial chr loc: {}", chr_loc);

    for n in 0..m {
        chr_loc = chr_loc + (n_Tile[n as usize] as i64) * 4;
    }

    //println!("Skip to new chr loc: {}", chr_loc);

    let mut nCnt: Vec<Vec<i32>> = Vec::new();
    for _ in 0..n_Tile.len() {
        nCnt.push(Vec::new());
    }

    let mut tIdx: Vec<Vec<i64>> = Vec::new();
    for _ in 0..n_Tile.len() {
        tIdx.push(Vec::new());
    }

    // TODO this block may be causing errors downstream when calculating overlaps

    for i in 0..m {
        let k = igd.nTile[i as usize];

        //println!("here is idx for i and k: {} {} ", i, k);
        let mut cnt = vec![0; k as usize]; //original code used calloc which does initialize arrays with 0's
        for kdx in 0..k {
            cnt[kdx as usize] = reader.read_i32::<LittleEndian>()?;
        }

        nCnt[i as usize] = cnt;

        let idx = vec![0; k as usize];

        tIdx[i as usize] = idx;
        tIdx[i as usize][0] = chr_loc;

        for j in 1..k {
            tIdx[i as usize][j as usize] = tIdx[i as usize][j as usize - 1]
                + (nCnt[i as usize][j as usize - 1] as i64) * gdsize as i64;
        }

        chr_loc = tIdx[i as usize][k as usize - 1]
            + nCnt[i as usize][k as usize - 1] as i64 * gdsize as i64;
        //println!("Skip to new chr loc after m_tile iteration: {}", chr_loc);
    }

    igd.nCnt = nCnt;
    igd.tIdx = tIdx;

    // Read cName

    let mut c_name = Vec::with_capacity(m as usize);
    for _ in 0..m {
        let mut buf = [0u8; 40];
        reader.read_exact(&mut buf)?;
        //println!("Raw bytes: {:x?}", buf);
        let name = String::from_utf8(buf.to_vec()).unwrap(); // TODO assumes utf 8, add handling for error later
        let name = name.trim_matches('\0');
        c_name.push(String::from(name)); // Maybe just have this be a String and not a vec<String>?
    }

    igd.cName = c_name.clone();

    // for name in c_name {
    //     println!("Retrieved chrom name (cName):  {}", name);
    // }

    // Place values in hash map
    for (i, name) in igd.cName.iter().enumerate() {
        hash_table.insert(name.to_string(), i as i32);
    }

    return Ok(igd);
}

pub fn get_file_info_tsv(tsv_path: PathBuf, igd: &mut igd_t_from_disk) -> Result<(), Error> {
    let path = Path::new(&tsv_path);

    let tsv_file = match OpenOptions::new()
        .create(true)
        .append(true)
        .read(true)
        .open(path)
    {
        Ok(temp_tile_file) => temp_tile_file,
        Err(err) => {
            println!("Error opening file: {}", err);
            return Err(err);
        }
    };

    let reader = BufReader::new(tsv_file);

    let mut lines = reader.lines();
    // Skip header
    lines.next();

    let mut infos: Vec<info_t> = Vec::new();

    let mut count = 0;

    for line in lines {
        //println!("Reading tsv lines...");
        count = count + 1;
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();

        let info = info_t {
            fileName: fields[1].to_string(),
            nr: fields[2].to_string().as_str().trim().parse().unwrap(),
            md: fields[3].to_string().as_str().trim().parse().unwrap(),
        };
        infos.push(info);
    }

    igd.nFiles = count;

    igd.file_info = infos;

    Ok(())
}

fn check_file_extension(path: &str, expected_extension: &str) -> Result<(), String> {
    let path = Path::new(path);
    let actual_extension = path
        .extension()
        .and_then(|ext| ext.to_str())
        .ok_or_else(|| format!("Invalid file path: {}", path.display()))?;

    if actual_extension != expected_extension {
        return Err(format!(
            "Incorrect file extension. Expected: {}, got: {}",
            expected_extension, actual_extension
        ));
    }

    Ok(())
}

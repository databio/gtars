use crate::common::consts::{BED_FILE_EXTENSION, IGD_FILE_EXTENSION};
use crate::igd::create::{gdata0_t, gdata_t, igd_t, MAX_CHROM_NAME_LEN};
use clap::ArgMatches;
use std::fs::{create_dir_all, DirEntry, File, OpenOptions};
use std::io::{BufRead, BufReader, Error, Read, Write};
use std::path::Path;
use byteorder::{LittleEndian,ReadBytesExt};

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
    pub file_info: info_t,
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
pub fn igd_search(database_path: &String, query_file_path: &String) -> Result<(), String> {
    // First check that BOTH the igd database and the query are the proper file types
    // else raise error

    let mode = 1;

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

    println!("\n {} \n {}", database_path, query_file_path);

    //Get file info from the associated TSV

    // Create IGD Struct from database
    let IGD: igd_t_from_disk = get_igd_info(database_path).expect("Could not open IGD");

    // If query "-q" used set to mode 1

    match mode {
        1 => {}
        _ => {
            println!("Invalid mode selected, exiting");
            return Ok(());
        }
    }

    println!("FINISHED");

    Ok(())
}
#[allow(unused_variables)]
pub fn get_igd_info(database_path: &String) -> Result<igd_t_from_disk, Error> {
    println!("hello from get_igd_info");

    let mut igd = igd_t_from_disk::new();

    // Open file

    let parent_path = database_path.clone();

    let dbpath = std::path::Path::new(&parent_path);

    let mut temp_tile_file = match OpenOptions::new()
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
    let mut buffer = [0u8; std::mem::size_of::<i32>()];

    reader.read_exact(&mut buffer)?;
    let nbp = i32::from_le_bytes(buffer);
    reader.read_exact(&mut buffer)?;
    let gType = i32::from_le_bytes(buffer);

    reader.read_exact(&mut buffer)?;
    let nCtg = i32::from_le_bytes(buffer);

    //println!("Found:\n nbp:{} gtype: {} nCtg: {}", nbp,gType,nCtg);

    igd.nbp = nbp;
    igd.gType = gType;
    igd.nCtg = nCtg;

    println!("Found:\n nbp:{} gtype: {} nCtg: {}", nbp,gType,nCtg);
    let gdsize = if gType == 0 {
        std::mem::size_of::<gdata0_t>()
    } else {
        std::mem::size_of::<gdata_t>()
    };

    let tileS = igd.nCtg;
    let m = igd.nCtg; //the idx of a tile in the chrom

    let mut n_Tile: Vec<i32> = Vec::with_capacity(m as usize);
    for _ in 0..m {
        n_Tile.push(reader.read_i32::<LittleEndian>()?);
    }

    igd.nTile = n_Tile.clone();
    // reader.read_exact(&mut buffer)?;
    // let nTile = i32::from_le_bytes(buffer);
    // igd.nTile = nTile;

    // This calculation is from og code.
    // TODO The above buffer size might throw it off and should be double checked
    let mut chr_loc = (12 + 44 * m) as i64; // originally this is the header size in bytes
    for n in 0..m {
        chr_loc = chr_loc + n as i64 * 4;
    }

    let mut nCnt: Vec<Vec<i32>> = Vec::with_capacity(n_Tile.len());
    let mut tIdx: Vec<Vec<i64>> = Vec::with_capacity(n_Tile.len());

    for (i, k) in n_Tile.iter().enumerate() {

        println!("\nFrom Enumeration, here is i: {},  k {}", i,k);
        println!("From Enumeration, here is chr_loc: {}", chr_loc);
        let mut cnt = vec![0; *k as usize];
        reader.read_exact(&mut cnt)?;

        // we read as u8 and then must convert back to i32. This seems like an unecessary step if we could just do everything as either u8 or i32...
        let i32_converted_cnt =  cnt.into_iter().map(|byte| byte as i32).collect();

        nCnt.push(i32_converted_cnt);


        let mut idx = vec![0; *k as usize];

        for j in 0..*k {
            if j > 0 {
                idx[j as usize] = idx[j as usize - 1] + (nCnt[i][j as usize - 1] as i64) * (gdsize as i64);
            }

            chr_loc = chr_loc + (nCnt[i][j as usize] as i64) * (gdsize as i64);
        }

        tIdx.push(idx);


    }

    igd.nCnt = nCnt;
    igd.tIdx = tIdx;

    // More of a direct port of the C code...
    // getting tile information

    // for i in 0..m {
    //     //k = iGD->nTile[i]
    //     let i_idx = i.clone() as usize;
    //     let k = igd.nTile[i_idx].clone();
    //     println!("\n k: {:?}, chrm_loc: {}", k, chr_loc);
    //     // og code, nCnt originally
    //     // k = iGD->nTile[i];
    //     // iGD->nCnt[i] = calloc(k, sizeof(int32_t));
    //     // ni = fread(iGD->nCnt[i], sizeof(int32_t)*k, 1, fp);
    //     reader.read_exact(&mut buffer)?;
    //     let current_nCnt = i32::from_le_bytes(buffer);
    //
    //     igd.nCnt.push(current_nCnt);
    //     //println!("\n k: {:?}, chrm_loc: {}", k, chr_loc);
    //
    //     // og code
    //     // iGD->tIdx[i] = calloc(k, sizeof(int64_t));
    //     // iGD->tIdx[i][0] = chr_loc;
    //
    //     //igd.tIdx.push(Vec::from(chr_loc.clone())); // vec of vecs
    //
    //     for j in 1..k {
    //         let idx = i as usize;
    //         let jdx = j as usize;
    //
    //         //igd.tIdx[idx][jdx];
    //     }
    // }

    // Read cName

    // Read cName data
    let mut c_name = Vec::with_capacity(m as usize);
    for _ in 0..m{

        let mut buf = [0u8; 40];
        reader.read_exact(&mut buf)?;
        println!("Raw bytes: {:x?}", buf);
        let name = String::from_utf8(buf.to_vec()).unwrap(); // TODO assumes utf 8, add handling for error later
        c_name.push(name); // Maybe just have this be a String and not a vec<String>?
    }

    igd.cName = c_name.clone();

    println!("Retrieved chrom name (cName):  {:?}", c_name);



    // Place values in hash map



    return Ok(igd);
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

use crate::common::consts::{BED_FILE_EXTENSION, IGD_FILE_EXTENSION};
use crate::igd::create::{gdata0_t, gdata_t, igd_t, MAX_CHROM_NAME_LEN};
use byteorder::{LittleEndian, ReadBytesExt};
use clap::ArgMatches;
use std::collections::HashMap;
use std::fs::{create_dir_all, DirEntry, File, OpenOptions};
use std::io::{BufRead, BufReader, Error, Read, Write};
use std::path::Path;

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
    let mut hash_table: HashMap<String, i32> = HashMap::new();

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

    //DEBUG
    //read_and_print_numbers(database_path.as_str());

    // Create IGD Struct from database
    let IGD: igd_t_from_disk =
        get_igd_info(database_path, &mut hash_table).expect("Could not open IGD");

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
fn read_and_print_numbers(filename: &str) -> std::io::Result<()> {
    // Just a debug function to determine what was actually written to a file.
    let file = File::open(filename)?;
    let mut reader = BufReader::new(file);

    let mut buffer = [0u8; 4];

    loop {
        match reader.read_exact(&mut buffer) {
            Ok(_) => {
                let number = u32::from_le_bytes(buffer);
                println!("{}", number);
            }
            Err(ref e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
            Err(e) => return Err(e),
        }
    }

    Ok(())
}
#[allow(unused_variables)]
pub fn get_igd_info(
    database_path: &String,
    hash_table: &mut HashMap<String, i32>,
) -> Result<igd_t_from_disk, Error> {
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
    //let mut buffer = [0u8; std::mem::size_of::<i32>()];
    let mut buffer = [0u8; 4];

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

    println!("Found:\n nbp:{} gtype: {} nCtg: {}", nbp, gType, nCtg);
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

    // This calculation is from og code.
    // TODO The above buffer size might throw it off and should be double checked
    let mut chr_loc = (12 + 44 * m) as i64; // originally this is the header size in bytes

    for n in 0..m {
        chr_loc = chr_loc + n as i64 * 4;
    }

    let mut nCnt: Vec<Vec<i32>> = Vec::with_capacity(n_Tile.len());
    let mut tIdx: Vec<Vec<i64>> = Vec::with_capacity(n_Tile.len());

    for (i, k) in n_Tile.iter().enumerate() {
        let mut cnt = Vec::with_capacity(*k as usize);
        for _ in 0..*k {
            cnt.push(reader.read_i32::<LittleEndian>()?);
        }
        nCnt.push(cnt);

        let mut idx = Vec::with_capacity(*k as usize);
        idx.push(chr_loc); // Assuming chr_loc is calculated outside this function
        for j in 1..*k {
            idx.push(
                idx[j as usize - 1] + (nCnt[i as usize][j as usize - 1] as i64) * gdsize as i64,
            );
        }
        tIdx.push(idx);
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

    for name in c_name {
        println!("Retrieved chrom name (cName):  {}", name);
    }

    // Place values in hash map
    for (i, name) in igd.cName.iter().enumerate() {
        hash_table.insert(name.to_string(), i as i32);
    }

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

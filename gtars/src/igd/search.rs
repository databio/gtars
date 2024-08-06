use clap::ArgMatches;
use std::path::Path;
use crate::common::consts::{BED_FILE_EXTENSION, IGD_FILE_EXTENSION};
use crate::igd::create::igd_t;
use std::fs::{create_dir_all, DirEntry, File, OpenOptions};
use std::io::{BufRead, BufReader, Error, Read, Write};

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
    pub nbp: i32, //data type: 0, 1, 2 etc; size differs
    pub gType: i32, //data type: 0, 1, 2 etc; size differs
    pub nCtg: i32, //data type: 0, 1, 2 etc; size differs
    // Original code uses pointer to pointers
    pub cName: String,
    pub nTile: i32,
    pub nCnt: i32,
    pub tIdx: i32,

}

impl igd_t_from_disk {
    /// Constructs new instance of IGD
    pub fn new() -> Self {
        Self::default()
    }
}

#[derive(Default)]
pub struct info_t {
    pub fileName: String, //dataset file
    pub nr: i32, // number of regions in dataset
    pub md: f64, // average width of the regions

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
        Err(e) => {;
            return Err(e);
        
        }
        ,
    }

    println!("\n {} \n {}", database_path,query_file_path);


    //Get file info from the associated TSV



    // Create IGD Struct from database
    let IGD: igd_t_from_disk = get_igd_info(database_path).expect("Could not open IGD");



    // If query "-q" used set to mode 1

    match mode {

        1 => {



        },
        _ => {
            println!("Invalid mode selected, exiting");
            return Ok(());
        },


    }


    println!("FINISHED");
    
    Ok(())

}

pub fn get_igd_info(database_path: &String) -> Result<igd_t_from_disk, Error>{

    println!("hello from get_igd_info");

    let igd = igd_t_from_disk::new();

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

    let mut buffer = [0u8; std::mem::size_of::<i32>()];

    reader.read_exact(&mut buffer)?;
    let nbp = i32::from_le_bytes(buffer);
    reader.read_exact(&mut buffer)?;
    let gType = i32::from_le_bytes(buffer);

    reader.read_exact(&mut buffer)?;
    let nctg = i32::from_le_bytes(buffer);

    println!("Found:\n nbp:{} gtype: {} nctg: {}", nbp,gType,nctg);



    return Ok(igd)


}

fn check_file_extension(path: &str, expected_extension: &str) -> Result<(), String> {
    let path = Path::new(path);
    let actual_extension = path
        .extension()
        .and_then(|ext| ext.to_str())
        .ok_or_else(|| format!("Invalid file path: {}", path.display()))?;

    if actual_extension != expected_extension {
        return Err(format!("Incorrect file extension. Expected: {}, got: {}", expected_extension, actual_extension));
    }

    Ok(())
}

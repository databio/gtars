use std::collections::HashMap;
use clap::ArgMatches;
use std::fs;
use std::fs::{DirEntry, File};
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use crate::common::consts::BED_FILE_EXTENSION;
//use clap::error::ContextValue::String;
//use polars::export::arrow::buffer::Buffer;
//use crate::vocab::consts;


pub const maxCount: i64 = 268435456;		//16* = 4GB memory  // original code had this as i32



#[derive(Default)]
pub struct gdata_t {
    pub idx: i32,
    pub start: i32,
    pub end: i32,
    pub value: i32,
}

#[derive(Default)]
pub struct tile_t {
    pub ncnts: i32, // batch counts
    pub nCnts: i32, // total (batch) counts
    pub mcnts: i32, //  max counts
    pub gList: gdata_t,
}
#[derive(Default)]
pub struct ctg_t {
    pub name: String,
    pub mTiles: i32,
    pub gTile: Vec<tile_t>,
}
impl ctg_t{

    /// Constructs new instance of IGD
    pub fn new() -> Self {Self::default()}

}

#[derive(Default)]
pub struct igd_t {
    // TODO create attributes for the IGD
    pub nbp: i32,
    pub gType: i32,
    pub nctg: i32,
    pub mctg: i32,
    pub total: i64,
    pub ctg: Vec<ctg_t>, // this might need to be a reference
}


// impl Default for igd_t{
//     pub fn default() -> Self {
//         todo!()
//     }
// }

impl igd_t{

    /// Constructs new instance of IGD
    pub fn new() -> Self {Self::default()}

}

impl tile_t{

    /// Constructs new instance of tile
    pub fn new() -> Self {Self::default()}

}

pub fn create_igd_f(matches: &ArgMatches){

    println!("HELLO FROM IGD SUBMODULE!");

    let output_path = matches
        .get_one::<String>("output")
        .expect("Output path is required");

    let filelist = matches
        .get_one::<String>("filelist")
        .expect("File list path is required");

    //Initialize IGD into Memory
    let mut igd = igd_t::new();

    igd.gType = 1;
    igd.nbp = 16384; // from og code tile_size = 16384;  -> this is the bin size (2^14) from the original paper
    igd.nctg = 0;
    igd.mctg = 32;
    igd.total=0;

    //Check that file path exists and get number of files
    let mut all_bed_files: Vec<PathBuf>  = Vec::new();
    //let mut all_bed_buffers = Vec::new();

    let mut ix = 0;
    let (mut start, mut end) = (0,0);

    ///--------------------
    /// Check each file and only keep the validated BED files
    ///
    /// -------------------

    for entry in fs::read_dir(filelist).unwrap() {

        // For now only take .bed files
        if let Some(extension) = entry.as_ref().unwrap().path().extension() {

            if extension != BED_FILE_EXTENSION.trim_start_matches('.') {
                continue;
            }
        } else {continue} // This will skip files that do not have an extension

        let entry = entry.unwrap();
        let file_type = entry.file_type().unwrap();

        if file_type.is_file() {

            // open bed file
            // TODO original code uses gzopen (I assume for .gz files?)
            let file = File::open(entry.path()).unwrap();

            let mut reader = BufReader::new(file);

            /// Read the very first line and see if it meets our criteria
            /// MUST USE by_ref() otherwise borrow checker won't let code compile
            /// ALSO bec careful to call by_ref() BEFORE .lines()
            ///
            let first_line = reader.by_ref().lines().next().unwrap().expect("expect");

            //TODO Need to do error handling to ensure we gracefully continue if there is no data in the file.
            let mut lines = reader.lines();

            // TODO Better name for og function?
            // TODO parse_bed -> parse_bed_file_line
            let ctg = parse_bed(&first_line, &mut start, &mut end);
            // if it parses, add it to collected lines, increment ix
            match ctg{

                Some(ctg) =>{
                    //all_bed_files.push(entry.path());
                    //all_bed_files.push(line);
                    //all_bed_buffers.push(lines);
                    all_bed_files.push(entry.path());
                    ix +=1;
                } ,
                None => continue,
            }

        }
    }

    //println!("ALL PARSED Lines from BED FILES:\n{:?}", all_bed_files);

    let n_files = ix;//all_bed_files.len();

    println!("Number of Bed Files found:\n{}", n_files);

    //Check that there is more than 0 files?

    //Prep memory allocation in a Rust-like manner
    // TODO original code checks that the bed file can be parsed BEFORE memory allocation
    // TODO but then re-parses the bed file again later.
    // TODO use something like avg.shrink_to_fit(); after we've collected all the files?
    // og C code:
    //    int32_t *nr = calloc(n_files, sizeof(int32_t));
    //     double *avg = calloc(n_files, sizeof(double));
    let mut avg: Vec<i32> = Vec::with_capacity(n_files); //Can we use arrays? Is this an array? no, can we put an array on files.
    avg.resize(n_files, 0);

    let mut nr: Vec<i32> = Vec::with_capacity(n_files);
    nr.resize(n_files, 0);

    ///--------------------
    /// READ VALIDATED FILES
    /// Note: this seems wasteful to load the file *again* using BufReader
    /// Is there a better way than below?
    /// -------------------
    // Initialize required variables
    let (mut i0, mut i1, mut L0, mut L1) = (0, 0, 0, 1);
    let (mut  va, mut i, mut j, mut k,
        mut ig, mut m, mut nL, mut nf10) =
        (0,0,0,0,0,0,0,n_files/10);


    while i0 < n_files {
        //from og code: 2.1 Start from (i0, L0): read till (i1, L1)
        ig = i0;
        m = 0;
        //from og code: 2.2 Read ~4GB data from files
        // og code skips first line (since its already in the vec but we need to reread the file.
        while m==0 && ig<n_files{ //og comment: m>0 defines breaks when reading maxCount

            // Have to take ref and then clone the PathBuf
            // TODO Is this the proper way to do it??
            let file_path_buf = &all_bed_files[ig]; // could not move all_bed_files, so using reference to the PathBuf
            let fp = file_path_buf.clone();

            let file = File::open(fp).unwrap();
            let mut reader = BufReader::new(file);

            nL=0;

            let mut buffer = String::new();

            while m==0 && reader.read_line(&mut buffer).unwrap() != 0{

                let ctg = parse_bed(&buffer, &mut start, &mut end);

                match ctg{

                    Some(ctg) =>{
                        // check that st>=0 and end <321000000   NOTE: these values taken from og code.
                        if start>=0 && end<321000000{
                            /// igd_add not yet implemented
                            igd_add(&mut igd, ctg, start, end, va, ig);
                            nr[ig] +=1;
                            avg[ig]+=end-start;
                            println!("DEBUG: after igd add");

                        }
                    } ,
                    None => continue,
                }

                nL+=1;

                if igd.total > maxCount{

                    m=1;
                    i1 =ig;
                    L1= nL;

                }

            }

            if m==0 {
                ig+=1;
            }

        }

        ///og: 2.3 save/append tiles to disc, add cnts to Cnts
        ///

        igd_saveT(&igd, output_path);

        i0 = ig;
        L0 = L1;
        L1 = 0;

    }

//TODO CODE TO save _index.tsv (part 3)

//TODO COde to sort tile data and save into single files per ctg (part 4)

}

pub fn igd_saveT(p0: &igd_t, p1: &String) {
    println!("HELLO from igd_saveT");
    //todo!()
}

pub fn igd_add(igd: &mut igd_t, chrm: String, start: i32, end: i32, v: i32, idx: usize) {
    ///Add an interval
    /// og code: layers: igd->ctg->gTile->gdata(list)

    if start>= end {

        println!("Start: {0} greater than End: {1}, returning from igd_add", start, end);
        return
    }
    let absent: i32;
    let i: i32;

    // Cloning chrm String because the hash table will own the key after insertion
    let mut key= chrm.clone();

    let n1 = start/igd.nbp;
    let n2 = (end-1)/igd.nbp;

    // create hash table
    let mut hash_table:HashMap<String, i32> = HashMap::new();

    let key_check = hash_table.contains_key(&key);


    if key_check == false{

        // Insert key and value (igd.nctg)
        hash_table.insert(key, igd.nctg);
        igd.nctg+=1;
        // initialize ctg
        let mut p = ctg_t::new();
        p.name = chrm;
        p.mTiles = 1 + n2;
        //p.gTile original code mallocs mTiles*sizeof title_t
        //p.gTile = Vec::with_capacity()

        for i in 0..p.mTiles{
            let mut new_tile: tile_t = tile_t::new();
            new_tile.ncnts = 0; //each batch
            new_tile.nCnts = 0; //total
            new_tile.mcnts =2 ;
            //new_tile.gList //tile->gList = malloc(tile->mcnts*sizeof(gdata_t));
            p.gTile.push(new_tile);

        }

        igd.ctg.push(p);

        // set key to name kh_key(h, k) = p->name;

    }

    println!("Here is hash map{:?}", hash_table);
    //let k = hash_table.insert()

    println!("HELLO from igd_add");


}

#[derive(PartialEq)] // So that we can do comparisons with equality operator
pub enum ParseBedResult {
    Str(String),
    Int(i32),
}

pub fn parse_bed(line: &String, start: &mut i32, end: &mut i32) -> Option<String> {

    println!("HERE IS THE LINE TO PARSE: {}", line);
    let mut fields = line.split('\t');
    // Get the first field which should be chromosome.
    let ctg = fields.next()?; // Why is ctg used as variable name in og code?
    println!("GOT CHR: {}", ctg);
    // Parse 2nd and 3rd string as integers or return -1 if failure
    let st = fields.next().and_then(|s| s.parse::<i32>().ok()).unwrap_or(-1);
    println!("GOT st: {}", st);
    let en = fields.next().and_then(|s| s.parse::<i32>().ok()).unwrap_or(-1);
    println!("GOT en: {}", en);

    // if fields.next().is_some() || !ctg.starts_with("chr") || ctg.len() >= 40 || en <= 0 {
    //     return None;
    // }
    if !ctg.starts_with("chr") || ctg.len() >= 40 || en <= 0 {
        println!("RETURNING NONE");
        return None;
    }


    *start = st;
    *end = en;

    println!("SUCCESSFULLY FINISHING PARSE");
    Some(ctg.parse().unwrap())

}

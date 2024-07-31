use std::collections::HashMap;
use clap::ArgMatches;
use std::{fs, io};
use std::fs::{create_dir_all, DirEntry, File, OpenOptions};
use std::io::{BufRead, BufReader, Read, Write, Error};
use std::path::{Path, PathBuf};
use std::mem;
use std::mem::size_of;
use crate::common::consts::BED_FILE_EXTENSION;
//use clap::error::ContextValue::String;
//use polars::export::arrow::buffer::Buffer;
//use crate::vocab::consts;
use anyhow::{Context, Result};

pub const maxCount: i64 = 268435456;		//16* = 4GB memory  // original code had this as i32



#[derive(Default)]
pub struct gdata_t {
    pub idx: usize, //genomic object--data set index
    pub start: i32, //region start
    pub end: i32, //region end
    pub value: i32,
}
impl gdata_t {

    /// Constructs new instance of a gdata_t
    pub fn new() -> Self {Self::default()}

}

#[derive(Default)]
pub struct tile_t {
    pub ncnts: i32, // batch counts
    pub nCnts: i32, // total (batch) counts
    pub mcnts: i32, //  max counts
    pub gList: Vec<gdata_t>, //genomic data
}
#[derive(Default)]
pub struct ctg_t {
    pub name: String, //name of the contig
    pub mTiles: i32,  //determined by the interval start and end
    pub gTile: Vec<tile_t>, //tile data
}
impl ctg_t{

    /// Constructs new instance of a ctg
    pub fn new() -> Self {Self::default()}

}

#[derive(Default)]
pub struct igd_t {
    // TODO create attributes for the IGD
    pub nbp: i32, //data type: 0, 1, 2 etc; size differs
    pub gType: i32, //data type: 0, 1, 2 etc; size differs
    pub nctg: i32, //data type: 0, 1, 2 etc; size differs
    pub mctg: i32, //data type: 0, 1, 2 etc; size differs
    pub total: i64, // total region in each ctg
    pub ctg: Vec<ctg_t>, // this is the list of contigs (of size n-ctg)  // this might need to be a reference
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

    let db_output_name = matches
        .get_one::<String>("dbname")
        .expect("File list path is required");

    //println!("{}",db_output_name);
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
    let nf10 = n_files/10;

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

            if nf10>1 {
                if ig % nf10 == 0 {
                    println!(".") // SHow progress for every 10 files
                }
            }

        }

        ///og: 2.3 save/append temp tiles to disc, add cnts to Cnts
        ///

        igd_saveT(&igd, output_path);

        i0 = ig;
        L0 = L1;
        L1 = 0;

    }

//TODO CODE TO save _index.tsv (part 3)

    //sprintf(idFile, "%s%s%s", oPath, igdName, "_index.tsv");
    let tsv_save_path = format!("{}{}{}",output_path,db_output_name,"_index.tsv");
    let tsv_parent_path = tsv_save_path.clone();
    let path = std::path::Path::new(&tsv_parent_path).parent().unwrap();
    let result = create_file_with_parents(path);

    match result {
        Ok(file) => println!("TSV File created or opened successfully!"),
        Err(err) => println!("Error creating file: {}", err),
    }
    let mut file = OpenOptions::new()
        .create(true)  // Create the file if it doesn't exist
        .append(true)  // Append data to the existing file if it does exist
        .open(tsv_save_path).unwrap();

    //fprintf(fpi, "Index\tFile\tNumber of regions\tAvg size\n");

    let initial_line = format!("Index\tFile\tNumber of Regions\t Avg size\n");
    let mut buffer = Vec::new();
    buffer.write_all((&initial_line).as_ref()).unwrap();

    let mut total_regions = 0;
    let mut total_avg_size = 0.0;

    for i in 0..n_files {

        let file_path = &all_bed_files[i].to_str().unwrap();

        // TODO this line isn't not grabbing the end name as desired
        let filename = file_path.rsplitn(1, '/',).next().unwrap_or(file_path);

        total_regions += nr[i];
        total_avg_size += avg[i] as f32;

        // Write file summary
        //writeln!(fpi, "{} \t {} \t {} \t {}", i, filename, nr[i], avg[i] / nr[i]).expect("Couldn't write to file");
        let current_line = format!("{} \t {} \t {} \t {}", i, filename, nr[i], avg[i] / nr[i]);
        buffer.write_all((&current_line).as_ref()).unwrap();
    }

    file.write_all(&buffer).unwrap();


//TODO Code to sort tile data and save into single files per ctg (part 4)

    // Sort tile data and save into single files per ctg
    igd_save_db(igd, output_path, db_output_name)

}

pub fn igd_save_db(igd: igd_t, output_path: &String, db_output_name: &String) {
    println!("HELLO from igd_save_db");
    // this is the igd_save func from the original c code

    // sprintf(idFile, "%s%s%s_%i", oPath, "data0/", ctg->name, j);
    let save_path = format!("{}{}{}",output_path,db_output_name,".igd");
    let parent_path = save_path.clone();

    let path = std::path::Path::new(&parent_path).parent().unwrap();
    let result = create_file_with_parents(path);

    match result {
        Ok(file) => println!("File created or opened successfully!"),
        Err(err) => println!("Error creating file: {}", err),
    }

    let mut file = OpenOptions::new()
        .create(true)  // Create the file if it doesn't exist
        .append(true)  // Append data to the existing file if it does exist
        .open(save_path).unwrap();

    let mut buffer = Vec::new();

    // for data in &current_tile.gList[..current_tile.ncnts as usize] {
    //     buffer.write_all(&data.idx.to_le_bytes()).unwrap();
    //     buffer.write_all(&data.start.to_le_bytes()).unwrap();
    //     buffer.write_all(&data.end.to_le_bytes()).unwrap();
    //     buffer.write_all(&data.value.to_le_bytes()).unwrap();
    // }
    //
    buffer.write_all(&igd.nbp.to_le_bytes()).unwrap();
    buffer.write_all(&igd.gType.to_le_bytes()).unwrap();
    buffer.write_all(&igd.nctg.to_le_bytes()).unwrap();


    for i in 0..igd.nctg{

        let idx = i.clone() as usize;
        let current_ctg = &igd.ctg[idx];


        buffer.write_all(&current_ctg.mTiles.to_le_bytes()).unwrap();

    }

    for i in 0..igd.nctg{
        let idx = i.clone() as usize;
        let current_ctg = &igd.ctg[idx];

        //let j = igd.nctg;

        let n = current_ctg.mTiles;

        for j in 0..n{
            let jdx = j.clone() as usize;

            buffer.write_all(&current_ctg.gTile[jdx].nCnts.to_le_bytes()).unwrap();
        }

    }

    for i in 0..igd.nctg{

        let idx = i.clone() as usize;
        let current_ctg = &igd.ctg[idx];

        buffer.write_all((&current_ctg.name).as_ref()).unwrap();

    }

    file.write_all(&buffer).unwrap();
    //2. SOrt and save tiles data

    let k: i32;

    for i in 0..igd.nctg{
        let idx = i.clone() as usize;

        let current_ctg = &igd.ctg[idx];
        let n = current_ctg.mTiles;

        for j in 0..n{
            let jdx = j.clone() as usize;

            let mut q = &current_ctg.gTile[jdx];
            let nrec = q.nCnts;

            if nrec>0{
                println!("nrec greater than 0");
                let save_path = format!("{}{}{}_{}{}",output_path,"data0/",current_ctg.name, j,".igd");
                let parent_path = save_path.clone();
                let path = std::path::Path::new(&parent_path).parent().unwrap();

                let mut file = OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(path);

                match file {
                    Ok(file) => {
                        println!("File created or opened successfully!");
                    }
                    Err(_) => {println!("Cannot open path!!!");
                    return;
                    }
                }

                // Read from Temp File
                //the next 4 lines are pulled from googling and are not quite right
                //let gdsize = nrec * std::mem::size_of::<gdata_t>() as i32;

                //let mut gdata = vec![gdata_t::default(); gdsize as usize];

                //let ni = file.read_exact(gdata.as_mut_slice().to_le_bytes());

                // Sort Data
                //gdata.sort_by_key(|d| d.start); // Sort by start value

                // Write to database after sorting
                //let _ = file.write_all(&gdata);

                // og code!!!!!!!!!!!!
                // gdsize = nrec*sizeof(gdata_t);
                // gdata_t *gdata = malloc(gdsize);
                // if(gdata==NULL){
                //     printf("Can't alloc mem %lld\n", (long long)gdsize);
                //     return;
                // }
                // ni = fread(gdata, gdsize, 1, fp0);
                // fclose(fp0);
                // //qsort(gdata, nrec, sizeof(gdata_t), compare_rstart);
                // radix_sort_intv(gdata, gdata+nrec);
                // fwrite(gdata, gdsize, 1, fp);
                // free(gdata);
                // remove(iname);


            }

            // todo set to zero but it claims that this is immutable
            //q.nCnts = 0;


        }

    }


    //file.write_all(&buffer).unwrap();


}

pub fn igd_saveT(igd: &igd_t, output_file_path: &String) {
    println!("HELLO from igd_saveT");

    // From OG COde:
    // TEMPORARILY save/append tiles to disc, add cnts to Cnts; reset tile.gList

    let mut nt =0;

    for i in 0..igd.nctg{

        let idx = i.clone() as usize;
        let idx_2 = idx;
        let current_ctg = &igd.ctg[idx_2];
        nt = nt + current_ctg.mTiles;

        for j in 0..current_ctg.mTiles{

            let jdx = j.clone() as usize;
            let jdx_2 = jdx;

            let current_tile = &current_ctg.gTile[jdx_2];

            if current_tile.ncnts>0{

                // Construct specific temp file on disk using this information

                // OG code
                // sprintf(idFile, "%s%s%s_%i", oPath, "data0/", ctg->name, j);
                let save_path = format!("{}{}{}_{}{}",output_file_path,"data0/",current_ctg.name, j,".igd");
                let parent_path = save_path.clone();

                println!("{}",save_path);

                //todo this needs to create the path if it does not already exist!!!

                let path = std::path::Path::new(&parent_path).parent().unwrap();
                let result = create_file_with_parents(path);

                match result {
                    Ok(file) => println!("File created or opened successfully!"),
                    Err(err) => println!("Error creating file: {}", err),
                }

                //let _ = create_dir_all(save_path.clone());
                //if let Ok(ret) = create_dir_all(save_path.clone());
                //
                // match result {
                //     Ok(_) => println!("Directory created successfully!"), // Optional: Print a success message
                //     Err(ref error) if error.kind() == fs:: => {
                //         println!("Directory already exists. Ignoring error.");
                //     },
                //     Err(error) => println!("Error creating directory: {}", error), // Handle other errors
                // }
                // let path = std::path::Path::new(&save_path);
                //
                // if let Some(parent) = path.parent() {
                //     std::fs::create_dir_all(parent).unwrap();
                // } else {
                //     anyhow::Error("Failed to create parent directories for gtok file!")
                // }

                let mut file = OpenOptions::new()
                    .create(true)  // Create the file if it doesn't exist
                    .append(true)  // Append data to the existing file if it does exist
                    .open(save_path).unwrap();

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


            }



        }

    }




}

fn create_file_with_parents(path: &Path) -> Result<File, Error> {
    // Create all parent directories if they don't exist (ignore errors)
    let _ = create_dir_all(path);  // Discard the result (success or error)

    // Open the file for creation or append, ignoring errors if it exists
    let file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(path);

    match file {
        Ok(file) => {
            println!("File created or opened successfully!");
            Ok(file)
        }
        Err(_) => Ok(File::open(path).unwrap_or_else(|_| File::create(path).unwrap())) // Handle existing file or create new one
    }

}

// fn create_file_with_parents(path: &Path) -> Result<File, Error> {
//     // Create all parent directories if they don't exist
//     let result = create_dir_all(path).unwrap();
//
//     match result {
//         Ok(file) => println!("File created or opened successfully!"),
//         Err(err) => println!("Error creating file: {}", err),
//     }
//
//
//     // Open the file for creation or append, ignoring errors if it exists
//     Ok(OpenOptions::new()
//         .create(true)
//         .append(true) // Optional: Append to existing file
//         .open(path)?)
// }


pub fn igd_add(igd: &mut igd_t, chrm: String, start: i32, end: i32, v: i32, idx: usize) {
    ///Add an interval
    /// og code: layers: igd->ctg->gTile->gdata(list)

    println!("HELLO from igd_add");

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

        println!("Key does not exist in hash map, creating for {}", key.clone());

        // Insert key and value (igd.nctg)
        hash_table.insert(key.clone(), igd.nctg);
        igd.nctg+=1;

        // initialize ctg
        let mut p = ctg_t::new();
        p.name = chrm;
        p.mTiles = 1 + n2;
        //p.gTile original code mallocs mTiles*sizeof title_t
        // however in Rust, structs have 0 size: https://doc.rust-lang.org/nomicon/exotic-sizes.html#zero-sized-types-zsts
        //p.gTile = Vec::with_capacity((p.mTiles as usize)*size_of(tile_t()));
        p.gTile = Vec::with_capacity((p.mTiles as usize));

        for i in 0..p.mTiles{

            let mut new_tile: tile_t = tile_t::new();

            new_tile.ncnts = 0; //each batch
            new_tile.nCnts = 0; //total
            new_tile.mcnts =2 ;
            //new_tile.gList //tile->gList = malloc(tile->mcnts*sizeof(gdata_t));
            //new_tile.gList = Vec::with_capacity((new_tile.mcnts as usize));

            for j in 0..new_tile.mcnts{
                new_tile.gList.push(gdata_t::new());
            }
            // for element in new_tile.gList.iter_mut() {
            //     //*element = gdata_t::new();  // Add new_value to each element
            //     //element.push(gdata_t::new());
            //     let element = &mut gdata_t::new();
            // }

            p.gTile.push(new_tile);

        }

        igd.ctg.push(p);

        // set key to name kh_key(h, k) = p->name;

    }

    // Retrieve values from Hash Map
    // println!("Here is hash map{:?}", hash_table);
    //let k = hash_table.insert()

    let keycloned = key.clone();

    let index = hash_table.get(&keycloned).unwrap();
    let cloned_index = index.clone();


    let p = &mut igd.ctg[cloned_index as usize];

    if (n2+1>=p.mTiles){

        println!("TRUE:{} vs {}", (n2+1), p.mTiles.clone());
        let tt = p.mTiles;

        p.mTiles = n2+1;
        // original code: p->gTile = realloc(p->gTile, p->mTiles*sizeof(tile_t));
        // Supposedly we may not need to do this ...  p.gTile = Vec::resize()   ???

        for i in tt..p.mTiles{

            let idx = i.clone() as usize;
            let idx_2 = idx as usize;

            let existing_tile: &mut tile_t = &mut p.gTile[idx_2];

            existing_tile.ncnts = 0;
            existing_tile.nCnts = 0;
            existing_tile.mcnts = 2;
            // og: tile->gList = malloc(tile->mcnts*sizeof(gdata_t));
            //existing_tile.gList = gdata_t::new(); // TODO Double check this, do we actually want to create a new struct?
            //existing_tile.gList = Vec::with_capacity((existing_tile.mcnts as usize));
            // for element in existing_tile.gList.iter_mut() {
            //     //*element = gdata_t::new();  // Add new_value to each element
            //     //element.push(gdata_t::new());
            //     let element = gdata_t::new();
            // }
            // existing_tile.gList = Vec::with_capacity(existing_tile.mcnts as usize)
            //     .iter_mut()  // Iterate over mutable references (not needed here)
            //     .map(|gdata_t: &mut gdata_t| gdata_t::new())  // Create new gdata_t for each element
            //     .collect();
            for j in 0..existing_tile.mcnts{
                existing_tile.gList.push(gdata_t::new());
            }

        }

    }

    for i in n1..=n2{ //this is inclusive of n1 and n2
        // Get index as usize
        let idx_1 = i.clone() as usize;
        let idx_2 = idx_1 as usize;
        // get the tile for the contig
        let existing_tile: &mut tile_t = &mut p.gTile[idx_2];
        // og code, not necessary in Rust?		if(tile->ncnts == tile->mcnts)
        // 			EXPAND(tile->gList, tile->mcnts);

        let tile_idx = existing_tile.ncnts.clone() as usize;
        let gdata = &mut existing_tile.gList[tile_idx];
        existing_tile.ncnts = existing_tile.ncnts+ 1;

        gdata.start = start;
        gdata.end = end;
        gdata.value = v;
        gdata.idx = idx;

    }

    println!("Finished from igd_add");
    return

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

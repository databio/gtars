use clap::ArgMatches;
use std::fs;
use std::fs::{DirEntry, File};
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
//use clap::error::ContextValue::String;
use polars::export::arrow::buffer::Buffer;
use crate::vocab::consts;


pub const maxCount: i32 = 268435456;		//16* = 4GB memory

#[derive(Default)]
pub struct IGD {
    // TODO create attributes for the IGD
    pub placeholder: String,
    pub total: i32,
}

impl IGD{

    /// Constructs new instance of IGD
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
    let mut igd = IGD::new();

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

            if extension != consts::FILE_EXTENSION.trim_start_matches('.') {
                continue;
            }
        }

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
            let mut lines = reader.lines();

            // TODO Better name for og function?
            // TODO parse_bed -> parse_bed_file_line
            let ctg = parse_bed(&first_line, start, end);
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
    let mut avg: Vec<i32> = Vec::with_capacity(n_files);
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

                let ctg = parse_bed(&buffer, start, end);

                match ctg{

                    Some(ctg) =>{
                        // check that st>=0 and end <321000000   NOTE: these values taken from og code.
                        if start>=0 && end<321000000{
                            /// igd_add not yet implemented
                            igd_add(&igd, ctg, start, end, va, ig);
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
            // if ig%nf10 == 0{
            //     println!(".") // og code: appears to be a debug line
            // }


            //
            // let first_line = reader.by_ref().lines().next().unwrap().expect("expect");
            // println!("Confirm reading first line: {}",first_line);
            // Get file from vec via index
            // read file
            ig +=1

        }

        ///og: 2.3 save/append tiles to disc, add cnts to cnts
        ///

        igd_saveT(&igd, output_path);

        i0 = ig;
        L0 = L1;
        L1 = 0;

    }

    // for path in all_bed_files{
    //
    //     // let file_path = path.unwrap()?;
    //
    //     println!("FIle path: {:?}", path);
    //
    // }
    // /// Debug check if first line is consumed...
    // for mut buf in all_bed_buffers{
    //     // CHECK IF first line consumed...
    //     for line in buf{
    //         println!("{:?}", line);
    //     }
    //
    // }
    // while i0 < n_files{
    //     //from og code: 2.1 Start from (i0, L0): read till (i1, L1)
    //     ig = i0;
    //     m = 0;
    //     //from og code: 2.2 Read ~4GB data from files
    //
    //
    //
    //
    // }



}

fn igd_saveT(p0: &IGD, p1: &String) {
    println!("HELLO from igd_saveT");
    //todo!()
}

fn igd_add(p0: &IGD, p1: String, p2: i32, p3: i32, p4: i32, p5: usize) {
    println!("HELLO from igd_add");
    //todo!()

}

#[derive(PartialEq)] // So that we can do comparisons with equality operator
pub enum ParseBedResult {
    Str(String),
    Int(i32),
}

pub fn parse_bed(line: &String, mut start: i32, mut end: i32) -> Option<String> {

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

    //*start = st; //Compiler said no.
    start = st;
    end = en;

    println!("SUCCESSFULLY FINISHING PARSE");
    Some(ctg.parse().unwrap())

}

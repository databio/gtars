use clap::ArgMatches;
use std::fs;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
//use clap::error::ContextValue::String;
use polars::export::arrow::buffer::Buffer;
use crate::vocab::consts;

#[derive(Default)]
pub struct IGD {
    // TODO create attributes for the IGD
    pub placeholder: String,
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


    // println!("Collected the following:");
    // println!("{0} \n {1} ",output_path, filelist)

    //Initialize IGD into Memory
    let mut igd = IGD::new();

    //Check that file path exists and get number of files
    let mut all_bed_files  = Vec::new();

    let mut ix = 0;
    let (mut start, mut end) = (0,0);

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

            let reader = BufReader::new(file);

            let mut buf = String::new();
            reader.buffer().read_to_string(&mut buf).expect("Cannot read buf string");
            // attempt to parse
            let ctg = parse_bed(buf, start, end);
            // if it parses, add it, increment ix


            match Some(ctg){

                Some(ctg) =>{
                    all_bed_files.push(entry.path());
                    ix +=1;
                } ,
                None => continue,
            }

        }
    }
    println!("ALL BED FILES:\n{:?}", all_bed_files);

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
    let mut avg: Vec<f64> = Vec::with_capacity(n_files);
    avg.resize(n_files, 0.0);

    let mut nr: Vec<i32> = Vec::with_capacity(n_files);
    nr.resize(n_files, 0);

    // READ FILES

    // Initialize required variables
    let (mut i0, mut i1, mut L0, mut L1) = (0, 0, 0, 1);
    let (mut  va, mut i, mut j, mut k, mut ig, mut m, mut nL, mut nf10) =
        (0,0,0,0,0,0,0,n_files/10);

    while i0 < n_files{

        println!("{}", i0);
        i0+=1;


    }

    for path in all_bed_files{

        println!("PATH: {:?}",path);



    }
    // Get file ids

    //Open files
    //Parse bed files
    //Close files

    // set number_of_files to the number of successfully opened and parsed files.




}

#[derive(PartialEq)] // So that we can do comparisons with equality operator
pub enum ParseBedResult {
    Str(String),
    Int(i32),
}

pub fn parse_bed<R: BufRead>(buf: String, start: i32, end: i32) -> Option<&str> {

    let str = String::from("Hello");

    let mut fields = buf.split('\t');

    let ctg = fields.next()?;

    let st = fields.next().and_then(|s| s.parse().ok())?;
    let en = fields.next().and_then(|s| s.parse().ok())?;

    if fields.next().is_some() || !ctg.starts_with("chr") || ctg.len() >= 40 || en <= 0 {
        return None;
    }

    *start = st;
    *end = en;

    Some(ctg)

}
// pub fn parse_bed(content: &[u8], start: i32, end: i32) -> ParseBedResult {
//
//     let str = String::from("Hello");
//
//
//
//     if !str.is_empty() {
//         ParseBedResult::Str(str)
//     }else{
//         ParseBedResult::Int(0)
//     }
//
// }
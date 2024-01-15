use clap::ArgMatches;
use std::fs;
use std::path::Path;
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
            all_bed_files.push(entry.path());

        }
    }
    println!("ALL BED FILES:\n{:?}", all_bed_files);

    //Check that there is more than 0 files
    // Get file ids

    //Open files
    //Parse bed files
    //Close files

    // set number_of_files to the number of successfully opened and parsed files.




}
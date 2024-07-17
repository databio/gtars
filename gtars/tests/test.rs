use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::fs::{File};

use rstest::*;
use tempfile::tempdir;

use gtars::uniwig::{parse_bed_file};

#[fixture]
fn path_to_data() -> &'static str {
    "tests/data"
}

#[fixture]
fn path_to_bed_file() -> &'static str {
    "tests/data/peaks.bed"
}

#[fixture]
fn path_to_sorted_small_bed_file() -> &'static str {
    "tests/data/test_sorted_small.bed"
}

#[fixture]
fn path_to_bed_file_gzipped() -> &'static str {
    "tests/data/peaks.bed.gz"
}

mod tests {
    use std::env::temp_dir;
    use gtars::uniwig::{Chromosome, read_bed_vec, uniwig_main};
    use gtars::igd::create::{parse_bed,create_igd_f,igd_add,igd_saveT};

    use super::*;

    // IGD TESTS

    #[rstest]
    fn test_igd_parse_bed_file() {

        // Given some random line from a  bed file...
        let bed_file_string = String::from("chr1	32481	32787	SRX4150706.05_peak_1	92	.	7.69231	13.22648	9.25988	155");

        //Placeholder start and end values
        let mut start = 0;
        let mut end = 0;

        let result = parse_bed(&bed_file_string, &mut start, &mut end).unwrap(); // this will return

        let unwrapped_result = result.as_str();

        assert_eq!(unwrapped_result, "chr1");

        // Ensure start and end is modified via parse_bed
        assert_eq!(start, 32481);
        assert_eq!(end, 32787);

    }



    // UNIWIG TESTS
    #[rstest]
    fn test_uniwig_parsed_bed_file(path_to_bed_file: &str) {
        let path = Path::new(path_to_bed_file);
        let file = File::open(path).unwrap();

        let mut reader = BufReader::new(file);
        let first_line = reader.by_ref().lines().next().unwrap().expect("expect");
        println!("{:?}", first_line);

        let result = parse_bed_file(&first_line);

        if let Some((ctg, st, en)) = result {

            println!("ctg: {}", ctg);
            println!("st: {}", st);
            println!("en: {}", en);
            assert_eq!(st, 7915738);
        } else {
            println!("Failed to parse BED record");
        }

    }

    #[rstest]
    fn test_uniwig_read_bed_vec(path_to_bed_file: &str, path_to_bed_file_gzipped: &str) {

        read_bed_vec(path_to_bed_file);
        read_bed_vec(path_to_bed_file_gzipped);

    }

    #[rstest]
    fn test_uniwig_read_bed_vec_length(path_to_sorted_small_bed_file: &str) {

        let chromosomes: Vec<Chromosome>  = read_bed_vec(path_to_sorted_small_bed_file);
        let num_chromosomes = chromosomes.len();

        assert_eq!(num_chromosomes, 5);

    }
    #[rstest]
    fn test_run_uniwig_main_wig_type(path_to_bed_file: &str) {

        let path_to_crate= env!("CARGO_MANIFEST_DIR");

        let tempbedpath = format!("{} {}",path_to_crate, "/tests/data/test5.bed");
        let combinedbedpath = tempbedpath.as_str();

        let chromsizerefpath: String = format!("{} {}",path_to_crate, "/tests/hg38.chrom.sizes");

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let bwfileheader_path = path.into_os_string().into_string().unwrap();
        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 5;
        let output_type ="wig";

        uniwig_main(smoothsize, combinedbedpath, &chromsizerefpath, bwfileheader, output_type)

    }

    #[rstest]
    fn test_run_uniwig_main_npy_type(path_to_bed_file: &str) {

        let path_to_crate= env!("CARGO_MANIFEST_DIR");

        let tempbedpath = format!("{} {}",path_to_crate, "/tests/data/test5.bed");
        let combinedbedpath = tempbedpath.as_str();

        let chromsizerefpath: String = format!("{} {}",path_to_crate, "/tests/hg38.chrom.sizes");

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let bwfileheader_path = path.into_os_string().into_string().unwrap();
        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 5;
        let output_type ="npy";

        uniwig_main(smoothsize, combinedbedpath, &chromsizerefpath, bwfileheader, output_type)

    }
}

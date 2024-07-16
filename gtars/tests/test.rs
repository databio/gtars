use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::fs::{File};

use rstest::*;
use tempfile::NamedTempFile;

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
    use gtars::common::utils::extract_regions_from_bed_file;
    use gtars::uniwig::{Chromosome, read_bed_vec, run_uniwig, uniwig_main};

    use super::*;

    #[rstest]
    fn test_parsed_bed_file(path_to_bed_file: &str) {
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
    fn test_read_bed_vec(path_to_bed_file: &str, path_to_bed_file_gzipped: &str) {

        read_bed_vec(path_to_bed_file);
        read_bed_vec(path_to_bed_file_gzipped);

    }

    #[rstest]
    fn test_read_bed_vec_length(path_to_sorted_small_bed_file: &str) {

        let mut chromosomes: Vec<Chromosome>  = read_bed_vec(path_to_sorted_small_bed_file);
        let num_chromosomes = chromosomes.len();

        assert_eq!(num_chromosomes, 5);

    }
    #[rstest]
    fn test_run_uniwig_main_wig_type(path_to_bed_file: &str) {

        let smoothsize: i32 = 5;
        let combinedbedpath: &str = "/home/drc/GITHUB/genimtools/genimtools/tests/data/test5.bed";
        let chromsizerefpath: String = "/home/drc/GITHUB/genimtools/genimtools/tests/hg38.chrom.sizes".to_string();
        let bwfileheader: &str = "/home/drc/Downloads/test_rust_wig/";
        let output_type ="wig";

        uniwig_main(smoothsize, combinedbedpath, &chromsizerefpath, bwfileheader, output_type)

    }

    #[rstest]
    fn test_run_uniwig_main_npy_type(path_to_bed_file: &str) {

        let smoothsize: i32 = 5;
        let combinedbedpath: &str = "/home/drc/GITHUB/genimtools/genimtools/tests/data/test5.bed";
        let chromsizerefpath: String = "/home/drc/GITHUB/genimtools/genimtools/tests/hg38.chrom.sizes".to_string();
        let bwfileheader: &str = "/home/drc/Downloads/test_rust_wig/";
        let output_type ="npy";

        uniwig_main(smoothsize, combinedbedpath, &chromsizerefpath, bwfileheader, output_type)

    }
    // #[rstest]
    // fn test_count_coordinate_reads() {
    //     // example input, marking read alignment locations
    //     let query: Vec<i32> = vec![2,2,2,3,3,7,10,12,12,12,12,15];
    //     let res = count_coordinate_reads(&query);
    //     // example output, counting number of reads at each position
    //     let answer = vec![0,3,2,0,0,0,1,0,0,1,0,4,0,0,1];
    //     assert_eq!(res, answer);
    //
    // }

    // #[rstest]
    // fn test_count_coordinate_reads_start_end() {
    //     // example input, marking read alignment locations
    //     let starts: Vec<i32> = vec![1,4,4,7,9,9];
    //     let ends: Vec<i32> = vec![3,6,6,9,10,11];
    //     let res = count_coordinate_reads_start_end(&starts, &ends);
    //
    // }
}

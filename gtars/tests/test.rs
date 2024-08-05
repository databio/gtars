use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use rstest::*;
use tempfile::tempdir;

use gtars::uniwig::parse_bed_file;

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
    use super::*;
    use gtars::uniwig::{read_bed_vec, read_chromosome_sizes, uniwig_main, Chromosome};
    use std::env::temp_dir;
    use std::ptr::read;

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
            panic!("Failed to parse BED record");
        }
    }

    #[rstest]
    fn test_read_bed_vec(path_to_bed_file: &str, path_to_bed_file_gzipped: &str) {
        let result1 = read_bed_vec(path_to_bed_file);
        assert_eq!(result1.len(), 20);

        let result2 = read_bed_vec(path_to_bed_file_gzipped);
        assert_eq!(result2.len(), 20);
    }

    #[rstest]
    fn test_read_bed_vec_length(path_to_sorted_small_bed_file: &str) {
        let chromosomes: Vec<Chromosome> = read_bed_vec(path_to_sorted_small_bed_file);
        let num_chromosomes = chromosomes.len();

        assert_eq!(num_chromosomes, 5);
    }
    #[rstest]
    fn test_run_uniwig_main_wig_type(
        path_to_bed_file: &str,
    ) -> Result<(), Box<(dyn std::error::Error + 'static)>> {
        // This test uses the bed file to determine chromsizes for speed
        let path_to_crate = env!("CARGO_MANIFEST_DIR");

        let tempbedpath = format!("{}{}", path_to_crate, "/tests/data/test5.bed");
        let combinedbedpath = tempbedpath.as_str();

        let chromsizerefpath = combinedbedpath;

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let bwfileheader_path = path.into_os_string().into_string().unwrap();
        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 5;
        let output_type = "wig";

        uniwig_main(
            smoothsize,
            combinedbedpath,
            chromsizerefpath,
            bwfileheader,
            output_type,
        )
        .expect("Uniwig main failed!");

        Ok(())
    }

    #[rstest]
    fn test_run_uniwig_main_npy_type(
        path_to_bed_file: &str,
    ) -> Result<(), Box<(dyn std::error::Error + 'static)>> {
        // This test uses the bed file to determine chromsizes for speed
        let path_to_crate = env!("CARGO_MANIFEST_DIR");

        let tempbedpath = format!("{}{}", path_to_crate, "/tests/data/test5.bed");
        let combinedbedpath = tempbedpath.as_str();

        let chromsizerefpath = combinedbedpath;

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let bwfileheader_path = path.into_os_string().into_string().unwrap();
        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 5;
        let output_type = "npy";

        uniwig_main(
            smoothsize,
            combinedbedpath,
            chromsizerefpath,
            bwfileheader,
            output_type,
        )
        .expect("Uniwig main failed!");
        Ok(())
    }

    #[rstest]
    fn test_reading_chrom_sizes(path_to_bed_file: &str) {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");

        // Read from sizes file
        let chromsizerefpath: String = format!("{}{}", path_to_crate, "/tests/hg38.chrom.sizes");
        let chrom_sizes = read_chromosome_sizes(chromsizerefpath.as_str()).unwrap();
        let chrom_name = String::from("chr13");
        let current_chrom_size = chrom_sizes[&chrom_name.clone()] as i32;
        assert_eq!(current_chrom_size, 114364328);

        // Read from BED file
        let tempbedpath = format!("{}{}", path_to_crate, "/tests/data/test5.bed");
        let combinedbedpath = tempbedpath.as_str();
        let chrom_sizes = read_chromosome_sizes(combinedbedpath).unwrap();
        let chrom_name = String::from("chr1");
        let current_chrom_size = chrom_sizes[&chrom_name.clone()] as i32;
        assert_eq!(current_chrom_size, 32);
    }

    #[rstest]
    fn test_uniwig_mismatched_chrom_sizes(path_to_bed_file: &str) {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");

        // Read from sizes file
        let chromsizerefpath: String = format!("{}{}", path_to_crate, "/tests/hg38.chrom.sizes");

        // Read from BED file that contains chromosomes not in size file
        let tempbedpath = format!("{}{}", path_to_crate, "/tests/data/test_unknown_chrom.bed");
        let combinedbedpath = tempbedpath.as_str();
        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let bwfileheader_path = path.into_os_string().into_string().unwrap();
        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 5;
        let output_type = "npy";

        let result = uniwig_main(
            smoothsize,
            combinedbedpath,
            &chromsizerefpath,
            bwfileheader,
            output_type,
        );

        assert!(result.is_err());
    }
}

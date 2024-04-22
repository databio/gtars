use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::fs::{File};

use rstest::*;
use tempfile::NamedTempFile;

use genimtools::common::models::{Region, RegionSet};
use genimtools::tokenizers::{Tokenizer, TreeTokenizer};
use genimtools::uniwig::{parse_bed_file, count_coordinate_reads, count_coordinate_reads_start_end};

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

#[fixture]
fn path_to_tokenize_bed_file() -> &'static str {
    "tests/data/to_tokenize.bed"
}

mod tests {
    use genimtools::common::utils::extract_regions_from_bed_file;
    use genimtools::uniwig::{Chromosome, read_bed_vec, run_uniwig, uniwig_main};

    use super::*;

    #[rstest]
    fn test_region() {
        let region = Region {
            chr: "chr1".to_string(),
            start: 100,
            end: 200,
        };

        assert_eq!(region.chr, "chr1");
        assert_eq!(region.start, 100);
        assert_eq!(region.end, 200);
    }

    #[rstest]
    fn test_extract_regions_from_bed_file(path_to_bed_file: &str) {
        let path = Path::new(path_to_bed_file);
        let regions = extract_regions_from_bed_file(path);
        assert!(regions.is_ok(), "Failed to extract regions from BED file");
        let regions = regions.unwrap();
        assert!(regions.len() == 25);
    }

    #[rstest]
    fn test_extract_regions_from_bed_file_gzipped(path_to_bed_file_gzipped: &str) {
        let path = Path::new(path_to_bed_file_gzipped);
        let regions = extract_regions_from_bed_file(path);
        assert!(regions.is_ok(), "Failed to extract regions from BED file");
        let regions = regions.unwrap();
        assert!(regions.len() == 25);
    }

    #[rstest]
    fn test_region_set_from_bed(path_to_bed_file: &str) {
        let path = Path::new(path_to_bed_file);
        let rs = RegionSet::try_from(path).unwrap();

        assert!(rs.regions.height() == 25);
    }

    #[rstest]
    fn test_region_set_to_bed(path_to_bed_file: &str) {
        let path = Path::new(path_to_bed_file);
        let rs = RegionSet::try_from(path).unwrap();

        // create a temporary file
        let tmp_file = NamedTempFile::new().unwrap();
        let tmp_path = tmp_file.into_temp_path();
        let tmp_path = Path::new(tmp_path.to_str().unwrap());

        // write the region set to the temporary file
        rs.to_bed(tmp_path).unwrap();

        // read the temporary file back in as a region set
        let rs2 = RegionSet::try_from(tmp_path).unwrap();

        assert!(rs2.regions.height() == 25);
    }

    #[rstest]
    fn test_create_tokenizer(path_to_bed_file: &str) {
        let tokenizer = TreeTokenizer::from(Path::new(path_to_bed_file));
        println!("{}", tokenizer.universe.len());
        assert!(tokenizer.universe.len() == 27); // 25 regions + 2 special tokens
    }

    #[rstest]
    fn test_tokenize_bed_file(path_to_bed_file: &str, path_to_tokenize_bed_file: &str) {
        let tokenizer = TreeTokenizer::from(Path::new(path_to_bed_file));
        let rs = RegionSet::try_from(Path::new(path_to_tokenize_bed_file)).unwrap();
        let tokenized_regions = tokenizer.tokenize_region_set(&rs).unwrap();

        println!("{}", tokenized_regions.len());
        assert!(tokenized_regions.len() == 4);

        // last should be the unknown token
        let unknown_token = tokenized_regions.regions[3].clone();
        assert!(unknown_token.chr == "chrUNK");
    }

    #[rstest]
    fn test_pretokenization_folder(path_to_data: &str, path_to_bed_file: &str) {
        let tokenizer = TreeTokenizer::from(Path::new(path_to_bed_file));
        let path_to_data = Path::new(path_to_data);
        let outdir = "tests/data/out";

        let res = genimtools::tools::pre_tokenize_data(path_to_data, outdir, &tokenizer);
        assert!(res.is_ok());
    }

    #[rstest]
    fn test_pretokenization_file(path_to_tokenize_bed_file: &str, path_to_bed_file: &str) {
        let tokenizer = TreeTokenizer::from(Path::new(path_to_bed_file));
        let path_to_data = Path::new(path_to_tokenize_bed_file);
        let outdir = "tests/data/out";

        let res = genimtools::tools::pre_tokenize_data(path_to_data, outdir, &tokenizer);
        assert!(res.is_ok());
    }

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
    fn test_run_uniwig_main(path_to_bed_file: &str) {

        let sorted: bool = true;
        let smoothsize: i32 = 5;
        let writesize: i32 = 1;
        let combinedbedpath: &str = "/home/drc/GITHUB/genimtools/genimtools/tests/data/test5.bed";
        let chromsizerefpath: String = "/home/drc/GITHUB/genimtools/genimtools/tests/hg38.chrom.sizes".to_string();
        let bwfileheader: &str = "/home/drc/Downloads/test_rust_wig/";
        let output_type ="wig";

        uniwig_main(sorted, smoothsize, combinedbedpath, &chromsizerefpath, bwfileheader, output_type)

    }

    #[rstest]
    fn test_count_coordinate_reads() {
        // example input, marking read alignment locations
        let query: Vec<i32> = vec![2,2,2,3,3,7,10,12,12,12,12,15];
        let res = count_coordinate_reads(&query);
        // example output, counting number of reads at each position
        let answer = vec![0,3,2,0,0,0,1,0,0,1,0,4,0,0,1];
        assert_eq!(res, answer);

    }

    #[rstest]
    fn test_count_coordinate_reads_start_end() {
        // example input, marking read alignment locations
        let starts: Vec<i32> = vec![1,4,4,7,9,9];
        let ends: Vec<i32> = vec![3,6,6,9,10,11];
        let res = count_coordinate_reads_start_end(&starts, &ends);

        // example output, counting number of reads at each position
        // let answer = vec![0,3,2,0,0,0,1,0,0,1,0,4,0,0,1];
        // assert_eq!(res, answer);

    }
}

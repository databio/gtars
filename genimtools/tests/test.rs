use std::path::Path;

use rstest::*;
use tempfile::NamedTempFile;

use genimtools::common::models::{Region, RegionSet};
use genimtools::tokenizers::{Tokenizer, TreeTokenizer};

#[fixture]
fn path_to_data() -> &'static str {
    "tests/data"
}

#[fixture]
fn path_to_bed_file() -> &'static str {
    "tests/data/peaks.bed"
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
}

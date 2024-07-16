//!
//! # Common, core utilities for `gtars`
//! This module contains core utilities across the `gtars` crate. While possible, it's usually not interfaced with directly
//! unless interacting with any of the [models].
//!
//! ## Examples
//! ### Create region set
//! ```rust
//! use std::path::Path;
//! use gtars::common::models::RegionSet;
//!
//! let path_to_tokenize_bed_file = "tests/data/to_tokenize.bed";
//! let rs = RegionSet::try_from(Path::new(path_to_tokenize_bed_file)).unwrap();
//!
//! println!("{:?}", rs.regions);
//! ```
//!

pub mod consts;
pub mod models;
pub mod utils;

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use tempfile::NamedTempFile;

    use super::models::{Region, RegionSet};
    use super::utils::extract_regions_from_bed_file;
    use std::io::Read;
    use std::path::Path;

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
    fn path_to_anndata_file() -> &'static str {
        "tests/data/pbmc_hg38.h5ad"
    }

    #[fixture]
    fn path_to_r2v_repo() -> &'static str {
        "databio/r2v-luecken2021-hg38-v2"
    }

    #[fixture]
    fn bb_bed_id() -> &'static str {
        "fa09672b962809b408b356728d81640e"
    }

    #[fixture]
    fn path_to_gtok_file() -> &'static str {
        "tests/data/out/tokens.gtok"
    }

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
        assert_eq!(regions.len(), 25);
    }

    #[rstest]
    fn test_region_set_from_bed(path_to_bed_file: &str) {
        let path = Path::new(path_to_bed_file);
        let rs = RegionSet::try_from(path).unwrap();

        assert_eq!(rs.len(), 25);
    }

    #[rstest]
    fn test_region_set_from_bytes(path_to_bed_file: &str) {
        let path = Path::new(path_to_bed_file);
        let rs = RegionSet::try_from(path).unwrap();

        let mut bytes: Vec<u8> = Vec::new();

        std::fs::File::open(path)
            .unwrap()
            .read_to_end(&mut bytes)
            .unwrap();

        let rs2 = RegionSet::from(bytes.as_slice());

        assert_eq!(rs2.len(), rs.len());
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

        assert_eq!(rs2.len(), 25);
    }
}

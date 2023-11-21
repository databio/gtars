use genimtools::common::models::{Region, RegionSet};
use rstest::*;

use std::path::Path;

#[fixture]
fn path_to_bed_file() -> &'static str {
    "tests/data/peaks.bed"
}

mod tests {
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
    fn test_region_set_from_bed(path_to_bed_file: &str) {
        let path = Path::new(path_to_bed_file);
        let rs = RegionSet::try_from(path).unwrap();
        assert!(rs.regions.len() > 0);
    }
}
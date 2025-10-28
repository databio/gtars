use std::path::Path;
use std::fs::File;

use std::io::{BufWriter, Write};
use flate2::write::GzEncoder;
use flate2::Compression;

use gtars_core::models::RegionSet;


pub trait BedWrite {
    ///
    /// Write data to disk as bed file
    ///
    /// # Arguments
    /// - path: the path to the file to dump to
    fn write_bed<T: AsRef<Path>>(&self, path: T) -> std::io::Result<()>;

    ///
    /// Write data to disk as bed.gz file
    ///
    /// # Arguments
    /// - path: the path to the file to dump to
    fn write_bed_gz<T: AsRef<Path>>(&self, path: T) -> std::io::Result<()>;
}

impl BedWrite for RegionSet {
    fn write_bed<T: AsRef<Path>>(&self, path: T) -> std::io::Result<()> {
        let path = path.as_ref();

        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        }

        let mut file = File::create(path)?;

        for region in &self.regions {
            writeln!(file, "{}", region.as_string())?;
        }
        Ok(())
    }

    fn write_bed_gz<T: AsRef<Path>>(&self, path: T) -> std::io::Result<()> {
        let path = path.as_ref();

        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        }

        let file = File::create(path)?;
        let mut encoder = GzEncoder::new(BufWriter::new(file), Compression::best());

        for region in &self.regions {
            writeln!(encoder, "{}", region.as_string())?;
        }

        encoder.finish()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    use std::path::PathBuf;

    use rstest::*;

    fn get_test_path(file_name: &str) -> Result<PathBuf, Box<dyn std::error::Error>> {
        let file_path: PathBuf = std::env::current_dir()
            .unwrap()
            .join("../tests/data/regionset")
            .join(file_name);
        Ok(file_path)
    }

    #[rstest]
    fn test_save_bed() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        let tempdir = tempfile::tempdir().unwrap();

        let mut new_file_path = tempdir.keep();
        new_file_path.push("new_bedfile.bed");

        assert!(region_set.write_bed(new_file_path.as_path()).is_ok());

        let new_region = RegionSet::try_from(new_file_path.as_path());
        assert!(new_region.is_ok());
        assert_eq!(new_region.unwrap().identifier(), region_set.identifier())
    }

    #[rstest]
    fn test_save_bed_gz() {
        let file_path = get_test_path("dummy.narrowPeak.bed.gz").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        let tempdir = tempfile::tempdir().unwrap();

        let mut new_file_path = tempdir.keep();
        new_file_path.push("new_file.bed.gz");

        assert!(region_set.write_bed_gz(new_file_path.as_path()).is_ok());

        let new_region = RegionSet::try_from(new_file_path.as_path());
        assert!(new_region.is_ok());
        assert_eq!(new_region.unwrap().identifier(), region_set.identifier())
    }

}
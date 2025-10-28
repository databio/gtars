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

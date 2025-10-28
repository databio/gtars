use std::collections::HashMap;
use std::path::Path;

use gtars_core::models::RegionSet;
use gtars_core::utils::get_chrom_sizes;

use bigtools::beddata::BedParserStreamingIterator;
use bigtools::{BedEntry, BigBedWrite as BigBedWriter};
use tokio::runtime;

pub trait BigBedWrite {
    ///
    /// Write data as bigBed (binary version of bed file)
    ///
    /// # Arguments
    /// - out_path: the path to the bigbed file which should be created
    /// - chrom_size: the path to chrom sizes file
    ///
    fn write_bigbed<P: AsRef<Path>>(&self, out_path: P, chrom_size: P) -> Result<(), std::io::Error>;
}


impl BigBedWrite for RegionSet {
    fn write_bigbed<P: AsRef<Path>>(&self, out_path: P, chrom_size: P) -> Result<(), std::io::Error> {
        let out_path = out_path.as_ref();
        if out_path.exists() {
            println!("Bed file already exists. Overwriting existing file")
        }

        if let Some(parent) = out_path.parent() {
            std::fs::create_dir_all(parent)?;
        }
        let chrom_sizes: HashMap<String, u32> = get_chrom_sizes(chrom_size);

        let mut warnings_count: i32 = 0;
        let region_vector = self.regions.iter().map(|i| {
            // This if it is removing regions that are not in chrom sizes file.
            if !chrom_sizes.contains_key(&i.chr) {
                eprintln!(
                    "Warning:: Chromosome is not found in Chrom sizes file. Chr: '{}'",
                    i.chr
                );
                if warnings_count > 40 {
                    panic!("Incorrect chrom sizes provided. Unable to create bigBed file!");
                }

                warnings_count += 1;
                return None;
            }
            Some(Ok::<_, std::io::Error>((
                i.chr.clone(),
                BedEntry {
                    start: i.start,
                    end: i.end,
                    rest: i
                        .rest
                        .as_deref()
                        .map_or(String::new(), |s| format!("\t{}", s)),
                },
            )))
        });

        #[allow(clippy::option_filter_map)]
        // I like this because its more readable and clear whats going on
        let region_vector = region_vector.filter(|e| e.is_some()).map(|e| e.unwrap());

        let runtime = runtime::Builder::new_multi_thread()
            .worker_threads(
                std::thread::available_parallelism()
                    .map(|c| c.into())
                    .unwrap_or(1),
            )
            .build()
            .expect("Unable to create thread pool.");

        let mut bb_out = BigBedWriter::create_file(out_path, chrom_sizes.clone())
            .expect("Failed to create bigBed file.");

        bb_out.options.max_zooms = 8;

        let data = BedParserStreamingIterator::wrap_iter(region_vector.into_iter(), true);
        bb_out
            .write(data, runtime)
            .map_err(|e| std::io::Error::other(e.to_string()))?;
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
    fn test_save_bigbed() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        let chrom_sizes_path: PathBuf = std::env::current_dir()
            .unwrap()
            .join("../tests/data/regionset/dummy_chrom_sizes");

        let tempdir = tempfile::tempdir().unwrap();
        let mut new_file_path = tempdir.keep();
        new_file_path.push("new.bigbed");

        assert!(region_set
            .write_bigbed(new_file_path.as_path(), chrom_sizes_path.as_path())
            .is_ok());
    }
}
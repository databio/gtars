use anyhow::Result;
use std::fs::File;
use std::io::BufRead;
use std::path::{Path, PathBuf};

use md5::{Digest, Md5};

use flate2::write::GzEncoder;
use flate2::Compression;
use std::collections::HashMap;
use std::fmt::Debug;
use std::fmt::{self, Display};
use std::io::{BufWriter, Error, ErrorKind, Write};
use tokio::runtime;

use bigtools::beddata::BedParserStreamingIterator;
use bigtools::{BedEntry, BigBedWrite};

use crate::common::models::Region;
use crate::common::utils::{get_chrom_sizes, get_dynamic_reader};

///
/// RegionSet struct, the representation of the interval region set file,
/// such as bed file.
///
#[derive(Clone, Debug)]
pub struct RegionSet {
    pub regions: Vec<Region>,
    pub header: Option<String>,
    pub path: Option<PathBuf>,
}

pub struct RegionSetIterator<'a> {
    region_set: &'a RegionSet,
    index: usize,
}

impl TryFrom<&Path> for RegionSet {
    type Error = anyhow::Error;

    ///
    /// Create a new [RegionSet] from a bed file.
    ///
    /// # Arguments:
    /// - value: path to bed file on disk.
    fn try_from(value: &Path) -> Result<Self> {
        let path = value;

        let mut new_regions: Vec<Region> = Vec::new();

        let reader = match path.is_file() {
            true => get_dynamic_reader(path).expect("!Can't read file"),
            false => anyhow::bail!("File not found!")
        };

        let mut header: String = String::new();

        for line in reader.lines() {
            let string_line = line?;

            let parts: Vec<String> = string_line.split('\t').map(|s| s.to_string()).collect();

            if parts.len() < 3 {
                if string_line.starts_with("browser")
                    | string_line.starts_with("track")
                    | string_line.starts_with("#")
                {
                    header.push_str(&string_line);
                }
                continue;
            }

            new_regions.push(Region {
                chr: parts[0].to_string(),

                // To ensure that lines are regions, and we can parse it, we are using Result matching
                // And it helps to skip lines that are headers.
                start: parts[1].parse()?,
                end: parts[2].parse()?,
                rest: Some(parts[3..].join("\t")).filter(|s| !s.is_empty()),
            });
        }
        if new_regions.len() <= 1 {
            let new_error = Error::new(
                ErrorKind::Other,
                format!(
                    "Corrupted file. 0 regions found in the file: {}",
                    path.display()
                ),
            );
            return Err(new_error.into());
        }

        Ok(RegionSet {
            regions: new_regions,
            header: match header.is_empty() {
                true => None,
                false => Some(header),
            },
            path: Some(value.to_owned()),
        })
    }
}

impl TryFrom<&str> for RegionSet {
    type Error = anyhow::Error;

    fn try_from(value: &str) -> Result<Self> {
        RegionSet::try_from(Path::new(value))
    }
}

impl TryFrom<String> for RegionSet {
    type Error = anyhow::Error;

    fn try_from(value: String) -> Result<Self> {
        RegionSet::try_from(Path::new(&value))
    }
}

impl From<Vec<Region>> for RegionSet {
    fn from(regions: Vec<Region>) -> Self {
        let path = None;

        RegionSet {
            regions,
            header: None,
            path,
        }
    }
}

impl From<&[u8]> for RegionSet {
    fn from(value: &[u8]) -> Self {
        let region_str = String::from_utf8_lossy(value);
        let regions: Vec<Region> = region_str
            .split('\n')
            .map(|line| {
                let parts = line.split('\t').collect::<Vec<&str>>();

                let chr = parts[0].to_string();
                let start = parts[1].parse::<u32>().unwrap();
                let end = parts[2].parse::<u32>().unwrap();
                let rest = Some(parts[3..].join("\t")).filter(|s| !s.is_empty());

                Region {
                    chr,
                    start,
                    end,
                    rest,
                }
            })
            .collect();

        RegionSet {
            regions,
            header: None,
            path: None,
        }
    }
}

impl<'a> Iterator for RegionSetIterator<'a> {
    type Item = &'a Region;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.region_set.regions.len() {
            let region = &self.region_set.regions[self.index];
            self.index += 1;
            Some(region)
        } else {
            None
        }
    }
}

impl<'a> IntoIterator for &'a RegionSet {
    type Item = &'a Region;
    type IntoIter = RegionSetIterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        RegionSetIterator {
            region_set: self,
            index: 0,
        }
    }
}

impl RegionSet {
    ///
    /// Save a regionset to disk as bed file
    ///
    /// # Arguments
    /// - path: the path to the file to dump to
    pub fn to_bed(&self, path: &Path) -> std::io::Result<()> {
        if path.exists() {
            println!("Bed file already exists. Overwriting existing file")
        }

        let mut file = File::create(path).unwrap();

        for region in &self.regions {
            writeln!(
                file,
                "{}\t{}\t{}\t{}",
                region.chr,
                region.start,
                region.end,
                region.rest.as_deref().unwrap_or("")
            )?;
        }
        Ok(())
    }

    ///
    /// Save a regionset to disk as bed.gz file
    ///
    /// # Arguments
    /// - path: the path to the file to dump to
    pub fn to_bed_gz(&self, path: &Path) -> std::io::Result<()> {
        if path.exists() {
            println!("Bed file already exists. Overwriting existing file")
        }

        let file = File::create(path)?;
        let mut buffer: String = String::new();

        for region in &self.regions {
            buffer.push_str(&format!(
                "{}\t{}\t{}\t{}\n",
                region.chr,
                region.start,
                region.end,
                region.rest.as_deref().unwrap_or("")
            ));
        }

        let mut encoder = GzEncoder::new(BufWriter::new(file), Compression::fast());
        encoder.write_all(buffer.as_bytes())?;

        Ok(())
    }

    ///
    /// Calculate identifier for RegionSet
    ///
    /// This function doesn't sort file, and identifer is based on
    /// unsorted first 3 columns.
    ///
    /// # Returns
    /// String containing RegionSet identifier
    pub fn identifier(&self) -> String {
        let mut chrs = String::new();
        let mut starts = String::new();
        let mut ends = String::new();

        let mut first = true;
        for region in &self.regions {
            if !first {
                chrs.push(',');
                starts.push(',');
                ends.push(',');
            }
            chrs.push_str(&region.chr);
            starts.push_str(&region.start.to_string());
            ends.push_str(&region.end.to_string());

            first = false;
        }

        // Update hasher with input data

        let mut hasher = Md5::new();
        hasher.update(chrs);
        let chrom_hash = hasher.finalize();

        let mut hasher = Md5::new();
        hasher.update(starts);
        let start_hash = hasher.finalize();

        let mut hasher = Md5::new();
        hasher.update(ends);
        let end_hash = hasher.finalize();

        let combined = format!("{:x},{:x},{:x}", chrom_hash, start_hash, end_hash);

        let mut hasher = Md5::new();
        hasher.update(combined);
        let hash = hasher.finalize();
        let bed_digest: String = format!("{:x}", hash);

        bed_digest
    }

    ///
    /// Save RegionSet as bigBed (binary version of bed file)
    ///
    /// # Arguments
    /// - out_path: the path to the bigbed file which should be created
    /// - chrom_size: the path to chrom sizes file
    ///
    pub fn to_bigbed<T: AsRef<Path>>(&self, out_path: &Path, chrom_size: T) -> Result<()> {
        let chrom_sizes: HashMap<String, u32> = get_chrom_sizes(chrom_size);

        let region_vector = self.regions.iter().map(|i| {
            // This if is removing regions that are not in chrom sizes file.
            if !chrom_sizes.contains_key(&i.chr) {
                return None;
            }
            Some(Ok::<_, Error>((
                i.chr.clone(),
                BedEntry {
                    start: i.start,
                    end: i.end,
                    rest: i.rest.as_deref().unwrap_or("").to_string(),
                },
            )))
        });

        let region_vector = region_vector.filter(|e| e.is_some()).map(|e| e.unwrap());

        let runtime = runtime::Builder::new_multi_thread()
            .worker_threads(
                std::thread::available_parallelism()
                    .map(|c| c.into())
                    .unwrap_or(1),
            )
            .build()
            .expect("Unable to create thread pool.");

        let mut bb_out = BigBedWrite::create_file(out_path, chrom_sizes.clone())
            .expect("Failed to create bigBed file.");

        bb_out.options.max_zooms = 8;

        let data = BedParserStreamingIterator::wrap_iter(region_vector.into_iter(), true);
        bb_out.write(data, runtime)?;
        Ok(())
    }

    ///
    /// Sort bed file based on first 3 columns
    ///
    pub fn sort(&mut self) -> () {
        self.regions
            .sort_by(|a, b| a.chr.cmp(&b.chr).then_with(|| a.start.cmp(&b.start)));
    }

    ///
    /// Get number of regions in RegionSet
    ///
    /// Returns:
    /// number of regions
    pub fn len(&self) -> usize {
        self.regions.len()
    }
}

impl Display for RegionSet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "RegionSet with {} regions.", self.len())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn get_test_path(file_name: &str) -> Result<PathBuf, Error> {
        let file_path: PathBuf = std::env::current_dir()
            .unwrap()
            .join("tests/data/regionset")
            .join(file_name);
        Ok(file_path)
    }

    #[test]
    fn test_open_from_path() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        assert!(RegionSet::try_from(file_path.as_path()).is_ok());
    }

    #[test]
    fn test_open_from_string() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        assert!(RegionSet::try_from(file_path.to_str().unwrap()).is_ok());
    }

    #[test]
    fn test_open_bed_gz() {
        let file_path = get_test_path("dummy.narrowPeak.bed.gz").unwrap();
        assert!(RegionSet::try_from(file_path.to_str().unwrap()).is_ok());
    }

    #[test]
    fn test_calculate_identifier() {
        let file_path = get_test_path("dummy.narrowPeak.bed.gz").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        assert_eq!("f0b2cf73383b53bd97ff525a0380f200", region_set.identifier());
    }

    #[test]
    fn test_save_bed_gz() {
        let file_path = get_test_path("dummy.narrowPeak.bed.gz").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        let tempdir = tempfile::tempdir().unwrap();

        let mut new_file_path = PathBuf::try_from(tempdir.into_path()).unwrap();
        new_file_path.push("new_file.bed.gz");

        assert!(region_set.to_bed_gz(new_file_path.as_path()).is_ok());

        let new_region = RegionSet::try_from(new_file_path.as_path());
        assert!(new_region.is_ok());
        assert_eq!(new_region.unwrap().identifier(), region_set.identifier())
    }

    #[test]
    fn test_save_bed() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        let tempdir = tempfile::tempdir().unwrap();

        let mut new_file_path = PathBuf::try_from(tempdir.into_path()).unwrap();
        new_file_path.push("new_bedfile.bed");

        assert!(region_set.to_bed(new_file_path.as_path()).is_ok());

        let new_region = RegionSet::try_from(new_file_path.as_path());
        assert!(new_region.is_ok());
        assert_eq!(new_region.unwrap().identifier(), region_set.identifier())
    }

    #[test]
    fn test_save_bigbed() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        let chrom_sizes_path: PathBuf = std::env::current_dir()
            .unwrap()
            .join("tests/data/regionset/dummy_chrom_sizes");

        let tempdir = tempfile::tempdir().unwrap();
        let mut new_file_path = PathBuf::try_from(tempdir.into_path()).unwrap();
        new_file_path.push("new.bigbed");

        assert!(region_set
            .to_bigbed(new_file_path.as_path(), chrom_sizes_path.as_path())
            .is_ok());
    }

    #[test]
    fn test_read_headers() {
        let file_path = get_test_path("dummy_headers.bed").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        assert!(!region_set.header.is_none());
        assert_eq!(region_set.path.unwrap(), file_path);
    }
}

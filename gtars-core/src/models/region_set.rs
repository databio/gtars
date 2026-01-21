use anyhow::Result;
use std::fs::File;
use std::io::BufRead;
use std::path::{Path, PathBuf};

use md5::{Digest, Md5};

use flate2::write::GzEncoder;
use flate2::Compression;
use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::fmt::{self, Display};
use std::io::{BufWriter, Error, Write};
#[cfg(feature = "bigbed")]
use tokio::runtime;

#[cfg(feature = "bigbed")]
use bigtools::beddata::BedParserStreamingIterator;
#[cfg(feature = "bigbed")]
use bigtools::{BedEntry, BigBedWrite};

use crate::models::Region;
#[cfg(feature = "http")]
use crate::utils::get_dynamic_reader_from_url;
use crate::utils::get_dynamic_reader;
#[cfg(feature = "bigbed")]
use crate::utils::get_chrom_sizes;

#[cfg(feature = "dataframe")]
use polars::prelude::*;
#[cfg(feature = "dataframe")]
use std::io::Cursor;

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
            #[cfg(feature = "http")]
            false => {
                match get_dynamic_reader_from_url(path) {
                    Ok(reader) => reader,
                    Err(_) => {
                        // Extract bbid from the path (e.g., the file stem)
                        let bbid = path.to_str().ok_or_else(|| {
                            anyhow::anyhow!("BEDbase identifier is not valid UTF-8: {:?}", path)
                        })?;

                        let fallback_url = format!(
                            "https://api.bedbase.org/v1/files/files/{}/{}/{}.bed.gz",
                            &bbid[0..1],
                            &bbid[1..2],
                            bbid
                        );

                        let fallback_path = PathBuf::from(fallback_url);

                        get_dynamic_reader_from_url(&fallback_path)
                            .expect("!Can't get file from path, url, or BEDbase identifier")
                    }
                }
            }
            #[cfg(not(feature = "http"))]
            false => {
                return Err(anyhow::anyhow!(
                    "File not found and HTTP feature not enabled: {}",
                    path.display()
                ));
            }
        };

        let mut header: String = String::new();

        let mut first_line: bool = true;

        for line in reader.lines() {
            let string_line = line?;

            let parts: Vec<String> = string_line.split('\t').map(|s| s.to_string()).collect();

            if string_line.starts_with("browser")
                | string_line.starts_with("track")
                | string_line.starts_with("#")
            {
                header.push_str(&string_line);
                first_line = false;
                continue;
            }

            // Handling column headers like `chr start end etc` without #
            if first_line {
                if parts.len() >= 3 {
                    let is_header: bool = match parts[1].parse::<u32>() {
                        Ok(_num) => false,
                        Err(_) => true,
                    };
                    if is_header {
                        header.push_str(&string_line);
                        first_line = false;
                        continue;
                    }
                }
                first_line = false;
            }

            new_regions.push(Region {
                chr: parts[0].to_owned(),

                // To ensure that lines are regions, and we can parse it, we are using Result matching
                start: match parts[1].parse() {
                    Ok(start) => start,
                    Err(_err) => {
                        return Err(Error::other(format!(
                            "Error in parsing start position: {:?}",
                            parts
                        ))
                        .into())
                    }
                },
                end: match parts[2].parse() {
                    Ok(end) => end,
                    Err(_err) => {
                        return Err(anyhow::Error::from(Error::other(format!(
                            "Error in parsing end position: {:?}",
                            parts
                        ))))
                    }
                },
                rest: Some(parts[3..].join("\t")).filter(|s| !s.is_empty()),
            });
        }
        if new_regions.is_empty() {
            let new_error = Error::other(format!(
                "Corrupted file. 0 regions found in the file: {}",
                path.display()
            ));
            return Err(new_error.into());
        }

        let mut rs = RegionSet {
            regions: new_regions,
            header: match header.is_empty() {
                true => None,
                false => Some(header),
            },
            path: Some(value.to_owned()),
        };
        // This line needed for correct calculate identifier
        rs.sort();

        Ok(rs)
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
        // println!("Converting String to Path: {}", value);
        RegionSet::try_from(Path::new(&value))
    }
}

impl TryFrom<PathBuf> for RegionSet {
    type Error = anyhow::Error;

    fn try_from(value: PathBuf) -> Result<Self> {
        RegionSet::try_from(value.as_path())
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
    pub fn to_bed<T: AsRef<Path>>(&self, path: T) -> std::io::Result<()> {
        let path = path.as_ref();
        if path.exists() {
            println!("Bed file already exists. Overwriting existing file")
        }

        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        }

        let mut file = File::create(path).unwrap();

        for region in &self.regions {
            writeln!(file, "{}", region.as_string())?;
        }
        Ok(())
    }

    ///
    /// Save a regionset to disk as bed.gz file
    ///
    /// # Arguments
    /// - path: the path to the file to dump to
    pub fn to_bed_gz<T: AsRef<Path>>(&self, path: T) -> std::io::Result<()> {
        let path = path.as_ref();
        if path.exists() {
            println!("Bed file already exists. Overwriting existing file")
        }

        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        }

        let file = File::create(path)?;
        let mut buffer: String = String::new();

        for region in &self.regions {
            buffer.push_str(&format!("{}\n", region.as_string(),));
        }

        let mut encoder = GzEncoder::new(BufWriter::new(file), Compression::best());
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

    pub fn file_digest(&self) -> String {
        let mut buffer: String = String::new();
        for region in &self.regions {
            buffer.push_str(&format!("{}\n", region.as_string(),));
        }

        let mut hasher = Md5::new();

        hasher.update(buffer);
        let hash = hasher.finalize();
        let file_digest: String = format!("{:x}", hash);

        file_digest
    }

    ///
    /// Iterate unique chromosomes located in RegionSet
    ///
    pub fn iter_chroms(&self) -> impl Iterator<Item = &String> {
        let unique_chroms: HashSet<&String> = self.regions.iter().map(|r| &r.chr).collect();
        unique_chroms.into_iter()
    }

    ///
    /// Iterate through regions located on specific Chromosome in RegionSet
    ///
    /// # Arguments
    /// - chr: chromosome name
    ///
    pub fn iter_chr_regions<'a>(&'a self, chr: &'a str) -> impl Iterator<Item = &'a Region> {
        self.regions.iter().filter(move |r| r.chr == chr)
    }

    ///
    /// Save RegionSet as bigBed (binary version of bed file)
    ///
    /// # Arguments
    /// - out_path: the path to the bigbed file which should be created
    /// - chrom_size: the path to chrom sizes file
    ///
    #[cfg(feature = "bigbed")]
    pub fn to_bigbed<T: AsRef<Path>>(&self, out_path: T, chrom_size: T) -> Result<()> {
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
            Some(Ok::<_, Error>((
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

        let mut bb_out = BigBedWrite::create_file(out_path, chrom_sizes.clone())
            .expect("Failed to create bigBed file.");

        bb_out.options.max_zooms = 8;

        let data = BedParserStreamingIterator::wrap_iter(region_vector.into_iter(), true);
        bb_out.write(data, runtime)?;
        Ok(())
    }

    ///
    /// Sort bed file based on first 3 columns.
    /// Sorting is happening inside the object,
    /// where original order will be overwritten
    ///
    pub fn sort(&mut self) {
        self.regions
            .sort_by(|a, b| a.chr.cmp(&b.chr).then_with(|| a.start.cmp(&b.start)));
    }

    ///
    /// Is regionSet empty?
    ///
    pub fn is_empty(&self) -> bool {
        if self.regions.len() == 0 {
            return true;
        }
        false
    }

    ///
    /// Calculate all regions width
    ///
    pub fn region_widths(&self) -> Result<Vec<u32>> {
        let mut widths: Vec<u32> = Vec::new();

        for region in &self.regions {
            widths.push(region.width())
        }

        Ok(widths)
    }

    ///
    /// Calculate mean region width for whole RegionSet
    ///
    pub fn mean_region_width(&self) -> f64 {
        let sum: u32 = self
            .regions
            .iter()
            .map(|region| region.end - region.start)
            .sum();
        let count: u32 = self.regions.len() as u32;

        // must be f64 because python doesn't understand f32
        ((sum as f64 / count as f64) * 100.0).round() / 100.0
    }

    ///
    /// Calculate middle point for each region, and return hashmap with midpoints for each chromosome
    ///
    pub fn calc_mid_points(&self) -> HashMap<String, Vec<u32>> {
        let mut mid_points: HashMap<String, Vec<u32>> = HashMap::new();
        for chromosome in self.iter_chroms() {
            let mut chr_mid_points: Vec<u32> = Vec::new();
            for region in self.iter_chr_regions(chromosome) {
                chr_mid_points.push(region.mid_point());
            }
            mid_points.insert(chromosome.clone(), chr_mid_points);
        }
        mid_points
    }

    ///
    /// Get number of regions in RegionSet
    ///
    /// Returns:
    /// number of regions
    pub fn len(&self) -> usize {
        self.regions.len()
    }

    ///
    /// Get the furthest region location for each region
    ///
    pub fn get_max_end_per_chr(&self) -> HashMap<String, u32> {
        let mut result: HashMap<String, u32> = HashMap::new();

        let mut current_chr: &String = &self.regions[0].chr;
        let mut max_end: u32 = self.regions[0].end;

        for r in &self.regions[1..] {
            if &r.chr == current_chr {
                // Same chromosome → update max end
                max_end = max_end.max(r.end);
            } else {
                // Chromosome changed → store previous one
                result.insert(current_chr.clone(), max_end);
                current_chr = &r.chr;
                max_end = r.end;
            }
        }

        // Store the last chromosome
        result.insert(current_chr.clone(), max_end);

        result
    }

    ///
    /// Get total nucleotide count
    ///
    pub fn nucleotides_length(&self) -> u32 {
        let mut total_count: u32 = 0;
        for r in &self.regions {
            total_count += r.width();
        }
        total_count
    }

    ///
    /// Create Polars DataFrame
    ///
    #[cfg(feature = "dataframe")]
    pub fn to_polars(&self) -> PolarsResult<DataFrame> {
        // Convert regions to tab-separated string format
        let data: String = self
            .regions
            .iter()
            .map(|region| {
                if let Some(rest) = region.rest.as_deref() {
                    format!("{}\t{}\t{}\t{}", region.chr, region.start, region.end, rest,)
                } else {
                    format!("{}\t{}\t{}", region.chr, region.start, region.end,)
                }
            })
            .collect::<Vec<_>>()
            .join("\n");

        let cursor = Cursor::new(data);

        let df = CsvReadOptions::default()
            .with_has_header(false)
            .map_parse_options(|parse_options| parse_options.with_separator(b'\t'))
            .with_infer_schema_length(Some(10000))
            .into_reader_with_file_handle(cursor)
            .finish()?;

        Ok(df)
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

    use pretty_assertions::assert_eq;
    use rstest::*;

    fn get_test_path(file_name: &str) -> Result<PathBuf, Error> {
        let file_path: PathBuf = std::env::current_dir()
            .unwrap()
            .join("../tests/data/regionset")
            .join(file_name);
        Ok(file_path)
    }

    #[rstest]
    fn test_open_from_path() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        assert!(RegionSet::try_from(file_path.as_path()).is_ok());
    }

    #[rstest]
    fn test_open_from_string() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        assert!(RegionSet::try_from(file_path.to_str().unwrap()).is_ok());
    }

    #[rstest]
    fn test_open_from_url() {
        let file_path = String::from(
            "https://www.encodeproject.org/files/ENCFF321QPN/@@download/ENCFF321QPN.bed.gz",
        );
        assert!(RegionSet::try_from(file_path).is_ok());
    }

    #[rstest]
    #[ignore = "Avoid BEDbase dependency in CI"]
    fn test_open_from_bedbase() {
        let bbid = String::from("6b2e163a1d4319d99bd465c6c78a9741");
        let region_set = RegionSet::try_from(bbid);
        assert_eq!(region_set.is_ok(), true);
        assert_eq!(
            region_set.unwrap().identifier(),
            "6b2e163a1d4319d99bd465c6c78a9741"
        );
    }

    #[rstest]
    fn test_open_bed_gz() {
        let file_path = get_test_path("dummy.narrowPeak.bed.gz").unwrap();
        assert!(RegionSet::try_from(file_path.to_str().unwrap()).is_ok());
    }

    #[rstest]
    fn test_calculate_identifier() {
        let file_path = get_test_path("dummy.narrowPeak.bed.gz").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        assert_eq!("f0b2cf73383b53bd97ff525a0380f200", region_set.identifier());
    }

    #[rstest]
    fn test_save_bed_gz() {
        let file_path = get_test_path("dummy.narrowPeak.bed.gz").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        let tempdir = tempfile::tempdir().unwrap();

        let mut new_file_path = tempdir.keep();
        new_file_path.push("new_file.bed.gz");

        assert!(region_set.to_bed_gz(new_file_path.as_path()).is_ok());

        let new_region = RegionSet::try_from(new_file_path.as_path());
        assert!(new_region.is_ok());
        assert_eq!(new_region.unwrap().identifier(), region_set.identifier())
    }

    #[rstest]
    fn test_save_bed() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        let tempdir = tempfile::tempdir().unwrap();

        let mut new_file_path = tempdir.keep();
        new_file_path.push("new_bedfile.bed");

        assert!(region_set.to_bed(new_file_path.as_path()).is_ok());

        let new_region = RegionSet::try_from(new_file_path.as_path());
        assert!(new_region.is_ok());
        assert_eq!(new_region.unwrap().identifier(), region_set.identifier())
    }

    #[cfg(feature = "bigbed")]
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
            .to_bigbed(new_file_path.as_path(), chrom_sizes_path.as_path())
            .is_ok());
    }

    #[rstest]
    fn test_read_headers() {
        let file_path = get_test_path("dummy_headers.bed").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        assert!(region_set.header.is_some());
        assert_eq!(region_set.path.unwrap(), file_path);
    }

    #[rstest]
    fn test_is_empty() {
        let file_path = get_test_path("dummy_headers.bed").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        assert!(!region_set.is_empty());
    }

    #[rstest]
    fn test_file_digest() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        assert_eq!(region_set.file_digest(), "6224c4d40832b3e0889250f061e01120");
        assert_eq!(region_set.identifier(), "f0b2cf73383b53bd97ff525a0380f200")
    }

    #[rstest]
    fn test_mean_region_width() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        assert_eq!(region_set.mean_region_width(), 4.22)
    }
    #[rstest]
    fn test_open_file_with_incorrect_headers() {
        let file_path = get_test_path("dummy_incorrect_headers.bed").unwrap();
        let _region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();
    }

    #[rstest]
    fn test_chr_length() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();
        assert_eq!(*region_set.get_max_end_per_chr().get("chr1").unwrap(), 36);
        assert_eq!(region_set.get_max_end_per_chr().len(), 1)
    }

    #[rstest]
    fn test_total_nucleotides_function() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        assert_eq!(region_set.nucleotides_length(), 38)
    }

    #[rstest]
    fn test_iter_chroms() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        assert_eq!(region_set.iter_chroms().collect::<Vec<_>>().len(), 1)
    }

    #[cfg(feature = "dataframe")]
    #[rstest]
    fn test_polars() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();
        let rs_polars = region_set.to_polars().unwrap();
        println!("Number of columns: {:?}", rs_polars.get_columns().len());
        assert_eq!(rs_polars.get_columns().len(), 10);
    }

    #[rstest]
    fn test_calc_mid_points() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        let mid_points = region_set.calc_mid_points();
        assert_eq!(mid_points.get("chr1").unwrap().len(), 9);
        assert_eq!(mid_points.len(), 1);
        assert_eq!(
            mid_points
                .get("chr1")
                .unwrap()
                .iter()
                .min()
                .copied()
                .unwrap(),
            6u32
        );
    }
}

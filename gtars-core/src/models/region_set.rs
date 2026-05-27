use std::fs::File;
use std::io::BufRead;
use std::path::{Path, PathBuf};

use md5::{Digest, Md5};

use flate2::write::GzEncoder;
use flate2::Compression;
use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::fmt::{self, Display};
use std::io::{BufWriter, Write};
#[cfg(feature = "bigbed")]
use tokio::runtime;

#[cfg(feature = "bigbed")]
use bigtools::beddata::BedParserStreamingIterator;
#[cfg(feature = "bigbed")]
use bigtools::{BedEntry, BigBedWrite};

use crate::errors::RegionSetError;
use crate::models::Region;
#[cfg(feature = "bigbed")]
use crate::utils::get_chrom_sizes;
use crate::utils::get_dynamic_reader;
#[cfg(feature = "http")]
use crate::utils::get_dynamic_reader_from_url;

#[cfg(feature = "dataframe")]
use polars::prelude::*;
#[cfg(feature = "dataframe")]
use std::io::Cursor;

///
/// RegionSet struct, the representation of the interval region set file,
/// such as bed file.
///
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct RegionSet {
    pub regions: Vec<Region>,
    pub header: Option<String>,
    #[cfg_attr(feature = "serde", serde(skip))]
    pub path: Option<PathBuf>,
}

pub struct RegionSetIterator<'a> {
    region_set: &'a RegionSet,
    index: usize,
}

impl TryFrom<&Path> for RegionSet {
    type Error = RegionSetError;

    ///
    /// Create a new [RegionSet] from a bed file.
    ///
    /// # Arguments:
    /// - value: path to bed file on disk.
    fn try_from(value: &Path) -> Result<Self, RegionSetError> {
        let path = value;

        let mut new_regions: Vec<Region> = Vec::new();

        let reader = match path.is_file() {
            true => get_dynamic_reader(path)
                .map_err(|e| RegionSetError::FileReadError(e.to_string()))?,
            #[cfg(feature = "http")]
            false => {
                match get_dynamic_reader_from_url(path) {
                    Ok(reader) => reader,
                    Err(_) => {
                        return Err(RegionSetError::InvalidPathOrUrl(format!("{:?}", path)));

                        // // This code should be disabled, because it potentially breaks bedboss pipeline
                        // // BEDbase identifiers are 32-character MD5 hashes
                        // if bbid.len() != 32 {
                        //     return Err(RegionSetError::InvalidPathOrUrl(format!("{:?}", path)));
                        // }
                        //
                        // let fallback_url = format!(
                        //     "https://api.bedbase.org/v1/files/files/{}/{}/{}.bed.gz",
                        //     &bbid[0..1],
                        //     &bbid[1..2],
                        //     bbid
                        // );
                        //
                        // let fallback_path = PathBuf::from(fallback_url);
                        //
                        // get_dynamic_reader_from_url(&fallback_path)
                        //     .map_err(|e| RegionSetError::BedbaseFetchError(e.to_string()))?
                    }
                }
            }
            #[cfg(not(feature = "http"))]
            false => {
                return Err(RegionSetError::HttpFeatureDisabled(
                    path.display().to_string(),
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

            if parts.len() < 3 {
                return Err(RegionSetError::RegionParseError(format!(
                    "Error in parsing start position: {:?}",
                    parts
                )));
            }

            new_regions.push(Region {
                chr: parts[0].to_owned(),

                // To ensure that lines are regions, and we can parse it, we are using Result matching
                start: match parts[1].parse() {
                    Ok(start) => start,
                    Err(_err) => {
                        return Err(RegionSetError::RegionParseError(format!(
                            "Error in parsing start position: {:?}",
                            parts
                        )));
                    }
                },
                end: match parts[2].parse() {
                    Ok(end) => end,
                    Err(_err) => {
                        return Err(RegionSetError::RegionParseError(format!(
                            "Error in parsing end position: {:?}",
                            parts
                        )));
                    }
                },
                rest: Some(parts[3..].join("\t")).filter(|s| !s.is_empty()),
            });
        }
        if new_regions.is_empty() {
            return Err(RegionSetError::EmptyRegionSet(path.display().to_string()));
        }

        let mut rs = RegionSet {
            regions: new_regions,
            header: match header.is_empty() {
                true => None,
                false => Some(header),
            },
            path: Some(value.to_owned()),
        };
        // This line needed for correct calculate identifier and to bigbed function
        rs.sort();

        Ok(rs)
    }
}

impl TryFrom<&str> for RegionSet {
    type Error = RegionSetError;

    fn try_from(value: &str) -> Result<Self, RegionSetError> {
        RegionSet::try_from(Path::new(value))
    }
}

impl TryFrom<String> for RegionSet {
    type Error = RegionSetError;

    fn try_from(value: String) -> Result<Self, RegionSetError> {
        RegionSet::try_from(Path::new(&value))
    }
}

impl TryFrom<PathBuf> for RegionSet {
    type Error = RegionSetError;

    fn try_from(value: PathBuf) -> Result<Self, RegionSetError> {
        RegionSet::try_from(value.as_path())
    }
}

impl From<Vec<Region>> for RegionSet {
    fn from(regions: Vec<Region>) -> Self {
        RegionSet {
            regions,
            header: None,
            path: None,
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
        let mut seen = HashSet::new();
        let mut unique_chroms = Vec::new();
        for r in &self.regions {
            if seen.insert(&r.chr) {
                unique_chroms.push(&r.chr);
            }
        }
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
    pub fn to_bigbed<T: AsRef<Path>>(
        &self,
        out_path: T,
        chrom_size: T,
    ) -> Result<(), RegionSetError> {
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

        let mut bb_out = BigBedWrite::create_file(out_path, chrom_sizes.clone())
            .expect("Failed to create bigBed file.");

        bb_out.options.max_zooms = 8;

        let data = BedParserStreamingIterator::wrap_iter(region_vector.into_iter(), true);
        bb_out
            .write(data, runtime)
            .map_err(|e| RegionSetError::BigBedError(e.to_string()))?;
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
    pub fn region_widths(&self) -> Vec<u32> {
        self.regions.iter().map(|region| region.width()).collect()
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

    /// Calculate midpoints using the specified coordinate convention.
    ///
    /// See [`Region::mid_point_with_mode`] for details on how each mode computes the midpoint.
    pub fn calc_mid_points_with_mode(
        &self,
        mode: super::coords::CoordinateMode,
    ) -> HashMap<String, Vec<u32>> {
        let mut mid_points: HashMap<String, Vec<u32>> = HashMap::new();
        for chromosome in self.iter_chroms() {
            let mut chr_mid_points: Vec<u32> = Vec::new();
            for region in self.iter_chr_regions(chromosome) {
                chr_mid_points.push(region.mid_point_with_mode(mode));
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

// ── SortedRegionSet ─────────────────────────────────────────────────────
//
// A newtype wrapper that guarantees the inner RegionSet is sorted by (chr, start).

/// A RegionSet whose regions are sorted by (chr, start).
///
/// Created via `SortedRegionSet::new(rs)`, which sorts in place.
pub struct SortedRegionSet(pub RegionSet);

impl SortedRegionSet {
    /// Consume a RegionSet and sort it in place.
    pub fn new(mut rs: RegionSet) -> Self {
        rs.sort();
        SortedRegionSet(rs)
    }
}

// ── Structural interval operations (inherent methods) ───────────────────

impl RegionSet {
    /// Merge overlapping and adjacent intervals per chromosome.
    ///
    /// Sorts by (chr, start), then sweeps to merge intervals where
    /// `next.start <= current.end`. Returns a minimal set of non-overlapping regions.
    pub fn reduce(&self) -> RegionSet {
        if self.regions.is_empty() {
            return RegionSet::from(Vec::<Region>::new());
        }

        let sorted = SortedRegionSet::new(RegionSet::from(self.regions.clone()));
        let regions = &sorted.0.regions;

        let mut merged: Vec<Region> = Vec::new();
        let mut current = regions[0].clone();

        for r in &regions[1..] {
            if r.chr == current.chr && r.start <= current.end {
                current.end = current.end.max(r.end);
            } else {
                merged.push(Region {
                    chr: current.chr.clone(),
                    start: current.start,
                    end: current.end,
                    rest: None,
                });
                current = r.clone();
            }
        }
        merged.push(Region {
            chr: current.chr,
            start: current.start,
            end: current.end,
            rest: None,
        });

        RegionSet::from(merged)
    }

    /// Combine two region sets without merging overlapping intervals.
    pub fn concat(&self, other: &RegionSet) -> RegionSet {
        let mut regions = self.regions.clone();
        regions.extend(other.regions.iter().cloned());
        RegionSet::from(regions)
    }

    /// Merge two region sets into a minimal non-overlapping set.
    ///
    /// Equivalent to `self.concat(other).reduce()`.
    pub fn union(&self, other: &RegionSet) -> RegionSet {
        self.concat(other).reduce()
    }

    /// Clamp regions to chromosome boundaries.
    pub fn trim(&self, chrom_sizes: &HashMap<String, u32>) -> RegionSet {
        let regions: Vec<Region> = self
            .regions
            .iter()
            .filter_map(|r| {
                let chrom_size = chrom_sizes.get(&r.chr)?;
                let start = r.start.min(*chrom_size);
                let end = r.end.min(*chrom_size);
                if start > end {
                    None
                } else {
                    Some(Region {
                        chr: r.chr.clone(),
                        start,
                        end,
                        rest: None,
                    })
                }
            })
            .collect();
        RegionSet::from(regions)
    }

    /// Shift all regions by a fixed offset.
    pub fn shift(&self, offset: i64) -> RegionSet {
        let regions: Vec<Region> = self
            .regions
            .iter()
            .map(|r| {
                let start = (r.start as i64 + offset).max(0) as u32;
                let end = (r.end as i64 + offset).max(start as i64) as u32;
                Region {
                    chr: r.chr.clone(),
                    start,
                    end,
                    rest: None,
                }
            })
            .collect();
        RegionSet::from(regions)
    }

    /// Generate flanking regions.
    pub fn flank(&self, width: u32, use_start: bool, both: bool) -> RegionSet {
        let regions: Vec<Region> = self
            .regions
            .iter()
            .map(|r| {
                if both {
                    let anchor = if use_start { r.start } else { r.end };
                    Region {
                        chr: r.chr.clone(),
                        start: anchor.saturating_sub(width),
                        end: anchor.saturating_add(width),
                        rest: None,
                    }
                } else if use_start {
                    Region {
                        chr: r.chr.clone(),
                        start: r.start.saturating_sub(width),
                        end: r.start,
                        rest: None,
                    }
                } else {
                    Region {
                        chr: r.chr.clone(),
                        start: r.end,
                        end: r.end.saturating_add(width),
                        rest: None,
                    }
                }
            })
            .collect();
        RegionSet::from(regions)
    }

    /// Resize regions to a fixed width, anchored at start, end, or center.
    pub fn resize(&self, width: u32, fix: &str) -> RegionSet {
        let regions: Vec<Region> = self
            .regions
            .iter()
            .map(|r| match fix {
                "end" => Region {
                    chr: r.chr.clone(),
                    start: r.end.saturating_sub(width),
                    end: r.end,
                    rest: None,
                },
                "center" => {
                    let mid = r.start + (r.end - r.start) / 2;
                    let half = width / 2;
                    Region {
                        chr: r.chr.clone(),
                        start: mid.saturating_sub(half),
                        end: mid.saturating_sub(half).saturating_add(width),
                        rest: None,
                    }
                }
                _ => Region {
                    chr: r.chr.clone(),
                    start: r.start,
                    end: r.start.saturating_add(width),
                    rest: None,
                },
            })
            .collect();
        RegionSet::from(regions)
    }

    /// Narrow each region by specifying a relative sub-range within it.
    pub fn narrow(&self, start: Option<u32>, end: Option<u32>, width: Option<u32>) -> RegionSet {
        let regions: Vec<Region> = self
            .regions
            .iter()
            .map(|r| {
                let region_width = r.end - r.start;
                let (rel_start, rel_end) = match (start, end, width) {
                    (Some(s), Some(e), None) => (s.saturating_sub(1), e),
                    (Some(s), None, Some(w)) => (s.saturating_sub(1), s.saturating_sub(1) + w),
                    (None, Some(e), Some(w)) => (e.saturating_sub(w), e),
                    _ => (0, region_width),
                };
                let new_start = r.start + rel_start.min(region_width);
                let new_end = r.start + rel_end.min(region_width);
                Region {
                    chr: r.chr.clone(),
                    start: new_start.min(new_end),
                    end: new_end.max(new_start),
                    rest: None,
                }
            })
            .collect();
        RegionSet::from(regions)
    }

    /// Generate promoter regions relative to each region's start position.
    pub fn promoters(&self, upstream: u32, downstream: u32) -> RegionSet {
        let regions: Vec<Region> = self
            .regions
            .iter()
            .map(|r| Region {
                chr: r.chr.clone(),
                start: r.start.saturating_sub(upstream),
                end: r.start.saturating_add(downstream),
                rest: None,
            })
            .collect();
        RegionSet::from(regions)
    }

    /// Pairwise intersection of two region sets by index position.
    pub fn pintersect(&self, other: &RegionSet) -> RegionSet {
        let regions: Vec<Region> = self
            .regions
            .iter()
            .zip(other.regions.iter())
            .map(|(a, b)| {
                if a.chr != b.chr {
                    return Region {
                        chr: a.chr.clone(),
                        start: a.start,
                        end: a.start,
                        rest: None,
                    };
                }
                let start = a.start.max(b.start);
                let end = a.end.min(b.end);
                if start >= end {
                    Region {
                        chr: a.chr.clone(),
                        start,
                        end: start,
                        rest: None,
                    }
                } else {
                    Region {
                        chr: a.chr.clone(),
                        start,
                        end,
                        rest: None,
                    }
                }
            })
            .collect();
        RegionSet::from(regions)
    }

    /// Break all regions into non-overlapping disjoint pieces.
    pub fn disjoin(&self) -> RegionSet {
        let mut by_chr: HashMap<String, Vec<u32>> = HashMap::new();
        for r in &self.regions {
            by_chr.entry(r.chr.clone()).or_default().push(r.start);
            by_chr.entry(r.chr.clone()).or_default().push(r.end);
        }

        let mut result: Vec<Region> = Vec::new();
        for (chr, mut boundaries) in by_chr {
            boundaries.sort();
            boundaries.dedup();
            for window in boundaries.windows(2) {
                result.push(Region {
                    chr: chr.clone(),
                    start: window[0],
                    end: window[1],
                    rest: None,
                });
            }
        }
        result.sort_by(|a, b| (&a.chr, a.start).cmp(&(&b.chr, b.start)));
        RegionSet::from(result)
    }

    /// Cluster nearby regions.
    pub fn cluster(&self, max_gap: u32) -> Vec<u32> {
        if self.regions.is_empty() {
            return vec![];
        }

        let n = self.regions.len();
        let mut result = vec![0u32; n];

        // Create sorted indices to preserve original order mapping
        let mut sorted_indices: Vec<usize> = (0..n).collect();
        sorted_indices.sort_by(|&i, &j| {
            self.regions[i]
                .chr
                .cmp(&self.regions[j].chr)
                .then(self.regions[i].start.cmp(&self.regions[j].start))
                .then(self.regions[i].end.cmp(&self.regions[j].end))
        });

        let mut cluster_id: u32 = 0;
        let mut cluster_end = self.regions[sorted_indices[0]].end;
        let mut current_chr = &self.regions[sorted_indices[0]].chr;
        result[sorted_indices[0]] = cluster_id;

        for &idx in &sorted_indices[1..] {
            let r = &self.regions[idx];
            if r.chr != *current_chr || r.start > cluster_end.saturating_add(max_gap) {
                cluster_id += 1;
                cluster_end = r.end;
                current_chr = &r.chr;
            } else {
                cluster_end = cluster_end.max(r.end);
            }
            result[idx] = cluster_id;
        }

        result
    }

    /// Find the nearest region in `other` for each region in `self`.
    pub fn closest(&self, other: &RegionSet) -> Vec<(usize, usize, i64)> {
        let sorted_other = SortedRegionSet::new(RegionSet::from(other.regions.clone()));

        let mut by_chr: HashMap<&str, Vec<(usize, &Region)>> = HashMap::new();
        for (i, r) in sorted_other.0.regions.iter().enumerate() {
            by_chr.entry(&r.chr).or_default().push((i, r));
        }

        let mut results = Vec::new();

        for (self_idx, self_region) in self.regions.iter().enumerate() {
            if let Some(chr_regions) = by_chr.get(self_region.chr.as_str()) {
                let mut best_dist = i64::MAX;
                let mut best_idx = 0usize;

                for &(other_idx, other_region) in chr_regions {
                    let dist = if self_region.end <= other_region.start {
                        (other_region.start - self_region.end) as i64
                    } else if other_region.end <= self_region.start {
                        (self_region.start - other_region.end) as i64
                    } else {
                        0
                    };

                    if dist < best_dist {
                        best_dist = dist;
                        best_idx = other_idx;
                    }
                }

                results.push((self_idx, best_idx, best_dist));
            }
        }

        results
    }
}

// ── Sweep-line helpers ──────────────────────────────────────────────────

/// Per-chromosome set difference using a sweep-line algorithm.
pub fn sweep_setdiff_chr(chr: &str, a: &[Region], b: &[Region]) -> Vec<Region> {
    let mut result = Vec::new();
    let mut b_idx = 0;

    for a_region in a {
        while b_idx < b.len() && b[b_idx].end <= a_region.start {
            b_idx += 1;
        }

        let mut pos = a_region.start;
        let mut j = b_idx;

        while j < b.len() && b[j].start < a_region.end && pos < a_region.end {
            if b[j].start > pos {
                result.push(Region {
                    chr: chr.to_string(),
                    start: pos,
                    end: b[j].start,
                    rest: None,
                });
            }
            pos = pos.max(b[j].end);
            j += 1;
        }

        if pos < a_region.end {
            result.push(Region {
                chr: chr.to_string(),
                start: pos,
                end: a_region.end,
                rest: None,
            });
        }
    }

    result
}

/// Per-chromosome intersection using a sweep-line algorithm.
pub fn sweep_intersect_chr(chr: &str, a: &[Region], b: &[Region]) -> Vec<Region> {
    let mut result = Vec::new();
    let mut b_idx = 0;

    for a_region in a {
        while b_idx < b.len() && b[b_idx].end <= a_region.start {
            b_idx += 1;
        }
        let mut j = b_idx;
        while j < b.len() && b[j].start < a_region.end {
            let start = a_region.start.max(b[j].start);
            let end = a_region.end.min(b[j].end);
            if start < end {
                result.push(Region {
                    chr: chr.to_string(),
                    start,
                    end,
                    rest: None,
                });
            }
            j += 1;
        }
    }

    result
}

// ── IntervalSetOps trait ────────────────────────────────────────────────

/// Two-set interval operations on genomic region sets.
///
/// Provides set algebra (setdiff, intersect) and similarity metrics
/// (jaccard, coverage, overlap_coefficient). Implementations may use
/// sweep-line or index-based algorithms.
pub trait IntervalSetOps {
    /// Set difference: remove portions of `self` that overlap with `other`.
    fn setdiff(&self, other: &RegionSet) -> RegionSet;

    /// Range-level intersection: positions covered by *both* sets.
    fn intersect(&self, other: &RegionSet) -> RegionSet;

    /// Nucleotide-level Jaccard similarity: `|intersection| / |union|`.
    fn jaccard(&self, other: &RegionSet) -> f64;

    /// Fraction of self's base pairs covered by other.
    fn coverage(&self, other: &RegionSet) -> f64;

    /// Overlap coefficient: `|intersection| / min(|self|, |other|)`.
    fn overlap_coefficient(&self, other: &RegionSet) -> f64;
}

impl IntervalSetOps for RegionSet {
    fn setdiff(&self, other: &RegionSet) -> RegionSet {
        let a = self.reduce();
        let b = other.reduce();

        let mut b_by_chr: HashMap<String, Vec<Region>> = HashMap::new();
        for r in &b.regions {
            b_by_chr.entry(r.chr.clone()).or_default().push(r.clone());
        }

        let mut result: Vec<Region> = Vec::new();

        let mut a_chr_start = 0;
        while a_chr_start < a.regions.len() {
            let chr = &a.regions[a_chr_start].chr;
            let mut a_chr_end = a_chr_start;
            while a_chr_end < a.regions.len() && a.regions[a_chr_end].chr == *chr {
                a_chr_end += 1;
            }

            let empty_vec = vec![];
            let b_chr = b_by_chr.get(chr.as_str()).unwrap_or(&empty_vec);
            result.extend(sweep_setdiff_chr(chr, &a.regions[a_chr_start..a_chr_end], b_chr));

            a_chr_start = a_chr_end;
        }

        RegionSet::from(result)
    }

    fn intersect(&self, other: &RegionSet) -> RegionSet {
        let a = self.reduce();
        let b = other.reduce();

        let mut b_by_chr: HashMap<String, Vec<Region>> = HashMap::new();
        for r in &b.regions {
            b_by_chr.entry(r.chr.clone()).or_default().push(r.clone());
        }

        let mut result: Vec<Region> = Vec::new();

        let mut a_i = 0;
        while a_i < a.regions.len() {
            let chr = &a.regions[a_i].chr;
            let mut a_end = a_i;
            while a_end < a.regions.len() && a.regions[a_end].chr == *chr {
                a_end += 1;
            }

            if let Some(b_chr) = b_by_chr.get(chr.as_str()) {
                result.extend(sweep_intersect_chr(chr, &a.regions[a_i..a_end], b_chr));
            }

            a_i = a_end;
        }

        RegionSet::from(result)
    }

    fn jaccard(&self, other: &RegionSet) -> f64 {
        let a_bp = self.reduce().nucleotides_length();
        let b_bp = other.reduce().nucleotides_length();
        let union_bp = self.union(other).nucleotides_length();
        if union_bp == 0 {
            return 0.0;
        }
        let intersection_bp = a_bp + b_bp - union_bp;
        intersection_bp as f64 / union_bp as f64
    }

    fn coverage(&self, other: &RegionSet) -> f64 {
        let self_reduced = self.reduce();
        let self_bp = self_reduced.nucleotides_length();
        if self_bp == 0 {
            return 0.0;
        }
        let diff = self_reduced.setdiff(other);
        let diff_bp = diff.nucleotides_length();
        1.0 - (diff_bp as f64 / self_bp as f64)
    }

    fn overlap_coefficient(&self, other: &RegionSet) -> f64 {
        let a_bp = self.reduce().nucleotides_length();
        let b_bp = other.reduce().nucleotides_length();
        let min_bp = a_bp.min(b_bp);
        if min_bp == 0 {
            return 0.0;
        }
        let union_bp = self.union(other).nucleotides_length();
        let intersection_bp = a_bp + b_bp - union_bp;
        intersection_bp as f64 / min_bp as f64
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

    fn get_test_path(file_name: &str) -> PathBuf {
        let file_path: PathBuf = std::env::current_dir()
            .unwrap()
            .join("../tests/data/regionset")
            .join(file_name);
        file_path
    }

    #[rstest]
    fn test_open_from_path() {
        let file_path = get_test_path("dummy.narrowPeak");
        assert!(RegionSet::try_from(file_path.as_path()).is_ok());
    }

    #[rstest]
    fn test_open_from_string() {
        let file_path = get_test_path("dummy.narrowPeak");
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
        let file_path = get_test_path("dummy.narrowPeak.bed.gz");
        assert!(RegionSet::try_from(file_path.to_str().unwrap()).is_ok());
    }

    #[rstest]
    fn test_calculate_identifier() {
        let file_path = get_test_path("dummy.narrowPeak.bed.gz");
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        assert_eq!("f0b2cf73383b53bd97ff525a0380f200", region_set.identifier());
    }

    #[rstest]
    fn test_save_bed_gz() {
        let file_path = get_test_path("dummy.narrowPeak.bed.gz");
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
        let file_path = get_test_path("dummy.narrowPeak");
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
        let file_path = get_test_path("dummy.narrowPeak");
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
        let file_path = get_test_path("dummy_headers.bed");
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        assert!(region_set.header.is_some());
        assert_eq!(region_set.path.unwrap(), file_path);
    }

    #[rstest]
    fn test_is_empty() {
        let file_path = get_test_path("dummy_headers.bed");
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        assert!(!region_set.is_empty());
    }

    #[rstest]
    fn test_file_digest() {
        let file_path = get_test_path("dummy.narrowPeak");
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        assert_eq!(region_set.file_digest(), "6224c4d40832b3e0889250f061e01120");
        assert_eq!(region_set.identifier(), "f0b2cf73383b53bd97ff525a0380f200")
    }

    #[rstest]
    fn test_mean_region_width() {
        let file_path = get_test_path("dummy.narrowPeak");
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        assert_eq!(region_set.mean_region_width(), 4.22)
    }
    #[rstest]
    fn test_open_file_with_incorrect_headers() {
        let file_path = get_test_path("dummy_incorrect_headers.bed");
        let _region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();
    }

    #[rstest]
    fn test_chr_length() {
        let file_path = get_test_path("dummy.narrowPeak");
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();
        assert_eq!(*region_set.get_max_end_per_chr().get("chr1").unwrap(), 36);
        assert_eq!(region_set.get_max_end_per_chr().len(), 1)
    }

    #[rstest]
    fn test_total_nucleotides_function() {
        let file_path = get_test_path("dummy.narrowPeak");
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        assert_eq!(region_set.nucleotides_length(), 38)
    }

    #[rstest]
    fn test_iter_chroms() {
        let file_path = get_test_path("dummy.narrowPeak");
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        assert_eq!(region_set.iter_chroms().collect::<Vec<_>>().len(), 1)
    }

    #[cfg(feature = "dataframe")]
    #[rstest]
    fn test_polars() {
        let file_path = get_test_path("dummy.narrowPeak");
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();
        let rs_polars = region_set.to_polars().unwrap();
        println!("Number of columns: {:?}", rs_polars.get_columns().len());
        assert_eq!(rs_polars.get_columns().len(), 10);
    }

    #[rstest]
    fn test_calc_mid_points() {
        let file_path = get_test_path("dummy.narrowPeak");
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

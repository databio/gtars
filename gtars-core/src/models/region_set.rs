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

/// A `RegionSet` that is guaranteed to be sorted by (chr, start).
///
/// Created by moving a `RegionSet` into `SortedRegionSet::new()`, which
/// sorts in place (no clone). Functions that require sorted input can
/// accept this type instead of re-sorting every time.
pub struct SortedRegionSet(pub RegionSet);

impl SortedRegionSet {
    /// Sort a RegionSet in place and wrap it. This is a move, not a clone.
    pub fn new(mut rs: RegionSet) -> Self {
        rs.sort();
        Self(rs)
    }
}

// ── Structural interval operations (pure coordinate math) ───────────────
//
// These are inherent methods on RegionSet — no trait import needed.
// All operations return new RegionSet instances (immutable pattern) with
// 0-based half-open coordinates. Metadata (`rest` field) is not preserved.

impl RegionSet {
    /// Clamp regions to chromosome boundaries.
    ///
    /// Regions extending past chromosome ends are trimmed to `[0, chrom_size)`.
    /// Regions on chromosomes not present in `chrom_sizes` are dropped.
    /// Empty regions (start > end after clamping) are dropped; zero-width
    /// regions where start == end are kept.
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

    /// Generate promoter regions relative to each region's start position.
    ///
    /// For each region, produces `[start - upstream, start + downstream)`.
    /// Uses saturating subtraction at coordinate 0.
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
    ///
    /// Clones regions from both sets into a single `RegionSet`. No sorting,
    /// deduplication, or merging is performed.
    pub fn concat(&self, other: &RegionSet) -> RegionSet {
        let mut regions = self.regions.clone();
        regions.extend(other.regions.iter().cloned());
        RegionSet::from(regions)
    }

    /// Shift all regions by a fixed offset.
    ///
    /// Adds `offset` to both start and end of every region. Negative offsets
    /// use saturating subtraction at coordinate 0.
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
    ///
    /// - `use_start = true`: flank upstream of start -> `[start - width, start)`
    /// - `use_start = false`: flank downstream of end -> `[end, end + width)`
    /// - `both = true`: flank on both sides -> `[anchor - width, anchor + width)`
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
    ///
    /// - `"start"`: keep start, set end = start + width
    /// - `"end"`: keep end, set start = end - width
    /// - `"center"`: keep midpoint, expand symmetrically
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
                    let mid = (r.start + r.end) / 2;
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
    ///
    /// Parameters are 1-based relative positions within each region (matching
    /// GenomicRanges convention). Exactly two of the three must be provided.
    pub fn narrow(&self, start: Option<u32>, end: Option<u32>, width: Option<u32>) -> RegionSet {
        let regions: Vec<Region> = self
            .regions
            .iter()
            .map(|r| {
                let region_width = r.end - r.start;
                let (rel_start, rel_end) = match (start, end, width) {
                    (Some(s), Some(e), None) => (s - 1, e),
                    (Some(s), None, Some(w)) => (s - 1, (s - 1) + w),
                    (None, Some(e), Some(w)) => (e.saturating_sub(w), e),
                    (None, None, Some(w)) => (0, w),
                    _ => (0, region_width),
                };
                let abs_start = r.start + rel_start.min(region_width);
                let abs_end = r.start + rel_end.min(region_width);
                Region {
                    chr: r.chr.clone(),
                    start: abs_start,
                    end: abs_end.max(abs_start),
                    rest: None,
                }
            })
            .collect();
        RegionSet::from(regions)
    }

    /// Break all regions into non-overlapping disjoint pieces.
    ///
    /// Every boundary in the input becomes a boundary in the output. The result
    /// tiles the covered positions exactly, with no overlaps.
    pub fn disjoin(&self) -> RegionSet {
        if self.regions.is_empty() {
            return RegionSet::from(Vec::<Region>::new());
        }

        let sorted = SortedRegionSet::new(RegionSet::from(self.regions.clone()));
        let regions = &sorted.0.regions;

        let mut result: Vec<Region> = Vec::new();

        let mut i = 0;
        while i < regions.len() {
            let chr = &regions[i].chr;
            let mut chr_end = i;
            while chr_end < regions.len() && regions[chr_end].chr == *chr {
                chr_end += 1;
            }

            let mut events: Vec<u32> = Vec::with_capacity((chr_end - i) * 2);
            for r in &regions[i..chr_end] {
                events.push(r.start);
                events.push(r.end);
            }
            events.sort_unstable();
            events.dedup();

            for w in events.windows(2) {
                let seg_start = w[0];
                let seg_end = w[1];
                let covered = regions[i..chr_end]
                    .iter()
                    .any(|r| r.start <= seg_start && r.end >= seg_end);
                if covered {
                    result.push(Region {
                        chr: chr.clone(),
                        start: seg_start,
                        end: seg_end,
                        rest: None,
                    });
                }
            }

            i = chr_end;
        }

        RegionSet::from(result)
    }

    /// Return the gaps between regions per chromosome.
    ///
    /// Reduces the input first, then emits intervals between consecutive
    /// reduced regions on each chromosome.
    pub fn gaps(&self) -> RegionSet {
        let reduced = self.reduce();
        if reduced.regions.is_empty() {
            return RegionSet::from(Vec::<Region>::new());
        }

        let mut result: Vec<Region> = Vec::new();

        let mut i = 0;
        while i < reduced.regions.len() {
            let chr = &reduced.regions[i].chr;
            let mut j = i + 1;
            while j < reduced.regions.len() && reduced.regions[j].chr == *chr {
                let gap_start = reduced.regions[j - 1].end;
                let gap_end = reduced.regions[j].start;
                if gap_start < gap_end {
                    result.push(Region {
                        chr: chr.clone(),
                        start: gap_start,
                        end: gap_end,
                        rest: None,
                    });
                }
                j += 1;
            }
            i = j.max(i + 1);
        }

        RegionSet::from(result)
    }

    /// Pairwise intersection of two region sets by index position.
    ///
    /// For each pair at the same index, computes `[max(a.start, b.start), min(a.end, b.end))`.
    /// Pairs with no overlap or mismatched chromosomes produce zero-width regions.
    /// If the region sets differ in length, intersections are computed only up
    /// to the length of the shorter set.
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

    /// Merge two region sets into a minimal non-overlapping set.
    pub fn union(&self, other: &RegionSet) -> RegionSet {
        self.concat(other).reduce()
    }

    /// Alias for `setdiff`.
    pub fn subtract(&self, other: &RegionSet) -> RegionSet {
        self.setdiff(other)
    }

    /// Find the nearest region in `other` for each region in `self`.
    ///
    /// Returns a vec of `(self_index, other_index, signed_distance)` tuples,
    /// one per region in `self` that shares a chromosome with at least one
    /// region in `other`.
    ///
    /// The algorithm sorts candidates by `start`, builds a prefix-max array
    /// of `end` values, then scans rightward and leftward from the binary
    /// search insertion point with early termination.
    pub fn closest(&self, other: &RegionSet) -> Vec<(usize, usize, i64)> {
        if other.regions.is_empty() {
            return Vec::new();
        }

        // Group other's regions by chromosome, sorted by start.
        // Each entry also carries a prefix-max of end values for pruning.
        let mut other_by_chr: HashMap<String, (Vec<(usize, &Region)>, Vec<u32>)> = HashMap::new();
        for (idx, r) in other.regions.iter().enumerate() {
            other_by_chr
                .entry(r.chr.clone())
                .or_insert_with(|| (Vec::new(), Vec::new()))
                .0
                .push((idx, r));
        }
        for (candidates, max_ends) in other_by_chr.values_mut() {
            candidates.sort_by_key(|(_, r)| r.start);
            // Build prefix-max of end values: max_ends[i] = max(candidates[0..=i].end)
            let mut running_max = 0u32;
            max_ends.clear();
            max_ends.reserve(candidates.len());
            for &(_, r) in candidates.iter() {
                running_max = running_max.max(r.end);
                max_ends.push(running_max);
            }
        }

        let mut result: Vec<(usize, usize, i64)> = Vec::new();

        for (self_idx, a) in self.regions.iter().enumerate() {
            let Some((candidates, max_ends)) = other_by_chr.get(&a.chr) else {
                continue;
            };

            let ins = candidates
                .binary_search_by_key(&a.start, |(_, r)| r.start)
                .unwrap_or_else(|x| x);

            let mut best_other_idx = 0usize;
            let mut best_dist = i64::MAX;

            // Scan rightward from insertion point
            for j in ins..candidates.len() {
                let (other_idx, b) = candidates[j];
                let dist = a.distance_to(b);
                if dist.abs() < best_dist.abs() {
                    best_dist = dist;
                    best_other_idx = other_idx;
                }
                if best_dist == 0 {
                    break;
                }
                // All further candidates start even farther right
                if b.start >= a.end && (b.start as i64 - a.end as i64) >= best_dist.abs() {
                    break;
                }
            }

            // Scan leftward from insertion point
            if best_dist != 0 {
                for j in (0..ins).rev() {
                    let (other_idx, b) = candidates[j];
                    let dist = a.distance_to(b);
                    if dist.abs() < best_dist.abs() {
                        best_dist = dist;
                        best_other_idx = other_idx;
                    }
                    if best_dist == 0 {
                        break;
                    }
                    // Prune using prefix-max: if max_ends[j] <= a.start, then no
                    // candidate at index 0..=j overlaps the query. Check if the
                    // gap from query start to the farthest-reaching end among
                    // 0..=j is already worse than best.
                    if max_ends[j] <= a.start
                        && (a.start as i64 - max_ends[j] as i64) >= best_dist.abs()
                    {
                        break;
                    }
                }
            }

            result.push((self_idx, best_other_idx, best_dist));
        }

        result
    }

    /// Cluster nearby regions within `max_gap` distance.
    pub fn cluster(&self, max_gap: u32) -> Vec<u32> {
        if self.regions.is_empty() {
            return Vec::new();
        }

        let n = self.regions.len();
        let mut result = vec![0u32; n];

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
            if r.chr != *current_chr
                || r.start > cluster_end.saturating_add(max_gap)
            {
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
}

// ── Shared sweep-line helpers ────────────────────────────────────────────
//
// Free functions that operate on per-chromosome slices of *already-reduced*
// (sorted, non-overlapping) regions.  Both RegionSet and MCO call these so
// the sweep-line logic lives in exactly one place.

/// Sweep-line set-difference on a single chromosome.
///
/// Both `a` and `b` must be sorted by start and non-overlapping (reduced).
/// Returns the portions of `a` that do **not** overlap any region in `b`.
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

/// Sweep-line intersection on a single chromosome.
///
/// Both `a` and `b` must be sorted by start and non-overlapping (reduced).
/// Returns regions covered by **both** `a` and `b`.
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
//
// Two-set operations defined as a trait so both RegionSet (sweep-line) and
// MultiChromOverlapper (index-based) can implement them.

/// Two-set interval operations on genomic region sets.
///
/// Provides set algebra (setdiff, intersect) and similarity metrics
/// (jaccard, coverage, overlap_coefficient). Implementations may use
/// sweep-line or index-based algorithms.
///
/// Note: `union`, `subtract`, `closest`, and `cluster` are inherent methods
/// on `RegionSet` rather than trait methods.
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

// ── Sweep-line implementation of IntervalSetOps for RegionSet ───────────

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

    // ── Helpers for structural / IntervalSetOps tests ────────────────

    fn make_region(chr: &str, start: u32, end: u32) -> Region {
        Region {
            chr: chr.to_string(),
            start,
            end,
            rest: None,
        }
    }

    fn make_regionset(regions: Vec<(&str, u32, u32)>) -> RegionSet {
        let regions: Vec<Region> = regions
            .into_iter()
            .map(|(chr, start, end)| make_region(chr, start, end))
            .collect();
        RegionSet::from(regions)
    }

    // ── SortedRegionSet tests ───────────────────────────────────────

    #[test]
    fn test_sorted_regionset_sorts_in_place() {
        let regions = vec![
            Region { chr: "chr1".into(), start: 100, end: 200, rest: None },
            Region { chr: "chr1".into(), start: 10, end: 20, rest: None },
            Region { chr: "chr2".into(), start: 5, end: 15, rest: None },
        ];
        let sorted = SortedRegionSet::new(RegionSet::from(regions));
        let starts: Vec<u32> = sorted.0.regions.iter().map(|r| r.start).collect();
        assert_eq!(starts, vec![10, 100, 5]);
    }

    // ── trim tests ──────────────────────────────────────────────────

    #[rstest]
    fn test_trim_clamps_past_chrom_end() {
        let rs = make_regionset(vec![("chr1", 90, 150)]);
        let chrom_sizes: HashMap<String, u32> =
            [("chr1".to_string(), 100)].into_iter().collect();
        let trimmed = rs.trim(&chrom_sizes);
        assert_eq!(trimmed.regions.len(), 1);
        assert_eq!(trimmed.regions[0].start, 90);
        assert_eq!(trimmed.regions[0].end, 100);
    }

    #[rstest]
    fn test_trim_drops_unknown_chrom() {
        let rs = make_regionset(vec![("chrX", 0, 50), ("chr1", 10, 20)]);
        let chrom_sizes: HashMap<String, u32> =
            [("chr1".to_string(), 100)].into_iter().collect();
        let trimmed = rs.trim(&chrom_sizes);
        assert_eq!(trimmed.regions.len(), 1);
        assert_eq!(trimmed.regions[0].chr, "chr1");
    }

    #[rstest]
    fn test_trim_keeps_zero_width_after_clamp() {
        let rs = make_regionset(vec![("chr1", 100, 200)]);
        let chrom_sizes: HashMap<String, u32> =
            [("chr1".to_string(), 100)].into_iter().collect();
        let trimmed = rs.trim(&chrom_sizes);
        assert_eq!(trimmed.regions.len(), 1);
        assert_eq!(trimmed.regions[0].start, 100);
        assert_eq!(trimmed.regions[0].end, 100);
    }

    #[rstest]
    fn test_trim_drops_inverted_after_clamp() {
        let rs = make_regionset(vec![("chr1", 150, 100)]);
        let chrom_sizes: HashMap<String, u32> =
            [("chr1".to_string(), 200)].into_iter().collect();
        let trimmed = rs.trim(&chrom_sizes);
        assert_eq!(trimmed.regions.len(), 0);
    }

    #[rstest]
    fn test_trim_no_change_when_within_bounds() {
        let rs = make_regionset(vec![("chr1", 10, 50)]);
        let chrom_sizes: HashMap<String, u32> =
            [("chr1".to_string(), 100)].into_iter().collect();
        let trimmed = rs.trim(&chrom_sizes);
        assert_eq!(trimmed.regions.len(), 1);
        assert_eq!(trimmed.regions[0].start, 10);
        assert_eq!(trimmed.regions[0].end, 50);
    }

    // ── promoters tests ─────────────────────────────────────────────

    #[rstest]
    fn test_promoters_standard() {
        let rs = make_regionset(vec![("chr1", 1000, 2000)]);
        let result = rs.promoters(500, 200);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0].start, 500);
        assert_eq!(result.regions[0].end, 1200);
    }

    #[rstest]
    fn test_promoters_saturating_sub_at_zero() {
        let rs = make_regionset(vec![("chr1", 100, 500)]);
        let result = rs.promoters(200, 50);
        assert_eq!(result.regions[0].start, 0);
        assert_eq!(result.regions[0].end, 150);
    }

    #[rstest]
    fn test_promoters_at_origin() {
        let rs = make_regionset(vec![("chr1", 0, 100)]);
        let result = rs.promoters(500, 200);
        assert_eq!(result.regions[0].start, 0);
        assert_eq!(result.regions[0].end, 200);
    }

    // ── reduce tests ────────────────────────────────────────────────

    #[rstest]
    fn test_reduce_dummy_bed() {
        let path = get_test_path("dummy.bed");
        let rs = RegionSet::try_from(path.to_str().unwrap()).unwrap();
        let reduced = rs.reduce();
        assert_eq!(reduced.regions.len(), 1);
        assert_eq!(reduced.regions[0].chr, "chr1");
        assert_eq!(reduced.regions[0].start, 2);
        assert_eq!(reduced.regions[0].end, 12);
    }

    #[rstest]
    fn test_reduce_non_overlapping() {
        let rs = make_regionset(vec![("chr1", 0, 5), ("chr1", 10, 15), ("chr1", 20, 25)]);
        let reduced = rs.reduce();
        assert_eq!(reduced.regions.len(), 3);
    }

    #[rstest]
    fn test_reduce_adjacent_merged() {
        let rs = make_regionset(vec![("chr1", 0, 10), ("chr1", 10, 20)]);
        let reduced = rs.reduce();
        assert_eq!(reduced.regions.len(), 1);
        assert_eq!(reduced.regions[0].start, 0);
        assert_eq!(reduced.regions[0].end, 20);
    }

    #[rstest]
    fn test_reduce_multi_chrom() {
        let rs = make_regionset(vec![
            ("chr1", 0, 10),
            ("chr1", 5, 15),
            ("chr2", 0, 10),
            ("chr2", 20, 30),
        ]);
        let reduced = rs.reduce();
        assert_eq!(reduced.regions.len(), 3);
        assert_eq!(reduced.regions[0], make_region("chr1", 0, 15));
        assert_eq!(reduced.regions[1], make_region("chr2", 0, 10));
        assert_eq!(reduced.regions[2], make_region("chr2", 20, 30));
    }

    #[rstest]
    fn test_reduce_lexicographic_chrom_order() {
        let rs = make_regionset(vec![
            ("chr10", 0, 10),
            ("chr2", 0, 10),
            ("chr1", 0, 10),
            ("chrX", 0, 10),
            ("chrM", 0, 10),
            ("chrY", 0, 10),
        ]);
        let reduced = rs.reduce();
        assert_eq!(reduced.regions.len(), 6);
        assert_eq!(reduced.regions[0].chr, "chr1");
        assert_eq!(reduced.regions[1].chr, "chr10");
        assert_eq!(reduced.regions[2].chr, "chr2");
        assert_eq!(reduced.regions[3].chr, "chrM");
        assert_eq!(reduced.regions[4].chr, "chrX");
        assert_eq!(reduced.regions[5].chr, "chrY");
    }

    #[rstest]
    fn test_reduce_empty() {
        let rs = RegionSet::from(Vec::<Region>::new());
        let reduced = rs.reduce();
        assert_eq!(reduced.regions.len(), 0);
    }

    // ── concat tests ────────────────────────────────────────────────

    #[rstest]
    fn test_concat() {
        let a = make_regionset(vec![("chr1", 0, 10), ("chr1", 20, 30)]);
        let b = make_regionset(vec![("chr1", 5, 15)]);
        let result = a.concat(&b);
        assert_eq!(result.regions.len(), 3);
    }

    #[rstest]
    fn test_concat_preserves_regions_and_order() {
        let a = make_regionset(vec![("chr1", 0, 10), ("chr1", 20, 30)]);
        let b = make_regionset(vec![("chr1", 5, 15)]);
        let result = a.concat(&b);
        assert_eq!(result.regions[0], make_region("chr1", 0, 10));
        assert_eq!(result.regions[1], make_region("chr1", 20, 30));
        assert_eq!(result.regions[2], make_region("chr1", 5, 15));
    }

    #[rstest]
    fn test_concat_both_empty() {
        let a = RegionSet::from(Vec::<Region>::new());
        let b = RegionSet::from(Vec::<Region>::new());
        let result = a.concat(&b);
        assert_eq!(result.regions.len(), 0);
    }

    // ── pintersect tests ────────────────────────────────────────────

    #[rstest]
    fn test_pintersect_overlapping_pair() {
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr1", 5, 15)]);
        let result = a.pintersect(&b);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0], make_region("chr1", 5, 10));
    }

    #[rstest]
    fn test_pintersect_no_overlap_zero_width() {
        let a = make_regionset(vec![("chr1", 0, 5)]);
        let b = make_regionset(vec![("chr1", 10, 20)]);
        let result = a.pintersect(&b);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0].start, 10);
        assert_eq!(result.regions[0].end, 10);
    }

    #[rstest]
    fn test_pintersect_chrom_mismatch_zero_width() {
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr2", 0, 10)]);
        let result = a.pintersect(&b);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0].chr, "chr1");
        assert_eq!(result.regions[0].start, 0);
        assert_eq!(result.regions[0].end, 0);
    }

    // ── disjoin tests ───────────────────────────────────────────────

    #[rstest]
    fn test_disjoin_overlapping() {
        let rs = make_regionset(vec![("chr1", 0, 10), ("chr1", 5, 15)]);
        let result = rs.disjoin();
        assert_eq!(result.regions.len(), 3);
        assert_eq!(result.regions[0], make_region("chr1", 0, 5));
        assert_eq!(result.regions[1], make_region("chr1", 5, 10));
        assert_eq!(result.regions[2], make_region("chr1", 10, 15));
    }

    #[rstest]
    fn test_disjoin_empty() {
        let rs = RegionSet::from(Vec::<Region>::new());
        assert_eq!(rs.disjoin().regions.len(), 0);
    }

    // ── gaps tests ──────────────────────────────────────────────────

    #[rstest]
    fn test_gaps_basic() {
        let rs = make_regionset(vec![("chr1", 0, 10), ("chr1", 20, 30)]);
        let result = rs.gaps();
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0], make_region("chr1", 10, 20));
    }

    #[rstest]
    fn test_gaps_empty() {
        let rs = RegionSet::from(Vec::<Region>::new());
        assert_eq!(rs.gaps().regions.len(), 0);
    }

    #[rstest]
    fn test_gaps_adjacent_no_gap() {
        let rs = make_regionset(vec![("chr1", 0, 10), ("chr1", 10, 20)]);
        let result = rs.gaps();
        assert_eq!(result.regions.len(), 0);
    }

    // ── shift tests ─────────────────────────────────────────────────

    #[rstest]
    fn test_shift_positive() {
        let rs = make_regionset(vec![("chr1", 10, 20)]);
        let result = rs.shift(5);
        assert_eq!(result.regions[0], make_region("chr1", 15, 25));
    }

    #[rstest]
    fn test_shift_negative_saturates() {
        let rs = make_regionset(vec![("chr1", 3, 10)]);
        let result = rs.shift(-5);
        assert_eq!(result.regions[0].start, 0);
        assert_eq!(result.regions[0].end, 5);
    }

    // ── flank tests ─────────────────────────────────────────────────

    #[rstest]
    fn test_flank_upstream_of_start() {
        let rs = make_regionset(vec![("chr1", 100, 200)]);
        let result = rs.flank(50, true, false);
        assert_eq!(result.regions[0], make_region("chr1", 50, 100));
    }

    #[rstest]
    fn test_flank_downstream_of_end() {
        let rs = make_regionset(vec![("chr1", 100, 200)]);
        let result = rs.flank(50, false, false);
        assert_eq!(result.regions[0], make_region("chr1", 200, 250));
    }

    #[rstest]
    fn test_flank_both_around_start() {
        let rs = make_regionset(vec![("chr1", 100, 200)]);
        let result = rs.flank(50, true, true);
        assert_eq!(result.regions[0], make_region("chr1", 50, 150));
    }

    // ── resize tests ────────────────────────────────────────────────

    #[rstest]
    fn test_resize_from_start() {
        let rs = make_regionset(vec![("chr1", 100, 200)]);
        let result = rs.resize(50, "start");
        assert_eq!(result.regions[0], make_region("chr1", 100, 150));
    }

    #[rstest]
    fn test_resize_from_end() {
        let rs = make_regionset(vec![("chr1", 100, 200)]);
        let result = rs.resize(50, "end");
        assert_eq!(result.regions[0], make_region("chr1", 150, 200));
    }

    #[rstest]
    fn test_resize_from_center() {
        let rs = make_regionset(vec![("chr1", 100, 200)]);
        let result = rs.resize(40, "center");
        assert_eq!(result.regions[0], make_region("chr1", 130, 170));
    }

    // ── narrow tests ────────────────────────────────────────────────

    #[rstest]
    fn test_narrow_start_and_end() {
        let rs = make_regionset(vec![("chr1", 100, 200)]);
        let result = rs.narrow(Some(1), Some(50), None);
        assert_eq!(result.regions[0], make_region("chr1", 100, 150));
    }

    #[rstest]
    fn test_narrow_start_and_width() {
        let rs = make_regionset(vec![("chr1", 100, 200)]);
        let result = rs.narrow(Some(1), None, Some(30));
        assert_eq!(result.regions[0], make_region("chr1", 100, 130));
    }

    // ── IntervalSetOps: setdiff tests ───────────────────────────────

    #[rstest]
    fn test_setdiff_middle_subtraction() {
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr1", 3, 7)]);
        let result = a.setdiff(&b);
        assert_eq!(result.regions.len(), 2);
        assert_eq!(result.regions[0], make_region("chr1", 0, 3));
        assert_eq!(result.regions[1], make_region("chr1", 7, 10));
    }

    #[rstest]
    fn test_setdiff_complete_subtraction() {
        let a = make_regionset(vec![("chr1", 3, 7)]);
        let b = make_regionset(vec![("chr1", 0, 10)]);
        let result = a.setdiff(&b);
        assert_eq!(result.regions.len(), 0);
    }

    #[rstest]
    fn test_setdiff_no_overlap() {
        let a = make_regionset(vec![("chr1", 0, 5)]);
        let b = make_regionset(vec![("chr1", 10, 20)]);
        let result = a.setdiff(&b);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0], make_region("chr1", 0, 5));
    }

    #[rstest]
    fn test_setdiff_multi_chrom() {
        let a = make_regionset(vec![("chr1", 0, 10), ("chr2", 0, 10)]);
        let b = make_regionset(vec![("chr1", 5, 15)]);
        let result = a.setdiff(&b);
        assert_eq!(result.regions.len(), 2);
        assert_eq!(result.regions[0], make_region("chr1", 0, 5));
        assert_eq!(result.regions[1], make_region("chr2", 0, 10));
    }

    #[rstest]
    fn test_setdiff_from_bed_files() {
        let path_a = get_test_path("dummy.bed");
        let path_b = get_test_path("dummy_b.bed");
        let a = RegionSet::try_from(path_a.to_str().unwrap()).unwrap();
        let b = RegionSet::try_from(path_b.to_str().unwrap()).unwrap();
        let result = a.setdiff(&b);
        assert_eq!(result.regions.len(), 3);
        assert_eq!(result.regions[0], make_region("chr1", 2, 3));
        assert_eq!(result.regions[1], make_region("chr1", 5, 8));
        assert_eq!(result.regions[2], make_region("chr1", 10, 12));
    }

    // ── IntervalSetOps: intersect tests ─────────────────────────────

    #[rstest]
    fn test_intersect_overlapping() {
        let a = make_regionset(vec![("chr1", 100, 200)]);
        let b = make_regionset(vec![("chr1", 150, 250)]);
        let result = a.intersect(&b);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0], make_region("chr1", 150, 200));
    }

    #[rstest]
    fn test_intersect_no_overlap() {
        let a = make_regionset(vec![("chr1", 100, 200)]);
        let b = make_regionset(vec![("chr1", 300, 400)]);
        let result = a.intersect(&b);
        assert_eq!(result.regions.len(), 0);
    }

    #[rstest]
    fn test_intersect_empty_inputs() {
        let empty = RegionSet::from(Vec::<Region>::new());
        let a = make_regionset(vec![("chr1", 100, 200)]);
        assert_eq!(a.intersect(&empty).regions.len(), 0);
        assert_eq!(empty.intersect(&a).regions.len(), 0);
    }

    // ── IntervalSetOps: union tests ─────────────────────────────────

    #[rstest]
    fn test_union_overlapping() {
        let a = make_regionset(vec![("chr1", 0, 15)]);
        let b = make_regionset(vec![("chr1", 10, 25)]);
        let result = a.union(&b);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0], make_region("chr1", 0, 25));
    }

    #[rstest]
    fn test_union_disjoint() {
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr1", 20, 30)]);
        let result = a.union(&b);
        assert_eq!(result.regions.len(), 2);
    }

    #[rstest]
    fn test_union_both_empty() {
        let a = RegionSet::from(Vec::<Region>::new());
        let b = RegionSet::from(Vec::<Region>::new());
        let result = a.union(&b);
        assert_eq!(result.regions.len(), 0);
    }

    // ── IntervalSetOps: jaccard tests ───────────────────────────────

    #[rstest]
    fn test_jaccard() {
        let path_a = get_test_path("dummy.bed");
        let path_b = get_test_path("dummy_b.bed");
        let a = RegionSet::try_from(path_a.to_str().unwrap()).unwrap();
        let b = RegionSet::try_from(path_b.to_str().unwrap()).unwrap();
        let j = a.jaccard(&b);
        assert!((j - 0.4).abs() < 1e-10);
    }

    #[rstest]
    fn test_jaccard_identical() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        assert!((a.jaccard(&a) - 1.0).abs() < 1e-10);
    }

    #[rstest]
    fn test_jaccard_disjoint() {
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr1", 20, 30)]);
        assert!((a.jaccard(&b)).abs() < 1e-10);
    }

    #[rstest]
    fn test_jaccard_empty() {
        let a = RegionSet::from(Vec::<Region>::new());
        let b = RegionSet::from(Vec::<Region>::new());
        assert!((a.jaccard(&b)).abs() < 1e-10);
    }

    // ── IntervalSetOps: coverage tests ──────────────────────────────

    #[rstest]
    fn test_coverage_identical() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        assert!((a.coverage(&a) - 1.0).abs() < 1e-10);
    }

    #[rstest]
    fn test_coverage_disjoint() {
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr1", 20, 30)]);
        assert!(a.coverage(&b).abs() < 1e-10);
    }

    #[rstest]
    fn test_coverage_partial_overlap() {
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr1", 5, 15)]);
        assert!((a.coverage(&b) - 0.5).abs() < 1e-10);
    }

    // ── IntervalSetOps: overlap_coefficient tests ───────────────────

    #[rstest]
    fn test_overlap_coefficient_identical() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        assert!((a.overlap_coefficient(&a) - 1.0).abs() < 1e-10);
    }

    #[rstest]
    fn test_overlap_coefficient_disjoint() {
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr1", 20, 30)]);
        assert!(a.overlap_coefficient(&b).abs() < 1e-10);
    }

    // ── IntervalSetOps: subtract tests ──────────────────────────────

    #[rstest]
    fn test_subtract_same_as_setdiff() {
        let a = make_regionset(vec![("chr1", 0, 20)]);
        let b = make_regionset(vec![("chr1", 5, 15)]);
        let setdiff_result = a.setdiff(&b);
        // subtract is now an inherent method (alias for setdiff)
        let subtract_result = a.subtract(&b);
        assert_eq!(setdiff_result.regions, subtract_result.regions);
    }

    // ── IntervalSetOps: closest tests ───────────────────────────────

    #[rstest]
    fn test_closest_overlapping() {
        let a = make_regionset(vec![("chr1", 100, 200)]);
        let b = make_regionset(vec![("chr1", 150, 250)]);
        let result = a.closest(&b);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], (0, 0, 0));
    }

    #[rstest]
    fn test_closest_upstream() {
        let a = make_regionset(vec![("chr1", 200, 300)]);
        let b = make_regionset(vec![("chr1", 100, 150)]);
        let result = a.closest(&b);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], (0, 0, 50));
    }

    #[rstest]
    fn test_closest_absent_chrom_omitted() {
        let a = make_regionset(vec![("chr2", 100, 200)]);
        let b = make_regionset(vec![("chr1", 100, 200)]);
        let result = a.closest(&b);
        assert_eq!(result.len(), 0);
    }

    /// Regression: a long interval far to the left overlaps the query but the
    /// old fixed-window algorithm never examined it.
    #[rstest]
    fn test_closest_long_interval_overlaps_query() {
        let other = make_regionset(vec![
            ("chr1", 0, 10000),
            ("chr1", 100, 101),
            ("chr1", 200, 201),
            ("chr1", 300, 301),
            ("chr1", 400, 401),
        ]);
        let query = make_regionset(vec![("chr1", 500, 501)]);
        let result = query.closest(&other);
        assert_eq!(result.len(), 1);
        // The long interval (0, 10000) overlaps (500, 501), so distance is 0
        assert_eq!(result[0].2, 0, "expected overlap (distance 0)");
        assert_eq!(result[0].1, 0, "expected match to index 0 in other");
    }

    /// Regression: a long interval far to the left is the nearest but does
    /// not overlap the query. The old algorithm would miss it.
    #[rstest]
    fn test_closest_long_interval_nearest_no_overlap() {
        let other = make_regionset(vec![
            ("chr1", 0, 490),
            ("chr1", 100, 101),
            ("chr1", 200, 201),
            ("chr1", 300, 301),
            ("chr1", 400, 401),
        ]);
        let query = make_regionset(vec![("chr1", 500, 501)]);
        let result = query.closest(&other);
        assert_eq!(result.len(), 1);
        // (0, 490) has distance 10 from (500, 501); (400, 401) has distance 99
        assert_eq!(result[0].2, 10, "expected distance 10 to (0,490)");
        assert_eq!(result[0].1, 0, "expected match to index 0 in other");
    }

    /// Stress test: many small intervals with one long interval far left
    /// that is the true nearest neighbor.
    #[rstest]
    fn test_closest_many_small_intervals_with_long_left() {
        let mut regions: Vec<(&str, u32, u32)> = Vec::new();
        // Long interval at the start that reaches close to query
        regions.push(("chr1", 0, 995));
        // 100 small intervals between positions 100-600
        for i in 0..100 {
            regions.push(("chr1", 100 + i * 5, 101 + i * 5));
        }
        let other = make_regionset(regions);
        let query = make_regionset(vec![("chr1", 1000, 1001)]);
        let result = query.closest(&other);
        assert_eq!(result.len(), 1);
        // (0, 995) has distance 5; the closest small interval is (595, 596) with distance 404
        assert_eq!(result[0].2, 5, "expected distance 5 to (0,995)");
        assert_eq!(result[0].1, 0, "expected match to index 0 in other");
    }

    /// Symmetry: a.closest(b) and b.closest(a) should give matching distances.
    #[rstest]
    fn test_closest_symmetric_distances() {
        let a = make_regionset(vec![
            ("chr1", 100, 200),
            ("chr1", 500, 600),
        ]);
        let b = make_regionset(vec![
            ("chr1", 0, 5000),
            ("chr1", 250, 260),
            ("chr1", 700, 800),
        ]);
        let ab = a.closest(&b);
        let ba = b.closest(&a);
        // For each pair found in ab, the distance should appear in ba as well
        for &(ai, bi, dist_ab) in &ab {
            // Find the ba entry for bi
            if let Some(&(_, matched_a, dist_ba)) = ba.iter().find(|&&(b_idx, _, _)| b_idx == bi) {
                // If b[bi]'s closest in a is a[ai], distances must match
                if matched_a == ai {
                    assert_eq!(
                        dist_ab.abs(),
                        dist_ba.abs(),
                        "symmetric distance mismatch for pair ({}, {})",
                        ai,
                        bi
                    );
                }
            }
        }
    }

    // ── IntervalSetOps: cluster tests ───────────────────────────────

    #[rstest]
    fn test_cluster_non_overlapping_separate() {
        let a = make_regionset(vec![("chr1", 0, 10), ("chr1", 20, 30), ("chr1", 40, 50)]);
        let ids = a.cluster(0);
        assert_eq!(ids.len(), 3);
        assert_ne!(ids[0], ids[1]);
        assert_ne!(ids[1], ids[2]);
    }

    #[rstest]
    fn test_cluster_overlapping_same() {
        let a = make_regionset(vec![("chr1", 0, 15), ("chr1", 10, 25)]);
        let ids = a.cluster(0);
        assert_eq!(ids[0], ids[1]);
    }

    #[rstest]
    fn test_cluster_empty() {
        let a = RegionSet::from(Vec::<Region>::new());
        let ids = a.cluster(0);
        assert!(ids.is_empty());
    }

    #[rstest]
    fn test_cluster_original_order() {
        let a = make_regionset(vec![
            ("chr1", 100, 110),
            ("chr1", 0, 10),
        ]);
        let ids = a.cluster(0);
        assert_eq!(ids.len(), 2);
        assert_eq!(ids[0], 1);
        assert_eq!(ids[1], 0);
    }
}

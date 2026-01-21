use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufRead, BufReader};
#[cfg(feature = "http")]
use std::io::Cursor;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use anyhow::{Context, Result};
#[cfg(feature = "http")]
use flate2::read::GzDecoder;
use flate2::read::MultiGzDecoder;
use std::error::Error;
#[cfg(feature = "http")]
use ureq::{get, Error as UreqError};

use crate::models::region::Region;

#[derive(Debug, Clone)]
#[allow(clippy::upper_case_acronyms)]
pub enum FileType {
    BED,
    BAM,
    NARROWPEAK,
    UNKNOWN, // Add an UNKNOWN variant for unhandled types
}

impl FromStr for FileType {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "bed" => Ok(FileType::BED),
            "bam" => Ok(FileType::BAM),
            "narrowpeak" => Ok(FileType::NARROWPEAK),
            _ => Ok(FileType::UNKNOWN), // Return UNKNOWN for unhandled types
                                        //_ => Err(format!("Invalid file type: {}", s)),
        }
    }
}

pub struct FileInfo {
    pub file_type: FileType,
    pub is_gzipped: bool,
}

pub fn get_file_info(path: &Path) -> FileInfo {
    let mut file_type = FileType::UNKNOWN;
    let mut is_gzipped = false;

    if let Some(os_str_filename) = path.file_name() {
        if let Some(filename) = os_str_filename.to_str() {
            // Check for .gz first
            if filename.ends_with(".gz") {
                is_gzipped = true;
                if let Some(base_filename) = filename.strip_suffix(".gz") {
                    // Try to get the extension before .gz
                    if let Some(ext) = PathBuf::from(base_filename)
                        .extension()
                        .and_then(|e| e.to_str())
                    {
                        file_type = FileType::from_str(ext).unwrap_or(FileType::UNKNOWN);
                    } else {
                        // If there's no extension before .gz (e.g., "my_data.gz"),
                        // you might want to handle this specifically or leave as UNKNOWN.
                        // For now, we'll try to parse the whole base_filename as a type
                        file_type = FileType::from_str(base_filename).unwrap_or(FileType::UNKNOWN);
                    }
                }
            } else {
                // Not gzipped, just get the direct extension
                if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
                    file_type = FileType::from_str(ext).unwrap_or(FileType::UNKNOWN);
                }
            }
        }
    }

    FileInfo {
        file_type,
        is_gzipped,
    }
}

/// Parses each line of given bed like file into a contig (chromosome), starts and ends
/// This ignores any other columns beyond start and ends.
pub fn parse_bedlike_file(line: &str) -> Option<(String, i32, i32)> {
    let mut fields = line.split('\t');
    // Get the first field which should be chromosome.
    let ctg = fields.next()?;
    // Parse 2nd and 3rd string as integers or return -1 if failure
    let st = fields
        .next()
        .and_then(|s| s.parse::<i32>().ok())
        .unwrap_or(-1);
    let en = fields
        .next()
        .and_then(|s| s.parse::<i32>().ok())
        .unwrap_or(-1);

    // Original code had a remainder of the line, r, but it does not appear to have been used
    // in any way

    Some((ctg.parse().unwrap(), st, en))
}

///
/// Get a reader for either a gzip'd or non-gzip'd file.
///
/// # Arguments
///
/// - path: path to the file to read
///
pub fn get_dynamic_reader(path: &Path) -> Result<BufReader<Box<dyn Read>>> {
    let is_gzipped = path.extension() == Some(OsStr::new("gz"));
    let file = File::open(path).with_context(|| format!("Failed to open file: {:?}", path))?;
    let file: Box<dyn Read> = match is_gzipped {
        true => Box::new(MultiGzDecoder::new(file)),
        false => Box::new(file),
    };

    let reader = BufReader::new(file);

    Ok(reader)
}

///
/// Get a reader for url ling. Either for gzip'd or non-gzip'd file
///
/// # Arguments
///
/// - path: path to the file to read
///
#[cfg(feature = "http")]
pub fn get_dynamic_reader_from_url(
    url: &Path,
) -> Result<BufReader<Box<dyn std::io::Read>>, Box<dyn Error>> {
    let mut url_str = url
        .to_str()
        .ok_or_else(|| "URL path is not valid UTF-8")?
        .to_string();

    let is_ftp = url_str.starts_with("ftp://");
    if is_ftp {
        println!("ftp is not fully implemented. Bugs could appear");
        url_str = url_str.replacen("ftp://", "http://", 1);
    }

    // Perform request
    let response = match get(&url_str).call() {
        Ok(resp) => resp,
        Err(UreqError::StatusCode(code)) => {
            return Err(format!("HTTP status {} when fetching {}", code, url_str).into())
        }
        Err(e) => return Err(format!("Request error when fetching {}: {}", url_str, e).into()),
    };

    // Read the entire HTTP response body into memory as a Vec<u8>
    let mut bytes = Vec::new();
    response
        .into_body()
        .into_reader()
        .read_to_end(&mut bytes)
        .map_err(|e| format!("Failed reading response body from {}: {}", url_str, e))?;

    let cursor = Cursor::new(bytes);

    let is_gzipped = url_str.ends_with(".gz");

    let reader: Box<dyn std::io::Read> = match is_gzipped {
        true => Box::new(GzDecoder::new(cursor)),
        false => Box::new(cursor),
    };

    Ok(BufReader::new(reader))
}

/// Get a reader for either a gzipped, non-gzipped file, or stdin
///
/// # Arguments
///
/// - file_path: path to the file to read, or '-' for stdin
///
/// # Returns
///
/// A `BufReader` object for a given file path or stdin.
pub fn get_dynamic_reader_w_stdin(file_path_str: &str) -> Result<BufReader<Box<dyn Read>>> {
    if file_path_str == "-" {
        Ok(BufReader::new(Box::new(std::io::stdin()) as Box<dyn Read>))
    } else {
        let file_path = Path::new(file_path_str);
        get_dynamic_reader(file_path)
    }
}

///
/// Create a region-to-id hash-map from a list of regions
///
/// # Arguments:
/// - regions: vec![] of [Region] structs
pub fn generate_region_to_id_map(regions: &[Region]) -> HashMap<Region, u32> {
    let mut current_id = 0;
    let mut region_to_id: HashMap<Region, u32> = HashMap::new();
    for region in regions.iter() {
        region_to_id.entry(region.to_owned()).or_insert_with(|| {
            let old_id = current_id;
            current_id += 1;
            old_id
        });
    }

    region_to_id
}

///
/// Generate an id-to-region hash-map from a list of regions
///
/// # Arguments:
/// - regions: vec![] of [Region] structs
pub fn generate_id_to_region_map(regions: &[Region]) -> HashMap<u32, Region> {
    let mut current_id = 0;
    let mut id_to_region: HashMap<u32, Region> = HashMap::new();

    for region in regions.iter() {
        id_to_region.entry(current_id).or_insert_with(|| {
            current_id += 1;
            region.clone()
        });
    }

    id_to_region
}

///
/// Create a region-to-id hash-map from a list of region strings
///
/// # Arguments:
/// - regions: vec![] of region strings in the form `chr:start-end`
pub fn generate_region_string_to_id_map(regions: &[String]) -> HashMap<String, u32> {
    let mut current_id = 0;
    let mut region_to_id: HashMap<String, u32> = HashMap::new();
    for region in regions.iter() {
        region_to_id.entry(region.to_owned()).or_insert_with(|| {
            let old_id = current_id;
            current_id += 1;
            old_id
        });
    }

    region_to_id
}

///
/// Generate an id-to-region string hash-map from a list of region strings
///
/// # Arguments:
/// - regions: vec![] of region strings in the form `chr:start-end`
pub fn generate_id_to_region_string_map(regions: &[String]) -> HashMap<u32, String> {
    let mut current_id = 0;
    let mut id_to_region: HashMap<u32, String> = HashMap::new();

    for region in regions.iter() {
        id_to_region.entry(current_id).or_insert_with(|| {
            current_id += 1;
            region.clone()
        });
    }

    id_to_region
}

pub fn get_chrom_sizes<T: AsRef<Path>>(path: T) -> HashMap<String, u32> {
    let chrom_sizes_file = File::open(path.as_ref())
        .with_context(|| "Failed to open chrom sizes file.")
        .unwrap();

    let mut chrom_sizes: HashMap<String, u32> = HashMap::new();

    let file_buf = BufReader::new(chrom_sizes_file);

    for line in file_buf.lines() {
        let line_string: String = match line {
            Ok(value) => value,
            Err(_) => panic!("Error while reading chrom sizes file"),
        };

        let line_parts: Vec<String> = line_string
            .split_whitespace()
            .map(|s| s.to_string())
            .collect();

        chrom_sizes.insert(line_parts[0].clone(), line_parts[1].parse::<u32>().unwrap());
    }

    chrom_sizes
}

///
/// Gen
pub fn generate_ordering_map_for_universe_regions<T: AsRef<Path>>(
    path: T,
) -> Result<HashMap<Region, f64>> {
    let mut map = HashMap::new();

    let reader = get_dynamic_reader(path.as_ref())?;

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();

        if parts.len() < 5 {
            anyhow::bail!("BED file line does not have at least 5 fields: {}. It needs to have chr, start, end, name, and score.", line);
        }

        // parse the fields
        let chr = parts[0];
        let start = parts[1].parse::<u32>().with_context(|| {
            format!("Failed to parse start position in BED file line: {}", line)
        })?;

        let end = parts[2]
            .parse::<u32>()
            .with_context(|| format!("Failed to parse end position in BED file line: {}", line))?;

        let score = parts[4]
            .parse::<f64>()
            .with_context(|| format!("Failed to parse score in BED file line: {}", line))?;

        let rest = Some(parts[3..].join("\t")).filter(|s| !s.is_empty());

        let region = Region {
            chr: chr.to_owned(),
            start,
            end,
            rest,
        };

        map.insert(region, score);
    }

    Ok(map)
}

pub fn read_bedset_file<P: AsRef<Path>>(file_path: P) -> Result<Vec<String>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    let bed_identifiers = reader
        .lines()
        .map(|line| line.map(|s| s.trim().to_string()))
        .collect::<Result<Vec<_>, _>>()?;

    Ok(bed_identifiers)
}

pub fn remove_all_extensions(path: &Path) -> String {
    let mut stem = path.file_stem().unwrap().to_string_lossy().to_string();

    let mut parent_path = path.with_file_name(stem.clone());
    while let Some(_extension) = parent_path.extension() {
        // Remove the extension by recreating the path without it
        parent_path = parent_path.with_extension("");
        stem = parent_path
            .file_stem()
            .unwrap()
            .to_string_lossy()
            .to_string();
    }

    stem
}

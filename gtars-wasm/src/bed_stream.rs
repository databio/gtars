//! Streaming BED parser for WASM - constant memory usage, gzip auto-detection.
//!
//! This module provides a streaming BED file parser that can process data in chunks,
//! maintaining constant memory usage regardless of file size. Ideal for WASM
//! environments where large files are fetched in chunks.
//!
//! # Streaming API
//!
//! ```javascript
//! const handle = bedAnalyzerNew();
//! try {
//!     const response = await fetch(url);
//!     const reader = response.body.getReader();
//!     while (true) {
//!         const { done, value } = await reader.read();
//!         if (done) break;
//!         bedAnalyzerUpdate(handle, value);
//!     }
//!     const result = bedAnalyzerFinish(handle);
//! } catch (err) {
//!     bedAnalyzerFree(handle);
//!     throw err;
//! }
//! ```

use std::collections::HashMap;
use std::io::Write;
use std::sync::Mutex;

use flate2::write::GzDecoder;
use serde::Serialize;
use wasm_bindgen::prelude::*;

use gtars_core::models::{Region, RegionSet};
use gtars_genomicdist::bed_classifier::classify_bed;
use gtars_genomicdist::statistics::GenomicIntervalSetStatistics;

use crate::regionset::{JsBedClassificationOutput, JsChromosomeStatistics, JsRegionDistribution};

// ============================================================================
// Global storage for streaming parser instances
// ============================================================================

/// Global storage for parser instances, protected by a mutex.
static BED_PARSER_STORAGE: Mutex<Option<BedParserStorage>> = Mutex::new(None);

struct BedParserStorage {
    parsers: HashMap<u32, BedStreamParser>,
    next_id: u32,
}

impl BedParserStorage {
    fn new() -> Self {
        Self {
            parsers: HashMap::new(),
            next_id: 1,
        }
    }

    fn insert(&mut self, parser: BedStreamParser) -> u32 {
        // Find an unused ID (handles wrap-around after ~4 billion allocations)
        let mut id = self.next_id;
        while self.parsers.contains_key(&id) || id == 0 {
            id = id.wrapping_add(1);
            if id == 0 {
                id = 1;
            }
        }
        self.next_id = id.wrapping_add(1);
        if self.next_id == 0 {
            self.next_id = 1;
        }
        self.parsers.insert(id, parser);
        id
    }

    fn get_mut(&mut self, id: u32) -> Option<&mut BedStreamParser> {
        self.parsers.get_mut(&id)
    }

    fn remove(&mut self, id: u32) -> Option<BedStreamParser> {
        self.parsers.remove(&id)
    }
}

fn with_storage<F, R>(f: F) -> R
where
    F: FnOnce(&mut BedParserStorage) -> R,
{
    let mut guard = BED_PARSER_STORAGE
        .lock()
        .expect("BED_PARSER_STORAGE mutex poisoned");
    if guard.is_none() {
        *guard = Some(BedParserStorage::new());
    }
    f(guard.as_mut().unwrap())
}

// ============================================================================
// Core streaming BED parser
// ============================================================================

/// Inner BED processor that implements Write.
/// Decompressed data flows into this via the Write trait.
struct BedProcessor {
    regions: Vec<Region>,
    line_buffer: String,
    lines_processed: usize,
    invalid_lines_skipped: usize,
    header: Option<String>,
    first_line: bool,
}

impl BedProcessor {
    fn new() -> Self {
        Self {
            regions: Vec::new(),
            line_buffer: String::with_capacity(8192),
            lines_processed: 0,
            invalid_lines_skipped: 0,
            header: None,
            first_line: true,
        }
    }

    fn process_byte(&mut self, byte: u8) {
        if byte == b'\n' || byte == b'\r' {
            if !self.line_buffer.is_empty() {
                self.process_line();
                self.line_buffer.clear();
            }
        } else {
            self.line_buffer.push(byte as char);
        }
    }

    fn process_line(&mut self) {
        let line = self.line_buffer.trim();
        if line.is_empty() {
            return;
        }

        self.lines_processed += 1;

        // Handle header lines
        if line.starts_with('#') || line.starts_with("browser") || line.starts_with("track") {
            if self.header.is_none() {
                self.header = Some(line.to_string());
            } else {
                self.header.as_mut().unwrap().push('\n');
                self.header.as_mut().unwrap().push_str(line);
            }
            self.first_line = false;
            return;
        }

        let parts: Vec<&str> = line.split('\t').collect();

        // Handle column headers like "chr start end" without #
        if self.first_line && parts.len() >= 3 {
            if parts[1].parse::<u32>().is_err() {
                // This is a header line
                if self.header.is_none() {
                    self.header = Some(line.to_string());
                }
                self.first_line = false;
                return;
            }
        }
        self.first_line = false;

        // Parse BED line
        if parts.len() < 3 {
            self.invalid_lines_skipped += 1;
            return;
        }

        let start = match parts[1].parse::<u32>() {
            Ok(v) => v,
            Err(_) => {
                self.invalid_lines_skipped += 1;
                return;
            }
        };

        let end = match parts[2].parse::<u32>() {
            Ok(v) => v,
            Err(_) => {
                self.invalid_lines_skipped += 1;
                return;
            }
        };

        let rest = if parts.len() > 3 {
            Some(parts[3..].join("\t"))
        } else {
            None
        };

        self.regions.push(Region {
            chr: parts[0].to_string(),
            start,
            end,
            rest,
        });
    }

    fn finish(mut self) -> RegionSet {
        // Process any remaining data in buffer
        if !self.line_buffer.is_empty() {
            self.process_line();
        }

        let mut rs = RegionSet {
            regions: self.regions,
            header: self.header,
            path: None,
        };
        rs.sort();
        rs
    }

    fn regions_count(&self) -> usize {
        self.regions.len()
    }

    fn lines_processed(&self) -> usize {
        self.lines_processed
    }

    fn invalid_lines_skipped(&self) -> usize {
        self.invalid_lines_skipped
    }
}

impl Write for BedProcessor {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        for &byte in buf {
            self.process_byte(byte);
        }
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// State for format detection and processing
enum ProcessorState {
    /// Haven't detected format yet
    Detecting,
    /// Plain BED (uncompressed)
    Plain(BedProcessor),
    /// Gzipped BED
    Gzipped(GzDecoder<BedProcessor>),
}

/// A streaming BED parser that processes data chunk-by-chunk.
///
/// This is designed for WASM environments where files are fetched in chunks.
/// Memory usage is constant regardless of file size (except for storing regions).
pub struct BedStreamParser {
    state: ProcessorState,
    bytes_processed: usize,
}

impl BedStreamParser {
    /// Create a new streaming BED parser.
    pub fn new() -> Self {
        Self {
            state: ProcessorState::Detecting,
            bytes_processed: 0,
        }
    }

    /// Process a chunk of BED data.
    ///
    /// Handles both plain text and gzip-compressed BED files.
    pub fn update(&mut self, chunk: &[u8]) -> Result<(), String> {
        if chunk.is_empty() {
            return Ok(());
        }

        self.bytes_processed += chunk.len();

        // Detect format on first non-empty chunk
        if matches!(self.state, ProcessorState::Detecting) {
            let is_gzipped = chunk.len() >= 2 && chunk[0] == 0x1f && chunk[1] == 0x8b;

            if is_gzipped {
                let processor = BedProcessor::new();
                let decoder = GzDecoder::new(processor);
                self.state = ProcessorState::Gzipped(decoder);
            } else {
                self.state = ProcessorState::Plain(BedProcessor::new());
            }
        }

        // Process the chunk
        match &mut self.state {
            ProcessorState::Detecting => unreachable!(),
            ProcessorState::Plain(processor) => {
                processor
                    .write_all(chunk)
                    .map_err(|e| format!("Write error: {}", e))?;
            }
            ProcessorState::Gzipped(decoder) => {
                decoder
                    .write_all(chunk)
                    .map_err(|e| format!("Gzip decompression error: {}", e))?;
            }
        }

        Ok(())
    }

    /// Finalize processing and return the RegionSet.
    pub fn finish(self) -> Result<RegionSet, String> {
        match self.state {
            ProcessorState::Detecting => {
                // No data was ever provided - return empty RegionSet
                Ok(RegionSet {
                    regions: Vec::new(),
                    header: None,
                    path: None,
                })
            }
            ProcessorState::Plain(processor) => Ok(processor.finish()),
            ProcessorState::Gzipped(decoder) => {
                let processor = decoder
                    .finish()
                    .map_err(|e| format!("Gzip finalization error: {}", e))?;
                Ok(processor.finish())
            }
        }
    }

    /// Get progress information.
    pub fn progress(&self) -> BedParserProgress {
        let (regions_found, lines_processed, invalid_lines_skipped) = match &self.state {
            ProcessorState::Detecting => (0, 0, 0),
            ProcessorState::Plain(p) => {
                (p.regions_count(), p.lines_processed(), p.invalid_lines_skipped())
            }
            ProcessorState::Gzipped(d) => {
                let p = d.get_ref();
                (p.regions_count(), p.lines_processed(), p.invalid_lines_skipped())
            }
        };

        BedParserProgress {
            bytes_processed: self.bytes_processed,
            regions_found,
            lines_processed,
            invalid_lines_skipped,
        }
    }
}

impl Default for BedStreamParser {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// WASM bindings
// ============================================================================

/// Progress information for streaming parser
#[derive(Serialize)]
pub struct BedParserProgress {
    pub bytes_processed: usize,
    pub regions_found: usize,
    pub lines_processed: usize,
    pub invalid_lines_skipped: usize,
}

/// Result from finishing a BED analysis
#[derive(Serialize)]
pub struct BedAnalysisResult {
    pub number_of_regions: i32,
    pub mean_region_width: f64,
    pub nucleotides_length: i32,
    pub identifier: String,
    pub chromosome_statistics: HashMap<String, JsChromosomeStatistics>,
    pub region_distribution: Vec<JsRegionDistribution>,
    pub classify: JsBedClassificationOutput,
}

/// Create a new streaming BED parser.
///
/// Returns a handle (u32) that must be used in subsequent calls.
#[wasm_bindgen(js_name = "bedAnalyzerNew")]
pub fn bed_analyzer_new() -> u32 {
    with_storage(|storage| storage.insert(BedStreamParser::new()))
}

/// Process a chunk of BED data.
///
/// This can be called multiple times with successive chunks.
/// Handles both plain text and gzip-compressed BED files.
#[wasm_bindgen(js_name = "bedAnalyzerUpdate")]
pub fn bed_analyzer_update(handle: u32, chunk: &[u8]) -> Result<(), JsError> {
    with_storage(|storage| {
        if let Some(parser) = storage.get_mut(handle) {
            parser.update(chunk).map_err(|e| JsError::new(&e))
        } else {
            Err(JsError::new("Invalid parser handle"))
        }
    })
}

/// Finalize the parser and return the analysis results.
///
/// This consumes the parser and frees its resources.
/// Returns a JavaScript object with all computed statistics.
#[wasm_bindgen(js_name = "bedAnalyzerFinish")]
pub fn bed_analyzer_finish(handle: u32) -> Result<JsValue, JsError> {
    let parser = with_storage(|storage| storage.remove(handle))
        .ok_or_else(|| JsError::new("Invalid parser handle"))?;

    let region_set = parser.finish().map_err(|e| JsError::new(&e))?;

    // Compute all statistics
    let chromosome_statistics = compute_chromosome_statistics(&region_set);
    let region_distribution = compute_region_distribution(&region_set, 300);
    let classify = compute_classification(&region_set);

    let result = BedAnalysisResult {
        number_of_regions: region_set.len() as i32,
        mean_region_width: region_set.mean_region_width(),
        nucleotides_length: region_set.nucleotides_length() as i32,
        identifier: region_set.identifier(),
        chromosome_statistics,
        region_distribution,
        classify,
    };

    serde_wasm_bindgen::to_value(&result)
        .map_err(|e| JsError::new(&format!("Serialization error: {}", e)))
}

/// Free a parser without getting results.
///
/// Use this if you need to cancel processing early.
#[wasm_bindgen(js_name = "bedAnalyzerFree")]
pub fn bed_analyzer_free(handle: u32) -> bool {
    with_storage(|storage| storage.remove(handle).is_some())
}

/// Get the current progress of a streaming parser.
#[wasm_bindgen(js_name = "bedAnalyzerProgress")]
pub fn bed_analyzer_progress(handle: u32) -> Result<JsValue, JsError> {
    let progress = with_storage(|storage| storage.get_mut(handle).map(|parser| parser.progress()));

    match progress {
        Some(p) => serde_wasm_bindgen::to_value(&p)
            .map_err(|e| JsError::new(&format!("Serialization error: {}", e))),
        None => Err(JsError::new("Invalid parser handle")),
    }
}

// ============================================================================
// Helper functions for computing statistics
// ============================================================================

fn compute_chromosome_statistics(rs: &RegionSet) -> HashMap<String, JsChromosomeStatistics> {
    let stats = rs.chromosome_statistics();
    let mut result: HashMap<String, JsChromosomeStatistics> = HashMap::new();

    for (key, value) in stats {
        result.insert(
            key.clone(),
            JsChromosomeStatistics {
                chromosome: value.chromosome.clone(),
                number_of_regions: value.number_of_regions,
                minimum_region_length: value.minimum_region_length,
                maximum_region_length: value.maximum_region_length,
                mean_region_length: value.mean_region_length,
                median_region_length: value.median_region_length,
                start_nucleotide_position: value.start_nucleotide_position,
                end_nucleotide_position: value.end_nucleotide_position,
            },
        );
    }

    result
}

fn compute_region_distribution(rs: &RegionSet, n_bins: u32) -> Vec<JsRegionDistribution> {
    let distribution = rs.region_distribution_with_bins(n_bins);
    let mut result: Vec<JsRegionDistribution> = Vec::new();

    for value in distribution.values() {
        result.push(JsRegionDistribution {
            chr: value.chr.clone(),
            start: value.start,
            end: value.end,
            n: value.n,
            rid: value.rid,
        });
    }

    result
}

fn compute_classification(rs: &RegionSet) -> JsBedClassificationOutput {
    match classify_bed(rs) {
        Ok(output) => JsBedClassificationOutput {
            bed_compliance: output.bed_compliance.clone(),
            data_format: format!("{:#?}", output.data_format),
            compliant_columns: output.compliant_columns,
            non_compliant_columns: output.non_compliant_columns,
        },
        Err(_) => JsBedClassificationOutput {
            bed_compliance: "unknown".to_string(),
            data_format: "unknown".to_string(),
            compliant_columns: 0,
            non_compliant_columns: 0,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_streaming_basic() {
        let mut parser = BedStreamParser::new();
        parser
            .update(b"chr1\t100\t200\nchr1\t300\t400\n")
            .expect("update");
        let rs = parser.finish().expect("finish");

        assert_eq!(rs.regions.len(), 2);
        assert_eq!(rs.regions[0].chr, "chr1");
        assert_eq!(rs.regions[0].start, 100);
        assert_eq!(rs.regions[0].end, 200);
    }

    #[test]
    fn test_streaming_chunked() {
        let mut parser = BedStreamParser::new();

        parser.update(b"chr1\t100\t200\n").expect("chunk 1");
        parser.update(b"chr1\t300\t400\nchr2").expect("chunk 2");
        parser.update(b"\t500\t600\n").expect("chunk 3");

        let rs = parser.finish().expect("finish");

        assert_eq!(rs.regions.len(), 3);
    }

    #[test]
    fn test_streaming_with_header() {
        let mut parser = BedStreamParser::new();
        parser
            .update(b"#header line\nchr1\t100\t200\n")
            .expect("update");
        let rs = parser.finish().expect("finish");

        assert_eq!(rs.regions.len(), 1);
        assert!(rs.header.is_some());
    }

    #[test]
    fn test_streaming_with_rest_columns() {
        let mut parser = BedStreamParser::new();
        parser
            .update(b"chr1\t100\t200\tname1\t500\t+\n")
            .expect("update");
        let rs = parser.finish().expect("finish");

        assert_eq!(rs.regions.len(), 1);
        assert_eq!(rs.regions[0].rest, Some("name1\t500\t+".to_string()));
    }

    #[test]
    fn test_streaming_invalid_lines() {
        let mut parser = BedStreamParser::new();
        parser
            .update(b"chr1\t100\t200\ninvalid\nchr1\t300\t400\n")
            .expect("update");
        let progress = parser.progress();
        let rs = parser.finish().expect("finish");

        assert_eq!(rs.regions.len(), 2);
        assert_eq!(progress.invalid_lines_skipped, 1);
    }

    #[test]
    fn test_streaming_gzipped() {
        use flate2::write::GzEncoder;
        use flate2::Compression;

        let bed_data = b"chr1\t100\t200\nchr1\t300\t400\n";
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(bed_data).expect("compress");
        let compressed = encoder.finish().expect("finish compression");

        let mut parser = BedStreamParser::new();
        parser.update(&compressed).expect("update");
        let rs = parser.finish().expect("finish");

        assert_eq!(rs.regions.len(), 2);
    }

    #[test]
    fn test_streaming_progress() {
        let mut parser = BedStreamParser::new();

        let progress = parser.progress();
        assert_eq!(progress.bytes_processed, 0);
        assert_eq!(progress.regions_found, 0);

        parser.update(b"chr1\t100\t200\n").expect("update");
        let progress = parser.progress();
        assert_eq!(progress.bytes_processed, 15);
        assert_eq!(progress.regions_found, 1);
    }
}

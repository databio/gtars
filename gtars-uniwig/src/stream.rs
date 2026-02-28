//! Streaming uniwig processor - constant memory usage regardless of chromosome size.
//!
//! This module provides a streaming BED processor that computes coverage counts
//! using a sliding window of O(smooth_size) memory, instead of O(chromosome_size).
//! Input is processed line-by-line via BufRead for efficient I/O; output is
//! buffered through BufWriter for batched syscalls.

use std::collections::HashMap;
use std::collections::VecDeque;
use std::io::{self, BufRead, BufReader, Read, Write};

// ──────────────────────────────────────────────
// Public types
// ──────────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CountType {
    Start,
    End,
    Core,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum OutputFormat {
    Wig,
    BedGraph,
}

#[derive(Debug, Clone, PartialEq)]
pub struct CountRecord {
    pub chrom: String,
    pub position: i32, // 1-based genomic position
    pub count: u32,
}

// ──────────────────────────────────────────────
// Internal types
// ──────────────────────────────────────────────

struct BedRecord {
    chrom: String,
    start: i32, // 0-based
    end: i32,   // 0-based exclusive
    score: i32, // default 1
}

#[derive(Clone, Debug, PartialEq)]
enum ChromState {
    AwaitingFirstRecord,
    InChromosome { name: String },
}

// ──────────────────────────────────────────────
// BED line parser
// ──────────────────────────────────────────────

fn parse_bed_line(line: &[u8]) -> Result<Option<BedRecord>, String> {
    if line.is_empty() {
        return Ok(None);
    }

    // Skip comment lines
    if line[0] == b'#' {
        return Ok(None);
    }

    let line_str = std::str::from_utf8(line).map_err(|e| format!("Invalid UTF-8: {}", e))?;
    let trimmed = line_str.trim();

    if trimmed.is_empty() {
        return Ok(None);
    }

    // Skip track/browser header lines
    if trimmed.starts_with("track") || trimmed.starts_with("browser") {
        return Ok(None);
    }

    // Split on whitespace (handles both tabs and spaces)
    let fields: Vec<&str> = trimmed.split_whitespace().collect();
    if fields.len() < 3 {
        return Err(format!(
            "BED line has fewer than 3 fields: '{}'",
            trimmed
        ));
    }

    let chrom = fields[0].to_string();
    let start: i32 = fields[1]
        .parse()
        .map_err(|e| format!("Cannot parse start '{}': {}", fields[1], e))?;
    let end: i32 = fields[2]
        .parse()
        .map_err(|e| format!("Cannot parse end '{}': {}", fields[2], e))?;

    // Score from column 4 (BED5+) or default to 1
    let score = if fields.len() >= 5 {
        fields[4]
            .parse::<i32>()
            .unwrap_or(1)
            .max(0) // Clamp negative scores to 0
    } else {
        1
    };

    Ok(Some(BedRecord {
        chrom,
        start,
        end,
        score,
    }))
}

// ──────────────────────────────────────────────
// Streaming processor
// ──────────────────────────────────────────────

/// Streaming processor for BED-like input.
/// Memory usage: O(smooth_size) for the count buffer.
/// Does NOT grow with chromosome size or file size.
pub struct UniwigStreamProcessor {
    state: ChromState,
    count_buffer: VecDeque<u32>, // sliding window
    buffer_start_pos: i32,       // genomic position of count_buffer[0]
    smooth_size: i32,
    step_size: i32,
    count_type: CountType,
    chrom_sizes: HashMap<String, u32>,
    output_records: Vec<CountRecord>,
    max_gap: i64, // gap-fill threshold: 0=sparse, N>0=fill gaps<=N, -1=fully dense
}

impl UniwigStreamProcessor {
    pub fn new(
        smooth_size: i32,
        step_size: i32,
        count_type: CountType,
        chrom_sizes: HashMap<String, u32>,
    ) -> Self {
        Self {
            state: ChromState::AwaitingFirstRecord,
            count_buffer: VecDeque::new(),
            buffer_start_pos: 0,
            smooth_size,
            step_size,
            count_type,
            chrom_sizes,
            output_records: Vec::new(),
            max_gap: 0,
        }
    }

    /// Set the maximum gap width to fill with zeros.
    /// 0 = fully sparse (default), N > 0 = fill gaps <= N, -1 = fully dense.
    pub fn set_max_gap(&mut self, max_gap: i64) {
        self.max_gap = max_gap;
    }

    /// Process a single BED line provided as a byte slice.
    /// This is the primary input method.
    pub fn process_line_bytes(&mut self, line: &[u8]) -> Result<(), String> {
        match parse_bed_line(line) {
            Ok(Some(record)) => {
                self.process_record(record);
                Ok(())
            }
            Ok(None) => Ok(()), // comment/empty/header line
            Err(e) => Err(e),
        }
    }

    /// Process a single BED record through the streaming algorithm.
    fn process_record(&mut self, record: BedRecord) {
        // Determine the count window for this record
        let (window_start, window_end) = match self.count_type {
            CountType::Start => {
                // Center on start+1 (convert 0-based BED to 1-based)
                let center = record.start + 1;
                let ws = (center - self.smooth_size).max(1);
                let we = center + self.smooth_size;
                (ws, we)
            }
            CountType::End => {
                // Center on end (already 1-based since BED end is exclusive)
                let center = record.end;
                let ws = (center - self.smooth_size).max(1);
                let we = center + self.smooth_size;
                (ws, we)
            }
            CountType::Core => {
                // Window spans [start+1, end-1] (no smoothing)
                // start+1 converts 0-based BED start to 1-based
                // end-1 because BED end is exclusive
                let ws = record.start + 1;
                let we = record.end - 1;
                if we < ws {
                    return; // degenerate interval, skip
                }
                (ws, we)
            }
        };

        // Handle chromosome transition
        match &self.state {
            ChromState::AwaitingFirstRecord => {
                self.state = ChromState::InChromosome {
                    name: record.chrom.clone(),
                };
                if self.max_gap < 0 {
                    // Fully dense: emit zeros from position 1 up to window_start
                    self.buffer_start_pos = 1;
                    self.emit_up_to(window_start);
                } else {
                    self.buffer_start_pos = window_start;
                }
            }
            ChromState::InChromosome { name } => {
                if *name != record.chrom {
                    // Chromosome transition: finalize the previous one
                    self.finalize_current_chromosome();
                    self.reset_counting_state();
                    self.state = ChromState::InChromosome {
                        name: record.chrom.clone(),
                    };
                    if self.max_gap < 0 {
                        self.buffer_start_pos = 1;
                        self.emit_up_to(window_start);
                    } else {
                        self.buffer_start_pos = window_start;
                    }
                }
            }
        }

        // Emit all finalized positions before window_start
        // (since input is sorted, no future record can affect these)
        self.emit_up_to(window_start);

        // Extend buffer to cover [window_start, window_end]
        self.ensure_buffer_covers(window_start, window_end);

        // Increment counts in the window by the record's score
        if record.score > 0 {
            let score = record.score as u32;
            for pos in window_start..=window_end {
                let idx = (pos - self.buffer_start_pos) as usize;
                if idx < self.count_buffer.len() {
                    self.count_buffer[idx] += score;
                }
            }
        }
    }

    /// Emit all positions from buffer_start_pos up to (but not including) up_to_pos.
    /// Only emits positions that have data in the buffer. Gaps between windows
    /// are NOT filled with zeros — the WIG writer handles them by emitting new
    /// fixedStep headers. This keeps memory and output size proportional to
    /// actual data, not chromosome size.
    fn emit_up_to(&mut self, up_to_pos: i32) {
        let chrom_name = match &self.state {
            ChromState::InChromosome { name } => name.clone(),
            ChromState::AwaitingFirstRecord => return,
        };

        while self.buffer_start_pos < up_to_pos && !self.count_buffer.is_empty() {
            let pos = self.buffer_start_pos;
            let count = self.count_buffer.pop_front().unwrap();

            // Only emit positions that are on step_size boundaries
            // relative to position 1
            if self.step_size <= 1 || ((pos - 1) % self.step_size == 0) {
                self.output_records.push(CountRecord {
                    chrom: chrom_name.clone(),
                    position: pos,
                    count,
                });
            }

            self.buffer_start_pos += 1;
        }

        if self.count_buffer.is_empty() && self.buffer_start_pos < up_to_pos {
            let gap_width = (up_to_pos - self.buffer_start_pos) as i64;
            let should_fill = self.max_gap < 0 || gap_width <= self.max_gap;

            if should_fill {
                // Fill gap with zeros (keeps contiguous fixedStep block)
                while self.buffer_start_pos < up_to_pos {
                    if self.step_size <= 1 || ((self.buffer_start_pos - 1) % self.step_size == 0) {
                        self.output_records.push(CountRecord {
                            chrom: chrom_name.clone(),
                            position: self.buffer_start_pos,
                            count: 0,
                        });
                    }
                    self.buffer_start_pos += 1;
                }
            } else {
                // Skip past the gap (new fixedStep header will be emitted by WigWriter)
                self.buffer_start_pos = up_to_pos;
            }
        }
    }

    /// Ensure the count buffer extends to cover positions [start, end].
    fn ensure_buffer_covers(&mut self, start: i32, end: i32) {
        // If buffer is empty, set the start position
        if self.count_buffer.is_empty() {
            self.buffer_start_pos = start;
        }

        // Extend the back of the buffer if needed
        let buffer_end_pos = self.buffer_start_pos + self.count_buffer.len() as i32 - 1;
        if end > buffer_end_pos {
            let needed = (end - buffer_end_pos) as usize;
            self.count_buffer.extend(std::iter::repeat(0).take(needed));
        }
    }

    /// Finalize the current chromosome: emit remaining buffer.
    fn finalize_current_chromosome(&mut self) {
        let chrom_name = match &self.state {
            ChromState::InChromosome { name } => name.clone(),
            ChromState::AwaitingFirstRecord => return,
        };

        // Emit everything remaining in the buffer
        while !self.count_buffer.is_empty() {
            let pos = self.buffer_start_pos;
            let count = self.count_buffer.pop_front().unwrap();

            if self.step_size <= 1 || ((pos - 1) % self.step_size == 0) {
                self.output_records.push(CountRecord {
                    chrom: chrom_name.clone(),
                    position: pos,
                    count,
                });
            }

            self.buffer_start_pos += 1;
        }

        // Fully dense mode (max_gap < 0): pad trailing zeros to chrom_size
        if self.max_gap < 0 {
            if let Some(&chrom_size) = self.chrom_sizes.get(&chrom_name) {
                let end_pos = chrom_size as i32 + 1; // 1-based, inclusive
                while self.buffer_start_pos < end_pos {
                    if self.step_size <= 1
                        || ((self.buffer_start_pos - 1) % self.step_size == 0)
                    {
                        self.output_records.push(CountRecord {
                            chrom: chrom_name.clone(),
                            position: self.buffer_start_pos,
                            count: 0,
                        });
                    }
                    self.buffer_start_pos += 1;
                }
            }
        }
    }

    /// Reset counting state for a new chromosome.
    fn reset_counting_state(&mut self) {
        self.count_buffer.clear();
        self.buffer_start_pos = 0;
    }

    /// Takes accumulated output records for periodic flushing.
    pub fn drain_output(&mut self) -> Vec<CountRecord> {
        std::mem::take(&mut self.output_records)
    }

    /// Finalizes last chromosome, returns all un-drained records. Consumes self.
    pub fn finish(mut self) -> Result<Vec<CountRecord>, String> {
        self.finalize_current_chromosome();
        Ok(self.output_records)
    }
}

// ──────────────────────────────────────────────
// Wiggle output formatting
// ──────────────────────────────────────────────

/// Stateful WIG writer that tracks position across multiple write calls.
pub struct WigWriter {
    current_chrom: Option<String>,
    last_pos: Option<i32>,
}

impl WigWriter {
    pub fn new() -> Self {
        Self {
            current_chrom: None,
            last_pos: None,
        }
    }

    /// Write CountRecords as fixedStep wiggle format, maintaining state across calls.
    pub fn write_records<W: Write>(
        &mut self,
        writer: &mut W,
        records: &[CountRecord],
    ) -> io::Result<()> {
        for record in records {
            let need_header = match &self.current_chrom {
                None => true,
                Some(chrom) => {
                    if chrom != &record.chrom {
                        true
                    } else {
                        match self.last_pos {
                            Some(lp) => record.position != lp + 1,
                            None => true,
                        }
                    }
                }
            };

            if need_header {
                writeln!(
                    writer,
                    "fixedStep chrom={} start={} step=1",
                    record.chrom, record.position
                )?;
                self.current_chrom = Some(record.chrom.clone());
            }

            writeln!(writer, "{}", record.count)?;
            self.last_pos = Some(record.position);
        }

        Ok(())
    }
}

/// Write CountRecords as fixedStep wiggle format to any Write implementor.
/// Convenience function for single-shot writing.
pub fn write_records_as_wig<W: Write>(
    writer: &mut W,
    records: &[CountRecord],
) -> io::Result<()> {
    WigWriter::new().write_records(writer, records)
}

/// Write CountRecords as bedGraph format to any Write implementor.
pub fn write_records_as_bedgraph<W: Write>(
    writer: &mut W,
    records: &[CountRecord],
) -> io::Result<()> {
    for record in records {
        // bedGraph is 0-based, half-open: chrom start end value
        writeln!(
            writer,
            "{}\t{}\t{}\t{}",
            record.chrom,
            record.position - 1,
            record.position,
            record.count
        )?;
    }
    Ok(())
}

// ──────────────────────────────────────────────
// Chrom sizes reader
// ──────────────────────────────────────────────

/// Parse a chrom.sizes file (tab-separated: name\tsize) from any reader.
pub fn read_chrom_sizes<R: BufRead>(reader: R) -> Result<HashMap<String, u32>, String> {
    let mut sizes = HashMap::new();
    for line in reader.lines() {
        let line = line.map_err(|e| format!("IO error reading chrom.sizes: {}", e))?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 2 {
            return Err(format!(
                "chrom.sizes line has fewer than 2 fields: '{}'",
                trimmed
            ));
        }
        let name = fields[0].to_string();
        let size: u32 = fields[1]
            .parse()
            .map_err(|e| format!("Cannot parse size '{}': {}", fields[1], e))?;
        sizes.insert(name, size);
    }
    Ok(sizes)
}

// ──────────────────────────────────────────────
// Helpers for the streaming pipeline
// ──────────────────────────────────────────────

/// Dispatch formatted output based on format type.
fn write_output<W: Write>(
    wig_writer: &mut WigWriter,
    output: &mut W,
    records: &[CountRecord],
    format: OutputFormat,
) -> io::Result<()> {
    match format {
        OutputFormat::Wig => wig_writer.write_records(output, records),
        OutputFormat::BedGraph => write_records_as_bedgraph(output, records),
    }
}

/// Iterate lines from a BufRead source and feed them to the processor,
/// periodically draining output to keep memory bounded.
fn process_lines<R: BufRead, W: Write>(
    mut reader: R,
    processor: &mut UniwigStreamProcessor,
    wig_writer: &mut WigWriter,
    output: &mut W,
    output_format: OutputFormat,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut line_buf = String::with_capacity(256);
    loop {
        line_buf.clear();
        let bytes_read = reader.read_line(&mut line_buf)?;
        if bytes_read == 0 {
            break;
        }
        let trimmed = line_buf.trim_end_matches(|c| c == '\n' || c == '\r');
        processor.process_line_bytes(trimmed.as_bytes())?;

        let records = processor.drain_output();
        if !records.is_empty() {
            write_output(wig_writer, output, &records, output_format)?;
        }
    }
    Ok(())
}

// ──────────────────────────────────────────────
// Top-level entry point
// ──────────────────────────────────────────────

/// Run the streaming uniwig pipeline.
/// Reads BED data from `input`, writes formatted output to `output`.
/// Automatically detects gzip-compressed input via magic number peek.
pub fn uniwig_streaming<R: Read, W: Write>(
    input: R,
    output: &mut W,
    chrom_sizes: HashMap<String, u32>,
    smooth_size: i32,
    step_size: i32,
    count_type: CountType,
    output_format: OutputFormat,
    max_gap: i64,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut processor = UniwigStreamProcessor::new(
        smooth_size, step_size, count_type, chrom_sizes,
    );
    processor.set_max_gap(max_gap);

    // Wrap output in BufWriter for batched syscalls
    let mut buf_output = io::BufWriter::with_capacity(65536, output);
    let mut wig_writer = WigWriter::new();

    // Buffer the input to allow peeking for gzip detection
    let mut buf_input = BufReader::with_capacity(65536, input);

    // Peek at first 2 bytes to detect gzip magic number
    let is_gzipped = {
        let peek = buf_input.fill_buf()?;
        peek.len() >= 2 && peek[0] == 0x1f && peek[1] == 0x8b
    };

    if is_gzipped {
        // Wrap in gzip decoder, then BufReader for line iteration
        let gz = flate2::bufread::GzDecoder::new(buf_input);
        let gz_reader = BufReader::with_capacity(65536, gz);
        process_lines(gz_reader, &mut processor, &mut wig_writer,
                      &mut buf_output, output_format)?;
    } else {
        // Plain text: buf_input already implements BufRead
        process_lines(buf_input, &mut processor, &mut wig_writer,
                      &mut buf_output, output_format)?;
    }

    // Finalize
    let records = processor.finish()?;
    if !records.is_empty() {
        write_output(&mut wig_writer, &mut buf_output, &records, output_format)?;
    }

    buf_output.flush()?;
    Ok(())
}

// ──────────────────────────────────────────────
// Unit tests
// ──────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    /// Helper: run the processor on inline BED data, return records
    fn run_processor(
        bed_data: &[u8],
        smooth_size: i32,
        step_size: i32,
        count_type: CountType,
        chrom_sizes: HashMap<String, u32>,
    ) -> Vec<CountRecord> {
        let mut processor =
            UniwigStreamProcessor::new(smooth_size, step_size, count_type, chrom_sizes);
        for line in bed_data.split(|&b| b == b'\n') {
            if !line.is_empty() {
                processor.process_line_bytes(line).expect("process_line_bytes failed");
            }
        }
        // Collect all records: drained during processing + remaining from finish
        let mut all_records = processor.drain_output();
        let final_records = processor.finish().expect("finish failed");
        all_records.extend(final_records);
        all_records
    }

    // Test 1: Single chrom, single read, smoothsize=0
    #[test]
    fn test_single_read_smooth0() {
        let bed = b"chr1\t100\t200\n";
        let records = run_processor(bed, 0, 1, CountType::Start, HashMap::new());
        // Start count type: center on start+1 = 101, smooth=0 => only pos 101
        let non_zero: Vec<_> = records.iter().filter(|r| r.count > 0).collect();
        assert_eq!(non_zero.len(), 1);
        assert_eq!(non_zero[0].position, 101);
        assert_eq!(non_zero[0].count, 1);
    }

    // Test 2: Single chrom, single read, smoothsize=5
    #[test]
    fn test_single_read_smooth5() {
        let bed = b"chr1\t10\t200\n";
        let records = run_processor(bed, 5, 1, CountType::Start, HashMap::new());
        // Center on 11, window = [6, 16]
        let non_zero: Vec<_> = records.iter().filter(|r| r.count > 0).collect();
        assert_eq!(non_zero.len(), 11); // positions 6 through 16
        assert_eq!(non_zero[0].position, 6);
        assert_eq!(non_zero[10].position, 16);
        for r in &non_zero {
            assert_eq!(r.count, 1);
        }
    }

    // Test 3: Overlapping reads
    #[test]
    fn test_overlapping_reads() {
        let bed = b"chr1\t10\t20\nchr1\t12\t25\n";
        let records = run_processor(bed, 2, 1, CountType::Start, HashMap::new());
        // Start 10 -> center 11, window [9, 13]
        // Start 12 -> center 13, window [11, 15]
        // Overlap at positions 11, 12, 13
        let r11 = records.iter().find(|r| r.position == 11).unwrap();
        assert_eq!(r11.count, 2); // both reads contribute
        let r13 = records.iter().find(|r| r.position == 13).unwrap();
        assert_eq!(r13.count, 2);
        let r9 = records.iter().find(|r| r.position == 9).unwrap();
        assert_eq!(r9.count, 1); // only first read
        let r15 = records.iter().find(|r| r.position == 15).unwrap();
        assert_eq!(r15.count, 1); // only second read
    }

    // Test 4: Chromosome transition
    #[test]
    fn test_chromosome_transition() {
        let bed = b"chr1\t10\t20\nchr2\t10\t20\n";
        let records = run_processor(bed, 0, 1, CountType::Start, HashMap::new());
        let chr1: Vec<_> = records.iter().filter(|r| r.chrom == "chr1").collect();
        let chr2: Vec<_> = records.iter().filter(|r| r.chrom == "chr2").collect();
        assert_eq!(chr1.len(), 1);
        assert_eq!(chr2.len(), 1);
        assert_eq!(chr1[0].position, 11);
        assert_eq!(chr2[0].position, 11);
        // Verify chr1 is fully emitted before chr2
        let chr1_idx = records
            .iter()
            .position(|r| r.chrom == "chr1")
            .unwrap();
        let chr2_idx = records
            .iter()
            .position(|r| r.chrom == "chr2")
            .unwrap();
        assert!(chr1_idx < chr2_idx);
    }

    // Test 5: Step size > 1
    #[test]
    fn test_step_size() {
        let bed = b"chr1\t0\t10\n";
        let records = run_processor(bed, 3, 2, CountType::Start, HashMap::new());
        // Center on 1, window [-2, 4] clamped to [1, 4]
        // Step=2 relative to 1: positions 1, 3
        for r in &records {
            assert!((r.position - 1) % 2 == 0, "Position {} not on step boundary", r.position);
        }
    }

    // Test 6: Edge case: read near position 0
    #[test]
    fn test_clamp_to_position_1() {
        let bed = b"chr1\t0\t10\n";
        let records = run_processor(bed, 5, 1, CountType::Start, HashMap::new());
        // Center on 1, window [-4, 6] clamped to [1, 6]
        assert_eq!(records[0].position, 1);
        assert!(records.iter().all(|r| r.position >= 1));
    }

    // Test 7: Core count type
    #[test]
    fn test_core_count_type() {
        let bed = b"chr1\t10\t15\n";
        let records = run_processor(bed, 0, 1, CountType::Core, HashMap::new());
        // Core: window = [start+1, end-1] = [11, 14]
        let non_zero: Vec<_> = records.iter().filter(|r| r.count > 0).collect();
        assert_eq!(non_zero.len(), 4); // positions 11, 12, 13, 14
        assert_eq!(non_zero[0].position, 11);
        assert_eq!(non_zero[3].position, 14);
    }

    // Test 8: End count type
    #[test]
    fn test_end_count_type() {
        let bed = b"chr1\t10\t15\n";
        let records = run_processor(bed, 2, 1, CountType::End, HashMap::new());
        // End: center on end=15, window = [13, 17]
        let non_zero: Vec<_> = records.iter().filter(|r| r.count > 0).collect();
        assert_eq!(non_zero.len(), 5);
        assert_eq!(non_zero[0].position, 13);
        assert_eq!(non_zero[4].position, 17);
    }

    // Test 9: Score support
    #[test]
    fn test_score_support() {
        // BED5 format: chrom start end name score
        let bed = b"chr1\t10\t20\tpeak1\t3\n";
        let records = run_processor(bed, 0, 1, CountType::Start, HashMap::new());
        let r = records.iter().find(|r| r.count > 0).unwrap();
        assert_eq!(r.count, 3);
    }

    // Test 10: Processor records match pipeline WIG output
    #[test]
    fn test_processor_matches_pipeline_output() {
        let bed = b"chr1\t10\t20\nchr1\t15\t25\n";

        // Line-by-line via processor
        let records = run_processor(bed, 1, 1, CountType::Start, HashMap::new());

        // Full pipeline via uniwig_streaming
        let mut output = Vec::new();
        uniwig_streaming(
            Cursor::new(bed), &mut output, HashMap::new(),
            1, 1, CountType::Start, OutputFormat::Wig, 0i64,
        ).unwrap();

        // Parse WIG output and compare against processor records
        let output_str = String::from_utf8(output).unwrap();
        let mut parsed_records: Vec<(i32, u32)> = Vec::new();
        let mut current_pos: i32 = 0;
        for line in output_str.lines() {
            if line.starts_with("fixedStep") {
                // Extract start position
                let start_str = line.split("start=").nth(1).unwrap()
                    .split_whitespace().next().unwrap();
                current_pos = start_str.parse().unwrap();
            } else {
                let count: u32 = line.parse().unwrap();
                parsed_records.push((current_pos, count));
                current_pos += 1;
            }
        }

        assert_eq!(records.len(), parsed_records.len(),
            "Processor produced {} records but pipeline WIG has {}",
            records.len(), parsed_records.len());
        for (rec, (pos, count)) in records.iter().zip(parsed_records.iter()) {
            assert_eq!(rec.position, *pos, "Position mismatch");
            assert_eq!(rec.count, *count, "Count mismatch at position {}", pos);
        }
    }

    // Test 11: Empty input
    #[test]
    fn test_empty_input() {
        let records = run_processor(b"", 1, 1, CountType::Start, HashMap::new());
        assert!(records.is_empty());
    }

    // Test 12: Comment/header lines
    #[test]
    fn test_comments_and_headers() {
        let bed = b"# comment line\ntrack name=test\nbrowser position chr1:1-100\nchr1\t10\t20\n";
        let records = run_processor(bed, 0, 1, CountType::Start, HashMap::new());
        let non_zero: Vec<_> = records.iter().filter(|r| r.count > 0).collect();
        assert_eq!(non_zero.len(), 1);
        assert_eq!(non_zero[0].position, 11);
    }

    // Test 13: Gzip input via uniwig_streaming
    #[test]
    fn test_gzip_input() {
        use flate2::write::GzEncoder;
        use flate2::Compression;

        let bed = b"chr1\t10\t20\nchr1\t15\t25\n";

        // Compress it
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(bed).unwrap();
        let compressed = encoder.finish().unwrap();

        // Plain via uniwig_streaming
        let mut plain_output = Vec::new();
        uniwig_streaming(
            Cursor::new(bed), &mut plain_output, HashMap::new(),
            1, 1, CountType::Start, OutputFormat::Wig, 0i64,
        ).unwrap();

        // Gzipped via uniwig_streaming
        let mut gz_output = Vec::new();
        uniwig_streaming(
            Cursor::new(&compressed), &mut gz_output, HashMap::new(),
            1, 1, CountType::Start, OutputFormat::Wig, 0i64,
        ).unwrap();

        assert_eq!(plain_output, gz_output);
    }

    // Test 14: Reference fixture: _start.wig
    #[test]
    fn test_reference_start_wig() {
        // dummy.bed content: chr1 2-6, 4-7, 5-9, 7-12
        let bed = b"chr1\t2\t6\nchr1\t4\t7\nchr1\t5\t9\nchr1\t7\t12\n";
        let mut chrom_sizes = HashMap::new();
        chrom_sizes.insert("chr1".to_string(), 20);

        // The reference fixture was generated with smoothsize=1
        let records = run_processor(bed, 1, 1, CountType::Start, chrom_sizes);

        // No trailing zero padding — only data positions are emitted.
        // Start positions (1-based): 3, 5, 6, 8. With smooth=1:
        //   3: window [2,4], 5: [4,6], 6: [5,7], 8: [7,9]
        // Contiguous range covers positions 2-9
        let expected_counts: Vec<u32> = vec![1, 1, 2, 2, 2, 2, 1, 1];

        assert_eq!(records.len(), 8, "Expected 8 records (positions 2-9)");
        assert_eq!(records[0].position, 2);

        for (i, expected) in expected_counts.iter().enumerate() {
            assert_eq!(
                records[i].count, *expected,
                "Mismatch at position {}: expected {}, got {}",
                records[i].position, expected, records[i].count
            );
        }
    }

    // Test 15: Reference fixture: _end.wig
    #[test]
    fn test_reference_end_wig() {
        let bed = b"chr1\t2\t6\nchr1\t4\t7\nchr1\t5\t9\nchr1\t7\t12\n";
        let mut chrom_sizes = HashMap::new();
        chrom_sizes.insert("chr1".to_string(), 20);

        let records = run_processor(bed, 1, 1, CountType::End, chrom_sizes);

        // End positions (BED exclusive): 6, 7, 9, 12. With smooth=1:
        //   6: [5,7], 7: [6,8], 9: [8,10], 12: [11,13]
        // Contiguous range covers positions 5-13 (no trailing zero padding)
        let expected_counts: Vec<u32> = vec![1, 2, 2, 2, 1, 1, 1, 1, 1];

        assert_eq!(records[0].position, 5);
        assert_eq!(records.len(), 9, "Expected 9 records (positions 5-13)");

        for (i, expected) in expected_counts.iter().enumerate() {
            assert_eq!(
                records[i].count, *expected,
                "End wig mismatch at position {}: expected {}, got {}",
                records[i].position, expected, records[i].count
            );
        }
    }

    // Test 16: Reference fixture: _core.wig
    #[test]
    fn test_reference_core_wig() {
        let bed = b"chr1\t2\t6\nchr1\t4\t7\nchr1\t5\t9\nchr1\t7\t12\n";
        let mut chrom_sizes = HashMap::new();
        chrom_sizes.insert("chr1".to_string(), 20);

        let records = run_processor(bed, 0, 1, CountType::Core, chrom_sizes);

        // Core ranges: [3,5], [5,6], [6,8], [8,11] => contiguous positions 3-11
        // No trailing zero padding
        let expected_counts: Vec<u32> = vec![1, 1, 2, 2, 1, 2, 1, 1, 1];

        assert_eq!(records[0].position, 3);
        assert_eq!(records.len(), 9, "Expected 9 records (positions 3-11)");

        for (i, expected) in expected_counts.iter().enumerate() {
            assert_eq!(
                records[i].count, *expected,
                "Core wig mismatch at position {}: expected {}, got {}",
                records[i].position, expected, records[i].count
            );
        }
    }

    // Test 17: Wig output formatting
    #[test]
    fn test_wig_output_formatting() {
        let records = vec![
            CountRecord { chrom: "chr1".to_string(), position: 5, count: 3 },
            CountRecord { chrom: "chr1".to_string(), position: 6, count: 2 },
            CountRecord { chrom: "chr1".to_string(), position: 7, count: 1 },
            CountRecord { chrom: "chr2".to_string(), position: 10, count: 5 },
        ];

        let mut output = Vec::new();
        write_records_as_wig(&mut output, &records).unwrap();
        let output_str = String::from_utf8(output).unwrap();

        assert!(output_str.contains("fixedStep chrom=chr1 start=5 step=1"));
        assert!(output_str.contains("fixedStep chrom=chr2 start=10 step=1"));
        let lines: Vec<&str> = output_str.lines().collect();
        // chr1 header, 3 values, chr2 header, 1 value
        assert_eq!(lines.len(), 6);
        assert_eq!(lines[1], "3");
        assert_eq!(lines[2], "2");
        assert_eq!(lines[3], "1");
        assert_eq!(lines[5], "5");
    }

    // Test 18: Multi-chrom BED
    #[test]
    fn test_multi_chrom_bed() {
        // test_sorted_small.bed content (with mixed tab/space delimiters)
        let bed = b"chr11 \t10\t50\nchr11\t20\t76\nchr12\t769\t2395\nchr13\t771\t3000\nchr14\t800\t2900\nchr21\t1\t30\nchr21\t2\t19\nchr21\t16\t31\n";
        let records = run_processor(bed, 0, 1, CountType::Start, HashMap::new());

        let chroms: Vec<String> = records
            .iter()
            .map(|r| r.chrom.clone())
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();

        assert!(chroms.contains(&"chr11".to_string()));
        assert!(chroms.contains(&"chr12".to_string()));
        assert!(chroms.contains(&"chr13".to_string()));
        assert!(chroms.contains(&"chr14".to_string()));
        assert!(chroms.contains(&"chr21".to_string()));
        assert_eq!(chroms.len(), 5);
    }

    // Test: chrom sizes reader
    #[test]
    fn test_read_chrom_sizes() {
        let data = b"chr1\t248956422\nchr2\t242193529\n";
        let sizes = read_chrom_sizes(Cursor::new(data)).unwrap();
        assert_eq!(sizes.get("chr1"), Some(&248956422));
        assert_eq!(sizes.get("chr2"), Some(&242193529));
    }

    // Helper: run processor with explicit max_gap setting
    fn run_processor_with_gap(
        bed_data: &[u8],
        smooth_size: i32,
        step_size: i32,
        count_type: CountType,
        chrom_sizes: HashMap<String, u32>,
        max_gap: i64,
    ) -> Vec<CountRecord> {
        let mut processor =
            UniwigStreamProcessor::new(smooth_size, step_size, count_type, chrom_sizes);
        processor.set_max_gap(max_gap);
        for line in bed_data.split(|&b| b == b'\n') {
            if !line.is_empty() {
                processor.process_line_bytes(line).expect("process_line_bytes failed");
            }
        }
        let mut all_records = processor.drain_output();
        let final_records = processor.finish().expect("finish failed");
        all_records.extend(final_records);
        all_records
    }

    // Test: max_gap=0 skips all gaps (fully sparse)
    #[test]
    fn test_max_gap_zero_skips_all_gaps() {
        // Two reads with a gap of 10 positions between their windows
        // Read1: start=10 -> center=11, window=[11,11] (smooth=0)
        // Read2: start=30 -> center=31, window=[31,31] (smooth=0)
        // Gap between pos 11 and 31 is 20 positions
        let bed = b"chr1\t10\t20\nchr1\t30\t40\n";
        let records = run_processor_with_gap(bed, 0, 1, CountType::Start, HashMap::new(), 0);

        // With max_gap=0, no zero-valued positions should appear between the two clusters
        let positions: Vec<i32> = records.iter().map(|r| r.position).collect();
        assert!(positions.contains(&11), "Should have pos 11");
        assert!(positions.contains(&31), "Should have pos 31");

        // No positions between 12 and 30 should be emitted
        for pos in 12..31 {
            assert!(!positions.contains(&pos), "Should NOT have pos {} with max_gap=0", pos);
        }
    }

    // Test: max_gap fills small gaps but not large ones
    #[test]
    fn test_max_gap_fills_small_gaps() {
        // Read1: start=10 -> center=11, window=[11,11] (smooth=0)
        // Read2: start=17 -> center=18, window=[18,18] (smooth=0)
        // Gap from 12 to 17 = width 6
        let bed = b"chr1\t10\t20\nchr1\t17\t25\n";

        // With max_gap=10, gap of 6 should be filled
        let records_filled = run_processor_with_gap(bed, 0, 1, CountType::Start, HashMap::new(), 10);
        let positions_filled: Vec<i32> = records_filled.iter().map(|r| r.position).collect();
        // Positions 12..17 should be emitted as zeros
        for pos in 12..18 {
            assert!(
                positions_filled.contains(&pos),
                "With max_gap=10, pos {} should be filled with zeros",
                pos
            );
        }

        // With max_gap=3, gap of 6 should NOT be filled
        let records_skipped = run_processor_with_gap(bed, 0, 1, CountType::Start, HashMap::new(), 3);
        let positions_skipped: Vec<i32> = records_skipped.iter().map(|r| r.position).collect();
        assert!(positions_skipped.contains(&11), "Should have pos 11");
        assert!(positions_skipped.contains(&18), "Should have pos 18");
        for pos in 12..18 {
            assert!(
                !positions_skipped.contains(&pos),
                "With max_gap=3, pos {} should be skipped",
                pos
            );
        }
    }

    // Test: max_gap=-1 is fully dense (matches old dense: true behavior)
    #[test]
    fn test_max_gap_negative_one_fully_dense() {
        let bed = b"chr1\t5\t10\n";
        let mut chrom_sizes = HashMap::new();
        chrom_sizes.insert("chr1".to_string(), 20u32);

        // With max_gap=-1: output starts at position 1 and extends to chrom_size=20
        let records = run_processor_with_gap(bed, 0, 1, CountType::Start, chrom_sizes, -1);
        let positions: Vec<i32> = records.iter().map(|r| r.position).collect();

        // Should start at position 1
        assert!(positions.contains(&1), "Fully dense should start at position 1");
        // Should end at position 20
        assert!(positions.contains(&20), "Fully dense should end at chrom_size=20");
        // All positions 1..=20 should be present
        for pos in 1i32..=20 {
            assert!(positions.contains(&pos), "Fully dense missing position {}", pos);
        }
        assert_eq!(records.len(), 20, "Should have exactly 20 records");
    }

    // Test: max_gap=i64::MAX fills all inter-data gaps but NOT leading/trailing zeros to chrom_size
    #[test]
    fn test_max_gap_large_value_fills_all_data_gaps() {
        // Read1: start=10 -> center=11, window=[11,11] (smooth=0)
        // Read2: start=30 -> center=31, window=[31,31] (smooth=0)
        // With max_gap=i64::MAX all gaps between data are filled, but no leading/trailing zeros
        let bed = b"chr1\t10\t20\nchr1\t30\t40\n";
        let mut chrom_sizes = HashMap::new();
        chrom_sizes.insert("chr1".to_string(), 100u32);

        let records = run_processor_with_gap(bed, 0, 1, CountType::Start, chrom_sizes, i64::MAX);
        let positions: Vec<i32> = records.iter().map(|r| r.position).collect();

        // Should NOT have leading zeros (positions 1-10)
        for pos in 1..11 {
            assert!(!positions.contains(&pos), "Should NOT have leading zero at pos {}", pos);
        }

        // Should have the gap filled between data clusters (12..31)
        for pos in 12..31 {
            assert!(
                positions.contains(&pos),
                "With max_gap=i64::MAX, pos {} should be filled",
                pos
            );
        }

        // Should NOT have trailing zeros (positions 32..=100)
        for pos in 32..=100 {
            assert!(!positions.contains(&pos), "Should NOT have trailing zero at pos {}", pos);
        }
    }

    // Test: max_gap boundary - gap of exactly N is filled, gap of N+1 is skipped
    #[test]
    fn test_max_gap_exact_boundary() {
        // Read1: start=10 -> center=11, window=[11,11] (smooth=0)
        // Read2: start=16 -> center=17, window=[17,17] (smooth=0)
        // Gap from pos 12 to 16 = width 5
        let bed = b"chr1\t10\t20\nchr1\t16\t25\n";

        // With max_gap=5: gap of exactly 5 should be filled (condition is <=)
        let records_filled = run_processor_with_gap(bed, 0, 1, CountType::Start, HashMap::new(), 5);
        let positions_filled: Vec<i32> = records_filled.iter().map(|r| r.position).collect();
        for pos in 12..17 {
            assert!(
                positions_filled.contains(&pos),
                "With max_gap=5, pos {} (gap=5) should be filled",
                pos
            );
        }

        // Read1: start=10 -> center=11, window=[11,11] (smooth=0)
        // Read2: start=17 -> center=18, window=[18,18] (smooth=0)
        // Gap from pos 12 to 17 = width 6
        let bed2 = b"chr1\t10\t20\nchr1\t17\t25\n";

        // With max_gap=5: gap of 6 should be skipped
        let records_skipped = run_processor_with_gap(bed2, 0, 1, CountType::Start, HashMap::new(), 5);
        let positions_skipped: Vec<i32> = records_skipped.iter().map(|r| r.position).collect();
        assert!(positions_skipped.contains(&11), "Should have data pos 11");
        assert!(positions_skipped.contains(&18), "Should have data pos 18");
        for pos in 12..18 {
            assert!(
                !positions_skipped.contains(&pos),
                "With max_gap=5, pos {} (gap=6) should be skipped",
                pos
            );
        }
    }
}

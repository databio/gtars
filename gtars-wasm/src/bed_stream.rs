use std::collections::HashMap;
use std::io::{Cursor, Read};
use std::sync::Mutex;

use flate2::read::GzDecoder;
use gtars_core::models::{Region, RegionSet};
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::models::BedEntries;

// ============================================================================
// BedStreamParser — streaming BED file parser with gzip auto-detection
// ============================================================================

struct BedStreamParser {
    regions: Vec<Region>,
    compressed_buf: Vec<u8>,
    text_buf: String,
    is_gzipped: Option<bool>,
    bytes_processed: usize,
}

impl BedStreamParser {
    fn new() -> Self {
        Self {
            regions: Vec::new(),
            compressed_buf: Vec::new(),
            text_buf: String::new(),
            is_gzipped: None,
            bytes_processed: 0,
        }
    }

    fn update(&mut self, chunk: &[u8]) -> Result<(), String> {
        if chunk.is_empty() {
            return Ok(());
        }

        self.bytes_processed += chunk.len();

        if self.is_gzipped.is_none() {
            self.is_gzipped = Some(chunk.len() >= 2 && chunk[0] == 0x1f && chunk[1] == 0x8b);
        }

        if self.is_gzipped == Some(true) {
            self.compressed_buf.extend_from_slice(chunk);
        } else {
            let text = String::from_utf8_lossy(chunk);
            self.text_buf.push_str(&text);
            self.parse_complete_lines();
        }

        Ok(())
    }

    fn finish(mut self) -> Result<RegionSet, String> {
        if self.is_gzipped == Some(true) {
            let cursor = Cursor::new(&self.compressed_buf);
            let mut decoder = GzDecoder::new(cursor);
            let mut decompressed = String::new();
            decoder
                .read_to_string(&mut decompressed)
                .map_err(|e| format!("Gzip decompression failed: {}", e))?;
            self.text_buf = decompressed;
            self.parse_complete_lines();
        }

        if !self.text_buf.is_empty() {
            let remaining = std::mem::take(&mut self.text_buf);
            self.parse_line(&remaining);
        }

        let mut rs = RegionSet::from(self.regions);
        rs.sort();
        Ok(rs)
    }

    fn region_count(&self) -> usize {
        self.regions.len()
    }

    fn bytes_processed(&self) -> usize {
        self.bytes_processed
    }

    fn parse_complete_lines(&mut self) {
        let last_newline = match self.text_buf.rfind('\n') {
            Some(pos) => pos,
            None => return,
        };

        let complete = self.text_buf[..last_newline].to_string();
        let remainder = self.text_buf[last_newline + 1..].to_string();
        self.text_buf = remainder;

        for line in complete.split('\n') {
            self.parse_line(line);
        }
    }

    fn parse_line(&mut self, line: &str) {
        let trimmed = line.trim();
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed.starts_with("browser")
            || trimmed.starts_with("track")
        {
            return;
        }

        let parts: Vec<&str> = trimmed.split('\t').collect();
        if parts.len() < 3 {
            return;
        }

        let start = match parts[1].parse::<u32>() {
            Ok(v) => v,
            Err(_) => return,
        };
        let end = match parts[2].parse::<u32>() {
            Ok(v) => v,
            Err(_) => return,
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
            rest: rest.filter(|s| !s.is_empty()),
        });
    }
}

// ============================================================================
// Global storage for streaming BED parser instances
// ============================================================================

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

fn with_bed_storage<F, R>(f: F) -> R
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
// WASM bindings
// ============================================================================

#[wasm_bindgen(js_name = "bedParserNew")]
pub fn bed_parser_new() -> u32 {
    with_bed_storage(|storage| storage.insert(BedStreamParser::new()))
}

#[wasm_bindgen(js_name = "bedParserUpdate")]
pub fn bed_parser_update(handle: u32, chunk: &[u8]) -> Result<(), JsError> {
    with_bed_storage(|storage| {
        if let Some(parser) = storage.get_mut(handle) {
            parser
                .update(chunk)
                .map_err(|e| JsError::new(&format!("Failed to process chunk: {}", e)))
        } else {
            Err(JsError::new("Invalid parser handle"))
        }
    })
}

/// Finalize parsing and return BedEntries (array of [chr, start, end, rest] tuples).
/// The result can be passed directly to `new RegionSet(entries)`.
#[wasm_bindgen(js_name = "bedParserFinish")]
pub fn bed_parser_finish(handle: u32) -> Result<JsValue, JsError> {
    let parser = with_bed_storage(|storage| storage.remove(handle))
        .ok_or_else(|| JsError::new("Invalid parser handle"))?;

    let region_set = parser
        .finish()
        .map_err(|e| JsError::new(&format!("Failed to finalize parser: {}", e)))?;

    let entries: Vec<(String, u32, u32, String)> = region_set
        .regions
        .into_iter()
        .map(|r| (r.chr, r.start, r.end, r.rest.unwrap_or_default()))
        .collect();

    serde_wasm_bindgen::to_value(&BedEntries(entries))
        .map_err(|e| JsError::new(&format!("Serialization error: {}", e)))
}

#[wasm_bindgen(js_name = "bedParserFree")]
pub fn bed_parser_free(handle: u32) -> bool {
    with_bed_storage(|storage| storage.remove(handle).is_some())
}

#[wasm_bindgen(js_name = "bedParserProgress")]
pub fn bed_parser_progress(handle: u32) -> Result<JsValue, JsError> {
    let progress = with_bed_storage(|storage| {
        storage.get_mut(handle).map(|parser| BedParserProgress {
            region_count: parser.region_count(),
            bytes_processed: parser.bytes_processed(),
        })
    });

    match progress {
        Some(p) => serde_wasm_bindgen::to_value(&p)
            .map_err(|e| JsError::new(&format!("Serialization error: {}", e))),
        None => Err(JsError::new("Invalid parser handle")),
    }
}

#[derive(Serialize, Deserialize)]
pub struct BedParserProgress {
    pub region_count: usize,
    pub bytes_processed: usize,
}

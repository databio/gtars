use std::io::{BufWriter, Write, Seek};
use std::fs::File;
use std::path::{Path, PathBuf};
use std::collections::HashMap;
use anyhow::{Result, Context};
use byteorder::{LittleEndian, WriteBytesExt};
use crate::models::{ScatrsFragment, ScatrsRegion};

const MAGIC_NUMBER: &[u8; 4] = b"SCFR";
const FORMAT_VERSION: u8 = 2;  // Version 2: simplified sequential IDs
const REGION_BUFFER_SIZE: usize = 100_000; // Buffer up to 100K fragments per region before flushing

/// Index entry for a region's fragments in the file
#[derive(Debug, Clone)]
struct RegionIndex {
    region_id: u32,
    region_start: u64,
    region_end: u64,
    file_offset: u64,
    fragment_count: u32,
}

/// Writer for chromosome-based fragment cache files
/// 
/// Stores all fragments for a single chromosome in one file with an embedded index.
/// This dramatically reduces file count from hundreds of thousands to ~25 files.
pub struct ChromosomeFragmentWriter {
    writer: BufWriter<File>,
    chromosome: String,
    current_offset: u64,
    total_fragments: u64,
    region_buffers: HashMap<u32, Vec<(u32, u32)>>, // region_id -> [(start, end)]
    region_indices: Vec<RegionIndex>,
    regions: HashMap<u32, (u64, u64)>, // region_id -> (start, end) coordinates
}

impl ChromosomeFragmentWriter {
    pub fn new(path: &Path, chromosome: String) -> Result<Self> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create chromosome fragment file: {:?}", path))?;
        
        let mut writer = BufWriter::with_capacity(1024 * 1024, file); // 1MB buffer
        
        // Write header (will be updated when closing)
        writer.write_all(MAGIC_NUMBER)?;
        writer.write_u8(FORMAT_VERSION)?;
        let chr_bytes = chromosome.as_bytes();
        writer.write_u8(chr_bytes.len() as u8)?;
        writer.write_all(chr_bytes)?;
        writer.write_u64::<LittleEndian>(0)?; // placeholder for total_fragment_count
        writer.write_u64::<LittleEndian>(0)?; // placeholder for index_offset
        
        let header_size = 4 + 1 + 1 + chr_bytes.len() + 8 + 8;
        
        Ok(Self {
            writer,
            chromosome,
            current_offset: header_size as u64,
            total_fragments: 0,
            region_buffers: HashMap::new(),
            region_indices: Vec::new(),
            regions: HashMap::new(),
        })
    }
    
    /// Add a fragment for a specific region
    pub fn add_fragment(&mut self, region_id: u32, region: &ScatrsRegion, fragment: &ScatrsFragment) -> Result<()> {
        // Verify fragment is on the correct chromosome
        if fragment.chrom != self.chromosome {
            return Err(anyhow::anyhow!(
                "Fragment chromosome {} doesn't match writer chromosome {}", 
                fragment.chrom, self.chromosome
            ));
        }
        
        // Store region coordinates if not seen before
        self.regions.entry(region_id)
            .or_insert((region.start, region.end));
        
        // Buffer the fragment
        self.region_buffers
            .entry(region_id)
            .or_insert_with(Vec::new)
            .push((fragment.start as u32, fragment.end as u32));
        
        self.total_fragments += 1;
        
        // Flush if buffer is getting large
        if let Some(buffer) = self.region_buffers.get(&region_id) {
            if buffer.len() >= REGION_BUFFER_SIZE {
                self.flush_region(region_id)?;
            }
        }
        
        Ok(())
    }
    
    /// Flush fragments for a specific region to disk
    fn flush_region(&mut self, region_id: u32) -> Result<()> {
        if let Some(fragments) = self.region_buffers.remove(&region_id) {
            if fragments.is_empty() {
                return Ok(());
            }
            
            let (region_start, region_end) = self.regions.get(&region_id)
                .copied()
                .unwrap_or((0, 0));
            
            // Record index entry
            let index_entry = RegionIndex {
                region_id,
                region_start,
                region_end,
                file_offset: self.current_offset,
                fragment_count: fragments.len() as u32,
            };
            self.region_indices.push(index_entry);
            
            // Write region data
            self.writer.write_u32::<LittleEndian>(region_id)?;
            self.writer.write_u32::<LittleEndian>(fragments.len() as u32)?;
            
            let fragment_count = fragments.len();
            for (start, end) in fragments {
                self.writer.write_u32::<LittleEndian>(start)?;
                self.writer.write_u32::<LittleEndian>(end)?;
            }
            
            // Update offset
            self.current_offset += 4 + 4 + (8 * fragment_count as u64);
        }
        
        Ok(())
    }
    
    /// Flush all remaining regions
    pub fn flush_all_regions(&mut self) -> Result<()> {
        let region_ids: Vec<u32> = self.region_buffers.keys().copied().collect();
        for region_id in region_ids {
            self.flush_region(region_id)?;
        }
        Ok(())
    }
    
    /// Finish writing and create the index
    pub fn finish(mut self) -> Result<()> {
        // Flush any remaining buffered regions
        self.flush_all_regions()?;
        
        // Write index section
        let index_offset = self.current_offset;
        
        // Write index header
        self.writer.write_u32::<LittleEndian>(self.region_indices.len() as u32)?;
        
        // Write index entries
        for entry in &self.region_indices {
            self.writer.write_u32::<LittleEndian>(entry.region_id)?;
            self.writer.write_u64::<LittleEndian>(entry.region_start)?;
            self.writer.write_u64::<LittleEndian>(entry.region_end)?;
            self.writer.write_u64::<LittleEndian>(entry.file_offset)?;
            self.writer.write_u32::<LittleEndian>(entry.fragment_count)?;
        }
        
        // Flush writer
        self.writer.flush()?;
        
        // Update header with actual counts
        let mut file = self.writer.into_inner()?;
        file.seek(std::io::SeekFrom::Start(4 + 1 + 1 + self.chromosome.len() as u64))?;
        file.write_u64::<LittleEndian>(self.total_fragments)?;
        file.write_u64::<LittleEndian>(index_offset)?;
        file.sync_all()?;
        
        Ok(())
    }
}

/// Manager for multiple chromosome writers
pub struct ChromosomeWriterManager {
    writers: HashMap<String, ChromosomeFragmentWriter>,
    base_dir: PathBuf,
    peak_regions: Vec<ScatrsRegion>,
    bg_regions: Vec<ScatrsRegion>,
}

impl ChromosomeWriterManager {
    pub fn new(base_dir: &Path, peak_regions: Vec<ScatrsRegion>, bg_regions: Vec<ScatrsRegion>) -> Self {
        Self {
            writers: HashMap::new(),
            base_dir: base_dir.to_path_buf(),
            peak_regions,
            bg_regions,
        }
    }
    
    /// Get or create a writer for the given chromosome
    fn get_writer(&mut self, chromosome: &str) -> Result<&mut ChromosomeFragmentWriter> {
        if !self.writers.contains_key(chromosome) {
            let filename = format!("{}_fragments.bin", chromosome);
            let path = self.base_dir.join(filename);
            let writer = ChromosomeFragmentWriter::new(&path, chromosome.to_string())?;
            self.writers.insert(chromosome.to_string(), writer);
        }
        
        Ok(self.writers.get_mut(chromosome).unwrap())
    }
    
    /// Add a fragment to the appropriate chromosome writer
    pub fn add_fragment(&mut self, region_idx: usize, is_peak: bool, fragment: &ScatrsFragment) -> Result<()> {
        // Get region and clone it to avoid borrow issues
        let region = if is_peak {
            self.peak_regions.get(region_idx)
                .ok_or_else(|| anyhow::anyhow!("Invalid peak region index: {}", region_idx))?
                .clone()
        } else {
            self.bg_regions.get(region_idx)
                .ok_or_else(|| anyhow::anyhow!("Invalid background region index: {}", region_idx))?
                .clone()
        };
        
        // SIMPLIFIED: Use sequential IDs - peaks use their index, background uses peak_count + index
        let region_id = if is_peak {
            region_idx as u32
        } else {
            self.peak_regions.len() as u32 + region_idx as u32
        };
        
        
        let writer = self.get_writer(&fragment.chrom)?;
        writer.add_fragment(region_id, &region, fragment)
    }
    
    /// Finish all writers in parallel
    pub fn finish_all(self) -> Result<()> {
        use rayon::prelude::*;
        
        let writers: Vec<ChromosomeFragmentWriter> = self.writers.into_iter()
            .map(|(_, w)| w)
            .collect();
        
        // Process writers in parallel
        writers.into_par_iter()
            .try_for_each(|writer| writer.finish())?;
        
        Ok(())
    }
}
use std::io::{BufWriter, Write, Seek, Read};
use std::fs::File;
use std::path::{Path, PathBuf};
use std::collections::HashMap;
use anyhow::{Result, Context};
use byteorder::{LittleEndian, WriteBytesExt};
use crate::models::ScatrsFragment;
use rayon::prelude::*;

const FORMAT_VERSION: u8 = 1;
const BUFFER_SIZE: usize = 10000; // Buffer fragments per chromosome before flushing
const MAX_OPEN_WRITERS: usize = 1000; // Maximum number of simultaneously open file writers

/// Binary format header for fragment cache files
#[derive(Debug)]
pub struct BinHeader {
    pub version: u8,
    pub bin_id: u32,
    pub total_fragment_count: u64,
    pub chromosome_count: u32,
}

/// Writer for binary fragment cache files
/// 
/// Streams fragments directly from BAM to binary files without accumulating
/// all fragments in memory. Organizes fragments by chromosome for efficient storage.
pub struct FragmentWriter {
    writer: BufWriter<File>,
    bin_id: u32,
    current_chr: Option<String>,
    chr_fragments: Vec<(u32, u32)>, // (start, end) pairs
    fragment_count: u64,
    chr_blocks: Vec<(String, u32)>, // Track (chromosome, fragment_count) for header
}

impl FragmentWriter {
    pub fn new(path: &Path, bin_id: u32) -> Result<Self> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create fragment cache file: {:?}", path))?;
        
        let mut writer = BufWriter::new(file);
        
        // Write placeholder header (will be updated when closing)
        writer.write_u8(FORMAT_VERSION)?;
        writer.write_u32::<LittleEndian>(bin_id)?;
        writer.write_u64::<LittleEndian>(0)?; // placeholder for total_fragment_count
        writer.write_u32::<LittleEndian>(0)?; // placeholder for chromosome_count
        
        Ok(Self {
            writer,
            bin_id,
            current_chr: None,
            chr_fragments: Vec::with_capacity(BUFFER_SIZE),
            fragment_count: 0,
            chr_blocks: Vec::new(),
        })
    }
    
    /// Open an existing file for appending more fragments
    pub fn append(path: &Path, bin_id: u32) -> Result<Self> {
        use std::fs::OpenOptions;
        use byteorder::ReadBytesExt;
        
        // Read existing header to get counts
        let mut file = File::open(path)?;
        let mut version = [0u8; 1];
        file.read_exact(&mut version)?;
        let _stored_bin_id = file.read_u32::<LittleEndian>()?;
        let existing_fragment_count = file.read_u64::<LittleEndian>()?;
        let _existing_chr_count = file.read_u32::<LittleEndian>()?;
        
        // TODO: Read existing chromosome blocks to populate chr_blocks
        let chr_blocks = Vec::new(); // Simplified for now
        
        // Open for append
        let file = OpenOptions::new()
            .write(true)
            .append(true)
            .open(path)?;
        
        let writer = BufWriter::new(file);
        
        Ok(Self {
            writer,
            bin_id,
            current_chr: None,
            chr_fragments: Vec::with_capacity(BUFFER_SIZE),
            fragment_count: existing_fragment_count,
            chr_blocks,
        })
    }
    
    /// Add a fragment to the writer
    /// 
    /// Fragments should be added in chromosome-sorted order (as they typically
    /// come from BAM files) for optimal efficiency.
    pub fn add_fragment(&mut self, fragment: &ScatrsFragment) -> Result<()> {
        // Check if we've switched to a new chromosome
        if self.current_chr.as_ref() != Some(&fragment.chrom) {
            // Flush previous chromosome if exists
            if self.current_chr.is_some() {
                self.flush_chromosome_block()?;
            }
            self.current_chr = Some(fragment.chrom.clone());
        }
        
        // Add fragment to buffer (store as u32 to save space)
        self.chr_fragments.push((fragment.start as u32, fragment.end as u32));
        self.fragment_count += 1;
        
        // Flush if buffer is full
        if self.chr_fragments.len() >= BUFFER_SIZE {
            self.flush_chromosome_block()?;
        }
        
        Ok(())
    }
    
    /// Flush the current chromosome block to disk
    fn flush_chromosome_block(&mut self) -> Result<()> {
        if self.chr_fragments.is_empty() {
            return Ok(());
        }
        
        let chr_name = self.current_chr.as_ref().unwrap();
        let fragment_count = self.chr_fragments.len() as u32;
        
        // Write chromosome block header
        let chr_name_bytes = chr_name.as_bytes();
        self.writer.write_u8(chr_name_bytes.len() as u8)?;
        self.writer.write_all(chr_name_bytes)?;
        self.writer.write_u32::<LittleEndian>(fragment_count)?;
        
        // Write fragments
        for (start, end) in &self.chr_fragments {
            self.writer.write_u32::<LittleEndian>(*start)?;
            self.writer.write_u32::<LittleEndian>(*end)?;
        }
        
        // Track this chromosome block
        self.chr_blocks.push((chr_name.clone(), fragment_count));
        
        // Clear buffer
        self.chr_fragments.clear();
        
        Ok(())
    }
    
    /// Finish writing and update header with final counts
    pub fn finish(mut self) -> Result<()> {
        // Flush any remaining fragments
        if !self.chr_fragments.is_empty() {
            self.flush_chromosome_block()?;
        }
        
        // Flush writer to ensure all data is written
        self.writer.flush()?;
        
        // Rewind and update header with actual counts
        let mut file = self.writer.into_inner()?;
        file.seek(std::io::SeekFrom::Start(0))?;
        
        file.write_u8(FORMAT_VERSION)?;
        file.write_u32::<LittleEndian>(self.bin_id)?;
        file.write_u64::<LittleEndian>(self.fragment_count)?;
        file.write_u32::<LittleEndian>(self.chr_blocks.len() as u32)?;
        
        file.sync_all()?;
        
        Ok(())
    }
}

/// Manager for multiple fragment writers (one per bin)
/// Uses LRU-style management to limit open file handles
pub struct FragmentWriterManager {
    writers: HashMap<u32, FragmentWriter>,
    pending_fragments: HashMap<u32, Vec<ScatrsFragment>>, // Buffer for bins without writers
    base_dir: PathBuf,
    is_peak: bool,
    access_order: Vec<u32>, // Track access order for LRU eviction
}

impl FragmentWriterManager {
    pub fn new(base_dir: &Path, is_peak: bool) -> Self {
        Self {
            writers: HashMap::new(),
            pending_fragments: HashMap::new(),
            base_dir: base_dir.to_path_buf(),
            is_peak,
            access_order: Vec::new(),
        }
    }
    
    /// Get or create a writer for the given bin, managing file handle limits
    fn get_writer(&mut self, bin_id: u32) -> Result<&mut FragmentWriter> {
        // Update access order for LRU
        if let Some(pos) = self.access_order.iter().position(|&id| id == bin_id) {
            self.access_order.remove(pos);
        }
        self.access_order.push(bin_id);
        
        // Check if we need to evict a writer to stay under the limit
        if !self.writers.contains_key(&bin_id) && self.writers.len() >= MAX_OPEN_WRITERS {
            // Evict least recently used writer
            if let Some(lru_id) = self.access_order.first().copied() {
                if lru_id != bin_id {
                    if let Some(writer) = self.writers.remove(&lru_id) {
                        // Finish and close the evicted writer
                        writer.finish()?;
                        self.access_order.retain(|&id| id != lru_id);
                    }
                }
            }
        }
        
        // Create writer if it doesn't exist
        if !self.writers.contains_key(&bin_id) {
            let prefix = if self.is_peak { "peak_bin_" } else { "bg_bin_" };
            let filename = format!("{}{}.bin", prefix, bin_id);
            let path = self.base_dir.join(filename);
            
            let mut writer = if path.exists() {
                // Append mode for previously evicted writers
                FragmentWriter::append(&path, bin_id)?
            } else {
                FragmentWriter::new(&path, bin_id)?
            };
            
            // Flush any pending fragments for this bin
            if let Some(pending) = self.pending_fragments.remove(&bin_id) {
                for fragment in pending {
                    writer.add_fragment(&fragment)?;
                }
            }
            
            self.writers.insert(bin_id, writer);
        }
        
        Ok(self.writers.get_mut(&bin_id).unwrap())
    }
    
    /// Add a fragment to the specified bin
    pub fn add_fragment(&mut self, bin_id: u32, fragment: &ScatrsFragment) -> Result<()> {
        // For bins that are accessed frequently, use direct writer
        // For rarely accessed bins, buffer in memory first
        if self.writers.contains_key(&bin_id) || self.writers.len() < MAX_OPEN_WRITERS {
            let writer = self.get_writer(bin_id)?;
            writer.add_fragment(fragment)
        } else {
            // Buffer fragment for later batch writing
            self.pending_fragments.entry(bin_id)
                .or_insert_with(Vec::new)
                .push(fragment.clone());
            
            // If buffer gets too large, flush it
            if let Some(pending) = self.pending_fragments.get(&bin_id) {
                if pending.len() >= BUFFER_SIZE {
                    let fragments = self.pending_fragments.remove(&bin_id).unwrap();
                    let writer = self.get_writer(bin_id)?;
                    for frag in fragments {
                        writer.add_fragment(&frag)?;
                    }
                }
            }
            Ok(())
        }
    }
    
    /// Finish all writers and close files in parallel
    pub fn finish_all(mut self) -> Result<()> {
        // First, flush all pending fragments to their writers
        let pending: Vec<(u32, Vec<ScatrsFragment>)> = self.pending_fragments.drain().collect();
        for (bin_id, fragments) in pending {
            if !fragments.is_empty() {
                let writer = self.get_writer(bin_id)?;
                for fragment in fragments {
                    writer.add_fragment(&fragment)?;
                }
            }
        }
        
        // Convert HashMap into a Vec for parallel processing
        let writers: Vec<(u32, FragmentWriter)> = self.writers.into_iter().collect();
        
        // Process writers in parallel using rayon
        writers.into_par_iter()
            .try_for_each(|(_, writer)| writer.finish())?;
        
        Ok(())
    }
}
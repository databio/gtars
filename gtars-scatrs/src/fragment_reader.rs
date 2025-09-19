use std::io::{BufReader, Read};
use std::fs::File;
use std::path::Path;
use anyhow::{Result, Context, anyhow};
use byteorder::{LittleEndian, ReadBytesExt};
use crate::models::ScatrsFragment;

const FORMAT_VERSION: u8 = 1;

/// Reader for binary fragment cache files
/// 
/// Provides streaming access to cached fragments, reconstructing full
/// ScatrsFragment objects with chromosome information.
pub struct FragmentReader {
    reader: BufReader<File>,
    bin_id: u32,
    total_fragment_count: u64,
    chromosome_count: u32,
    current_chr: Option<String>,
    current_chr_fragments: Vec<(u32, u32)>,
    current_fragment_idx: usize,
    chromosomes_read: u32,
}

impl FragmentReader {
    pub fn new(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open fragment cache file: {:?}", path))?;
        
        let mut reader = BufReader::new(file);
        
        // Read header
        let version = reader.read_u8()?;
        if version != FORMAT_VERSION {
            return Err(anyhow!("Unsupported fragment cache format version: {}", version));
        }
        
        let bin_id = reader.read_u32::<LittleEndian>()?;
        let total_fragment_count = reader.read_u64::<LittleEndian>()?;
        let chromosome_count = reader.read_u32::<LittleEndian>()?;
        
        Ok(Self {
            reader,
            bin_id,
            total_fragment_count,
            chromosome_count,
            current_chr: None,
            current_chr_fragments: Vec::new(),
            current_fragment_idx: 0,
            chromosomes_read: 0,
        })
    }
    
    /// Get total fragment count in this bin
    pub fn fragment_count(&self) -> u64 {
        self.total_fragment_count
    }
    
    /// Get bin ID
    pub fn bin_id(&self) -> u32 {
        self.bin_id
    }
    
    /// Read the next chromosome block from the file
    fn read_next_chromosome_block(&mut self) -> Result<bool> {
        if self.chromosomes_read >= self.chromosome_count {
            return Ok(false); // No more chromosomes
        }
        
        // Read chromosome name
        let chr_name_len = self.reader.read_u8()? as usize;
        let mut chr_name_bytes = vec![0u8; chr_name_len];
        self.reader.read_exact(&mut chr_name_bytes)?;
        let chr_name = String::from_utf8(chr_name_bytes)
            .context("Invalid UTF-8 in chromosome name")?;
        
        // Read fragment count for this chromosome
        let fragment_count = self.reader.read_u32::<LittleEndian>()? as usize;
        
        // Read all fragments for this chromosome
        let mut fragments = Vec::with_capacity(fragment_count);
        for _ in 0..fragment_count {
            let start = self.reader.read_u32::<LittleEndian>()?;
            let end = self.reader.read_u32::<LittleEndian>()?;
            fragments.push((start, end));
        }
        
        self.current_chr = Some(chr_name);
        self.current_chr_fragments = fragments;
        self.current_fragment_idx = 0;
        self.chromosomes_read += 1;
        
        Ok(true)
    }
    
    /// Get the next fragment from the reader
    pub fn next_fragment(&mut self) -> Result<Option<ScatrsFragment>> {
        // Check if we need to load the next chromosome block
        if self.current_fragment_idx >= self.current_chr_fragments.len() {
            if !self.read_next_chromosome_block()? {
                return Ok(None); // No more fragments
            }
        }
        
        // Get the current fragment
        let chr = self.current_chr.as_ref().unwrap().clone();
        let (start, end) = self.current_chr_fragments[self.current_fragment_idx];
        self.current_fragment_idx += 1;
        
        Ok(Some(ScatrsFragment::new(chr, start as u64, end as u64)))
    }
}

/// Iterator adapter for FragmentReader
pub struct FragmentIterator {
    reader: FragmentReader,
}

impl FragmentIterator {
    pub fn new(reader: FragmentReader) -> Self {
        Self { reader }
    }
}

impl Iterator for FragmentIterator {
    type Item = Result<ScatrsFragment>;
    
    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.next_fragment() {
            Ok(Some(fragment)) => Some(Ok(fragment)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

impl FragmentReader {
    /// Create an iterator over all fragments in the file
    pub fn iter_fragments(self) -> FragmentIterator {
        FragmentIterator::new(self)
    }
    
    /// Collect all fragments into a vector (use with caution on large files)
    pub fn collect_all(mut self) -> Result<Vec<ScatrsFragment>> {
        let mut fragments = Vec::with_capacity(self.total_fragment_count as usize);
        
        while let Some(fragment) = self.next_fragment()? {
            fragments.push(fragment);
        }
        
        Ok(fragments)
    }
}
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::fs::File;
use std::path::{Path, PathBuf};
use std::collections::HashMap;
use anyhow::{Result, Context, anyhow};
use byteorder::{LittleEndian, ReadBytesExt};
use crate::models::ScatrsFragment;

const MAGIC_NUMBER: &[u8; 4] = b"SCFR";
const FORMAT_VERSION: u8 = 2;  // Version 2: simplified sequential IDs

/// Index entry for fast region lookup
#[derive(Debug, Clone)]
pub struct RegionIndex {
    pub region_id: u32,
    pub region_start: u64,
    pub region_end: u64,
    pub file_offset: u64,
    pub fragment_count: u32,
}

/// Reader for chromosome-based fragment cache files
/// 
/// Provides efficient random access to fragments by region using an embedded index.
pub struct ChromosomeFragmentReader {
    reader: BufReader<File>,
    chromosome: String,
    total_fragments: u64,
    index: HashMap<u32, RegionIndex>,
}

impl ChromosomeFragmentReader {
    pub fn new(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open chromosome fragment file: {:?}", path))?;
        
        let mut reader = BufReader::with_capacity(1024 * 1024, file); // 1MB buffer
        
        // Read and validate header
        let mut magic = [0u8; 4];
        reader.read_exact(&mut magic)?;
        if &magic != MAGIC_NUMBER {
            return Err(anyhow!("Invalid file format: magic number mismatch"));
        }
        
        let version = reader.read_u8()?;
        if version != FORMAT_VERSION {
            return Err(anyhow!("Unsupported format version: {}", version));
        }
        
        // Read chromosome name
        let chr_len = reader.read_u8()? as usize;
        let mut chr_bytes = vec![0u8; chr_len];
        reader.read_exact(&mut chr_bytes)?;
        let chromosome = String::from_utf8(chr_bytes)
            .context("Invalid UTF-8 in chromosome name")?;
        
        // Read counts and index offset
        let total_fragments = reader.read_u64::<LittleEndian>()?;
        let index_offset = reader.read_u64::<LittleEndian>()?;
        
        // Read index from end of file
        reader.seek(SeekFrom::Start(index_offset))?;
        let index = Self::read_index(&mut reader)?;
        
        Ok(Self {
            reader,
            chromosome,
            total_fragments,
            index,
        })
    }
    
    /// Read the index section from the file
    fn read_index(reader: &mut BufReader<File>) -> Result<HashMap<u32, RegionIndex>> {
        let region_count = reader.read_u32::<LittleEndian>()?;
        let mut index = HashMap::with_capacity(region_count as usize);
        
        for _ in 0..region_count {
            let region_id = reader.read_u32::<LittleEndian>()?;
            let region_start = reader.read_u64::<LittleEndian>()?;
            let region_end = reader.read_u64::<LittleEndian>()?;
            let file_offset = reader.read_u64::<LittleEndian>()?;
            let fragment_count = reader.read_u32::<LittleEndian>()?;
            
            index.insert(region_id, RegionIndex {
                region_id,
                region_start,
                region_end,
                file_offset,
                fragment_count,
            });
}
        
        Ok(index)
    }
    
    /// Get the chromosome name
    pub fn chromosome(&self) -> &str {
        &self.chromosome
    }
    
    /// Get total fragment count
    pub fn fragment_count(&self) -> u64 {
        self.total_fragments
    }
    
    /// Get fragments for a specific region
    pub fn get_fragments_for_region(&mut self, region_id: u32) -> Result<Vec<ScatrsFragment>> {
        let index_entry = self.index.get(&region_id)
            .ok_or_else(|| {
                // Provide detailed error message with debugging information
                let total_regions = self.index.len();
                
                anyhow!(
                    "Region {} has no cached fragments on chromosome {}.\n\
                    This is EXPECTED for regions with zero coverage in the input BAM files.\n\
                    The cache only stores regions that contain fragments ({} regions with data).\n\
                    This region will be skipped in the simulation.",
                    region_id, self.chromosome, total_regions
                )
            })?
            .clone();
        
        // Seek to the region's data
        self.reader.seek(SeekFrom::Start(index_entry.file_offset))?;
        
        // Read region header
        let stored_region_id = self.reader.read_u32::<LittleEndian>()?;
        if stored_region_id != region_id {
            return Err(anyhow!("Region ID mismatch: expected {}, got {}", region_id, stored_region_id));
        }
        
        let fragment_count = self.reader.read_u32::<LittleEndian>()?;
        if fragment_count != index_entry.fragment_count {
            return Err(anyhow!("Fragment count mismatch for region {}", region_id));
        }
        
        // Read fragments
        let mut fragments = Vec::with_capacity(fragment_count as usize);
        for _ in 0..fragment_count {
            let start = self.reader.read_u32::<LittleEndian>()? as u64;
            let end = self.reader.read_u32::<LittleEndian>()? as u64;
            fragments.push(ScatrsFragment::new(
                self.chromosome.clone(),
                start,
                end,
            ));
        }
        
        Ok(fragments)
    }
    
    /// Get fragments for multiple regions efficiently
    pub fn get_fragments_for_regions(&mut self, region_ids: &[u32]) -> Result<HashMap<u32, Vec<ScatrsFragment>>> {
        let mut result = HashMap::new();
        
        // Sort regions by file offset for sequential reading
        let mut sorted_regions: Vec<(u32, u64)> = region_ids.iter()
            .filter_map(|&id| {
                self.index.get(&id).map(|entry| (id, entry.file_offset))
            })
            .collect();
        sorted_regions.sort_by_key(|&(_, offset)| offset);
        
        for (region_id, _) in sorted_regions {
            let fragments = self.get_fragments_for_region(region_id)?;
            result.insert(region_id, fragments);
        }
        
        Ok(result)
    }
    
    /// Iterate over all fragments in the file
    pub fn iter_all_fragments(&mut self) -> Result<impl Iterator<Item = Result<ScatrsFragment>> + '_> {
        // Get all region IDs sorted by file offset
        let mut sorted_regions: Vec<u32> = self.index.keys().copied().collect();
        sorted_regions.sort_by_key(|&id| self.index[&id].file_offset);
        
        Ok(ChromosomeFragmentIterator {
            reader: self,
            region_ids: sorted_regions,
            current_region: 0,
            current_fragments: Vec::new(),
            current_fragment_idx: 0,
        })
    }
}

/// Iterator over all fragments in a chromosome file
pub struct ChromosomeFragmentIterator<'a> {
    reader: &'a mut ChromosomeFragmentReader,
    region_ids: Vec<u32>,
    current_region: usize,
    current_fragments: Vec<ScatrsFragment>,
    current_fragment_idx: usize,
}

impl<'a> Iterator for ChromosomeFragmentIterator<'a> {
    type Item = Result<ScatrsFragment>;
    
    fn next(&mut self) -> Option<Self::Item> {
        // Check if we have fragments in current batch
        if self.current_fragment_idx < self.current_fragments.len() {
            let fragment = self.current_fragments[self.current_fragment_idx].clone();
            self.current_fragment_idx += 1;
            return Some(Ok(fragment));
        }
        
        // Load next region's fragments
        if self.current_region < self.region_ids.len() {
            let region_id = self.region_ids[self.current_region];
            self.current_region += 1;
            
            match self.reader.get_fragments_for_region(region_id) {
                Ok(fragments) => {
                    self.current_fragments = fragments;
                    self.current_fragment_idx = 0;
                    // Recursively call to get first fragment from new batch
                    self.next()
                }
                Err(e) => Some(Err(e)),
            }
        } else {
            None
        }
    }
}

/// Manager for reading fragments across multiple chromosomes
pub struct ChromosomeReaderManager {
    base_dir: PathBuf,
    readers: HashMap<String, ChromosomeFragmentReader>,
}

impl ChromosomeReaderManager {
    pub fn new(base_dir: &Path) -> Self {
        Self {
            base_dir: base_dir.to_path_buf(),
            readers: HashMap::new(),
        }
    }
    
    /// Get or open a reader for the given chromosome
    pub fn get_reader(&mut self, chromosome: &str) -> Result<&mut ChromosomeFragmentReader> {
        if !self.readers.contains_key(chromosome) {
            let filename = format!("{}_fragments.bin", chromosome);
            let path = self.base_dir.join(filename);
            
            if !path.exists() {
                return Err(anyhow!("Chromosome fragment file not found: {:?}", path));
            }
            
            let reader = ChromosomeFragmentReader::new(&path)?;
            self.readers.insert(chromosome.to_string(), reader);
        }
        
        Ok(self.readers.get_mut(chromosome).unwrap())
    }
    
    /// Check if fragment cache exists for any chromosome
    pub fn has_cache(&self) -> bool {
        // Check for any chromosome fragment files
        if let Ok(entries) = std::fs::read_dir(&self.base_dir) {
            for entry in entries.flatten() {
                if let Some(name) = entry.file_name().to_str() {
                    if name.ends_with("_fragments.bin") {
                        return true;
                    }
                }
            }
        }
        false
    }
    
    /// Get fragments for specific regions across chromosomes
    pub fn get_fragments_for_regions(
        &mut self,
        regions_by_chr: HashMap<String, Vec<u32>>,
    ) -> Result<HashMap<u32, Vec<ScatrsFragment>>> {
        let mut all_fragments = HashMap::new();
        
        for (chromosome, region_ids) in regions_by_chr {
            let reader = self.get_reader(&chromosome)?;
            let chr_fragments = reader.get_fragments_for_regions(&region_ids)?;
            all_fragments.extend(chr_fragments);
        }
        
        Ok(all_fragments)
    }
}
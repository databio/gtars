use crate::models::{ScatrsFragment, ScatrsRegion, Peak};
use crate::consts::VALID_CHROMOSOMES;
use crate::chr_frag_writer::{ChromosomeWriterManager};
use crate::chr_frag_reader::{ChromosomeReaderManager};
use noodles::bam;
use noodles::bgzf;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::Record as SamRecord;
use noodles::sam;
use std::path::Path;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::collections::HashMap;
use anyhow::{Result, Context, anyhow};
use indicatif::{ProgressBar, ProgressStyle};

// ============================================================================
// BAM Streaming
// ============================================================================

/// Streaming BAM reader that processes fragments one at a time
/// 
/// CORE EFFICIENCY PRINCIPLE: Never load all fragments into memory.
/// This struct provides true streaming capability for BAM files,
/// processing and caching fragments during the single pass through the BAM.
/// Reservoir sampling happens DURING the stream, not after loading.
pub struct BamStreamer {
    reader: bam::io::Reader<bgzf::Reader<std::fs::File>>,
    header: sam::Header,
    total_reads: u64,
    filtered_reads: u64,
    progress_bar: Option<ProgressBar>,
}

impl BamStreamer {
    pub fn new(bam_path: &Path) -> Result<Self> {
        let file = File::open(bam_path)
            .with_context(|| format!("Failed to open BAM file: {:?}", bam_path))?;
        let mut reader = bam::io::reader::Builder::default()
            .build_from_reader(file);
            
        let header = reader.read_header()
            .with_context(|| "Failed to read BAM header")?;
        
        let pb = ProgressBar::new_spinner();
        pb.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} [{elapsed_precise}] {msg} ({pos} fragments)")
                .unwrap()
        );
        pb.set_message(format!("Streaming BAM file: {:?}", bam_path.file_name().unwrap_or_default()));
        
        Ok(Self {
            reader,
            header,
            total_reads: 0,
            filtered_reads: 0,
            progress_bar: Some(pb),
        })
    }
    
    /// Read the next fragment from the BAM file
    /// Returns None when no more fragments are available
    pub fn next_fragment(&mut self) -> Result<Option<ScatrsFragment>> {
        let mut record = bam::Record::default();
        
        loop {
            match self.reader.read_record(&mut record) {
                Ok(0) => {
                    // End of file (no bytes read)
                    if let Some(ref pb) = self.progress_bar {
                        pb.finish_with_message(format!(
                            "Processed {} reads, kept {} fragments",
                            self.total_reads, self.filtered_reads
                        ));
                    }
                    return Ok(None);
                }
                Ok(_) => {
                    self.total_reads += 1;
                    
                    // Update progress bar much less frequently (every 1M reads)
                    if self.total_reads % 1000000 == 0 {
                        if let Some(ref pb) = self.progress_bar {
                            pb.set_position(self.filtered_reads);
                        }
                    }
                    
                    let flags = record.flags();
                    
                    // Filter for paired, mapped, non-duplicate reads
                    if !Self::should_keep_read(&flags) {
                        continue;
                    }
                    
                    // Process only the first read of the pair to avoid duplicates
                    if !flags.is_first_segment() {
                        continue;
                    }
                    
                    // Check if mate is mapped and on same chromosome
                    if flags.is_mate_unmapped() {
                        continue;
                    }
                    
                    // Extract chromosome information
                    let chrom = match record.reference_sequence_id() {
                        Some(id) => {
                            let idx = match id {
                                Ok(i) => i,
                                _ => continue,
                            };
                            let name = self.header.reference_sequences().get_index(idx)
                                .map(|(name, _)| name.to_string());
                            match name {
                                Some(n) => n,
                                None => continue,
                            }
                        },
                        None => continue,
                    };
                    
                    // Check if mate is on same chromosome
                    let mate_chrom_id = match record.mate_reference_sequence_id() {
                        Some(id) => match id {
                            Ok(i) => i,
                            _ => continue,
                        },
                        None => continue,
                    };
                    
                    let read_chrom_id = match record.reference_sequence_id() {
                        Some(id) => match id {
                            Ok(i) => i,
                            _ => continue,
                        },
                        None => continue,
                    };
                    
                    // Skip if mate is on different chromosome
                    if mate_chrom_id != read_chrom_id {
                        continue;
                    }
                    
                    // Get read position and length
                    let read_start = match record.alignment_start() {
                        Some(pos) => match pos {
                            Ok(p) => p.get() as u64 - 1, // Convert to 0-based
                            _ => continue,
                        },
                        _ => continue,
                    };
                    
                    let read_span = match record.alignment_span() {
                        Some(s) => match s {
                            Ok(sp) => sp as u64,
                            _ => continue,
                        },
                        _ => continue,
                    };
                    
                    // Get mate position
                    let mate_start = match record.mate_alignment_start() {
                        Some(pos) => match pos {
                            Ok(p) => p.get() as u64 - 1, // Convert to 0-based
                            _ => continue,
                        },
                        _ => continue,
                    };
                    
                    // For paired-end reads, the fragment spans from the leftmost position
                    // to the rightmost position of the pair
                    // Assume mate has similar read length as current read for end calculation
                    let fragment_start = read_start.min(mate_start);
                    let fragment_end = if read_start < mate_start {
                        // Current read is leftmost, mate is rightmost
                        mate_start + read_span // Approximate mate end with read length
                    } else {
                        // Mate is leftmost, current read is rightmost
                        read_start + read_span
                    };
                    
                    // Validate fragment size (typical ATAC-seq fragments are 50-2000bp)
                    let fragment_length = fragment_end - fragment_start;
                    if fragment_length < 50 || fragment_length > 2000 {
                        // Log warning for unusual fragment sizes if debug enabled
                        if std::env::var("SCATRS_DEBUG").is_ok() && self.filtered_reads < 10 {
                            eprintln!("Fragment length {} outside expected range [50-2000] for {}:{}-{}", 
                                     fragment_length, chrom, fragment_start, fragment_end);
                        }
                        // Still keep the fragment, just warn
                    }
                    
                    let strand = if flags.is_reverse_complemented() {
                        Some('-')
                    } else {
                        Some('+')
                    };
                    
                    let quality = record.mapping_quality()
                        .map(|q| q.get())
                        .unwrap_or(0);
                    
                    let fragment = ScatrsFragment::new(chrom, fragment_start, fragment_end)
                        .with_strand(strand.unwrap())
                        .with_quality(quality);
                    
                    // Debug logging for first few fragments
                    if std::env::var("SCATRS_DEBUG").is_ok() && self.filtered_reads < 10 {
                        eprintln!("Fragment {}: {}:{}-{} (length: {}bp)", 
                                 self.filtered_reads + 1, fragment.chrom, fragment.start, 
                                 fragment.end, fragment.end - fragment.start);
                    }
                    
                    self.filtered_reads += 1;
                    return Ok(Some(fragment));
                }
                Err(e) => {
                    eprintln!("Warning: Error reading BAM record: {}", e);
                    continue;
                }
            }
        }
    }
    
    fn should_keep_read(flags: &Flags) -> bool {
        // Keep only properly paired, mapped, non-duplicate reads
        flags.is_segmented()  // equivalent to is_paired
            && flags.is_properly_segmented()  // equivalent to is_properly_aligned
            && !flags.is_unmapped()
            && !flags.is_duplicate()
            && !flags.is_secondary()
            && !flags.is_supplementary()
    }
    
    pub fn get_stats(&self) -> (u64, u64) {
        (self.total_reads, self.filtered_reads)
    }
}

// ============================================================================
// BAM Processing
// ============================================================================

/// Processes BAM files to extract ATAC-seq fragments
/// 
/// Filters for properly paired, non-duplicate reads and extracts
/// fragment coordinates representing accessible chromatin regions.
pub struct BamProcessor;

impl BamProcessor {
    pub fn filter_fragments(bam_path: &Path) -> Result<Vec<ScatrsFragment>> {
        let file = File::open(bam_path)
            .with_context(|| format!("Failed to open BAM file: {:?}", bam_path))?;
        let mut reader = bam::io::reader::Builder::default()
            .build_from_reader(file);
            
        let header = reader.read_header()
            .with_context(|| "Failed to read BAM header")?;
        
        let pb = ProgressBar::new_spinner();
        pb.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} [{elapsed_precise}] {msg} ({pos} fragments)")
                .unwrap()
        );
        pb.set_message(format!("Processing BAM file: {:?}", bam_path.file_name().unwrap_or_default()));
        
        let mut fragments = Vec::new();
        let mut total_reads = 0u64;
        let mut filtered_reads = 0u64;
        
        for result in reader.records() {
            match result {
                Ok(record) => {
                    total_reads += 1;
                    
                    if total_reads % 10000 == 0 {
                        pb.set_position(filtered_reads);
                    }
                    
                    let flags = record.flags();
                    
                    // Filter for paired, mapped, non-duplicate reads
                    if !Self::should_keep_read(&flags) {
                        continue;
                    }
                    
                    // Process only the first read of the pair to avoid duplicates
                    if !flags.is_first_segment() {
                        continue;
                    }
                    
                    // Check if mate is mapped
                    if flags.is_mate_unmapped() {
                        continue;
                    }
                    
                    // Extract chromosome information
                    let chrom = match record.reference_sequence_id() {
                        Some(id) => {
                            let idx = match id {
                                Ok(i) => i,
                                _ => continue,
                            };
                            let name = header.reference_sequences().get_index(idx)
                                .map(|(name, _)| name.to_string());
                            match name {
                                Some(n) => n,
                                None => continue,
                            }
                        },
                        None => continue,
                    };
                    
                    // Check if mate is on same chromosome
                    let mate_chrom_id = match record.mate_reference_sequence_id() {
                        Some(id) => match id {
                            Ok(i) => i,
                            _ => continue,
                        },
                        None => continue,
                    };
                    
                    let read_chrom_id = match record.reference_sequence_id() {
                        Some(id) => match id {
                            Ok(i) => i,
                            _ => continue,
                        },
                        None => continue,
                    };
                    
                    // Skip if mate is on different chromosome
                    if mate_chrom_id != read_chrom_id {
                        continue;
                    }
                    
                    // Get read position and length
                    let read_start = match record.alignment_start() {
                        Some(pos) => match pos {
                            Ok(p) => p.get() as u64 - 1, // Convert to 0-based
                            _ => continue,
                        },
                        _ => continue,
                    };
                    
                    let read_span = match record.alignment_span() {
                        Some(s) => match s {
                            Ok(sp) => sp as u64,
                            _ => continue,
                        },
                        _ => continue,
                    };
                    
                    // Get mate position
                    let mate_start = match record.mate_alignment_start() {
                        Some(pos) => match pos {
                            Ok(p) => p.get() as u64 - 1, // Convert to 0-based
                            _ => continue,
                        },
                        _ => continue,
                    };
                    
                    // For paired-end reads, the fragment spans from the leftmost position
                    // to the rightmost position of the pair
                    // Assume mate has similar read length as current read for end calculation
                    let fragment_start = read_start.min(mate_start);
                    let fragment_end = if read_start < mate_start {
                        // Current read is leftmost, mate is rightmost
                        mate_start + read_span // Approximate mate end with read length
                    } else {
                        // Mate is leftmost, current read is rightmost
                        read_start + read_span
                    };
                    
                    // Validate fragment size (typical ATAC-seq fragments are 50-2000bp)
                    let fragment_length = fragment_end - fragment_start;
                    if fragment_length < 50 || fragment_length > 2000 {
                        // Log warning for unusual fragment sizes if debug enabled
                        if std::env::var("SCATRS_DEBUG").is_ok() && filtered_reads < 10 {
                            eprintln!("Fragment length {} outside expected range [50-2000] for {}:{}-{}", 
                                     fragment_length, chrom, fragment_start, fragment_end);
                        }
                        // Still keep the fragment, just warn
                    }
                    
                    let strand = if flags.is_reverse_complemented() {
                        Some('-')
                    } else {
                        Some('+')
                    };
                    
                    let quality = record.mapping_quality()
                        .map(|q| q.get())
                        .unwrap_or(0);
                    
                    let fragment = ScatrsFragment::new(chrom, fragment_start, fragment_end)
                        .with_strand(strand.unwrap())
                        .with_quality(quality);
                    
                    // Debug logging for first few fragments
                    if std::env::var("SCATRS_DEBUG").is_ok() && filtered_reads < 10 {
                        eprintln!("Fragment {}: {}:{}-{} (length: {}bp)", 
                                 filtered_reads + 1, fragment.chrom, fragment.start, 
                                 fragment.end, fragment.end - fragment.start);
                    }
                    
                    fragments.push(fragment);
                    filtered_reads += 1;
                }
                Err(e) => {
                    eprintln!("Warning: Error reading BAM record: {}", e);
                    continue;
                }
            }
        }
        
        pb.finish_with_message(format!(
            "Processed {} reads, kept {} fragments",
            total_reads, filtered_reads
        ));
        
        Ok(fragments)
    }
    
    fn should_keep_read(flags: &Flags) -> bool {
        // Keep only properly paired, mapped, non-duplicate reads
        flags.is_segmented()  // equivalent to is_paired
            && flags.is_properly_segmented()  // equivalent to is_properly_aligned
            && !flags.is_unmapped()
            && !flags.is_duplicate()
            && !flags.is_secondary()
            && !flags.is_supplementary()
    }
    
    pub fn count_fragments_in_regions(
        fragments: &[ScatrsFragment],
        regions: &[ScatrsRegion],
    ) -> Vec<u32> {
        use rayon::prelude::*;
        
        regions.par_iter().map(|region| {
            fragments.iter().filter(|frag| {
                frag.chrom == region.chrom 
                    && frag.start < region.end 
                    && frag.end > region.start
            }).count() as u32
        }).collect()
    }
    
    /// Stream through BAM file and count fragments using chromosome-based caching
    /// This dramatically reduces file count from 500K+ to ~25 files
    pub fn count_fragments_streaming_with_chromosome_cache(
        bam_path: &Path,
        peaks: &[ScatrsRegion],
        background: &[ScatrsRegion],
        cache_dir: Option<&Path>,
    ) -> Result<(Vec<u32>, Vec<u32>)> {
        use std::time::Instant;
        
        let mut peak_counts = vec![0u32; peaks.len()];
        let mut bg_counts = vec![0u32; background.len()];
        
        // Initialize chromosome writer if caching is enabled
        let mut chr_writer = cache_dir.map(|dir| {
            ChromosomeWriterManager::new(dir, peaks.to_vec(), background.to_vec())
        });
        
        // Show concise initial message
        let cache_msg = if cache_dir.is_some() { " and caching fragments by chrom" } else { "" };
        println!("  Counting fragments in {} peaks and {} background regions{}...", 
            peaks.len(), background.len(), cache_msg);
        
        let start_time = Instant::now();
        let mut streamer = BamStreamer::new(bam_path)?;
        let mut fragment_count = 0u64;
        let mut last_report_time = Instant::now();
        let mut regions_with_fragments = std::collections::HashSet::new();
        
        // Create chromosome-indexed maps for efficient lookup
        let mut peaks_by_chrom: HashMap<String, Vec<(usize, &ScatrsRegion)>> = HashMap::new();
        for (i, peak) in peaks.iter().enumerate() {
            peaks_by_chrom.entry(peak.chrom.clone())
                .or_insert_with(Vec::new)
                .push((i, peak));
        }
        
        let mut bg_by_chrom: HashMap<String, Vec<(usize, &ScatrsRegion)>> = HashMap::new();
        for (i, bg) in background.iter().enumerate() {
            bg_by_chrom.entry(bg.chrom.clone())
                .or_insert_with(Vec::new)
                .push((i, bg));
        }
        
        // Sort regions within each chromosome for efficient search
        for regions in peaks_by_chrom.values_mut() {
            regions.sort_by_key(|(_, r)| r.start);
        }
        for regions in bg_by_chrom.values_mut() {
            regions.sort_by_key(|(_, r)| r.start);
        }
        
        // Process ONE fragment at a time
        while let Some(fragment) = streamer.next_fragment()? {
            fragment_count += 1;
            
            // Count this fragment in peaks using binary search optimization
            if let Some(chrom_peaks) = peaks_by_chrom.get(&fragment.chrom) {
                // Use binary search to find starting point
                let start_idx = match chrom_peaks.binary_search_by_key(&fragment.start, |(_, r)| r.start) {
                    Ok(i) => i,
                    Err(i) => i.saturating_sub(1),
                };
                
                // Only check regions that could possibly overlap
                for idx in start_idx..chrom_peaks.len() {
                    let (i, peak) = &chrom_peaks[idx];
                    if fragment.start >= peak.end {
                        continue;  // Skip this one
                    }
                    if peak.start >= fragment.end {
                        break;  // No more overlaps possible
                    }
                    if fragment.end > peak.start && fragment.start < peak.end {
                        peak_counts[*i] += 1;
                        regions_with_fragments.insert(format!("peak_{}", i));
                        // IMPORTANT: Stream fragments directly to cache without storing in memory
                        // This ensures we never load all fragments at once (core efficiency principle)
                        if let Some(ref mut writer) = chr_writer {
                            writer.add_fragment(*i, true, &fragment)?;
                        }
                    }
                }
            }
            
            // Count this fragment in background regions using binary search
            if let Some(chrom_bg) = bg_by_chrom.get(&fragment.chrom) {
                let start_idx = match chrom_bg.binary_search_by_key(&fragment.start, |(_, r)| r.start) {
                    Ok(i) => i,
                    Err(i) => i.saturating_sub(1),
                };
                
                for idx in start_idx..chrom_bg.len() {
                    let (i, bg) = &chrom_bg[idx];
                    if fragment.start >= bg.end {
                        continue;
                    }
                    if bg.start >= fragment.end {
                        break;
                    }
                    if fragment.end > bg.start && fragment.start < bg.end {
                        bg_counts[*i] += 1;
                        regions_with_fragments.insert(format!("bg_{}", i));
                        // IMPORTANT: Stream background fragments to cache in real-time
                        // Avoids memory bloat from accumulating fragments
                        if let Some(ref mut writer) = chr_writer {
                            writer.add_fragment(*i, false, &fragment)?;
                        }
                    }
                }
            }
            
            // Progress updates every 10M fragments or every 15 seconds
            if fragment_count % 10000000 == 0 {
                let elapsed = start_time.elapsed();
                let rate = fragment_count as f64 / elapsed.as_secs_f64();
                let regions_count = regions_with_fragments.len();
                println!("    {} fragments processed | {} regions with fragments | {:.0} fragments/sec", 
                    fragment_count, regions_count, rate);
                last_report_time = Instant::now();
            } else if last_report_time.elapsed().as_secs() >= 15 {
                // Report every 15 seconds if no fragment milestone
                let elapsed = start_time.elapsed();
                let rate = fragment_count as f64 / elapsed.as_secs_f64();
                let regions_count = regions_with_fragments.len();
                println!("    {} fragments processed | {} regions with fragments | {:.0} fragments/sec", 
                    fragment_count, regions_count, rate);
                last_report_time = Instant::now();
            }
        }
        
        let (_total_reads, filtered_reads) = streamer.get_stats();
        let total_elapsed = start_time.elapsed();
        
        // Finish and close fragment writers if caching is enabled
        if let Some(writer) = chr_writer {
            println!("    Writing chromosome fragment indices (much faster with ~25 files)...");
            let cache_start = Instant::now();
            writer.finish_all()?;
            let cache_elapsed = cache_start.elapsed();
            println!("    Chromosome fragment cache written in {:.1}s", cache_elapsed.as_secs_f64());
        }
        
        // Report summary statistics
        let peaks_with_counts: usize = peak_counts.iter().filter(|&&c| c > 0).count();
        let bg_with_counts: usize = bg_counts.iter().filter(|&&c| c > 0).count();
        let cache_msg = if cache_dir.is_some() { " and cached fragments (chromosome-based)" } else { "" };
        println!("    Counted {} fragments in {:.1}s ({} peaks, {} background regions with counts){}", 
            filtered_reads, total_elapsed.as_secs_f64(), peaks_with_counts, bg_with_counts, cache_msg);
        
        Ok((peak_counts, bg_counts))
    }
    
    /// Count total fragments in a BAM file without loading into memory
    pub fn count_total_fragments(bam_path: &Path) -> Result<u64> {
        let mut streamer = BamStreamer::new(bam_path)?;
        let mut total = 0u64;
        
        while let Some(_) = streamer.next_fragment()? {
            total += 1;
        }
        
        Ok(total)
    }
    
    /// Extract fragments that overlap with given regions
    pub fn extract_fragments_for_regions(
        bam_path: &Path,
        regions: &[ScatrsRegion],
    ) -> Result<std::collections::HashMap<String, Vec<ScatrsFragment>>> {
        use std::collections::HashMap;
        use std::time::Instant;
        
        let mut region_fragments: HashMap<String, Vec<ScatrsFragment>> = HashMap::new();
        
        // Initialize empty vectors for each region
        for region in regions {
            let key = format!("{}:{}-{}", region.chrom, region.start, region.end);
            region_fragments.insert(key, Vec::new());
        }
        
        let start_time = Instant::now();
        let mut streamer = BamStreamer::new(bam_path)?;
        let mut fragment_count = 0u64;
        let mut extracted_count = 0u64;
        let mut last_report_time = Instant::now();
        
        println!("Extracting fragments for {} regions from {:?}...", regions.len(), bam_path);
        
        // Create chromosome-indexed maps for efficient lookup
        let mut regions_by_chrom: HashMap<String, Vec<&ScatrsRegion>> = HashMap::new();
        for region in regions {
            regions_by_chrom.entry(region.chrom.clone())
                .or_insert_with(Vec::new)
                .push(region);
        }
        
        // Sort regions within each chromosome for efficient search
        for chrom_regions in regions_by_chrom.values_mut() {
            chrom_regions.sort_by_key(|r| r.start);
        }
        
        while let Some(fragment) = streamer.next_fragment()? {
            fragment_count += 1;
            
            // Check if fragment overlaps with any region on the same chromosome
            if let Some(chrom_regions) = regions_by_chrom.get(&fragment.chrom) {
                for region in chrom_regions {
                    if fragment.end > region.start && fragment.start < region.end {
                        let key = format!("{}:{}-{}", region.chrom, region.start, region.end);
                        if let Some(fragments) = region_fragments.get_mut(&key) {
                            fragments.push(fragment.clone());
                            extracted_count += 1;
                        }
                    }
                }
            }
            
            // Report progress every 10 seconds or every 1M fragments
            let now = Instant::now();
            if now.duration_since(last_report_time).as_secs() >= 10 || fragment_count % 1_000_000 == 0 {
                let elapsed = start_time.elapsed();
                let rate = fragment_count as f64 / elapsed.as_secs_f64();
                println!("  [{:.1}s] Processed {} fragments ({:.0}/s), extracted {} fragments", 
                    elapsed.as_secs_f64(), fragment_count, rate, extracted_count);
                last_report_time = now;
            }
        }
        
        let total_elapsed = start_time.elapsed();
        let final_rate = fragment_count as f64 / total_elapsed.as_secs_f64();
        println!("  Completed in {:.1}s: {} total fragments processed ({:.0}/s), {} fragments extracted", 
            total_elapsed.as_secs_f64(), fragment_count, final_rate, extracted_count);
        
        // Report summary of extracted fragments by region
        let regions_with_fragments = region_fragments.values().filter(|v| !v.is_empty()).count();
        println!("  {} of {} regions have extracted fragments", regions_with_fragments, regions.len());
        
        Ok(region_fragments)
    }
    
    /// Extract fragments from chromosome-based cache files
    pub fn extract_fragments_from_chromosome_cache(
        cache_dir: &Path,
        peak_indices: &[usize],
        bg_indices: &[usize],
        peaks: &[ScatrsRegion],
        background: &[ScatrsRegion],
    ) -> Result<HashMap<String, Vec<ScatrsFragment>>> {
        use std::time::Instant;
        
        println!("    Reading cached fragments from chromosome files...");
        let start_time = Instant::now();
        
        // Validate inputs
        if peak_indices.is_empty() && bg_indices.is_empty() {
            return Err(anyhow!("No region indices provided for fragment extraction"));
        }
        
        // Validate that indices are within bounds
        for &idx in peak_indices {
            if idx >= peaks.len() {
                return Err(anyhow!(
                    "Peak index {} out of bounds (have {} peaks). This usually indicates a mismatch between staging and simulation.",
                    idx, peaks.len()
                ));
            }
        }
        
        for &idx in bg_indices {
            if idx >= background.len() {
                return Err(anyhow!(
                    "Background index {} out of bounds (have {} background regions). This usually indicates a mismatch between staging and simulation.",
                    idx, background.len()
                ));
            }
        }
        
        let mut region_fragments: HashMap<String, Vec<ScatrsFragment>> = HashMap::new();
        let _missing_regions: Vec<String> = Vec::new();
        
        // Group regions by chromosome for efficient reading
        println!("    Grouping {} regions by chromosome...", peak_indices.len() + bg_indices.len());
        let mut regions_by_chr: HashMap<String, Vec<(u32, String)>> = HashMap::new();
        
        // Process peak regions - use simple sequential IDs matching the staging phase
        for &idx in peak_indices {
            if idx < peaks.len() {
                let region = &peaks[idx];
                // CRITICAL: Region ID consistency - peaks use their index (0, 1, 2, ...)
                // This MUST match how IDs are assigned during staging phase
                let region_id = idx as u32;
                let region_key = format!("{}:{}-{}", region.chrom, region.start, region.end);
                regions_by_chr.entry(region.chrom.clone())
                    .or_insert_with(Vec::new)
                    .push((region_id, region_key));
            }
        }
        
        // Process background regions - use offset IDs matching the staging phase
        let peak_count = peaks.len() as u32;
        for &idx in bg_indices {
            if idx < background.len() {
                let region = &background[idx];
                // CRITICAL: Region ID consistency - background regions use peak_count + index
                // This offset scheme ensures no ID collisions between peaks and background
                let region_id = peak_count + idx as u32;
                let region_key = format!("{}:{}-{}", region.chrom, region.start, region.end);
                regions_by_chr.entry(region.chrom.clone())
                    .or_insert_with(Vec::new)
                    .push((region_id, region_key));
            }
        }
        
        // Read fragments for each chromosome sequentially (parallel causes file I/O deadlock)
        println!("    Processing {} chromosomes sequentially...", regions_by_chr.len());
        let mut chr_results = Vec::new();
        let mut reader_manager = ChromosomeReaderManager::new(cache_dir);
        
        let total_chromosomes = regions_by_chr.len();
        let mut processed_chromosomes = 0;
        
        for (chromosome, regions) in regions_by_chr {
            processed_chromosomes += 1;
            println!("    [{}/{}] Processing chromosome {} ({} regions)...", 
                processed_chromosomes, total_chromosomes, chromosome, regions.len());
            
            let mut chr_fragments = HashMap::new();
            
            match reader_manager.get_reader(&chromosome) {
                Ok(reader) => {
                    let mut processed_regions = 0;
                    let mut missing_regions = 0;
                    for (region_id, region_key) in regions {
                        match reader.get_fragments_for_region(region_id) {
                            Ok(fragments) => {
                                chr_fragments.insert(region_key.clone(), fragments);
                                processed_regions += 1;
                            }
                            Err(_) => {
                                // This region has no fragments - this is expected for many regions
                                // Just skip it silently
                                missing_regions += 1;
                            }
                        }
                    }
                    if processed_regions > 0 || missing_regions > 0 {
                        println!("      Loaded fragments from {} regions on {} ({} regions had no fragments)", 
                            processed_regions, chromosome, missing_regions);
                    }
                }
                Err(e) => {
                    eprintln!("Warning: Could not read cache for chromosome {}: {}", chromosome, e);
                }
            }
            chr_results.push(chr_fragments);
        }
        
        // Merge results from all chromosomes
        for chr_fragments in chr_results {
            region_fragments.extend(chr_fragments);
        }
        
        // No need to report missing regions - they simply have no fragments which is expected
        
        let total_fragments: usize = region_fragments.values().map(|v| v.len()).sum();
        let elapsed = start_time.elapsed();
        println!("    Loaded {} cached fragments from {} regions in {:.1}s ({:.0} fragments/sec)", 
            total_fragments, region_fragments.len(), elapsed.as_secs_f64(), 
            total_fragments as f64 / elapsed.as_secs_f64());
        
        if region_fragments.is_empty() {
            return Err(anyhow!(
                "No fragments could be loaded from cache. This usually means the staged data is incomplete or corrupted. Please re-run the stage command."
            ));
        }
        
        Ok(region_fragments)
    }
    
    /// Check if chromosome-based fragment cache exists
    pub fn has_chromosome_cache(sample_dir: &Path) -> bool {
        // Check if any chromosome fragment files exist
        if let Ok(entries) = std::fs::read_dir(sample_dir) {
            for entry in entries.flatten() {
                if let Some(name) = entry.file_name().to_str() {
                    if name.ends_with("_fragments.bin") && name != "fragments.bin" {
                        return true;
                    }
                }
            }
        }
        false
    }
    
    /// Extract all fragments from a BAM file (for no-peaks mode)
    pub fn extract_all_fragments(
        bam_path: &Path,
        max_fragments: Option<usize>,
    ) -> Result<Vec<ScatrsFragment>> {
        use std::time::Instant;
        
        let mut streamer = BamStreamer::new(bam_path)?;
        let mut fragments = Vec::new();
        let mut fragment_count = 0usize;
        let start_time = Instant::now();
        let mut last_report_time = Instant::now();
        
        println!("Extracting all fragments from {:?}...", bam_path);
        if let Some(max) = max_fragments {
            println!("  Maximum fragment limit: {}", max);
        }
        
        while let Some(fragment) = streamer.next_fragment()? {
            fragments.push(fragment);
            fragment_count += 1;
            
            if let Some(max) = max_fragments {
                if fragment_count >= max {
                    println!("  Reached maximum fragment limit: {}", max);
                    break;
                }
            }
            
            // Report progress every 10 seconds or every 1M fragments
            let now = Instant::now();
            if now.duration_since(last_report_time).as_secs() >= 10 || fragment_count % 1_000_000 == 0 {
                let elapsed = start_time.elapsed();
                let rate = fragment_count as f64 / elapsed.as_secs_f64();
                println!("  [{:.1}s] Loaded {} fragments ({:.0}/s)", 
                    elapsed.as_secs_f64(), fragment_count, rate);
                last_report_time = now;
            }
        }
        
        let total_elapsed = start_time.elapsed();
        let final_rate = fragment_count as f64 / total_elapsed.as_secs_f64();
        println!("  Completed in {:.1}s: {} total fragments loaded ({:.0}/s)", 
            total_elapsed.as_secs_f64(), fragment_count, final_rate);
        
        Ok(fragments)
    }
    
    /// Stream through BAM file with optional early termination for testing
    pub fn count_fragments_streaming_with_limit(
        bam_path: &Path,
        peaks: &[ScatrsRegion],
        background: &[ScatrsRegion],
        max_fragments: Option<usize>,
        timeout_secs: Option<u64>,
    ) -> Result<(Vec<u32>, Vec<u32>)> {
        use std::time::{Duration, Instant};
        
        let mut peak_counts = vec![0u32; peaks.len()];
        let mut bg_counts = vec![0u32; background.len()];
        
        println!("Starting fragment counting (with limits) across {} peaks and {} background regions...", 
            peaks.len(), background.len());
        
        if let Some(max) = max_fragments {
            println!("Will process up to {} fragments", max);
        }
        if let Some(timeout) = timeout_secs {
            println!("Will timeout after {} seconds", timeout);
        }
        
        let start_time = Instant::now();
        let timeout_duration = timeout_secs.map(Duration::from_secs);
        let mut streamer = BamStreamer::new(bam_path)?;
        let mut fragment_count = 0u64;
        
        // Create chromosome-indexed maps (same as before)
        let mut peaks_by_chrom: HashMap<String, Vec<(usize, &ScatrsRegion)>> = HashMap::new();
        for (i, peak) in peaks.iter().enumerate() {
            peaks_by_chrom.entry(peak.chrom.clone())
                .or_insert_with(Vec::new)
                .push((i, peak));
        }
        
        let mut bg_by_chrom: HashMap<String, Vec<(usize, &ScatrsRegion)>> = HashMap::new();
        for (i, bg) in background.iter().enumerate() {
            bg_by_chrom.entry(bg.chrom.clone())
                .or_insert_with(Vec::new)
                .push((i, bg));
        }
        
        // Sort for binary search
        for regions in peaks_by_chrom.values_mut() {
            regions.sort_by_key(|(_, r)| r.start);
        }
        for regions in bg_by_chrom.values_mut() {
            regions.sort_by_key(|(_, r)| r.start);
        }
        
        // Process fragments with early termination checks
        while let Some(fragment) = streamer.next_fragment()? {
            fragment_count += 1;
            
            // Check for early termination conditions
            if let Some(max) = max_fragments {
                if fragment_count > max as u64 {
                    println!("Reached fragment limit of {}, stopping early", max);
                    break;
                }
            }
            
            if let Some(timeout) = timeout_duration {
                if start_time.elapsed() > timeout {
                    eprintln!("Warning: Fragment counting timed out after {} seconds", timeout_secs.unwrap());
                    eprintln!("Processed {} fragments before timeout", fragment_count);
                    break;
                }
            }
            
            // Count fragments (using binary search optimization from above)
            if let Some(chrom_peaks) = peaks_by_chrom.get(&fragment.chrom) {
                let start_idx = match chrom_peaks.binary_search_by_key(&fragment.start, |(_, r)| r.start) {
                    Ok(i) => i,
                    Err(i) => i.saturating_sub(1),
                };
                
                for idx in start_idx..chrom_peaks.len() {
                    let (i, peak) = &chrom_peaks[idx];
                    if peak.start >= fragment.end {
                        break;
                    }
                    if fragment.end > peak.start && fragment.start < peak.end {
                        peak_counts[*i] += 1;
                    }
                }
            }
            
            if let Some(chrom_bg) = bg_by_chrom.get(&fragment.chrom) {
                let start_idx = match chrom_bg.binary_search_by_key(&fragment.start, |(_, r)| r.start) {
                    Ok(i) => i,
                    Err(i) => i.saturating_sub(1),
                };
                
                for idx in start_idx..chrom_bg.len() {
                    let (i, bg) = &chrom_bg[idx];
                    if bg.start >= fragment.end {
                        break;
                    }
                    if fragment.end > bg.start && fragment.start < bg.end {
                        bg_counts[*i] += 1;
                    }
                }
            }
            
            // Progress reporting - every 100M fragments
            if fragment_count % 100000000 == 0 {
                let elapsed = start_time.elapsed();
                let rate = fragment_count as f64 / elapsed.as_secs_f64();
                println!("    {} fragments processed ({:.0} fragments/sec)", fragment_count, rate);
            }
        }
        
        println!("  Completed: {} fragments processed", fragment_count);
        Ok((peak_counts, bg_counts))
    }
    
    /// Sample fragments from BAM file using reservoir sampling during streaming
    /// Never loads all fragments into memory
    pub fn sample_fragments_streaming(
        bam_path: &Path,
        sample_size: usize,
        seed: Option<u64>,
    ) -> Result<Vec<ScatrsFragment>> {
        use crate::sampling::a_expj::AExpJSampler;
        
        let mut sampler = AExpJSampler::new(sample_size, seed);
        let mut streamer = BamStreamer::new(bam_path)?;
        
        println!("Sampling {} fragments using reservoir sampling...", sample_size);
        
        // Stream and sample - this is TRUE streaming with reservoir sampling
        while let Some(fragment) = streamer.next_fragment()? {
            // The reservoir sampler handles the sampling logic internally
            // It keeps only sample_size fragments in memory at any time
            sampler.add_weighted(fragment, 1.0);
        }
        
        let sampled = sampler.get_reservoir();
        println!("Sampled {} fragments from stream", sampled.len());
        
        Ok(sampled)
    }
}

// ============================================================================
// Peak Merging
// ============================================================================

pub struct PeakMerger;

impl PeakMerger {
    pub fn merge_peaks(
        peaks: Vec<Vec<Peak>>,
        merge_distance: u32,
    ) -> Result<Vec<ScatrsRegion>> {
        // Collect all regions and sort by chromosome and start
        let mut all_regions: Vec<ScatrsRegion> = peaks
            .into_iter()
            .flatten()
            .map(|p| p.region)
            .collect();
        
        println!("Total peaks before merging: {}", all_regions.len());
        
        // Sort regions for efficient merging
        all_regions.sort();
        
        let pb = ProgressBar::new(all_regions.len() as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} {msg}")
                .unwrap()
        );
        pb.set_message("Merging peaks");
        
        // Merge overlapping or nearby peaks
        let mut merged = Vec::new();
        let mut current: Option<ScatrsRegion> = None;
        
        for region in all_regions {
            pb.inc(1);
            
            match current {
                None => current = Some(region),
                Some(ref mut curr) => {
                    if curr.chrom == region.chrom && 
                       region.start <= curr.end + merge_distance as u64 {
                        // Extend current region
                        curr.end = curr.end.max(region.end);
                        // Keep the higher score if available
                        if let (Some(s1), Some(s2)) = (curr.score, region.score) {
                            curr.score = Some(s1.max(s2));
                        }
                    } else {
                        // Save current and start new
                        merged.push(current.take().unwrap());
                        current = Some(region);
                    }
                }
            }
        }
        
        // Don't forget the last region
        if let Some(curr) = current {
            merged.push(curr);
        }
        
        pb.finish_with_message("Peak merging complete");

        // Filter for standard chromosomes
        let valid_chroms: Vec<String> = VALID_CHROMOSOMES
            .iter()
            .map(|&c| c.to_string())
            .collect();

        let filtered: Vec<ScatrsRegion> = merged
            .into_iter()
            .filter(|r| valid_chroms.contains(&r.chrom))
            .collect();

        println!("Merged peaks: {} (filtered to standard chromosomes)", filtered.len());

        // Deduplicate merged peaks
        println!("Deduplicating merged peaks...");
        let deduplicated = Self::deduplicate_regions(filtered);
        println!("After deduplication: {} unique peaks", deduplicated.len());

        Ok(deduplicated)
    }
    
    pub fn extend_peaks(regions: &[ScatrsRegion], extend_bp: u32) -> Vec<ScatrsRegion> {
        regions.iter().map(|r| r.extend(extend_bp as u64)).collect()
    }
    
    pub fn deduplicate_regions(mut regions: Vec<ScatrsRegion>) -> Vec<ScatrsRegion> {
        use std::collections::HashSet;

        // Sort for efficient comparison
        regions.sort();

        // Track seen regions by their genomic coordinates
        let mut seen = HashSet::new();
        let mut unique_regions = Vec::new();

        for region in regions {
            let key = format!("{}:{}-{}", region.chrom, region.start, region.end);
            if seen.insert(key) {
                unique_regions.push(region);
            }
        }

        unique_regions
    }

    pub fn diagnose_peak_duplicates(peaks: &[ScatrsRegion]) -> (usize, usize) {
        use std::collections::HashSet;
        let mut seen = HashSet::new();
        let total = peaks.len();
        let unique = peaks.iter().filter(|p| {
            let key = format!("{}:{}-{}", p.chrom, p.start, p.end);
            seen.insert(key)
        }).count();
        (total, unique)
    }
}

// ============================================================================
// Background Generation
// ============================================================================

pub struct BackgroundGenerator;

impl BackgroundGenerator {
    pub fn deduplicate_regions(mut regions: Vec<ScatrsRegion>) -> Vec<ScatrsRegion> {
        use std::collections::HashSet;

        // Sort for efficient comparison
        regions.sort();

        // Track seen regions by their genomic coordinates
        let mut seen = HashSet::new();
        let mut unique_regions = Vec::new();

        for region in regions {
            let key = format!("{}:{}-{}", region.chrom, region.start, region.end);
            if seen.insert(key) {
                unique_regions.push(region);
            }
        }

        unique_regions
    }

    pub fn generate_background(
        chrom_sizes: &HashMap<String, u64>,
        extended_peaks: &[ScatrsRegion],
        blacklist: &[ScatrsRegion],
        bin_size: u32,
    ) -> Result<Vec<ScatrsRegion>> {
        use rayon::prelude::*;
        
        let total_bins_estimate = chrom_sizes.values().sum::<u64>() / bin_size as u64;
        let pb = ProgressBar::new(total_bins_estimate);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} {msg}")
                .unwrap()
        );
        pb.set_message("Generating background regions (parallel)");
        
        // Process chromosomes in parallel
        let all_background: Vec<Vec<ScatrsRegion>> = VALID_CHROMOSOMES
            .par_iter()
            .map(|&chrom| {
                let chrom_str = chrom.to_string();
                
                // Skip if chromosome not in sizes
                let size = match chrom_sizes.get(&chrom_str) {
                    Some(s) => *s,
                    None => return vec![],
                };
                
                // Get peaks and blacklist for this chromosome
                let mut chrom_peaks: Vec<_> = extended_peaks
                    .iter()
                    .filter(|r| r.chrom == chrom_str)
                    .collect();
                chrom_peaks.sort_by_key(|r| r.start);
                
                let mut chrom_blacklist: Vec<_> = blacklist
                    .iter()
                    .filter(|r| r.chrom == chrom_str)
                    .collect();
                chrom_blacklist.sort_by_key(|r| r.start);
                
                // Generate bins for this chromosome
                let mut chrom_background = Vec::new();
                for start in (0..size).step_by(bin_size as usize) {
                    pb.inc(1);
                    
                    let end = std::cmp::min(start + bin_size as u64, size);
                    let bin = ScatrsRegion::new(chrom_str.clone(), start, end)
                        .with_name(format!("bg_{}_{}", chrom, start));
                    
                    // Check if bin overlaps with peaks using binary search for efficiency
                    let overlaps_peak = Self::has_overlap(&bin, &chrom_peaks);
                    let overlaps_blacklist = Self::has_overlap(&bin, &chrom_blacklist);
                    
                    if !overlaps_peak && !overlaps_blacklist {
                        chrom_background.push(bin);
                    }
                }
                
                chrom_background
            })
            .collect();
        
        // Flatten the results
        let background_regions: Vec<ScatrsRegion> = all_background
            .into_iter()
            .flatten()
            .collect();
        
        pb.finish_with_message(format!("Generated {} background regions", background_regions.len()));
        
        Ok(background_regions)
    }
    
    fn has_overlap(query: &ScatrsRegion, sorted_regions: &[&ScatrsRegion]) -> bool {
        // Binary search for potential overlaps
        let idx = sorted_regions.binary_search_by_key(&query.start, |r| r.start);
        
        let start_idx = match idx {
            Ok(i) => i,
            Err(i) => i.saturating_sub(1),
        };
        
        // Check nearby regions for overlap
        for &region in sorted_regions.iter().skip(start_idx).take(5) {
            if region.start >= query.end {
                break;
            }
            if query.overlaps(region) {
                return true;
            }
        }
        
        false
    }
    
    pub fn read_blacklist(blacklist_file: &Path) -> Result<Vec<ScatrsRegion>> {
        let file = File::open(blacklist_file)
            .with_context(|| format!("Failed to open blacklist file: {:?}", blacklist_file))?;
        let reader = BufReader::new(file);
        let mut blacklist = Vec::new();
        
        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 3 {
                continue;
            }
            
            let region = ScatrsRegion::new(
                fields[0].to_string(),
                fields[1].parse()?,
                fields[2].parse()?,
            );
            
            blacklist.push(region);
        }
        
        Ok(blacklist)
    }
    
    pub fn read_chrom_sizes(chrom_sizes_file: &Path) -> Result<HashMap<String, u64>> {
        let file = File::open(chrom_sizes_file)
            .with_context(|| format!("Failed to open chromosome sizes file: {:?}", chrom_sizes_file))?;
        let reader = BufReader::new(file);
        let mut sizes = HashMap::new();
        
        for line in reader.lines() {
            let line = line?;
            if line.is_empty() {
                continue;
            }
            
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 2 {
                continue;
            }
            
            sizes.insert(
                fields[0].to_string(),
                fields[1].parse()?,
            );
        }
        
        Ok(sizes)
    }
}
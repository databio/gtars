use crate::scatrs::models::{ScatrsFragment, ScatrsRegion, Peak};
use crate::scatrs::consts::VALID_CHROMOSOMES;
use noodles::bam;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::Record;
use std::path::Path;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::collections::HashMap;
use anyhow::{Result, Context};
use indicatif::{ProgressBar, ProgressStyle};

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
                    
                    // Extract fragment information
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
                    
                    let start = match record.alignment_start() {
                        Some(pos) => match pos {
                            Ok(p) => p.get() as u64 - 1, // Convert to 0-based
                            _ => continue,
                        },
                        _ => continue,
                    };
                    
                    // Calculate end position using start + span
                    let span = match record.alignment_span() {
                        Some(s) => match s {
                            Ok(sp) => sp as u64,
                            _ => continue,
                        },
                        _ => continue,
                    };
                    let end = start + span;
                    
                    let strand = if flags.is_reverse_complemented() {
                        Some('-')
                    } else {
                        Some('+')
                    };
                    
                    let quality = record.mapping_quality()
                        .map(|q| q.get())
                        .unwrap_or(0);
                    
                    let fragment = ScatrsFragment::new(chrom, start, end)
                        .with_strand(strand.unwrap())
                        .with_quality(quality);
                    
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
        
        Ok(filtered)
    }
    
    pub fn extend_peaks(regions: &[ScatrsRegion], extend_bp: u32) -> Vec<ScatrsRegion> {
        regions.iter().map(|r| r.extend(extend_bp as u64)).collect()
    }
    
    pub fn deduplicate_regions(mut regions: Vec<ScatrsRegion>) -> Vec<ScatrsRegion> {
        regions.sort();
        regions.dedup();
        regions
    }
}

// ============================================================================
// Background Generation
// ============================================================================

pub struct BackgroundGenerator;

impl BackgroundGenerator {
    pub fn generate_background(
        chrom_sizes: &HashMap<String, u64>,
        extended_peaks: &[ScatrsRegion],
        blacklist: &[ScatrsRegion],
        bin_size: u32,
    ) -> Result<Vec<ScatrsRegion>> {
        let mut background_regions = Vec::new();
        
        let total_bins_estimate = chrom_sizes.values().sum::<u64>() / bin_size as u64;
        let pb = ProgressBar::new(total_bins_estimate);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} {msg}")
                .unwrap()
        );
        pb.set_message("Generating background regions");
        
        // Process each chromosome
        for chrom in VALID_CHROMOSOMES {
            let chrom_str = chrom.to_string();
            
            // Skip if chromosome not in sizes
            let size = match chrom_sizes.get(&chrom_str) {
                Some(s) => *s,
                None => continue,
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
            
            // Generate bins and check for overlaps
            for start in (0..size).step_by(bin_size as usize) {
                pb.inc(1);
                
                let end = std::cmp::min(start + bin_size as u64, size);
                let bin = ScatrsRegion::new(chrom_str.clone(), start, end)
                    .with_name(format!("bg_{}_{}", chrom, start));
                
                // Check if bin overlaps with peaks using binary search for efficiency
                let overlaps_peak = Self::has_overlap(&bin, &chrom_peaks);
                let overlaps_blacklist = Self::has_overlap(&bin, &chrom_blacklist);
                
                if !overlaps_peak && !overlaps_blacklist {
                    background_regions.push(bin);
                }
            }
        }
        
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
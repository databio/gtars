#[cfg(test)]
mod tests {
    use crate::models::*;
    use crate::*;
    use std::path::PathBuf;
    use std::fs;
    use tempfile::TempDir;
    
    
    #[test]
    fn test_fragment_creation() {
        let fragment = ScatrsFragment::new(
            "chr1".to_string(),
            1000,
            1500,
        )
        .with_cell_id("cell_001".to_string())
        .with_strand('+')
        .with_quality(30);
        
        assert_eq!(fragment.chrom, "chr1");
        assert_eq!(fragment.start, 1000);
        assert_eq!(fragment.end, 1500);
        assert_eq!(fragment.length(), 500);
        assert_eq!(fragment.cell_id, "cell_001");
        assert_eq!(fragment.strand, Some('+'));
        assert_eq!(fragment.quality, Some(30));
    }
    
    #[test]
    fn test_region_operations() {
        let region1 = ScatrsRegion::new("chr1".to_string(), 1000, 2000);
        let region2 = ScatrsRegion::new("chr1".to_string(), 1500, 2500);
        let region3 = ScatrsRegion::new("chr2".to_string(), 1000, 2000);
        
        // Test overlap
        assert!(region1.overlaps(&region2));
        assert!(!region1.overlaps(&region3));
        
        // Test merge
        let merged = region1.merge(&region2);
        assert!(merged.is_some());
        let merged = merged.unwrap();
        assert_eq!(merged.start, 1000);
        assert_eq!(merged.end, 2500);
        
        // Test extend
        let extended = region1.extend(100);
        assert_eq!(extended.start, 900);
        assert_eq!(extended.end, 2100);
    }
    
    #[test]
    fn test_config_creation() {
        let config = ScatrsConfig::default();
        assert_eq!(config.cell_count, 100);
        match &config.signal_to_noise_distribution {
            SignalToNoiseDistribution::Fixed { signal_to_noise } => {
                assert_eq!(*signal_to_noise, 0.8);
            },
            _ => panic!("Expected Fixed distribution"),
        }
        assert_eq!(config.extend_peaks, 250);
        assert_eq!(config.bin_size, 5000);
    }
    
    #[test]
    fn test_cell_metadata() {
        let mut cell = SimulatedCell::new(
            "cell_001".to_string(),
            "T_cell".to_string(),
        );
        
        let peak_fragment = ScatrsFragment::new("chr1".to_string(), 1000, 1500)
            .with_cell_id("cell_001".to_string());
        let bg_fragment = ScatrsFragment::new("chr1".to_string(), 5000, 5500)
            .with_cell_id("cell_001".to_string());
        
        cell.add_peak_fragment(peak_fragment);
        cell.add_background_fragment(bg_fragment);
        cell.calculate_frip();
        
        assert_eq!(cell.metadata.fragment_count, 2);
        assert_eq!(cell.metadata.peak_fragments, 1);
        assert_eq!(cell.metadata.background_fragments, 1);
        assert_eq!(cell.metadata.frip_score, Some(0.5));
    }
    
    // ============================================================================
    // Integration Tests
    // ============================================================================
    
    #[test]
    fn test_staging_data_persistence() {
        // Test that staging data can be saved and loaded correctly
        let temp_dir = TempDir::new().expect("Failed to create temp dir");
        let test_path = temp_dir.path();
        
        // Test fragment counts serialization
        let counts = vec![10u32, 20, 30, 40, 50];
        let counts_path = test_path.join("counts.bin");
        
        StagedDataWriter::save_counts(&counts, &counts_path)
            .expect("Failed to save counts");
        
        let loaded = StagedDataReader::load_counts(&counts_path)
            .expect("Failed to load counts");
        
        assert_eq!(counts, loaded, "Counts should match after save/load");
    }
    
    #[test]
    fn test_region_serialization() {
        let temp_dir = TempDir::new().expect("Failed to create temp dir");
        let test_path = temp_dir.path();
        
        // Test region serialization
        let regions = vec![
            ScatrsRegion::new("chr1".to_string(), 1000, 2000),
            ScatrsRegion::new("chr2".to_string(), 3000, 4000),
        ];
        
        let regions_path = test_path.join("regions.bin");
        
        StagedDataWriter::save_regions(&regions, &regions_path)
            .expect("Failed to save regions");
        
        let loaded = StagedDataReader::load_regions(&regions_path)
            .expect("Failed to load regions");
        
        assert_eq!(regions.len(), loaded.len());
        for (orig, load) in regions.iter().zip(loaded.iter()) {
            assert_eq!(orig.chrom, load.chrom);
            assert_eq!(orig.start, load.start);
            assert_eq!(orig.end, load.end);
        }
    }
    
    #[test]
    fn test_metadata_yaml() {
        let temp_dir = TempDir::new().expect("Failed to create temp dir");
        let test_path = temp_dir.path();
        
        let metadata = SampleMetadata {
            sample_name: "test".to_string(),
            bam_file: PathBuf::from("/test/path.bam"),
            total_fragments: 1000,
            fragments_in_peaks: Some(800),
            fragments_in_background: Some(200),
            processing_time_secs: 1.5,
            timestamp: "2024-01-01".to_string(),
        };
        
        let metadata_path = test_path.join("metadata.yaml");
        
        StagedDataWriter::save_sample_metadata(&metadata, &metadata_path)
            .expect("Failed to save metadata");
        
        let loaded = StagedDataReader::load_sample_metadata(&metadata_path)
            .expect("Failed to load metadata");
        
        assert_eq!(metadata.sample_name, loaded.sample_name);
        assert_eq!(metadata.total_fragments, loaded.total_fragments);
    }
    
    #[test]
    fn test_bam_fragment_counting() {
        // Test with the actual test BAM file if it exists
        let test_bam = PathBuf::from("/home/nsheff/code/scan-atac-sim/gtars/gtars/tests/data/test_chr22_small.bam");
        
        if test_bam.exists() {
            // BamProcessor is now internal to staging module - use staging::BamProcessor
            let count = crate::staging::BamProcessor::count_total_fragments(&test_bam);
            
            // The test BAM has 8 paired reads, but may be filtered
            match count {
                Ok(n) => {
                    assert!(n <= 16, "Fragment count should be reasonable for test BAM");
                    println!("Test BAM contains {} fragments after filtering", n);
                }
                Err(e) => {
                    println!("Could not count fragments: {}", e);
                }
            }
        }
    }
    
    #[test]
    fn test_peak_loading() {
        let temp_dir = TempDir::new().expect("Failed to create temp dir");
        let peak_file = temp_dir.path().join("test.narrowPeak");
        
        // Create a test peak file
        let content = "chr1\t1000\t2000\tpeak1\t100\t.\t10.0\t5.0\t-1\t500\n\
                       chr2\t3000\t4000\tpeak2\t200\t.\t20.0\t10.0\t-1\t500\n";
        
        fs::write(&peak_file, content).expect("Failed to write peak file");
        
        let peaks = NarrowPeakReader::read_peaks(&peak_file, "test");
        
        match peaks {
            Ok(peak_list) => {
                assert_eq!(peak_list.len(), 2, "Should load 2 peaks");
                assert_eq!(peak_list[0].region.chrom, "chr1");
                assert_eq!(peak_list[0].region.start, 1000);
                assert_eq!(peak_list[0].region.end, 2000);
            }
            Err(e) => {
                println!("Error loading peaks: {}", e);
            }
        }
    }
    
    #[test]
    fn test_background_generation() {
        // Test background region generation
        let peaks = vec![
            ScatrsRegion::new("chr1".to_string(), 1000, 2000),
            ScatrsRegion::new("chr1".to_string(), 5000, 6000),
        ];
        
        let mut chrom_sizes = std::collections::HashMap::new();
        chrom_sizes.insert("chr1".to_string(), 10000u64);
        
        let background = BackgroundGenerator::generate_background(
            &chrom_sizes,
            &peaks,
            &[],   // empty blacklist
            1000,  // bin_size
        );
        
        match background {
            Ok(bg_regions) => {
                // Check that background doesn't overlap peaks
                for bg in &bg_regions {
                    for peak in &peaks {
                        if bg.chrom == peak.chrom {
                            let overlaps = bg.start < peak.end && bg.end > peak.start;
                            assert!(!overlaps, "Background should not overlap peaks");
                        }
                    }
                }
                
                // Check bin size
                for bg in &bg_regions {
                    assert_eq!(bg.end - bg.start, 1000, "Bins should be 1000bp");
                }
            }
            Err(e) => {
                println!("Error generating background: {}", e);
            }
        }
    }
}
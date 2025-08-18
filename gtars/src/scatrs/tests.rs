#[cfg(test)]
mod tests {
    use crate::scatrs::models::*;
    
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
        assert_eq!(config.fragments_per_cell, 10000);
        assert_eq!(config.signal_to_noise, 0.8);
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
}
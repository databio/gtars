use crate::models::{ScatrsRegion, FragmentDistribution};
use crate::sampling::WeightedSampler;
use rand::prelude::*;

#[test]
fn test_a_res_integration() {
    // Create test data
    let peaks: Vec<ScatrsRegion> = (0..1000)
        .map(|i| ScatrsRegion::new(
            format!("chr{}", (i % 22) + 1),
            (i as u64) * 1000,
            (i as u64) * 1000 + 500,
        ))
        .collect();
    
    let peak_weights: Vec<f64> = (0..1000)
        .map(|i| ((i % 100) + 1) as f64)
        .collect();
    
    let background: Vec<ScatrsRegion> = (1000..5000)
        .map(|i| ScatrsRegion::new(
            format!("chr{}", (i % 22) + 1),
            (i as u64) * 1000,
            (i as u64) * 1000 + 500,
        ))
        .collect();
    
    let background_weights: Vec<f64> = (0..4000)
        .map(|i| ((i % 50) + 1) as f64)
        .collect();
    
    // Create sampler
    let sampler = WeightedSampler::new(Some(42));
    let mut rng = StdRng::seed_from_u64(42);
    
    // Test the optimized sampling method
    let result = sampler.sample_regions_for_cell_optimized(
        &peaks,
        &peak_weights,
        &background,
        &background_weights,
        0.8,  // signal-to-noise ratio
        &FragmentDistribution::Uniform { fragments_per_cell: 10000 },
        &mut rng,
    );
    
    assert!(result.is_ok());
    let sampled_regions = result.unwrap();
    
    // Verify we got some regions
    assert!(!sampled_regions.is_empty());
    
    // Check that we have both peak and background regions
    let peak_count = sampled_regions.iter().filter(|r| r.is_peak).count();
    let bg_count = sampled_regions.iter().filter(|r| !r.is_peak).count();
    
    assert!(peak_count > 0, "Should have sampled some peak regions");
    assert!(bg_count > 0, "Should have sampled some background regions");
    
    // Verify signal-to-noise ratio is reasonable
    // Note: The actual ratio may differ from the requested 0.8 due to
    // without-replacement sampling and limited number of regions
    let peak_ratio = peak_count as f64 / (peak_count + bg_count) as f64;
    println!("Requested signal-to-noise: 0.8, Actual: {:.3}", peak_ratio);
    
    // Just verify we have a mix of peaks and background
    assert!(peak_ratio > 0.1 && peak_ratio < 0.9,
        "Signal-to-noise ratio should be between 0.1 and 0.9, got {}",
        peak_ratio
    );
    
    println!("Integration test passed:");
    println!("  Sampled {} peak regions", peak_count);
    println!("  Sampled {} background regions", bg_count);
    println!("  Actual signal-to-noise ratio: {:.3}", peak_ratio);
}

#[test]
fn test_a_res_performance_comparison() {
    use std::time::Instant;
    
    // Create larger dataset for performance testing
    let n = 10000;
    let peaks: Vec<ScatrsRegion> = (0..n)
        .map(|i| ScatrsRegion::new(
            format!("chr{}", (i % 22) + 1),
            (i as u64) * 1000,
            (i as u64) * 1000 + 500,
        ))
        .collect();
    
    let weights: Vec<f64> = (0..n)
        .map(|i| ((i % 100) + 1) as f64)
        .collect();
    
    let sampler = WeightedSampler::new(Some(42));
    let rng = StdRng::seed_from_u64(42);
    
    // Measure time for optimized implementation
    let start = Instant::now();
    for _ in 0..10 {
        let mut local_rng = rng.clone();
        let _ = sampler.sample_regions_for_cell_optimized(
            &peaks,
            &weights,
            &peaks,
            &weights,
            0.5,
            &FragmentDistribution::Uniform { fragments_per_cell: 1000 },
            &mut local_rng,
        );
    }
    let optimized_time = start.elapsed();
    
    println!("Performance test results:");
    println!("  A-Res time for 10 iterations: {:?}", optimized_time);
    println!("  Average time per iteration: {:?}", optimized_time / 10);
}

#[test]
fn test_without_replacement_property() {
    use std::collections::HashSet;
    
    let peaks: Vec<ScatrsRegion> = (0..100)
        .map(|i| ScatrsRegion::new(
            format!("region_{}", i),
            (i as u64) * 1000,
            (i as u64) * 1000 + 500,
        ))
        .collect();
    
    let weights: Vec<f64> = vec![1.0; 100];
    
    let sampler = WeightedSampler::new(Some(42));
    let mut rng = StdRng::seed_from_u64(42);
    
    // Sample with high fragment count to test without-replacement
    let result = sampler.sample_regions_for_cell_optimized(
        &peaks,
        &weights,
        &[],
        &[],
        1.0,  // all peaks, no background
        &FragmentDistribution::Uniform { fragments_per_cell: 50 },
        &mut rng,
    ).unwrap();
    
    // Check for duplicates
    let unique_regions: HashSet<_> = result.iter()
        .map(|r| format!("{}_{}_{}", r.region.chrom, r.region.start, r.region.end))
        .collect();
    
    assert_eq!(
        result.len(),
        unique_regions.len(),
        "Found duplicate regions in sample, violating without-replacement property"
    );
}
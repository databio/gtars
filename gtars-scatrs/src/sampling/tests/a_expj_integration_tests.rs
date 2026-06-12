use crate::models::ScatrsRegion;
use crate::sampling::a_expj::{AExpJSampler, AExpJStrategy};
use std::collections::HashSet;

#[test]
fn test_a_expj_integration_basic() {
    // Create test data
    let regions: Vec<ScatrsRegion> = (0..1000)
        .map(|i| ScatrsRegion::new(
            format!("chr{}", (i % 22) + 1),
            (i as u64) * 1000,
            (i as u64) * 1000 + 500,
        ))
        .collect();
    
    let weights: Vec<f64> = (0..1000)
        .map(|i| ((i % 100) + 1) as f64)
        .collect();
    
    // Test A-ExpJ sampling
    let mut sampler = AExpJSampler::new(100, Some(42));
    let result = sampler.sample_stream(
        regions.iter().zip(weights.iter()),
        |(_, &weight)| weight,
    );
    
    assert_eq!(result.len(), 100);
    
    // Verify no duplicates
    let unique_regions: HashSet<_> = result.iter()
        .map(|(r, _)| format!("{}_{}_{}", r.chrom, r.start, r.end))
        .collect();
    assert_eq!(result.len(), unique_regions.len());
    
    // Check statistics
    let stats = sampler.get_statistics();
    println!("A-ExpJ Integration Test Stats:");
    println!("  Items processed: {}", stats.items_processed);
    println!("  Items skipped: {}", stats.items_skipped);
    println!("  Jump efficiency: {:.2}%", stats.jump_efficiency * 100.0);
}

#[test]
fn test_a_expj_strategies() {
    let n = 5000;
    let k = 50;
    let regions: Vec<ScatrsRegion> = (0..n)
        .map(|i| ScatrsRegion::new(
            format!("chr{}", (i % 22) + 1),
            (i as u64) * 1000,
            (i as u64) * 1000 + 500,
        ))
        .collect();
    
    let weights: Vec<f64> = vec![1.0; n];
    
    for strategy in [
        AExpJStrategy::Standard,
        AExpJStrategy::Conservative,
        AExpJStrategy::Aggressive,
        AExpJStrategy::Disabled,
    ] {
        let mut sampler = AExpJSampler::new(k, Some(42))
            .with_strategy(strategy);
        
        let result = sampler.sample_stream(
            regions.iter().zip(weights.iter()),
            |(_, &weight)| weight,
        );
        
        assert_eq!(result.len(), k, "Strategy {:?} should sample k items", strategy);
        
        let stats = sampler.get_statistics();
        println!("Strategy {:?}: efficiency={:.2}%", 
                 strategy, stats.jump_efficiency * 100.0);
        
        if matches!(strategy, AExpJStrategy::Disabled) {
            assert_eq!(stats.items_skipped, 0, "Disabled strategy should skip no items");
        }
    }
}

#[test]
fn test_a_expj_large_scale_performance() {
    use std::time::Instant;
    
    // Large-scale test
    let n = 100_000;
    let k = 1_000; // 1% sample ratio
    
    let regions: Vec<ScatrsRegion> = (0..n)
        .map(|i| ScatrsRegion::new(
            format!("chr{}", (i % 22) + 1),
            (i as u64) * 1000,
            (i as u64) * 1000 + 500,
        ))
        .collect();
    
    // Mixed weights to simulate realistic scenario
    let weights: Vec<f64> = (0..n)
        .map(|i| {
            if i % 1000 < 50 {
                10.0 + (i % 10) as f64 // Peak regions with higher weight
            } else {
                1.0 + (i % 3) as f64 // Background regions
            }
        })
        .collect();
    
    let start = Instant::now();
    let mut sampler = AExpJSampler::new(k, Some(42));
    let result = sampler.sample_stream(
        regions.iter().zip(weights.iter()),
        |(_, &weight)| weight,
    );
    let elapsed = start.elapsed();
    
    assert_eq!(result.len(), k);
    
    let stats = sampler.get_statistics();
    println!("Large-scale A-ExpJ performance:");
    println!("  Time: {:?}", elapsed);
    println!("  Jump efficiency: {:.2}%", stats.jump_efficiency * 100.0);
    println!("  Items skipped: {}", stats.items_skipped);
    
    // Should achieve reasonable jump efficiency for 1% sample ratio
    assert!(stats.jump_efficiency > 0.1, 
            "Should achieve >10% jump efficiency for 1% sample ratio");
}


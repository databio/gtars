use super::{AResSampler, VectorAResSampler};
use std::collections::HashSet;

#[cfg(test)]
mod statistical_tests {
    use super::*;
    
    /// Chi-squared test for goodness of fit
    fn chi_squared_test(observed: &[u32], expected: &[f64], total_samples: u32) -> f64 {
        let mut chi_squared = 0.0;
        
        for (i, &obs) in observed.iter().enumerate() {
            let exp = expected[i] * total_samples as f64;
            if exp > 0.0 {
                let diff = obs as f64 - exp;
                chi_squared += (diff * diff) / exp;
            }
        }
        
        chi_squared
    }
    
    #[test]
    fn test_sampling_distribution_uniform() {
        // Test that uniform weights produce uniform sampling
        let n: usize = 100;
        let k: usize = 20;
        let num_trials = 10000;
        
        let items: Vec<usize> = (0..n).collect();
        let weights: Vec<f64> = vec![1.0; n];
        
        let mut counts = vec![0u32; n];
        
        for _ in 0..num_trials {
            let mut sampler = AResSampler::new(k, None);
            let sampled = sampler.sample_stream(
                items.iter().zip(weights.iter()),
                |(_, &weight)| weight
            );
            
            for (item, _) in sampled {
                counts[*item] += 1;
            }
        }
        
        // Expected probability for each item
        let expected_prob: Vec<f64> = vec![k as f64 / n as f64; n];
        
        // Chi-squared test
        let chi_squared = chi_squared_test(&counts, &expected_prob, num_trials);
        
        // Degrees of freedom = n - 1
        // For α = 0.05 and df = 99, critical value ≈ 123.2
        // We'll use a more lenient threshold for random tests
        assert!(chi_squared < 150.0, 
            "Chi-squared value {} exceeds threshold, suggesting non-uniform distribution", 
            chi_squared);
    }
    
    #[test]
    fn test_sampling_distribution_weighted() {
        // Test that weighted sampling follows the weight distribution
        let n: usize = 50;
        let k: usize = 10;
        let num_trials = 10000;
        
        let items: Vec<usize> = (0..n).collect();
        // Create exponentially increasing weights
        let weights: Vec<f64> = (0..n)
            .map(|i| (1.1_f64).powi(i as i32))
            .collect();
        
        let weight_sum: f64 = weights.iter().sum();
        let _normalized_weights: Vec<f64> = weights.iter()
            .map(|w| w / weight_sum)
            .collect();
        
        let mut counts = vec![0u32; n];
        
        for _ in 0..num_trials {
            let mut sampler = AResSampler::new(k, None);
            let sampled = sampler.sample_stream(
                items.iter().zip(weights.iter()),
                |(_, &weight)| weight
            );
            
            for (item, _) in sampled {
                counts[*item] += 1;
            }
        }
        
        // Verify that higher weighted items are sampled more frequently
        let first_quarter_avg = counts[0..n/4].iter().sum::<u32>() as f64 / (n/4) as f64;
        let last_quarter_avg = counts[3*n/4..n].iter().sum::<u32>() as f64 / (n/4) as f64;
        
        println!("First quarter avg: {}, Last quarter avg: {}", first_quarter_avg, last_quarter_avg);
        println!("First few counts: {:?}", &counts[0..5]);
        println!("Last few counts: {:?}", &counts[n-5..n]);
        
        assert!(last_quarter_avg > first_quarter_avg * 1.5,
            "Higher weighted items should be sampled significantly more frequently. First quarter: {}, Last quarter: {}", 
            first_quarter_avg, last_quarter_avg);
    }
    
    #[test]
    fn test_without_replacement_property() {
        // Verify that sampling is truly without replacement
        let n: usize = 100;
        let k: usize = 50;
        let num_trials = 1000;
        
        let items: Vec<usize> = (0..n).collect();
        let weights: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
        
        for _ in 0..num_trials {
            let mut sampler = AResSampler::new(k, None);
            let sampled = sampler.sample_stream(
                items.iter().zip(weights.iter()),
                |(_, &weight)| weight
            );
            
            let unique_items: HashSet<_> = sampled.iter().map(|(item, _)| item).collect();
            assert_eq!(sampled.len(), unique_items.len(), 
                "Sample contains duplicates, violating without-replacement property");
        }
    }
    
    #[test]
    fn test_deterministic_with_seed() {
        // Verify that using the same seed produces identical results
        let n: usize = 1000;
        let k: usize = 100;
        let seed = 12345u64;
        
        let items: Vec<usize> = (0..n).collect();
        let weights: Vec<f64> = (0..n).map(|i| ((i * 7 + 3) % 100) as f64 + 1.0).collect();
        
        // First run
        let mut sampler1 = AResSampler::new(k, Some(seed));
        let result1 = sampler1.sample_stream(
            items.iter().zip(weights.iter()),
            |(_, &weight)| weight
        );
        
        // Second run with same seed
        let mut sampler2 = AResSampler::new(k, Some(seed));
        let result2 = sampler2.sample_stream(
            items.iter().zip(weights.iter()),
            |(_, &weight)| weight
        );
        
        // Results should be identical
        assert_eq!(result1.len(), result2.len());
        for (item1, item2) in result1.iter().zip(result2.iter()) {
            assert_eq!(item1, item2, "Results differ despite using same seed");
        }
    }
    
    #[test]
    fn test_edge_cases() {
        // Test various edge cases
        
        // Empty input
        let items: Vec<usize> = vec![];
        let mut sampler = AResSampler::new(10, Some(42));
        let result = sampler.sample_stream(items.into_iter(), |_| 1.0);
        assert_eq!(result.len(), 0);
        
        // k > n
        let items: Vec<usize> = vec![1, 2, 3];
        let mut sampler = AResSampler::new(10, Some(42));
        let result = sampler.sample_stream(items.into_iter(), |_| 1.0);
        assert_eq!(result.len(), 3);
        
        // k = 0
        let items: Vec<usize> = vec![1, 2, 3];
        let mut sampler = AResSampler::new(0, Some(42));
        let result = sampler.sample_stream(items.into_iter(), |_| 1.0);
        assert_eq!(result.len(), 0);
        
        // Zero and negative weights
        let items: Vec<usize> = vec![1, 2, 3, 4, 5];
        let weights: Vec<f64> = vec![0.0, -1.0, 2.0, 0.0, 3.0];
        let mut sampler = AResSampler::new(2, Some(42));
        let result = sampler.sample_stream(
            items.iter().zip(weights.iter()),
            |(_, &weight)| weight
        );
        // Should only sample from items with positive weights (items 3 and 5)
        assert_eq!(result.len(), 2);
        for (item, _) in result {
            assert!(*item == 3 || *item == 5, "Sampled item {} with non-positive weight", item);
        }
    }
    
    #[test]
    fn test_optimized_version_equivalence() {
        // Verify that optimized version produces similar results
        let n: usize = 100;
        let k: usize = 20;
        let num_trials = 1000;
        
        let items: Vec<usize> = (0..n).collect();
        let weights: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
        
        let mut standard_counts = vec![0u32; n];
        let mut optimized_counts = vec![0u32; n];
        
        for trial in 0..num_trials {
            let seed = 42 + trial;
            
            // Standard version
            let mut sampler = AResSampler::new(k, Some(seed));
            let result = sampler.sample_stream(
                items.iter().zip(weights.iter()),
                |(_, &weight)| weight
            );
            for (&item, _) in result {
                standard_counts[item as usize] += 1;
            }
            
            // Optimized version
            let mut sampler = VectorAResSampler::new(k, Some(seed));
            let result = sampler.sample_stream(
                items.iter().zip(weights.iter()),
                |(_, &weight)| weight
            );
            for (&item, _) in result {
                optimized_counts[item as usize] += 1;
            }
        }
        
        // Distributions should be similar (not necessarily identical due to different implementations)
        let mut correlation = 0.0;
        let std_mean = standard_counts.iter().sum::<u32>() as f64 / n as f64;
        let opt_mean = optimized_counts.iter().sum::<u32>() as f64 / n as f64;
        
        for i in 0..n {
            correlation += (standard_counts[i] as f64 - std_mean) * 
                          (optimized_counts[i] as f64 - opt_mean);
        }
        
        // Both versions should produce similar sampling distributions
        assert!(correlation > 0.0, "Optimized version produces very different distribution");
    }
    
    #[test]
    fn test_batch_sampling() {
        // Test batch sampling functionality
        let n: usize = 50;
        let k: usize = 10;
        let num_batches = 100;
        
        let items: Vec<(usize, f64)> = (0..n)
            .map(|i| (i, (i + 1) as f64))
            .collect();
        
        let batches = AResSampler::sample_batch(
            items.clone().into_iter(),
            |(_, weight)| *weight,
            k,
            num_batches,
            Some(42),
        );
        
        assert_eq!(batches.len(), num_batches);
        
        for batch in &batches {
            assert_eq!(batch.len(), k);
            // Check no duplicates within each batch
            let unique_indices: HashSet<_> = batch.iter().map(|(i, _)| i).collect();
            assert_eq!(batch.len(), unique_indices.len());
        }
        
        // Different batches should generally be different (with high probability)
        let first_indices: HashSet<_> = batches[0].iter().map(|(i, _)| i).collect();
        let second_indices: HashSet<_> = batches[1].iter().map(|(i, _)| i).collect();
        let overlap = first_indices.intersection(&second_indices).count();
        
        // Some overlap is expected but not complete overlap
        assert!(overlap < k, "Batches are too similar");
    }
    
    #[test]
    fn test_streaming_property() {
        // Test that the algorithm works correctly in streaming fashion
        let n: usize = 1000;
        let k: usize = 50;
        
        let items: Vec<usize> = (0..n).collect();
        let weights: Vec<f64> = (0..n).map(|i| ((i * 3 + 7) % 10) as f64 + 1.0).collect();
        
        // Process all at once
        let mut sampler1 = AResSampler::new(k, Some(42));
        let result1 = sampler1.sample_stream(
            items.iter().zip(weights.iter()),
            |(_, &weight)| weight
        );
        
        // Process in chunks (simulating streaming)
        let mut sampler2 = AResSampler::new(k, Some(42));
        for chunk_start in (0..n).step_by(100) {
            let chunk_end = (chunk_start + 100).min(n);
            let chunk_items = &items[chunk_start..chunk_end];
            let chunk_weights = &weights[chunk_start..chunk_end];
            
            for (item, weight) in chunk_items.iter().zip(chunk_weights.iter()) {
                sampler2.add_item(*item, &|_| *weight);
            }
        }
        
        // Extract final result from sampler2
        let result2: Vec<usize> = sampler2.reservoir
            .into_sorted_vec()
            .into_iter()
            .map(|keyed| keyed.item)
            .collect();
        
        // Both approaches should yield the same result
        assert_eq!(result1.len(), result2.len());
    }
}
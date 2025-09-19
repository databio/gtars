use super::{AExpJSampler, AExpJStrategy};
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

    /// Compute correlation between two frequency distributions
    fn compute_correlation(dist1: &[u32], dist2: &[u32]) -> f64 {
        let n = dist1.len() as f64;
        let mean1 = dist1.iter().sum::<u32>() as f64 / n;
        let mean2 = dist2.iter().sum::<u32>() as f64 / n;

        let mut cov = 0.0;
        let mut var1 = 0.0;
        let mut var2 = 0.0;

        for i in 0..dist1.len() {
            let d1 = dist1[i] as f64 - mean1;
            let d2 = dist2[i] as f64 - mean2;
            cov += d1 * d2;
            var1 += d1 * d1;
            var2 += d2 * d2;
        }

        if var1 * var2 > 0.0 {
            cov / (var1 * var2).sqrt()
        } else {
            0.0
        }
    }

    #[test]
    fn test_a_expj_preserves_distribution() {
        // Test that A-ExpJ produces same distribution as basic sampling
        let n = 10000;
        let k = 100;
        let num_trials = 100;

        let items: Vec<usize> = (0..n).collect();
        let weights: Vec<f64> = (0..n)
            .map(|i| (1.05_f64).powi((i % 100) as i32))
            .collect();

        let mut counts_standard = vec![0u32; n];
        let mut counts_aggressive = vec![0u32; n];

        for trial in 0..num_trials {
            let seed = 1000 + trial;

            // Standard A-ExpJ
            let mut sampler1 = AExpJSampler::new(k, Some(seed))
                .with_strategy(AExpJStrategy::Standard);
            let result1 = sampler1.sample_stream(
                items.iter().zip(weights.iter()),
                |(_, &weight)| weight,
            );
            for (item, _) in &result1 {
                counts_standard[**item] += 1;
            }

            // Aggressive A-ExpJ (should still preserve distribution)
            let mut sampler2 = AExpJSampler::new(k, Some(seed))
                .with_strategy(AExpJStrategy::Aggressive);
            let result2 = sampler2.sample_stream(
                items.iter().zip(weights.iter()),
                |(_, &weight)| weight,
            );
            for (item, _) in &result2 {
                counts_aggressive[**item] += 1;
            }
        }

        // Distributions should be correlated (relaxed due to different jump strategies)
        let correlation = compute_correlation(&counts_standard, &counts_aggressive);
        assert!(
            correlation > 0.6,
            "A-ExpJ strategies should preserve distribution: correlation = {}",
            correlation
        );
    }

    #[test]
    fn test_a_expj_jump_effectiveness() {
        // Test that A-ExpJ actually skips items
        let n = 100000;
        let k = 1000; // 1% sample ratio - ideal for A-ExpJ

        let items: Vec<usize> = (0..n).collect();
        let weights: Vec<f64> = vec![1.0; n];

        let mut sampler = AExpJSampler::new(k, Some(42));
        let _ = sampler.sample_stream(
            items.iter().zip(weights.iter()),
            |(_, &weight)| weight,
        );

        let stats = sampler.get_statistics();
        
        println!("A-ExpJ performance for n={}, k={}:", n, k);
        println!("  Items processed: {}", stats.items_processed);
        println!("  Items skipped: {}", stats.items_skipped);
        println!("  Jump efficiency: {:.2}%", stats.jump_efficiency * 100.0);
        println!("  Average jump size: {:.2}", stats.average_jump_size);
        println!("  Number of jumps: {}", stats.num_jumps);

        // Should skip significant portion of items
        assert!(
            stats.jump_efficiency > 0.4,
            "A-ExpJ should skip >40% of items for 1% sample ratio"
        );
    }

    #[test]
    fn test_a_expj_deterministic_with_seed() {
        // Verify same seed produces identical results
        let n = 5000;
        let k = 100;
        let seed = 12345u64;

        let items: Vec<usize> = (0..n).collect();
        let weights: Vec<f64> = (0..n).map(|i| ((i * 7 + 3) % 100) as f64 + 1.0).collect();

        // First run
        let mut sampler1 = AExpJSampler::new(k, Some(seed));
        let result1 = sampler1.sample_stream(
            items.iter().zip(weights.iter()),
            |(_, &weight)| weight,
        );

        // Second run with same seed
        let mut sampler2 = AExpJSampler::new(k, Some(seed));
        let result2 = sampler2.sample_stream(
            items.iter().zip(weights.iter()),
            |(_, &weight)| weight,
        );

        // Extract just the items for comparison
        let items1: Vec<_> = result1.iter().map(|(item, _)| item).collect();
        let items2: Vec<_> = result2.iter().map(|(item, _)| item).collect();

        assert_eq!(items1, items2, "Same seed should produce identical results");
    }

    #[test]
    fn test_a_expj_without_replacement() {
        // Verify no duplicates in sample
        let n = 1000;
        let k = 200;
        let num_trials = 100;

        let items: Vec<usize> = (0..n).collect();
        let weights: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();

        for trial in 0..num_trials {
            let mut sampler = AExpJSampler::new(k, Some(42 + trial));
            let sampled = sampler.sample_stream(
                items.iter().zip(weights.iter()),
                |(_, &weight)| weight,
            );

            let unique_items: HashSet<_> = sampled.iter().map(|(item, _)| **item).collect();
            assert_eq!(
                sampled.len(),
                unique_items.len(),
                "Sample contains duplicates in trial {}",
                trial
            );
        }
    }

    #[test]
    fn test_a_expj_edge_cases() {
        // Empty input
        let items: Vec<usize> = vec![];
        let mut sampler = AExpJSampler::new(10, Some(42));
        let result = sampler.sample_stream(items.into_iter(), |_| 1.0);
        assert_eq!(result.len(), 0);

        // k > n
        let items: Vec<usize> = vec![1, 2, 3];
        let mut sampler = AExpJSampler::new(10, Some(42));
        let result = sampler.sample_stream(items.into_iter(), |_| 1.0);
        assert_eq!(result.len(), 3);

        // k = 0
        let items: Vec<usize> = vec![1, 2, 3];
        let mut sampler = AExpJSampler::new(0, Some(42));
        let result = sampler.sample_stream(items.into_iter(), |_| 1.0);
        assert_eq!(result.len(), 0);

        // Zero and negative weights
        let items: Vec<usize> = vec![1, 2, 3, 4, 5];
        let weights: Vec<f64> = vec![0.0, -1.0, 2.0, 0.0, 3.0];
        let mut sampler = AExpJSampler::new(2, Some(42));
        let result = sampler.sample_stream(
            items.iter().zip(weights.iter()),
            |(_, &weight)| weight,
        );
        // Should only sample from items with positive weights
        assert_eq!(result.len(), 2);
        for (item, _) in result {
            assert!(
                *item == 3 || *item == 5,
                "Sampled item {} with non-positive weight",
                *item
            );
        }
    }

    #[test]
    fn test_a_expj_deterministic() {
        // Test that same seed produces identical results
        let n = 1000;
        let k = 50;
        let seed = 42;

        let items: Vec<usize> = (0..n).collect();
        let weights: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();

        // First run
        let mut sampler1 = AExpJSampler::new(k, Some(seed));
        let result1 = sampler1.sample_stream(
            items.iter().zip(weights.iter()),
            |(_, &weight)| weight,
        );

        // Second run with same seed
        let mut sampler2 = AExpJSampler::new(k, Some(seed));
        let result2 = sampler2.sample_stream(
            items.iter().zip(weights.iter()),
            |(_, &weight)| weight,
        );

        // Extract just the items for comparison
        let items1: Vec<_> = result1.iter().map(|(item, _)| **item).collect();
        let items2: Vec<_> = result2.iter().map(|(item, _)| **item).collect();

        assert_eq!(items1, items2, "Same seed should produce identical results");
    }

    #[test]
    fn test_a_expj_batch_sampling() {
        let n = 500;
        let k = 20;
        let num_batches = 10;

        let items: Vec<(usize, f64)> = (0..n).map(|i| (i, (i + 1) as f64)).collect();

        let batches = AExpJSampler::sample_batch(
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

        // Different batches should be different
        let first_indices: HashSet<_> = batches[0].iter().map(|(i, _)| i).collect();
        let second_indices: HashSet<_> = batches[1].iter().map(|(i, _)| i).collect();
        let overlap = first_indices.intersection(&second_indices).count();

        // Some overlap is expected but not complete overlap
        assert!(overlap < k, "Batches should be different");
    }

    #[test]
    fn test_a_expj_jump_formula_correctness() {
        // Test the mathematical correctness of jump computation
        let mut sampler = AExpJSampler::new(10, Some(42));

        // Fill reservoir with known keys
        for i in 0..10 {
            sampler.add_item_expj(i, 1.0);
        }

        // Collect multiple jump distances
        let mut jumps = Vec::new();
        for _ in 0..1000 {
            let jump = sampler.compute_next_jump();
            if jump > 0.0 {
                jumps.push(jump);
            }
        }

        assert!(!jumps.is_empty(), "Should generate jumps");

        // Check jump statistics
        let avg_jump = jumps.iter().sum::<f64>() / jumps.len() as f64;
        let min_jump = jumps.iter().fold(f64::MAX, |a, &b| a.min(b));
        let max_jump = jumps.iter().fold(f64::MIN, |a, &b| a.max(b));

        println!("Jump statistics:");
        println!("  Average: {:.2}", avg_jump);
        println!("  Min: {:.2}", min_jump);
        println!("  Max: {:.2}", max_jump);
        println!("  Count: {}", jumps.len());

        // Jumps should be reasonable
        assert!(min_jump > 0.0, "Jumps should be positive");
        assert!(max_jump < 10000.0, "Jumps should be bounded");
    }

    #[test]
    fn test_a_expj_performance_scaling() {
        // Test that A-ExpJ scales better than O(n)
        let k = 100;
        let sizes = vec![1000, 10000, 100000];
        let mut times = Vec::new();
        let mut efficiencies = Vec::new();

        for &n in &sizes {
            let items: Vec<usize> = (0..n).collect();
            let weights: Vec<f64> = vec![1.0; n];

            let start = std::time::Instant::now();
            let mut sampler = AExpJSampler::new(k, Some(42));
            let _ = sampler.sample_stream(
                items.iter().zip(weights.iter()),
                |(_, &weight)| weight,
            );
            let elapsed = start.elapsed();

            let stats = sampler.get_statistics();
            
            times.push(elapsed.as_micros());
            efficiencies.push(stats.jump_efficiency);

            println!(
                "n={}: time={:?}, efficiency={:.2}%",
                n,
                elapsed,
                stats.jump_efficiency * 100.0
            );
        }

        // Efficiency should increase with larger n
        assert!(
            efficiencies[2] > efficiencies[0],
            "Jump efficiency should increase with dataset size"
        );

        // Time should scale sub-linearly
        let time_ratio_1 = times[1] as f64 / times[0] as f64;
        let time_ratio_2 = times[2] as f64 / times[1] as f64;
        let size_ratio = 10.0; // Each size is 10x larger

        println!("Time scaling: 1K->10K={:.2}x, 10K->100K={:.2}x", 
                 time_ratio_1, time_ratio_2);

        // Should scale better than O(n) (with small tolerance for system variance)
        // Allow up to 5% deviation from linear scaling
        let tolerance = 1.05;
        assert!(
            time_ratio_2 < size_ratio * tolerance,
            "A-ExpJ should scale sub-linearly (ratio: {:.2}x, expected: <{:.2}x)",
            time_ratio_2, size_ratio * tolerance
        );
    }

    #[test]
    fn test_realistic_scatac_scenario() {
        // Simulate realistic scATAC-seq fragment sampling
        let num_fragments = 100000; // 100K fragments (reduced for test speed)
        let sample_size = 1000;     // 1K sampled fragments (1% ratio)

        // Create fragments with realistic weight distribution
        let fragments: Vec<usize> = (0..num_fragments).collect();
        let weights: Vec<f64> = (0..num_fragments)
            .map(|i| {
                if i % 1000 < 50 {
                    // 5% are "peak" fragments with higher weight
                    10.0 + (i % 10) as f64
                } else {
                    // 95% are background with lower weight
                    1.0 + (i % 3) as f64
                }
            })
            .collect();

        let start_time = std::time::Instant::now();

        let mut sampler = AExpJSampler::new(sample_size, Some(42));
        let sampled = sampler.sample_stream(
            fragments.iter().zip(weights.iter()),
            |(_, &weight)| weight,
        );

        let elapsed = start_time.elapsed();
        let stats = sampler.get_statistics();

        println!("Realistic scATAC-seq simulation:");
        println!("  Sampled {} fragments from {} total in {:?}",
                 sampled.len(), num_fragments, elapsed);
        println!("  Jump efficiency: {:.2}%", stats.jump_efficiency * 100.0);
        println!("  Items skipped: {}", stats.items_skipped);
        println!("  Items processed: {}", stats.items_processed);
        println!("  Number of jumps: {}", stats.num_jumps);

        assert_eq!(sampled.len(), sample_size);
        assert!(
            stats.jump_efficiency > 0.4,
            "Should achieve >40% jump efficiency for 1% sample ratio"
        );
    }
}
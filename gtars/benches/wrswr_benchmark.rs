use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};
use gtars::scatrs::models::ScatrsRegion;
use gtars::scatrs::sampling::{WeightedSampler, wrswr_skip::{WRSWRSkipSampler, OptimizedWRSWRSkip}};
use rand::prelude::*;
use std::time::Duration;

fn generate_test_regions(n: usize) -> Vec<ScatrsRegion> {
    (0..n)
        .map(|i| ScatrsRegion::new(
            format!("chr{}", (i % 22) + 1),
            (i as u64) * 1000,
            (i as u64) * 1000 + 500,
        ))
        .collect()
}

fn generate_test_weights(n: usize) -> Vec<f64> {
    let mut rng = StdRng::seed_from_u64(42);
    (0..n).map(|_| rng.gen::<f64>() * 100.0 + 1.0).collect()
}

fn benchmark_wrswr_skip(c: &mut Criterion) {
    let mut group = c.benchmark_group("wrswr_skip_sampling");
    group.measurement_time(Duration::from_secs(10));
    
    // Test different dataset sizes and sample ratios
    let test_cases = vec![
        (1_000, 100),      // 10% sample ratio
        (10_000, 500),     // 5% sample ratio
        (100_000, 1_000),  // 1% sample ratio
        (100_000, 10_000), // 10% sample ratio
    ];
    
    for (n, k) in test_cases {
        let regions = generate_test_regions(n);
        let weights = generate_test_weights(n);
        
        // Benchmark WRSWR-SKIP
        group.bench_with_input(
            BenchmarkId::new("wrswr_skip", format!("n={}_k={}", n, k)),
            &(&regions, &weights, k),
            |b, (regions, weights, k)| {
                b.iter(|| {
                    let mut sampler = WRSWRSkipSampler::new(*k, Some(42));
                    sampler.sample_stream(
                        regions.iter().cloned().zip(weights.iter().cloned()),
                        |(_, weight)| weight
                    )
                });
            },
        );
        
        // Benchmark Optimized WRSWR-SKIP
        group.bench_with_input(
            BenchmarkId::new("optimized_wrswr", format!("n={}_k={}", n, k)),
            &(&regions, &weights, k),
            |b, (regions, weights, k)| {
                b.iter(|| {
                    let mut sampler = OptimizedWRSWRSkip::new(*k, Some(42));
                    sampler.sample_stream(
                        regions.iter().cloned().zip(weights.iter().cloned()),
                        |(_, weight)| weight
                    )
                });
            },
        );
        
        // Benchmark current implementation (if k is small enough to be practical)
        if n <= 10_000 && k <= 1_000 {
            group.bench_with_input(
                BenchmarkId::new("current_impl", format!("n={}_k={}", n, k)),
                &(&regions, &weights, k),
                |b, (regions, weights, k)| {
                    b.iter(|| {
                        let sampler = WeightedSampler::new(Some(42));
                        let mut rng = StdRng::seed_from_u64(42);
                        // Call the private method through reflection or make it public for benchmarking
                        // For now, we'll skip this as it requires modifying the existing code
                    });
                },
            );
        }
    }
    
    group.finish();
}

fn benchmark_batch_sampling(c: &mut Criterion) {
    let mut group = c.benchmark_group("batch_sampling");
    group.measurement_time(Duration::from_secs(10));
    
    let n = 10_000;
    let k = 100;
    let regions: Vec<ScatrsRegion> = generate_test_regions(n);
    let weights: Vec<f64> = generate_test_weights(n);
    
    for num_batches in [1, 10, 100, 1000] {
        group.bench_with_input(
            BenchmarkId::new("batch_wrswr", format!("batches={}", num_batches)),
            &(&regions, &weights, k, num_batches),
            |b, (regions, weights, k, num_batches)| {
                b.iter(|| {
                    WRSWRSkipSampler::sample_batch(
                        regions.iter().zip(weights.iter()).map(|(r, &w)| (r.clone(), w)),
                        |(_, weight)| weight,
                        *k,
                        *num_batches,
                        Some(42),
                    )
                });
            },
        );
    }
    
    group.finish();
}

fn benchmark_memory_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("memory_scaling");
    
    // Test memory efficiency with different k values
    for k in [100, 1_000, 10_000] {
        let n = 100_000;
        let regions = generate_test_regions(n);
        let weights = generate_test_weights(n);
        
        group.bench_with_input(
            BenchmarkId::new("memory_k", k),
            &(&regions, &weights, k),
            |b, (regions, weights, k)| {
                b.iter(|| {
                    let mut sampler = OptimizedWRSWRSkip::new(*k, Some(42));
                    sampler.sample_stream(
                        regions.iter().cloned().zip(weights.iter().cloned()),
                        |(_, weight)| weight
                    )
                });
            },
        );
    }
    
    group.finish();
}

criterion_group!(
    benches, 
    benchmark_wrswr_skip,
    benchmark_batch_sampling,
    benchmark_memory_scaling
);
criterion_main!(benches);
use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use gtars::scatrs::models::ScatrsRegion;
use gtars::scatrs::sampling::a_expj::{AExpJSampler, AExpJStrategy};
use gtars::scatrs::sampling::a_res::{AResSampler, VectorAResSampler};
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

fn benchmark_a_expj_vs_a_res(c: &mut Criterion) {
    let mut group = c.benchmark_group("a_expj_vs_a_res");
    group.measurement_time(Duration::from_secs(10));
    
    // Test different dataset sizes and sample ratios
    let test_cases = vec![
        (10_000, 100),      // 1% sample ratio - ideal for A-ExpJ
        (10_000, 500),      // 5% sample ratio
        (100_000, 1_000),   // 1% sample ratio
        (100_000, 10_000),  // 10% sample ratio
    ];
    
    for (n, k) in test_cases {
        let regions = generate_test_regions(n);
        let weights = generate_test_weights(n);
        
        // Benchmark A-Res (no exponential jumps)
        group.bench_with_input(
            BenchmarkId::new("a_res", format!("n={}_k={}", n, k)),
            &(&regions, &weights, k),
            |b, (regions, weights, k)| {
                b.iter(|| {
                    let mut sampler = AResSampler::new(*k, Some(42));
                    sampler.sample_stream(
                        regions.iter().cloned().zip(weights.iter().cloned()),
                        |(_, weight)| *weight
                    )
                });
            },
        );
        
        // Benchmark A-ExpJ (with exponential jumps)
        group.bench_with_input(
            BenchmarkId::new("a_expj", format!("n={}_k={}", n, k)),
            &(&regions, &weights, k),
            |b, (regions, weights, k)| {
                b.iter(|| {
                    let mut sampler = AExpJSampler::new(*k, Some(42));
                    sampler.sample_stream(
                        regions.iter().cloned().zip(weights.iter().cloned()),
                        |(_, weight)| *weight
                    )
                });
            },
        );
        
    }
    
    group.finish();
}

fn benchmark_a_expj_strategies(c: &mut Criterion) {
    let mut group = c.benchmark_group("a_expj_strategies");
    
    let n = 100_000;
    let k = 1_000; // 1% sample ratio
    let regions = generate_test_regions(n);
    let weights = generate_test_weights(n);
    
    for strategy in [
        AExpJStrategy::Standard,
        AExpJStrategy::Conservative,
        AExpJStrategy::Aggressive,
        AExpJStrategy::Disabled,
    ] {
        group.bench_with_input(
            BenchmarkId::new("strategy", format!("{:?}", strategy)),
            &(&regions, &weights, k, strategy),
            |b, (regions, weights, k, strategy)| {
                b.iter(|| {
                    let mut sampler = AExpJSampler::new(*k, Some(42))
                        .with_strategy(*strategy);
                    sampler.sample_stream(
                        regions.iter().cloned().zip(weights.iter().cloned()),
                        |(_, weight)| *weight
                    )
                });
            },
        );
    }
    
    group.finish();
}

fn benchmark_weight_distributions(c: &mut Criterion) {
    let mut group = c.benchmark_group("weight_distributions");
    
    let n = 50_000;
    let k = 500; // 1% sample ratio
    
    // Uniform weights
    let uniform_regions = generate_test_regions(n);
    let uniform_weights: Vec<f64> = vec![1.0; n];
    
    // Exponential weights
    let exp_regions = generate_test_regions(n);
    let exp_weights: Vec<f64> = (0..n)
        .map(|i| (1.01_f64).powi((i % 1000) as i32))
        .collect();
    
    // Heavy-tailed weights (power law)
    let heavy_regions = generate_test_regions(n);
    let heavy_weights: Vec<f64> = (0..n)
        .map(|i| 1.0 / ((i % 1000) + 1) as f64)
        .collect();
    
    // Benchmark uniform weights
    group.bench_function("uniform_a_expj", |b| {
        b.iter(|| {
            let mut sampler = AExpJSampler::new(k, Some(42));
            black_box(sampler.sample_stream(
                uniform_regions.iter().cloned().zip(uniform_weights.iter().cloned()),
                |(_, weight)| *weight
            ))
        });
    });
    
    group.bench_function("uniform_a_res", |b| {
        b.iter(|| {
            let mut sampler = AResSampler::new(k, Some(42));
            black_box(sampler.sample_stream(
                uniform_regions.iter().cloned().zip(uniform_weights.iter().cloned()),
                |(_, weight)| *weight
            ))
        });
    });
    
    // Benchmark exponential weights
    group.bench_function("exponential_a_expj", |b| {
        b.iter(|| {
            let mut sampler = AExpJSampler::new(k, Some(42));
            black_box(sampler.sample_stream(
                exp_regions.iter().cloned().zip(exp_weights.iter().cloned()),
                |(_, weight)| *weight
            ))
        });
    });
    
    group.bench_function("exponential_a_res", |b| {
        b.iter(|| {
            let mut sampler = AResSampler::new(k, Some(42));
            black_box(sampler.sample_stream(
                exp_regions.iter().cloned().zip(exp_weights.iter().cloned()),
                |(_, weight)| *weight
            ))
        });
    });
    
    // Benchmark heavy-tailed weights
    group.bench_function("heavy_tail_a_expj", |b| {
        b.iter(|| {
            let mut sampler = AExpJSampler::new(k, Some(42));
            black_box(sampler.sample_stream(
                heavy_regions.iter().cloned().zip(heavy_weights.iter().cloned()),
                |(_, weight)| *weight
            ))
        });
    });
    
    group.bench_function("heavy_tail_a_res", |b| {
        b.iter(|| {
            let mut sampler = AResSampler::new(k, Some(42));
            black_box(sampler.sample_stream(
                heavy_regions.iter().cloned().zip(heavy_weights.iter().cloned()),
                |(_, weight)| *weight
            ))
        });
    });
    
    group.finish();
}

fn benchmark_batch_sampling(c: &mut Criterion) {
    let mut group = c.benchmark_group("batch_sampling_a_expj");
    group.measurement_time(Duration::from_secs(10));
    
    let n = 10_000;
    let k = 100;
    let regions: Vec<ScatrsRegion> = generate_test_regions(n);
    let weights: Vec<f64> = generate_test_weights(n);
    
    for num_batches in [1, 10, 100] {
        group.bench_with_input(
            BenchmarkId::new("batch_a_expj", format!("batches={}", num_batches)),
            &(&regions, &weights, k, num_batches),
            |b, (regions, weights, k, num_batches)| {
                b.iter(|| {
                    AExpJSampler::sample_batch(
                        regions.iter().zip(weights.iter()).map(|(r, &w)| (r.clone(), w)),
                        |(_, weight)| *weight,
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

fn benchmark_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("scaling");
    
    // Test how A-ExpJ scales with different n values for fixed k
    let k = 1000;
    for n in [10_000, 50_000, 100_000, 500_000] {
        let regions = generate_test_regions(n);
        let weights = vec![1.0; n]; // Uniform for consistent comparison
        
        group.bench_with_input(
            BenchmarkId::new("a_expj_scaling", n),
            &(&regions, &weights, k),
            |b, (regions, weights, k)| {
                b.iter(|| {
                    let mut sampler = AExpJSampler::new(*k, Some(42));
                    sampler.sample_stream(
                        regions.iter().cloned().zip(weights.iter().cloned()),
                        |(_, weight)| *weight
                    )
                });
            },
        );
    }
    
    group.finish();
}

criterion_group!(
    benches,
    benchmark_a_expj_vs_a_res,
    benchmark_a_expj_strategies,
    benchmark_weight_distributions,
    benchmark_batch_sampling,
    benchmark_scaling
);
criterion_main!(benches);
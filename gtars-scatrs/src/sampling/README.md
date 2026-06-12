# ScATRS Sampling Module

This module implements weighted reservoir sampling algorithms for single-cell ATAC-seq simulation.

## Structure

```
sampling/
├── mod.rs                 # Module exports and high-level sampling logic
├── a_res.rs              # A-Res: Basic Algorithm A (reservoir sampling)
├── a_res_tests.rs        # Unit tests for A-Res
├── a_expj.rs             # A-ExpJ: Algorithm A with Exponential Jumps
├── a_expj_tests.rs       # Unit tests for A-ExpJ
├── tests/                # Integration tests
│   ├── mod.rs
│   ├── integration_tests.rs       # Original A-Res integration tests
│   └── a_expj_integration_tests.rs # A-ExpJ integration tests
└── benches/              # Performance benchmarks
    ├── README.md
    ├── a_res_benchmark.rs
    └── a_expj_benchmark.rs
```

## Algorithms

### A-Res (Algorithm A with Reservoir)
- **Paper**: Efraimidis & Spirakis (2006)
- **Complexity**: O(n + k log k)
- **Description**: Basic weighted reservoir sampling without replacement
- **Use case**: General purpose weighted sampling

### A-ExpJ (Algorithm A with Exponential Jumps)
- **Paper**: Efraimidis & Spirakis (2006)
- **Complexity**: O(k log(n/k) + k log k)
- **Description**: Optimized version using exponential jumps to skip items
- **Use case**: Large datasets with small sample ratios (k/n < 10%)

## Key Features

1. **Without Replacement**: Each item can only be sampled once
2. **Weighted Sampling**: Items are sampled proportional to their weights
3. **Streaming**: Can process items one at a time without knowing total count
4. **Deterministic**: Same seed produces identical results
5. **Memory Efficient**: O(k) memory usage

## Usage

```rust
use crate::scatrs::sampling::a_expj::AExpJSampler;

// Sample 100 items from a weighted stream
let mut sampler = AExpJSampler::new(100, Some(42));
let result = sampler.sample_stream(
    items.iter().zip(weights.iter()),
    |(_, &weight)| weight,
);
```

## Performance

For small sample ratios (k/n < 10%), A-ExpJ provides significant speedup:
- 1% ratio: 5-10x faster
- 5% ratio: 2-5x faster
- 10% ratio: 1.5-2x faster

## Testing

```bash
# Run unit tests
cargo test sampling::a_res
cargo test sampling::a_expj

# Run integration tests
cargo test sampling::tests

# Run benchmarks
cargo bench --bench a_res_benchmark
cargo bench --bench a_expj_benchmark
```

## Implementation Notes

- Both algorithms use the key generation formula: `key = u^(1/weight)` where u ~ Uniform(0,1)
- The reservoir maintains the k items with the largest keys
- A-ExpJ adds jump computation: `Xw = log(random) / log(Tw)` to skip items
- Two implementations available: heap-based and vector-based (better cache locality)
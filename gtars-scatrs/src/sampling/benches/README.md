# ScATRS Sampling Benchmarks

This directory contains performance benchmarks for the weighted reservoir sampling algorithms used in the ScATRS (Single-cell ATAC-seq Simulation) module.

## Available Benchmarks

### 1. `a_res_benchmark`
Benchmarks for the A-Res (Algorithm A with Reservoir) implementation - the basic weighted reservoir sampling without replacement algorithm.

### 2. `a_expj_benchmark`
Benchmarks for the A-ExpJ (Algorithm A with Exponential Jumps) implementation - an optimized version that uses exponential jumps to skip items that cannot enter the reservoir.

## Running Benchmarks

### Run All Benchmarks
```bash
cargo bench
```

### Run Specific Benchmark
```bash
# Run only A-Res benchmarks
cargo bench --bench a_res_benchmark

# Run only A-ExpJ benchmarks
cargo bench --bench a_expj_benchmark
```

### Run with Baseline Comparison
To compare against a saved baseline:
```bash
# Save a baseline first
cargo bench --bench a_expj_benchmark -- --save-baseline main

# Compare against baseline
cargo bench --bench a_expj_benchmark -- --baseline main
```

### Filter Specific Tests
Run only specific benchmark groups:
```bash
# Run only the scaling benchmarks
cargo bench --bench a_expj_benchmark -- scaling

# Run only vector implementation benchmarks
cargo bench --bench a_res_benchmark -- vector
```

## Interpreting Results

The benchmarks use the Criterion.rs framework, which provides:
- **Median time**: The median execution time across iterations
- **Standard deviation**: Variance in execution times
- **Throughput**: Items processed per second (where applicable)
- **Statistical significance**: Automatic detection of performance regressions

Results are saved in `target/criterion/` with HTML reports for detailed analysis.

## Benchmark Scenarios

### A-ExpJ vs A-Res Comparison
Tests different dataset sizes and sample ratios:
- **Small datasets**: 10K items with various sample ratios (1%, 5%, 10%)
- **Large datasets**: 100K-1M items to test scaling behavior
- **Weight distributions**: Uniform, exponential, and heavy-tailed distributions

### Key Performance Metrics

#### Expected Performance (A-ExpJ vs A-Res)
| Sample Ratio | Expected Speedup | Use Case |
|--------------|------------------|----------|
| 1% | 5-10x | Typical scATAC-seq simulation |
| 5% | 2-5x | Moderate sampling |
| 10% | 1.5-2x | Higher sampling rates |
| >20% | 1-1.2x | A-ExpJ overhead may not be worth it |

### Memory Usage
Both algorithms use O(k) memory where k is the sample size. The benchmarks test:
- Heap-based implementation (better for dynamic k)
- Vector-based implementation (better cache locality)

## Adding New Benchmarks

To add a new benchmark:
1. Create a new file in `benches/` (e.g., `my_benchmark.rs`)
2. Add the benchmark configuration to `Cargo.toml`:
   ```toml
   [[bench]]
   name = "my_benchmark"
   harness = false
   ```
3. Use the Criterion API for consistent reporting

## Profiling

For detailed profiling beyond benchmarks:
```bash
# CPU profiling with flamegraph
cargo bench --bench a_expj_benchmark -- --profile-time=10

# Memory profiling with Valgrind (Linux)
cargo bench --bench a_expj_benchmark -- --profile-memory
```

## Continuous Benchmarking

For CI/CD integration, use the `--save-baseline` flag to track performance over time:
```bash
# In CI pipeline
cargo bench -- --save-baseline $COMMIT_SHA
```

## Understanding A-ExpJ Optimization

The A-ExpJ algorithm improves on A-Res by:
- **Skipping items**: Uses exponential jumps to skip items that cannot enter the reservoir
- **Complexity**: Reduces from O(n) to O(k log(n/k)) for small sample ratios
- **Best for**: Small k/n ratios (<10%) common in genomics applications

The benchmarks demonstrate these improvements across various scenarios relevant to scATAC-seq simulation.
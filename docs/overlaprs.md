# gtars-overlaprs

Core infrastructure for all genomic interval overlap computations.

## Purpose

**gtars-overlaprs** provides the foundational overlap computation primitives that all other modules build upon. This is the single source of truth for efficient interval operations in gtars.

## Features

- Efficient interval intersection algorithms using AIList data structure
- Coverage calculations with bit vectors
- Support for large-scale genomic datasets (millions of intervals)
- Optimized memory usage and logarithmic query time
- Common `Overlapper` trait for consistent interface

## Usage

### From Rust
```rust
use gtars_overlaprs::{AIList, Interval};

// Create an AIList from intervals
let mut ailist = AIList::new();
ailist.add(Interval::new(100, 200));
ailist.add(Interval::new(150, 250));

// Query overlaps
let overlaps = ailist.query(120, 180);
```

## Data Structures

- **AIList** - Augmented Interval List for fast overlap queries
- **Bits** - Bit vector implementation for coverage calculations
- **Overlapper** - Trait for implementing overlap operations

## Performance

Optimized for:

- Large BED files (millions of intervals)
- Fast overlap queries (logarithmic time)
- Memory efficiency through compact data structures
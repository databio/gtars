# gtars-scatrs

SCATRS (Single-Cell ATAC-seq Region Simulator) - High-performance single-cell ATAC-seq fragment simulation.

## Features

- Realistic fragment generation based on empirical distributions
- Cell type-specific simulation
- Quality score modeling
- Batch and barcode generation
- Parallel processing with rayon

## Usage

### Command Line
```bash
scatrs simulate \
  --peaks peaks.bed \
  --cells 10000 \
  --output fragments.tsv.gz
```

### From Rust
```rust
use gtars_scatrs::{Simulator, SimulationParams};

let params = SimulationParams {
    num_cells: 10000,
    mean_fragments_per_cell: 5000,
    ..Default::default()
};

let simulator = Simulator::new(params);
simulator.run("peaks.bed", "output.tsv.gz")?;
```

## Simulation Parameters

- Number of cells
- Fragment size distribution
- Coverage depth
- Cell type proportions
- Technical noise parameters

## Output Format

Standard 10x Genomics fragment file format:
```
chr1    1000    1500    AAACCCAAGAAACACT-1    1
chr1    2000    2300    AAACCCAAGAAACACT-1    1
```
# gtars-scoring

Wrapper around gtars-overlaprs for generating X-by-peak count matrices.

## Purpose

**gtars-scoring** wraps the core overlap infrastructure from gtars-overlaprs to produce count matrices for genomic analysis. This module handles all use cases that need X-by-peak matrices:

- **Cell-by-peak** matrices (single-cell ATAC-seq)
- **Sample-by-peak** matrices (bulk ATAC-seq/ChIP-seq)
- **Pseudobulk-by-peak** matrices (aggregated single-cell data)

## What it does

Generates count matrices from fragment-peak overlaps where:

- Rows = analysis units (cells, samples, or pseudobulk aggregations)
- Columns = consensus regions (peaks)
- Values = overlap counts

## CLI Usage

```bash
gtars fscoring <fragments> <consensus> /
  --mode <atac|chip> /
  --output <output.csv.gz>
```

**Arguments:**

- `fragments`: Glob pattern for fragment files (e.g., `"*.bed.gz"`)
- `consensus`: Path to consensus peak set (BED file)
- `--mode`: Scoring mode (`atac` or `chip`, default: `atac`)
  - `atac`: Counts 5' and 3' cut sites with shifts (+4/-5 bp)
  - `chip`: Counts full fragment overlaps
- `--output`: Output file path (default: `count_matrix.csv.gz`)

**Example:**
```bash
gtars fscoring "fragments/*.bed.gz" peaks.bed /
  --mode atac /
  --output counts.csv.gz
```

## Rust Library Usage

```rust
use gtars::scoring::{
    region_scoring_from_fragments, 
    ConsensusSet, 
    FragmentFileGlob, 
    ScoringMode
};

// Set up inputs
let mut fragments = FragmentFileGlob::new("path/to/*.bed.gz")?;
let consensus = ConsensusSet::new("peaks.bed".into())?;

// Generate count matrix
let count_matrix = region_scoring_from_fragments(
    &mut fragments, 
    &consensus, 
    ScoringMode::Atac
)?;

// Write to file
count_matrix.write_to_file("output.csv.gz")?;

// Or access counts directly
let count = count_matrix.get(row, col);
```

## Output Format

Gzipped CSV file with:

- Each row = one fragment file
- Each column = one consensus region
- Values = count of overlapping fragments

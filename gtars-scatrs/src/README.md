# Scatrs - Single-cell ATAC-seq Transformer in Rust

A high-performance module for simulating single-cell ATAC-seq data from bulk ATAC-seq data, with realistic artifacts like doublets and ambient contamination.

## Overview: Two-Step Process

Scatrs uses a two-step approach to efficiently simulate single-cell data:

1. **Stage** - Process bulk ATAC-seq BAM files once to extract fragments and peaks. This computationally intensive step creates reusable staged data and outputs a `stage_config.yaml` file.

2. **Simulate** - Generate single cells by sampling from the staged data. This fast step can be run multiple times with different parameters (cell counts, doublet rates, contamination levels) without re-processing BAMs.

### Unified Configuration System

Scatrs uses a unified YAML configuration file that serves both the staging and simulation steps:

- **Unified config (scatrs_config.yaml)** - Contains all parameters for both staging and simulation
- **CSV manifest (samples.csv)** - Lists samples with paths and weights, referenced by the unified config
- **Stage config (stage_config.yaml)** - Created automatically by `stage` command, contains staging metadata

This design allows you to:
- Stage data once, simulate many times
- Quickly iterate on simulation parameters
- Test different experimental conditions without reprocessing
- Share staged data across simulations

## Features

- **Weighted Cell Types**: Control cell type proportions via sample weights
- **Doublet Simulation**: Create multiplets with configurable rates and weighted mixing
- **Ambient DNA**: Add realistic background contamination
- **Pipeline Agnostic**: Works with any peak caller (MACS2, HOMER, etc.) 
- **High Performance**: Parallel processing with progress indicators
- **Flexible Modes**: Simple (BAM folder) or Advanced (CSV manifest)

## Quick Start

### Unified Configuration Mode (Recommended)

```bash
# Step 1: Generate configuration template
gtars scatrs config --output scatrs_config.yaml

# Step 2: Create your sample manifest (samples.csv)
cat > samples.csv << EOF
sample_name,bam_file,peaks_file,weight
CD4_T,/data/cd4t.bam,/data/cd4t_peaks.narrowPeak,3.0
CD8_T,/data/cd8t.bam,/data/cd8t_peaks.narrowPeak,2.0
B_cells,/data/bcells.bam,/data/bcells_peaks.narrowPeak,1.5
NK_cells,/data/nk.bam,,0.5
EOF

# Step 3: Edit scatrs_config.yaml to set your parameters
# The config contains both staging and simulation sections

# Step 4: Stage data (do this once)
gtars scatrs stage --config scatrs_config.yaml
# Updates the config with staging_config path

# Step 5: Simulate cells (can run multiple times)
gtars scatrs simulate --config scatrs_config.yaml
# Uses the same config file, now with staging_config filled in

# To change simulation parameters without re-staging:
# 1. Edit the simulation section in scatrs_config.yaml
# 2. Or override with CLI args:
gtars scatrs simulate \
    --config scatrs_config.yaml \
    --cells 5000 \
    --doublet-rate 0.15 \
    --output results_v2/
```

### Legacy Modes (Backward Compatible)

#### Simple Mode (Equal cell type distribution)

```bash
# Stage BAM files without a config
gtars scatrs stage \
    --bam-dir /path/to/bams/ \
    --chrom-sizes hg38.sizes \
    --output staged/

# Then create a config for simulation...
```

#### Advanced Mode (CSV manifest without unified config)

```bash
# Stage with manifest directly
gtars scatrs stage \
    --manifest samples.csv \
    --chrom-sizes hg38.sizes \
    --output staged/
```

## How It Works

### Unified Configuration Workflow

The unified YAML config file contains all parameters for both staging and simulation:

```yaml
# scatrs_config.yaml
sample_table: samples.csv      # CSV manifest with sample paths and weights
chrom_sizes: hg38.sizes        # Chromosome sizes file
blacklist: blacklist.bed       # Optional blacklist regions

staging:
  output_dir: staged/
  extend_peaks: 250
  bin_size: 5000
  merge_distance: 20

staging_config: null           # Filled automatically after staging

simulation:
  output_dir: results/
  cell_count: 1000
  fragments_per_cell: 10000
  signal_to_noise: 0.8
  doublet_rate: 0.1
  ambient_contamination: 0.02
  fragment_distribution:
    type: gaussian
    mean: 10000.0
    variance: 4000000.0
    min: 1000
    max: 50000
```

### Stage Step
During staging, Scatrs:
1. **Reads unified config** to get all parameters and sample table location
2. **Loads sample manifest** (CSV) to get sample names, BAM/peak paths, and weights
3. **Extracts fragments** from BAM files (paired, non-duplicate reads)
4. **Merges peaks** across samples (if provided) to create consensus peak set
5. **Generates background regions** from non-peak genomic bins
6. **Counts fragments** per sample and in peaks/background for weighted sampling
7. **Saves staged data** in an efficient format for reuse
8. **Creates stage_config.yaml** with staging metadata
9. **Updates unified config** by adding the staging_config path

This is the slowest step but only needs to be done once per dataset.

### Simulate Step  
During simulation, Scatrs:
1. **Reads unified config** including the staging_config path
2. **Loads stage config** to get staged data location and metadata
3. **Reads sample table** to get current weights (can be modified without re-staging!)
4. **Assigns cell types** based on weights from the CSV manifest
5. **Samples regions** for each cell (peaks vs background based on signal-to-noise)
6. **Samples fragments** uniformly from selected regions
7. **Creates doublets** by merging cells (if doublet-rate > 0)
8. **Adds ambient DNA** from the total fragment pool (if contamination > 0)
9. **Outputs files** in standard formats (10x fragments, BED, etc.)

This step is fast and can be repeated with different parameters.

## Key Parameters

### Unified Config Mode
- `--config`: Unified YAML configuration file containing all parameters

### Staging (`stage` command)
When using unified config, parameters come from the `staging` section:
- `output_dir`: Where to save staged data
- `extend_peaks`: Extend peaks by N bp (default: 250)
- `bin_size`: Background region bin size (default: 5000)
- `merge_distance`: Merge peaks within N bp (default: 20)

Legacy CLI options (when not using unified config):
- `--bam-dir`: Directory with BAM files (simple mode)
- `--manifest`: CSV with samples, paths, and weights (advanced mode)
- `--chrom-sizes`: Chromosome sizes file (required)
- `--blacklist`: Optional blacklist regions BED file

### Simulation (`simulate` command)
When using unified config, parameters come from the `simulation` section:
- `cell_count`: Number of cells to simulate
- `fragments_per_cell`: Average fragments per cell
- `signal_to_noise`: Peak vs background ratio (0-1, only with peaks)
- `doublet_rate`: Fraction of doublets (0-1)
- `ambient_contamination`: Fraction of ambient DNA (0-1)
- `fragment_distribution`: Distribution type and parameters
- `seed`: Random seed for reproducibility
- `compress`: Gzip output files

CLI overrides (work with unified config):
- `--cells`: Override cell count
- `--fragments`: Override fragments per cell
- `--doublet-rate`: Override doublet rate
- `--ambient-contamination`: Override contamination rate
- `--output`: Override output directory

### Helper Commands
- `config`: Generate a template unified YAML configuration file
- `manifest`: Generate an example CSV manifest file

## Cell Type Weights

The `weight` column in the manifest controls cell type proportions:
- Weights are relative: `3:2:1:0.5` → approximately `43%:29%:14%:7%`
- Simple mode uses equal weights (1.0 for all)
- Higher weight = more cells of that type

## Doublet Simulation

Doublets are created by merging fragments from two cells:
- Uses Beta distribution (α=2, β=2) for realistic mixing ratios
- Not just 50/50 splits - can be 70/30, 60/40, etc.
- Outputs `doublet_annotations.tsv` for benchmarking doublet detectors

## Ambient Contamination

Simulates background DNA from lysed cells:
- Fragments sampled uniformly from total pool across all samples
- Realistic contamination rates typically 1-5%
- Tracked in cell metadata for analysis

## Output Files

- `fragments.tsv[.gz]` - 10x Genomics fragment format
- `fragments.bed[.gz]` - BED format
- `cell_metadata.tsv` - Per-cell statistics including doublet status
- `barcodes.tsv` - Cell barcodes
- `doublet_annotations.tsv` - Ground truth doublet labels
- `summary.txt` - Overall simulation statistics

## Example: Benchmarking Doublet Detection

```python
import pandas as pd

# Load ground truth
doublets = pd.read_csv('results/doublet_annotations.tsv', sep='\t')
metadata = pd.read_csv('results/cell_metadata.tsv', sep='\t')

# Your doublet detector predictions
predictions = your_doublet_detector('results/fragments.tsv.gz')

# Compare against ground truth
true_positives = predictions[predictions.barcode.isin(doublets.cell_id)]
accuracy = calculate_metrics(predictions, doublets)
```

## Configuration Files

### Unified Configuration (scatrs_config.yaml)

The `config` command generates a unified YAML configuration template:

```bash
# Generate template config
gtars scatrs config --output scatrs_config.yaml
```

This single config file is used for both staging and simulation:

```yaml
# scatrs_config.yaml
sample_table: samples.csv      # CSV manifest with samples
chrom_sizes: hg38.sizes        # Chromosome sizes
blacklist: blacklist.bed       # Optional blacklist

staging:
  output_dir: staged/
  extend_peaks: 250
  bin_size: 5000
  merge_distance: 20

# This field is added automatically after staging
staging_config: staged/stage_config.yaml

simulation:
  output_dir: results/
  cell_count: 1000
  fragments_per_cell: 10000
  signal_to_noise: 0.8
  doublet_rate: 0.1
  ambient_contamination: 0.02
  fragment_distribution:
    type: gaussian
    mean: 10000.0
    variance: 4000000.0
    min: 1000
    max: 50000
  seed: 42
  compress: false
```

### Stage Configuration (stage_config.yaml)

Created automatically by the `stage` command, this file contains:
- Staging metadata and statistics
- Pointers to staged data location
- Sample information from processing

```yaml
version: "1.0.0"
timestamp: "2024-01-15T10:30:00Z"
staged_dir: /path/to/staged
chrom_sizes_file: hg38.sizes
manifest_file: samples.csv

samples:
  - name: CD4_T
    weight: 3.0              # Original weight from CSV
    fragment_count: 1500000  # Calculated during staging
    peak_count: 25000        # Calculated during staging
  - name: CD8_T
    weight: 2.0
    fragment_count: 1300000
    peak_count: 23000

parameters:
  extend_peaks: 250
  bin_size: 5000
  merge_distance: 20
  blacklist_file: blacklist.bed

stats:
  total_fragments: 5000000
  total_peaks: 35000
  total_background_regions: 450000
  samples_processed: 4
```

### Sample Manifest (samples.csv)

The CSV manifest lists samples with their files and weights:

```csv
sample_name,bam_file,peaks_file,weight
CD4_T,/data/cd4t.bam,/data/cd4t_peaks.narrowPeak,3.0
CD8_T,/data/cd8t.bam,/data/cd8t_peaks.narrowPeak,2.0
B_cells,/data/bcells.bam,/data/bcells_peaks.narrowPeak,1.5
NK_cells,/data/nk.bam,,0.5
```

**Important**: You can modify the `weight` column and re-run `simulate` without re-staging!

## Architecture

- `models.rs` - Core data structures
- `io.rs` - File I/O operations  
- `staging.rs` - BAM processing and peak operations
- `sampling.rs` - Fragment sampling algorithms
- `cli.rs` - Command-line interface
- `consts.rs` - Constants and defaults

## Performance Tips

- Use `--compress` for large simulations to save disk space
- Staging is the slowest step but only needs to be done once
- Parallel processing uses all available cores automatically
- For >10M fragments, ensure 8GB+ RAM available

## Testing

```bash
cargo test scatrs::tests
```
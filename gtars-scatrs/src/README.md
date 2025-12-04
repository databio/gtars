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
# Specify which simulation to run by name (from the simulations array in config)
gtars scatrs simulate --config scatrs_config.yaml default
# The positional argument "default" is the simulation name

# To change simulation parameters without re-staging:
# 1. Edit the simulation section in scatrs_config.yaml
# 2. Or override with CLI args:
gtars scatrs simulate \
    --config scatrs_config.yaml \
    default \
    --cells 5000 \
    --doublet-rate 0.15 \
    --output results_v2/

# You can define multiple simulations in the config and run them separately:
# scatrs_config.yaml can contain multiple simulations with different parameters
gtars scatrs simulate --config scatrs_config.yaml high_doublets
gtars scatrs simulate --config scatrs_config.yaml low_contamination
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

# Multiple simulation configurations can be defined
simulations:
  - name: default              # Name used to select this simulation
    output_dir: results/
    cell_count: 1000
    signal_to_noise_distribution:
      type: fixed
      signal_to_noise: 0.8
    doublet_rate: 0.1
    ambient_contamination: 0.02
    fragment_distribution:
      type: gaussian
      mean_fragments_per_cell: 10000.0
      sd: 2000.0
      min: 1000
      max: 50000

  - name: high_doublets        # Alternative simulation with different parameters
    output_dir: results_high_doublets/
    cell_count: 1000
    signal_to_noise_distribution:
      type: fixed
      signal_to_noise: 0.8
    doublet_rate: 0.3          # Higher doublet rate
    ambient_contamination: 0.02
    fragment_distribution:
      type: gaussian
      mean_fragments_per_cell: 10000.0
      sd: 2000.0
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

#### Specifying Which Simulation to Run
When using unified config, you must specify which simulation to run by providing its name as a positional argument:
```bash
gtars scatrs simulate --config scatrs_config.yaml SIMULATION_NAME
```

The config file's `simulations` array can contain multiple simulation configurations, each with a unique `name` field. This allows you to:
- Run different simulation scenarios from the same staged data
- Compare different parameter combinations (e.g., varying doublet rates)
- Test different fragment distributions without re-staging

#### Simulation Parameters
Each simulation in the `simulations` array supports these parameters:
- `name`: Unique identifier for this simulation (required)
- `output_dir`: Where to save simulation results
- `cell_count`: Number of cells to simulate
- `signal_to_noise_distribution`: Peak vs background ratio distribution (see below)
- `doublet_rate`: Fraction of doublets (0-1)
- `ambient_contamination`: Fraction of ambient DNA (0-1)
- `fragment_distribution`: Distribution type and parameters (see below)
- `seed`: Random seed for reproducibility
- `compress`: Gzip output files
- `weight_column`: Which weight column to use from CSV (optional)

#### Signal-to-Noise Distribution Types
The `signal_to_noise_distribution` parameter controls the ratio of peak to background fragments for each cell:

**Fixed** (constant for all cells):
```yaml
signal_to_noise_distribution:
  type: fixed
  signal_to_noise: 0.8  # All cells have 80% signal, 20% background
```

**Beta Distribution** (recommended for realistic variation):
```yaml
signal_to_noise_distribution:
  type: beta
  mean: 0.8    # Mean signal-to-noise across cells
  sd: 0.1      # Standard deviation (controls spread)
```
- Ideal for values bounded between 0 and 1
- Variance must satisfy: sd² < mean × (1 - mean)
- Creates realistic heterogeneity in cell quality

**Gaussian Distribution** (with bounds):
```yaml
signal_to_noise_distribution:
  type: gaussian
  mean: 0.8
  sd: 0.1
  min: 0.0     # Optional, defaults to 0.0
  max: 1.0     # Optional, defaults to 1.0
```
- Values clamped to [0, 1] range
- Use when you want normal variation with hard limits

#### Fragment Distribution Types
The `fragment_distribution` parameter controls how many fragments each cell receives:

**Gaussian Distribution** (recommended for realistic data):
```yaml
fragment_distribution:
  type: gaussian
  mean_fragments_per_cell: 10000.0  # Mean of the distribution
  sd: 2000.0                        # Standard deviation
  min: 1000                         # Optional: minimum fragments per cell
  max: 50000                        # Optional: maximum fragments per cell
```
- Fragments per cell vary naturally following a normal distribution
- Most cells get values near the mean, with some variation
- Use `min`/`max` to prevent extreme outliers

**Negative Binomial Distribution** (most realistic for scATAC-seq):
```yaml
fragment_distribution:
  type: negativebinomial
  mean: 10000.0    # Mean fragments per cell
  sd: 4000.0       # Standard deviation (sd² must be > mean for overdispersion)
  min: 1000        # Optional: minimum fragments per cell
  max: 50000       # Optional: maximum fragments per cell
```
- Models the overdispersion commonly seen in real scATAC-seq data
- sd² > mean indicates overdispersion (more variability than Poisson)
- Better captures the long tail of high-fragment cells
- Use when simulating realistic biological variability

**Uniform Distribution** (random values in a range):
```yaml
fragment_distribution:
  type: uniform
  min: 8000    # Minimum fragments per cell
  max: 12000   # Maximum fragments per cell
```
- Each cell gets a random fragment count uniformly distributed between min and max
- All values in the range are equally likely
- Mean = (min + max) / 2, SD = (max - min) / √12

**Fixed Distribution** (constant fragment count):
```yaml
fragment_distribution:
  type: fixed
  fragments_per_cell: 10000  # Exact count for every cell
```
- Every cell gets exactly the same number of fragments
- Zero variation (SD = 0)
- Useful for controlled testing where you want identical cells

CLI overrides (work with unified config):
- `--cells`: Override cell count
- `--fragments`: Override fragments per cell
- `--doublet-rate`: Override doublet rate
- `--ambient-contamination`: Override contamination rate
- `--output`: Override output directory

### Helper Commands
- `config`: Generate a template unified YAML configuration file
- `manifest`: Generate an example CSV manifest file

## Multiple Simulations

The unified config system allows you to define multiple simulation configurations in the same YAML file. This is useful for:
- **Parameter sweeps**: Test different doublet rates, contamination levels, or fragment counts
- **Distribution comparisons**: Compare Gaussian vs Uniform vs Fixed fragment distributions
- **Cell count variations**: Generate datasets with different numbers of cells
- **Reproducibility**: All simulation parameters in one file for easy sharing

### Example: Multiple Simulations in Config

```yaml
simulations:
  - name: baseline
    output_dir: results_baseline/
    cell_count: 1000
    doublet_rate: 0.05
    ambient_contamination: 0.01
    fragment_distribution:
      type: gaussian
      mean_fragments_per_cell: 10000.0
      sd: 2000.0

  - name: high_quality
    output_dir: results_high_quality/
    cell_count: 2000
    doublet_rate: 0.02      # Lower doublet rate
    ambient_contamination: 0.005
    fragment_distribution:
      type: gaussian
      mean_fragments_per_cell: 15000.0
      sd: 2500.0

  - name: low_quality
    output_dir: results_low_quality/
    cell_count: 1000
    doublet_rate: 0.15      # Higher doublet rate
    ambient_contamination: 0.05
    fragment_distribution:
      type: gaussian
      mean_fragments_per_cell: 5000.0
      sd: 1500.0

  - name: uniform_fragments
    output_dir: results_uniform/
    cell_count: 1000
    doublet_rate: 0.1
    ambient_contamination: 0.02
    fragment_distribution:
      type: uniform        # Random values in range
      min: 8000
      max: 12000

  - name: fixed_fragments
    output_dir: results_fixed/
    cell_count: 1000
    doublet_rate: 0.1
    ambient_contamination: 0.02
    fragment_distribution:
      type: fixed          # Exact count for all cells
      fragments_per_cell: 10000

  - name: realistic_nb
    output_dir: results_nb/
    cell_count: 1000
    doublet_rate: 0.1
    ambient_contamination: 0.02
    fragment_distribution:
      type: negativebinomial  # Overdispersed like real data
      mean: 10000.0
      sd: 4000.0              # Variance > mean for overdispersion
```

### Running Multiple Simulations

```bash
# Stage once
gtars scatrs stage --config scatrs_config.yaml

# Run each simulation
gtars scatrs simulate --config scatrs_config.yaml baseline
gtars scatrs simulate --config scatrs_config.yaml high_quality
gtars scatrs simulate --config scatrs_config.yaml low_quality
gtars scatrs simulate --config scatrs_config.yaml uniform_fragments
gtars scatrs simulate --config scatrs_config.yaml fixed_fragments
gtars scatrs simulate --config scatrs_config.yaml realistic_nb

# Or run them in parallel
gtars scatrs simulate --config scatrs_config.yaml baseline &
gtars scatrs simulate --config scatrs_config.yaml high_quality &
wait
```

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

# Define one or more simulation configurations
simulations:
  - name: default
    output_dir: results/
    cell_count: 1000
    signal_to_noise_distribution:
      type: fixed
      signal_to_noise: 0.8
    doublet_rate: 0.1
    ambient_contamination: 0.02
    fragment_distribution:
      type: gaussian
      mean_fragments_per_cell: 10000.0
      sd: 2000.0
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
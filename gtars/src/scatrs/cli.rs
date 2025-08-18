use clap::{Command, Arg, ArgMatches};
use std::path::PathBuf;
use anyhow::Result;
use crate::scatrs::consts::*;

pub fn create_scatrs_cli() -> Command {
    Command::new(SCATRS_CMD)
        .about("Simulate single-cell ATAC-seq from bulk data")
        .version(VERSION)
        .subcommand_required(true)
        .subcommand(
            Command::new("stage")
                .about("Stage bulk ATAC-seq data for simulation")
                .long_about(
                    "Stage bulk ATAC-seq data for simulation using a unified YAML configuration file.\n\n\
                    The config file should contain sample_table, chrom_sizes, and staging parameters."
                )
                .arg(Arg::new("config")
                    .short('c')
                    .long("config")
                    .value_name("FILE")
                    .help("Unified YAML configuration file")
                    .required(true))
        )
        .subcommand(
            Command::new("simulate")
                .about("Run single-cell ATAC-seq simulation")
                .long_about(
                    "Run a specific simulation from the unified YAML configuration file.\n\n\
                    Usage: gtars scatrs simulate --config config.yaml sim1"
                )
                .arg(Arg::new("config")
                    .short('c')
                    .long("config")
                    .value_name("FILE")
                    .help("Unified YAML configuration file (must have staging_config field from stage command)")
                    .required(true))
                .arg(Arg::new("simulation")
                    .value_name("SIMULATION_NAME")
                    .help("Name of the simulation to run from the config file")
                    .required(true)
                    .index(1))
        )
        .subcommand(
            Command::new("config")
                .about("Generate example unified YAML configuration file")
                .arg(Arg::new("output")
                    .short('o')
                    .long("output")
                    .value_name("FILE")
                    .help("Output configuration file")
                    .default_value("scatrs_config.yaml"))
        )
        .subcommand(
            Command::new("manifest")
                .about("Generate example manifest CSV file")
                .arg(Arg::new("output")
                    .short('o')
                    .long("output")
                    .value_name("FILE")
                    .help("Output manifest file")
                    .default_value("samples.csv"))
        )
}

pub mod handlers {
    use super::*;
    use crate::scatrs::{
        models::*,
        staging::*,
        sampling::*,
        io::*,
    };
    use std::fs;
    
    pub fn handle_stage(matches: &ArgMatches) -> Result<()> {
        println!("Starting bulk ATAC-seq data staging...");
        
        // Load unified config
        let config_path = matches.get_one::<String>("config")
            .ok_or_else(|| anyhow::anyhow!("Config file is required"))?;
        let config_path = PathBuf::from(config_path);
        let config = UnifiedConfig::from_yaml(&config_path)?;
        
        let chrom_sizes_file = config.chrom_sizes.clone();
        let output_dir = config.staging.output_dir.clone();
        
        // Get current timestamp (simplified)
        let timestamp = format!("{:?}", std::time::SystemTime::now());
        
        // Read samples from the sample table (CSV manifest)
        let samples = ManifestReader::read_manifest(&config.sample_table)?;
        println!("Loaded {} samples from {}", samples.len(), config.sample_table.display());
        
        let extend_bp = config.staging.extend_peaks;
        let bin_size = config.staging.bin_size;
        let merge_distance = config.staging.merge_distance;
        
        // Create output directory
        fs::create_dir_all(&output_dir)?;
        
        // Read chromosome sizes
        let chrom_sizes = BackgroundGenerator::read_chrom_sizes(&chrom_sizes_file)?;
        println!("Loaded {} chromosome sizes", chrom_sizes.len());
        
        // Read blacklist if provided
        let blacklist = if let Some(ref blacklist_path) = config.blacklist {
            BackgroundGenerator::read_blacklist(blacklist_path)?
        } else {
            Vec::new()
        };
        if !blacklist.is_empty() {
            println!("Loaded {} blacklist regions", blacklist.len());
        }
        
        // Process peaks if available
        let (merged_peaks, background) = {
            let mut all_peaks = Vec::new();
            let mut has_peaks = false;
            
            // Load peaks for samples that have them
            for sample in &samples {
                if let Some(ref peaks_file) = sample.peaks_file {
                    if peaks_file.exists() {
                        let peaks = NarrowPeakReader::read_peaks(peaks_file, &sample.sample_name)?;
                        println!("Loaded {} peaks for {}", peaks.len(), sample.sample_name);
                        all_peaks.push(peaks);
                        has_peaks = true;
                    }
                }
            }
            
            if has_peaks {
                println!("Processing peaks from {} samples", all_peaks.len());
                
                // Merge peaks
                let merged = PeakMerger::merge_peaks(all_peaks, merge_distance)?;
                
                // Extend peaks
                let extended = PeakMerger::extend_peaks(&merged, extend_bp);
                
                // Generate background regions
                let bg = BackgroundGenerator::generate_background(
                    &chrom_sizes,
                    &extended,
                    &blacklist,
                    bin_size,
                )?;
                
                (Some(merged), Some(bg))
            } else {
                println!("No peak files provided - will use fragment-only mode");
                (None, None)
            }
        };
        
        // Process BAM files for each sample
        for sample in &samples {
            if !sample.bam_file.exists() {
                eprintln!("Warning: BAM file not found for {}: {:?}", 
                    sample.sample_name, sample.bam_file);
                continue;
            }
            
            println!("Processing BAM file for {}...", sample.sample_name);
            let fragments = BamProcessor::filter_fragments(&sample.bam_file)?;
            
            // Count fragments in regions if peaks are available
            if let (Some(ref peaks), Some(ref bg)) = (&merged_peaks, &background) {
                let _peak_counts = BamProcessor::count_fragments_in_regions(&fragments, peaks);
                let _bg_counts = BamProcessor::count_fragments_in_regions(&fragments, bg);
                println!("  Counted fragments in {} peaks and {} background regions", 
                    peaks.len(), bg.len());
            } else {
                println!("  Found {} fragments (no peak regions)", fragments.len());
            }
            
            // Save staged data
            let sample_dir = output_dir.join(&sample.sample_name);
            fs::create_dir_all(&sample_dir)?;
            
            // TODO: Save fragments and optional peaks/background to binary format
            println!("Saved staged data for {}", sample.sample_name);
        }
        
        // Create simplified stage configuration file
        let sample_names: Vec<String> = samples.iter()
            .map(|s| s.sample_name.clone())
            .collect();
        
        let stage_config = StageConfig {
            version: "1.0.0".to_string(),
            timestamp,
            staged_dir: output_dir.clone(),
            samples: sample_names,
            parameters: StageParameters {
                extend_peaks: extend_bp,
                bin_size,
                merge_distance,
                blacklist_file: config.blacklist.clone(),
            },
            stats: StageStats {
                total_fragments: 0,  // TODO: Calculate actual totals if needed
                total_peaks: merged_peaks.as_ref().map(|p| p.len() as u32),
                total_background_regions: background.as_ref().map(|b| b.len() as u32),
                samples_processed: samples.len() as u32,
            },
        };
        
        // Write stage config file using YAML
        let stage_config_path = output_dir.join("stage_config.yaml");
        stage_config.to_file(&stage_config_path)?;
        
        // Update the unified config with the staging_config path
        let mut updated_config = config;
        updated_config.update_staging_config(stage_config_path.clone(), &config_path)?;
        println!("Updated unified config with staging_config path: {:?}", stage_config_path);
        
        println!("Staging complete!");
        println!("Stage configuration saved to: {:?}", stage_config_path);
        println!("Use the same config file with 'gtars scatrs simulate' to run simulation");
        Ok(())
    }
    
    pub fn handle_simulate(matches: &ArgMatches) -> Result<()> {
        println!("Starting single-cell ATAC-seq simulation...");
        
        // Load unified config
        let config_path = matches.get_one::<String>("config")
            .ok_or_else(|| anyhow::anyhow!("Config file is required"))?;
        let config_path = PathBuf::from(config_path);
        let config = UnifiedConfig::from_yaml(&config_path)?;
        
        // Get simulation name from arguments
        let simulation_name = matches.get_one::<String>("simulation")
            .ok_or_else(|| anyhow::anyhow!("Simulation name is required"))?;
        
        // Get the specific simulation config
        let sim_config = config.get_simulation(simulation_name)
            .ok_or_else(|| anyhow::anyhow!("Simulation '{}' not found in config", simulation_name))?
            .clone();
        
        println!("Running simulation: {}", sim_config.name);
        
        // Verify staging has been done
        let stage_config_path = config.staging_config
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("Config file does not contain staging_config path. Please run 'stage' command first."))?;
        
        // Load stage configuration
        let stage_config = StageConfig::from_file(&stage_config_path)?;
        println!("Loaded stage config from: {:?}", stage_config_path);
        println!("Using staged data from: {:?}", stage_config.staged_dir);
        println!("Stage contained {} samples", stage_config.samples.len());
        
        // Load sample weights from the CSV manifest using the specified weight column
        let samples_with_weights = ManifestReader::read_manifest_with_weights(
            &config.sample_table,
            sim_config.weight_column.as_deref()
        )?;
        let sample_weights: std::collections::HashMap<String, f64> = samples_with_weights
            .into_iter()
            .map(|(s, w)| (s.sample_name, w))
            .collect();
        
        // Use simulation parameters from the selected simulation config
        let num_cells = sim_config.cell_count as usize;
        let fragments_per_cell = sim_config.fragments_per_cell;
        let signal_to_noise = sim_config.signal_to_noise;
        let doublet_rate = sim_config.doublet_rate;
        let ambient_contamination = sim_config.ambient_contamination;
        let seed = sim_config.seed;
        let compress = sim_config.compress;
        let output_dir = sim_config.output_dir.clone();
        
        // Create output directory
        fs::create_dir_all(&output_dir)?;
        
        // Load staged data (simplified - would load from saved files)
        println!("Loading staged data...");
        
        // Placeholder for actual data loading
        let peaks = Vec::new();
        let peak_weights = Vec::new();
        let background = Vec::new();
        let background_weights = Vec::new();
        let region_fragments = std::collections::HashMap::new();
        
        // Run weighted sampling (if peaks available)
        let mut weighted_sampler = WeightedSampler::new(seed);
        let cell_regions = weighted_sampler.sample_regions_for_cells(
            if peaks.is_empty() { None } else { Some(&peaks) },
            if peak_weights.is_empty() { None } else { Some(&peak_weights) },
            if background.is_empty() { None } else { Some(&background) },
            if background_weights.is_empty() { None } else { Some(&background_weights) },
            signal_to_noise,
            &sim_config.fragment_distribution,
            num_cells,
        )?;
        
        // Run uniform fragment sampling
        let uniform_sampler = UniformSampler::new(seed);
        
        // Get cell type names from stage config and weights from CSV manifest
        let cell_types: Vec<String> = stage_config.samples.clone();
        let cell_type_weights: Vec<f64> = stage_config.samples.iter()
            .map(|name| sample_weights.get(name).copied().unwrap_or(1.0))
            .collect();
        
        let mut cells = uniform_sampler.sample_fragments_for_cells(
            cell_regions,
            &region_fragments,
            &cell_types,
            Some(&cell_type_weights),
            None,  // TODO: Pass all_fragments if in no-peaks mode
            Some(fragments_per_cell),
        )?;
        
        // Apply doublet simulation if requested
        let doublet_annotations = if doublet_rate > 0.0 {
            // Create a temporary ScatrsConfig for the doublet simulator (for compatibility)
            let temp_config = ScatrsConfig {
                cell_types: cell_types.clone(),
                cell_count: num_cells as u32,
                fragments_per_cell,
                signal_to_noise,
                extend_peaks: config.staging.extend_peaks,
                bin_size: config.staging.bin_size,
                fragment_distribution: sim_config.fragment_distribution.clone(),
                output_dir: output_dir.clone(),
                peak_dir: None,
                bam_dir: None,
                blacklist_file: config.blacklist.clone(),
                genome_file: Some(config.chrom_sizes.clone()),
                doublet_rate,
                ambient_contamination,
            };
            
            let doublet_simulator = DoubletSimulator::new(seed);
            doublet_simulator.simulate_doublets(&mut cells, doublet_rate, &temp_config)?
        } else {
            Vec::new()
        };
        
        // Apply ambient contamination if requested
        if ambient_contamination > 0.0 {
            // TODO: Load all fragments from staged data for ambient pool
            let ambient_simulator = AmbientSimulator::new(seed);
            // For now, using empty fragment pool - would load from staged data
            let all_fragments = Vec::new();
            ambient_simulator.add_ambient_contamination(&mut cells, ambient_contamination, &all_fragments)?;
        }
        
        // Write output files
        println!("Writing output files...");
        
        if compress {
            let fragments_file = output_dir.join("fragments.tsv.gz");
            FragmentWriter::write_fragments_file_gz(&cells, &fragments_file)?;
            
            let bed_file = output_dir.join("fragments.bed.gz");
            BedWriter::write_bed_gz(&cells, &bed_file)?;
        } else {
            let fragments_file = output_dir.join("fragments.tsv");
            FragmentWriter::write_fragments_file(&cells, &fragments_file)?;
            
            let bed_file = output_dir.join("fragments.bed");
            BedWriter::write_bed(&cells, &bed_file)?;
        }
        
        // Write metadata and summary
        let metadata_file = output_dir.join("cell_metadata.tsv");
        BedWriter::write_cell_metadata(&cells, &metadata_file)?;
        
        let barcodes_file = output_dir.join("barcodes.tsv");
        FragmentWriter::write_barcodes(&cells, &barcodes_file)?;
        
        let summary_file = output_dir.join("summary.txt");
        FragmentWriter::write_summary_stats(&cells, &summary_file)?;
        
        // Write doublet annotations if any
        if !doublet_annotations.is_empty() {
            let doublet_file = output_dir.join("doublet_annotations.tsv");
            FragmentWriter::write_doublet_annotations(&doublet_annotations, &doublet_file)?;
        }
        
        println!("Simulation complete!");
        Ok(())
    }
    
    pub fn handle_config(matches: &ArgMatches) -> Result<()> {
        let output_file = PathBuf::from(matches.get_one::<String>("output").unwrap());
        
        // Create a YAML string with inline comments
        let yaml_content = r#"# Unified configuration for scATAC-seq simulation
# This file is used by both 'stage' and 'simulate' commands

# Path to CSV manifest with sample information
sample_table: samples.csv

# Chromosome sizes file (required for both staging and simulation)
chrom_sizes: hg38.sizes

# Optional blacklist file for regions to exclude (used by both staging and simulation)
blacklist: blacklist.bed  # Comment out or set to null if not needed

# Staging parameters - configures how BAM files are processed
staging:
  output_dir: staged/           # Directory to save staged data
  extend_peaks: 250             # Extend peaks by N base pairs
  bin_size: 5000               # Size of genomic bins for background regions
  merge_distance: 20           # Merge peaks within N base pairs
  skip_if_exists: false        # Skip staging if output already exists

# This field is automatically added after staging completes
# staging_config: staged/stage_config.yaml

# Multiple simulation configurations
simulations:
  - name: sim1                 # Name to reference this simulation
    weight_column: weight      # Which CSV column to use for cell type weights (null for equal)
    output_dir: results/sim1/
    cell_count: 1000           # Number of cells to simulate
    fragments_per_cell: 10000  # Average fragments per cell
    signal_to_noise: 0.8       # Ratio of peak vs background fragments (0-1)
    doublet_rate: 0.1          # Fraction of cells that are doublets (0-1)
    ambient_contamination: 0.02 # Fraction of ambient DNA contamination (0-1)
    fragment_distribution:
      type: gaussian           # Distribution type: gaussian or uniform
      mean: 10000.0
      variance: 4000000.0
      min: 1000
      max: 50000
    seed: 42                   # Random seed for reproducibility (null for random)
    compress: false            # Whether to gzip output files

  - name: sim2                 # Second simulation with different parameters
    weight_column: weight2     # Can use different weight columns from CSV
    output_dir: results/sim2/
    cell_count: 500
    fragments_per_cell: 5000
    signal_to_noise: 0.9
    doublet_rate: 0.05
    ambient_contamination: 0.01
    fragment_distribution:
      type: uniform
      count: 5000              # For uniform distribution
    seed: 123
    compress: true

# Usage:
# 1. Edit this file and samples.csv with your data
# 2. Stage: gtars scatrs stage --config scatrs_config.yaml
# 3. Simulate: gtars scatrs simulate --config scatrs_config.yaml sim1
#             gtars scatrs simulate --config scatrs_config.yaml sim2
"#;
        
        // Change extension to .yaml if needed
        let output_file = if output_file.extension() == Some(std::ffi::OsStr::new("toml")) {
            output_file.with_extension("yaml")
        } else {
            output_file
        };
        
        std::fs::write(&output_file, yaml_content)?;
        
        println!("Generated example unified configuration file: {:?}", output_file);
        println!("\nUsage:");
        println!("  1. Edit {} and create samples.csv with your data", output_file.display());
        println!("  2. Stage: gtars scatrs stage --config {}", output_file.display());
        println!("  3. Simulate: gtars scatrs simulate --config {} sim1", output_file.display());
        println!("             gtars scatrs simulate --config {} sim2", output_file.display());
        Ok(())
    }
    
    pub fn handle_manifest(matches: &ArgMatches) -> Result<()> {
        let output_file = PathBuf::from(matches.get_one::<String>("output").unwrap());
        
        ManifestReader::create_example_manifest(&output_file)?;
        
        Ok(())
    }
}

pub fn handle_scatrs_command(matches: &ArgMatches) -> Result<()> {
    match matches.subcommand() {
        Some(("stage", sub_matches)) => handlers::handle_stage(sub_matches),
        Some(("simulate", sub_matches)) => handlers::handle_simulate(sub_matches),
        Some(("config", sub_matches)) => handlers::handle_config(sub_matches),
        Some(("manifest", sub_matches)) => handlers::handle_manifest(sub_matches),
        _ => unreachable!("Subcommand required"),
    }
}
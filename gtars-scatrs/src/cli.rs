use clap::{Command, Arg, ArgMatches};
use std::path::PathBuf;
use anyhow::Result;
use crate::consts::*;

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
                .arg(Arg::new("threads")
                    .short('t')
                    .long("threads")
                    .value_name("NUMBER")
                    .help("Number of threads for parallel region processing (default: 4)")
                    .value_parser(clap::value_parser!(usize))
                    .required(false))
                .arg(Arg::new("cache_fragments")
                    .long("cache-fragments")
                    .help("Cache fragments to binary files for faster simulation (default: true)")
                    .action(clap::ArgAction::SetTrue)
                    .default_value("true"))
                .arg(Arg::new("no_cache")
                    .long("no-cache")
                    .help("Disable fragment caching")
                    .action(clap::ArgAction::SetTrue)
                    .conflicts_with("cache_fragments"))
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
                .arg(Arg::new("threads")
                    .short('t')
                    .long("threads")
                    .value_name("NUMBER")
                    .help("Number of threads for parallel region processing (default: 4)")
                    .value_parser(clap::value_parser!(usize))
                    .required(false))
                .arg(Arg::new("debug_regions")
                    .long("debug-regions")
                    .help("Enable debug mode to track region sampling and verify cache integrity")
                    .action(clap::ArgAction::SetTrue))
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
    use crate::{
        models::*,
        staging::*,
        sampling::*,
        io::*,
    };
    use std::fs;
    use std::collections::{HashMap, HashSet};
    
    pub fn handle_stage(matches: &ArgMatches) -> Result<()> {
        println!("Starting bulk ATAC-seq data staging...");
        let overall_start_time = std::time::Instant::now();
        
        // Configure thread pool
        let thread_count = get_thread_count(matches.get_one::<usize>("threads").copied());
        println!("Using {} threads for parallel region processing", thread_count);
        
        rayon::ThreadPoolBuilder::new()
            .num_threads(thread_count)
            .build_global()
            .expect("Failed to configure thread pool");
        
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
        println!("Processing {} samples sequentially (optimal for I/O)...", samples.len());
        
        let extend_bp = config.staging.extend_peaks;
        let bin_size = config.staging.bin_size;
        let merge_distance = config.staging.merge_distance;
        
        // Create output directory
        fs::create_dir_all(&output_dir)
            .map_err(|e| anyhow::anyhow!("Failed to create staging output directory {:?}: {}", output_dir, e))?;
        
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
            
            let start_time = std::time::Instant::now();
            let sample_dir = output_dir.join(&sample.sample_name);
            
            // Clean old chromosome cache files to avoid ID scheme conflicts
            if sample_dir.exists() {
                for entry in fs::read_dir(&sample_dir)? {
                    if let Ok(entry) = entry {
                        let path = entry.path();
                        if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
                            if name.ends_with("_fragments.bin") {
                                println!("  Removing old cache file: {}", name);
                                fs::remove_file(path)?;
                            }
                        }
                    }
                }
            }
            
            fs::create_dir_all(&sample_dir)
                .map_err(|e| anyhow::anyhow!("Failed to create sample directory {:?}: {}", sample_dir, e))?;
            
            // Stream and count fragments in regions
            if let (Some(ref peaks), Some(ref bg)) = (&merged_peaks, &background) {
                // Determine if fragment caching should be enabled
                let cache_fragments = !matches.get_flag("no_cache");
                let cache_dir = if cache_fragments {
                    Some(sample_dir.as_path())
                } else {
                    None
                };
                
                // Use chromosome-based caching by default (dramatically reduces file count)
                let (peak_counts, bg_counts) = BamProcessor::count_fragments_streaming_with_chromosome_cache(
                    &sample.bam_file,
                    peaks,
                    bg,
                    cache_dir
                )?;
                
                // Save the counts
                StagedDataWriter::save_counts(
                    &peak_counts,
                    &sample_dir.join("peak_counts.bin")
                )?;
                StagedDataWriter::save_counts(
                    &bg_counts,
                    &sample_dir.join("bg_counts.bin")
                )?;
                
                // Calculate statistics
                let fragments_in_peaks: u64 = peak_counts.iter().map(|&c| c as u64).sum();
                let fragments_in_bg: u64 = bg_counts.iter().map(|&c| c as u64).sum();
                
                // Save metadata
                let metadata = SampleMetadata {
                    sample_name: sample.sample_name.clone(),
                    bam_file: sample.bam_file.clone(),
                    total_fragments: fragments_in_peaks + fragments_in_bg,
                    fragments_in_peaks: Some(fragments_in_peaks),
                    fragments_in_background: Some(fragments_in_bg),
                    processing_time_secs: start_time.elapsed().as_secs_f64(),
                    timestamp: timestamp.clone(),
                };
                StagedDataWriter::save_sample_metadata(
                    &metadata,
                    &sample_dir.join("metadata.yaml")
                )?;
                
                println!("  Saved {} peak counts and {} background counts",
                    peak_counts.len(), bg_counts.len());
                println!("  Fragments: {} in peaks, {} in background",
                    fragments_in_peaks, fragments_in_bg);
                
                // Validate staging results
                if fragments_in_peaks == 0 && fragments_in_bg == 0 {
                    eprintln!("  ⚠ Warning: No fragments found in any regions for {}!", sample.sample_name);
                    eprintln!("    This may indicate:");
                    eprintln!("    1. BAM file has no valid fragments");
                    eprintln!("    2. Regions don't overlap with fragments");
                    eprintln!("    3. Chromosome naming mismatch (e.g., 'chr1' vs '1')");
                }
                
                let elapsed = start_time.elapsed();
                println!("  Processing time: {:.1}s", elapsed.as_secs_f64());
                    
            } else {
                // Fragment-only mode (no peaks)
                let total = BamProcessor::count_total_fragments(&sample.bam_file)?;
                
                let metadata = SampleMetadata {
                    sample_name: sample.sample_name.clone(),
                    bam_file: sample.bam_file.clone(),
                    total_fragments: total,
                    fragments_in_peaks: None,
                    fragments_in_background: None,
                    processing_time_secs: start_time.elapsed().as_secs_f64(),
                    timestamp: timestamp.clone(),
                };
                StagedDataWriter::save_sample_metadata(
                    &metadata,
                    &sample_dir.join("metadata.yaml")
                )?;
                
                println!("  Found {} fragments (no peak regions)", total);
            }
            
            println!("  Saved staged data to {:?}", sample_dir);
        }
        
        // Save global regions
        if let Some(ref peaks) = merged_peaks {
            StagedDataWriter::save_regions(
                peaks,
                &output_dir.join("merged_peaks.bin")
            )?;
            println!("Saved {} merged peaks", peaks.len());
        }
        
        if let Some(ref bg) = background {
            StagedDataWriter::save_regions(
                bg,
                &output_dir.join("background_regions.bin")
            )?;
            println!("Saved {} background regions", bg.len());
        }
        
        // Create simplified stage configuration file
        let sample_names: Vec<String> = samples.iter()
            .map(|s| s.sample_name.clone())
            .collect();
        
        // Calculate total fragments from all samples
        let mut total_fragments = 0u64;
        for sample in &samples {
            let sample_dir = output_dir.join(&sample.sample_name);
            let metadata_path = sample_dir.join("metadata.yaml");
            if metadata_path.exists() {
                if let Ok(metadata) = StagedDataReader::load_sample_metadata(&metadata_path) {
                    total_fragments += metadata.total_fragments;
                }
            }
        }
        
        let stage_config = StageConfig {
            version: "1.0.0".to_string(),
            timestamp,
            samples: sample_names,
            parameters: StageParameters {
                extend_peaks: extend_bp,
                bin_size,
                merge_distance,
                blacklist_file: config.blacklist.clone(),
            },
            stats: StageStats {
                total_fragments,
                total_peaks: merged_peaks.as_ref().map(|p| p.len() as u32),
                total_background_regions: background.as_ref().map(|b| b.len() as u32),
                samples_processed: samples.len() as u32,
            },
        };
        
        // Write stage config file using YAML
        let stage_config_path = output_dir.join("stage_config.yaml");
        stage_config.to_file(&stage_config_path)?;
        
        // Update the unified config with the staging_config path relative to config file
        let mut updated_config = config;
        
        // Make the stage_config_path relative to the config file's directory
        let config_dir = config_path.parent()
            .ok_or_else(|| anyhow::anyhow!("Config file path has no parent directory"))?;
        let relative_stage_path = pathdiff::diff_paths(&stage_config_path, config_dir)
            .unwrap_or_else(|| stage_config_path.clone());
        
        updated_config.update_staging_config(relative_stage_path.clone(), &config_path)?;
        println!("Updated unified config with staging_config path: {:?}", relative_stage_path);
        
        // Validate overall staging results
        let mut staging_warnings = Vec::new();
        let mut cached_region_count = 0;
        
        for sample_name in &samples.iter().map(|s| &s.sample_name).collect::<Vec<_>>() {
            let sample_dir = output_dir.join(sample_name);
            
            // Check if chromosome cache was created
            if BamProcessor::has_chromosome_cache(&sample_dir) {
                cached_region_count += 1;
            } else if BamProcessor::has_fragment_cache(&sample_dir) {
                staging_warnings.push(format!("{}: Using old per-bin cache format (many files)", sample_name));
            } else {
                staging_warnings.push(format!("{}: No fragment cache created", sample_name));
            }
        }
        
        let total_elapsed = overall_start_time.elapsed();
        let minutes = total_elapsed.as_secs() / 60;
        let seconds = total_elapsed.as_secs() % 60;
        
        println!("\n══════════════════════════════════════════════════════");
        println!("Staging complete! Total time: {}m {}s", minutes, seconds);
        println!("══════════════════════════════════════════════════════");
        
        // Report staging summary
        if let Some(ref peaks) = merged_peaks {
            println!("Regions: {} peaks, {} background regions", 
                peaks.len(), background.as_ref().map_or(0, |b| b.len()));
        }
        println!("Samples: {} processed, {} with cached fragments", 
            samples.len(), cached_region_count);
        
        if !staging_warnings.is_empty() {
            println!("\n⚠ Staging warnings:");
            for warning in staging_warnings.iter().take(5) {
                println!("  - {}", warning);
            }
            if staging_warnings.len() > 5 {
                println!("  ... and {} more", staging_warnings.len() - 5);
            }
        }
        
        println!("\nStage configuration saved to: {:?}", stage_config_path);
        println!("Use the same config file with 'gtars scatrs simulate' to run simulation");
        Ok(())
    }
    
    fn collect_sampled_regions(cell_regions: &[Vec<SampledRegion>]) -> Vec<ScatrsRegion> {
        let mut regions = Vec::new();
        let mut seen = std::collections::HashSet::new();
        
        println!("  Collecting unique regions from {} cells...", cell_regions.len());
        let mut processed_cells = 0;
        
        for cell in cell_regions {
            for sampled_region in cell {
                let region_key = format!(
                    "{}:{}-{}", 
                    sampled_region.region.chrom,
                    sampled_region.region.start,
                    sampled_region.region.end
                );
                if seen.insert(region_key) {
                    // insert returns true if the value was not already present
                    regions.push(sampled_region.region.clone());
                }
            }
            
            processed_cells += 1;
            if processed_cells % 500 == 0 {
                println!("    Processed {}/{} cells, found {} unique regions so far...", 
                    processed_cells, cell_regions.len(), regions.len());
            }
        }
        
        println!("  Found {} unique regions total", regions.len());
        regions
    }
    
    pub fn handle_simulate(matches: &ArgMatches) -> Result<()> {
        println!("\n╔══════════════════════════════════════════════════════╗");
        println!("║     Starting single-cell ATAC-seq simulation...     ║");
        println!("╚══════════════════════════════════════════════════════╝");
        
        // Configure thread pool
        let thread_count = get_thread_count(matches.get_one::<usize>("threads").copied());
        println!("Using {} threads for parallel region processing", thread_count);
        
        rayon::ThreadPoolBuilder::new()
            .num_threads(thread_count)
            .build_global()
            .expect("Failed to configure thread pool");
        
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
            .ok_or_else(|| {
                let available_sims: Vec<String> = config.simulations.iter()
                    .map(|s| s.name.clone())
                    .collect();
                if available_sims.is_empty() {
                    anyhow::anyhow!(
                        "Simulation '{}' not found. No simulations are defined in the config file.\n\
                        Please add a simulation configuration to your YAML file.",
                        simulation_name
                    )
                } else {
                    anyhow::anyhow!(
                        "Simulation '{}' not found in config.\n\
                        Available simulations: {}\n\
                        Please use one of the available simulation names.",
                        simulation_name,
                        available_sims.join(", ")
                    )
                }
            })?
            .clone();
        
        println!("Running simulation: {}", sim_config.name);
        
        // Check if staging has been done
        let (staged_dir, stage_config) = if let Some(ref stage_config_path) = config.staging_config {
            // Load from staged data
            let stage_config = StageConfig::from_file(&stage_config_path)?;
            println!("Loaded stage config from: {:?}", stage_config_path);
            
            let staged_dir = stage_config_path.parent()
                .ok_or_else(|| anyhow::anyhow!("Stage config path has no parent directory"))?
                .to_path_buf();
            println!("Using staged data from: {:?}", staged_dir);
            println!("Stage contained {} samples", stage_config.samples.len());
            
            (Some(staged_dir), Some(stage_config))
        } else {
            // Staging config is missing - this will be handled below
            (None, None)
        };
        
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
        let signal_to_noise = sim_config.signal_to_noise;
        let doublet_rate = sim_config.doublet_rate;
        let ambient_contamination = sim_config.ambient_contamination;
        let seed = sim_config.seed;
        let compress = sim_config.compress;
        let output_dir = sim_config.output_dir.clone();
        
        // Create output directory
        fs::create_dir_all(&output_dir)
            .map_err(|e| anyhow::anyhow!("Failed to create simulation output directory {:?}: {}", output_dir, e))?;
        
        // Load or prepare data based on whether staging was done
        let (peaks, background, sample_peak_weights, sample_bg_weights, _bam_paths, samples_list) = 
            if let Some(staged_dir_path) = &staged_dir {
                // Load from staged data
                println!("\n📂 Loading staged data...");
                
                // Load global regions if available
                let peaks = if staged_dir_path.join("merged_peaks.bin").exists() {
                    let peaks = StagedDataReader::load_regions(&staged_dir_path.join("merged_peaks.bin"))?;
                    println!("  ✓ Loaded {} peak regions", peaks.len());
                    peaks
                } else {
                    println!("  ℹ No peak regions found");
                    Vec::new()
                };
                
                let background = if staged_dir_path.join("background_regions.bin").exists() {
                    let bg = StagedDataReader::load_regions(&staged_dir_path.join("background_regions.bin"))?;
                    println!("  ✓ Loaded {} background regions", bg.len());
                    bg
                } else {
                    println!("  ℹ No background regions found");
                    Vec::new()
                };
                
                // Load sample-specific data and store weighted counts per sample
                let mut sample_peak_weights: HashMap<String, Vec<f64>> = HashMap::new();
                let mut sample_bg_weights: HashMap<String, Vec<f64>> = HashMap::new();
                let mut bam_paths = HashMap::new();
                
                let stage_config = stage_config.as_ref().unwrap();
                for sample_name in &stage_config.samples {
                    let sample_dir = staged_dir_path.join(sample_name);
                    
                    // Load metadata
                    if sample_dir.join("metadata.yaml").exists() {
                        let metadata = StagedDataReader::load_sample_metadata(
                            &sample_dir.join("metadata.yaml")
                        )?;
                        
                        // Store BAM path for later fragment extraction
                        bam_paths.insert(sample_name.clone(), metadata.bam_file);
                        
                        // Load counts if available and store per-sample
                        if sample_dir.join("peak_counts.bin").exists() {
                            let peak_counts = StagedDataReader::load_counts(
                                &sample_dir.join("peak_counts.bin")
                            )?;
                            
                            // DO NOT apply sample weight to counts - weights are only for cell proportions!
                            // The counts already represent the biological signal
                            let weighted_peaks: Vec<f64> = peak_counts.iter()
                                .map(|&count| count as f64)
                                .collect();
                            sample_peak_weights.insert(sample_name.clone(), weighted_peaks);
                        }
                        
                        if sample_dir.join("bg_counts.bin").exists() {
                            let bg_counts = StagedDataReader::load_counts(
                                &sample_dir.join("bg_counts.bin")
                            )?;
                            
                            // DO NOT apply sample weight to counts - weights are only for cell proportions!
                            let weighted_bg: Vec<f64> = bg_counts.iter()
                                .map(|&count| count as f64)
                                .collect();
                            sample_bg_weights.insert(sample_name.clone(), weighted_bg);
                        }
                    }
                }
                
                (peaks, background, sample_peak_weights, sample_bg_weights, bam_paths, stage_config.samples.clone())
            } else {
                // Staging is mandatory - return error with clear instructions
                return Err(anyhow::anyhow!(
                    "Staging has not been completed!\n\n\
                    The SCATRS module requires data to be staged before simulation.\n\
                    This ensures efficient memory usage and fast fragment access.\n\n\
                    To fix this:\n\
                    1. First run the staging command:\n   \
                       gtars scatrs stage --config {}\n\n\
                    2. Then run your simulation:\n   \
                       gtars scatrs simulate --config {} {}\n\n\
                    The staging step will:\n\
                    - Process BAM files in a streaming fashion\n\
                    - Create efficient fragment caches\n\
                    - Generate region indices\n\
                    - Save the staging_config path to your config file\n\n\
                    After staging, the simulate command will use the staged data.",
                    config_path.display(),
                    config_path.display(), 
                    simulation_name
                ));
            };
        
        // region_fragments is no longer used - we use fragments_by_type instead
        
        // Fragment extraction will be done after region sampling to only extract needed regions
        
        // Pre-assign cell types based on sample weights
        let samples_with_weights: Vec<(String, f64)> = samples_list.iter()
            .map(|name| (name.clone(), sample_weights.get(name).copied().unwrap_or(1.0)))
            .collect();
        
        let cell_type_assignments = assign_cell_types_by_weight(
            &samples_with_weights,
            num_cells,
            seed
        );
        
        // Print cell type distribution
        let mut type_counts: HashMap<String, usize> = HashMap::new();
        for cell_type in &cell_type_assignments {
            *type_counts.entry(cell_type.clone()).or_insert(0) += 1;
        }
        println!("\n📊 Cell type distribution:");
        for (cell_type, count) in &type_counts {
            let percentage = (*count as f64 / num_cells as f64) * 100.0;
            println!("  {}: {} cells ({:.1}%)", cell_type, count, percentage);
        }
        
        // Run weighted sampling with per-sample weights
        println!("\n🎲 Sampling regions for {} cells with cell-type-specific patterns...", num_cells);
        
        // Display ambient contamination calculation if applicable
        if ambient_contamination > 0.0 {
            let target_fragments = match &sim_config.fragment_distribution {
                crate::models::FragmentDistribution::Uniform { fragments_per_cell } => *fragments_per_cell,
                crate::models::FragmentDistribution::Gaussian { mean_fragments_per_cell, .. } => *mean_fragments_per_cell as u32,
            };
            let source_fragments = ((target_fragments as f64) * (1.0 - ambient_contamination)) as u32;
            let ambient_fragments = ((target_fragments as f64) * ambient_contamination) as u32;
            println!("  With {:.1}% ambient contamination: {} source + {} ambient = {} total fragments per cell",
                ambient_contamination * 100.0, source_fragments, ambient_fragments, target_fragments);
        }
        
        let mut weighted_sampler = WeightedSampler::new(seed);
        let cell_regions = weighted_sampler.sample_regions_for_cells_by_type(
            if peaks.is_empty() { None } else { Some(&peaks) },
            if sample_peak_weights.is_empty() { None } else { Some(&sample_peak_weights) },
            if background.is_empty() { None } else { Some(&background) },
            if sample_bg_weights.is_empty() { None } else { Some(&sample_bg_weights) },
            signal_to_noise,
            &sim_config.fragment_distribution,
            &cell_type_assignments,
            ambient_contamination,
        )?;
        
        // Extract real fragments from staged data
        println!("\n🔍 Extracting real fragments from staged data...");
        
        // Get unique regions from sampled cells (with their actual region indices)
        let sampled_regions = collect_sampled_regions(&cell_regions);
        println!("  Need fragments for {} unique regions across {} cells", 
            sampled_regions.len(), cell_regions.len());
        
        // Debug mode: Show first few sampled regions
        let debug_regions = matches.get_flag("debug_regions");
        if debug_regions {
            println!("\n🔍 DEBUG: First 10 sampled regions:");
            for (i, region) in sampled_regions.iter().take(10).enumerate() {
                println!("    {}: {}:{}-{}", i + 1, region.chrom, region.start, region.end);
            }
        }
        
        // Build HashMaps for O(1) lookups instead of O(n) linear searches
        println!("  Building region lookup maps for fast indexing...");
        let peak_lookup: std::collections::HashMap<String, usize> = peaks
            .iter()
            .enumerate()
            .map(|(idx, r)| (format!("{}:{}-{}", r.chrom, r.start, r.end), idx))
            .collect();
        
        let bg_lookup: std::collections::HashMap<String, usize> = background
            .iter()
            .enumerate()
            .map(|(idx, r)| (format!("{}:{}-{}", r.chrom, r.start, r.end), idx))
            .collect();
        
        println!("  Built lookup maps: {} peaks, {} background regions", 
            peak_lookup.len(), bg_lookup.len());

        // Group cells by their assigned type and collect region indices per type
        // CRITICAL: This grouping ensures we load fragments separately for each cell type
        // preventing cross-contamination between different cell types' fragments
        println!("  Grouping cells by type and collecting region indices...");
        let mut cells_by_type: HashMap<String, Vec<usize>> = HashMap::new();
        let mut regions_by_type: HashMap<String, (HashSet<usize>, HashSet<usize>)> = HashMap::new();
        
        for (cell_idx, cell_regions_item) in cell_regions.iter().enumerate() {
            let cell_type = &cell_type_assignments[cell_idx];
            
            // Track which cells belong to each type
            cells_by_type.entry(cell_type.clone())
                .or_insert_with(Vec::new)
                .push(cell_idx);
            
            // Collect regions for this cell type
            let (peak_indices, bg_indices) = regions_by_type.entry(cell_type.clone())
                .or_insert_with(|| (HashSet::new(), HashSet::new()));
            
            for sampled_region in cell_regions_item {
                let region_key = format!("{}:{}-{}", 
                    sampled_region.region.chrom,
                    sampled_region.region.start,
                    sampled_region.region.end
                );
                
                let region_idx = if sampled_region.is_peak {
                    peak_lookup.get(&region_key).copied()
                } else {
                    bg_lookup.get(&region_key).copied()
                };
                
                if let Some(idx) = region_idx {
                    if sampled_region.is_peak {
                        peak_indices.insert(idx);
                    } else {
                        bg_indices.insert(idx);
                    }
                }
            }
        }
        
        // Print cell type grouping info
        for (cell_type, cell_indices) in &cells_by_type {
            let (peak_set, bg_set) = &regions_by_type[cell_type];
            println!("    {}: {} cells, {} peak regions, {} background regions",
                cell_type, cell_indices.len(), peak_set.len(), bg_set.len());
        }
        
        // Load fragments separately for each cell type
        // FIX: This is the critical fix for cell-type assignment bug
        // Previously, all fragments were loaded into a single pool, causing cross-contamination
        // Now we maintain separate fragment pools for each cell type
        println!("\n📂 Loading fragments per cell type (preserving cell type specificity)...");
        let mut fragments_by_type: HashMap<String, HashMap<String, Vec<crate::models::ScatrsFragment>>> = HashMap::new();
        
        for (cell_type, _) in &cells_by_type {
            println!("  Processing fragments for {} cells...", cell_type);
            
            let (type_peak_indices, type_bg_indices) = &regions_by_type[cell_type];
            let peak_indices_vec: Vec<usize> = type_peak_indices.iter().copied().collect();
            let bg_indices_vec: Vec<usize> = type_bg_indices.iter().copied().collect();
            
            println!("    Need fragments from {} peaks and {} background regions", 
                peak_indices_vec.len(), bg_indices_vec.len());
            
            // Debug mode: Show region indices that will be loaded
            if debug_regions {
                println!("\n🔍 DEBUG: Region indices for {}:", cell_type);
                println!("    Peak indices (first 10): {:?}", peak_indices_vec.iter().take(10).collect::<Vec<_>>());
                println!("    Background indices (first 10): {:?}", bg_indices_vec.iter().take(10).collect::<Vec<_>>());
            }
            
            // Load fragments ONLY from the matching sample
            // KEY: Each cell type loads fragments from its own BAM file only
            // This prevents cross-contamination between different cell types
            let sample_dir = staged_dir.as_ref().expect("staged_dir should be set").join(cell_type);
            
            let mut type_fragments: HashMap<String, Vec<crate::models::ScatrsFragment>> = HashMap::new();
            
            if BamProcessor::has_chromosome_cache(&sample_dir) {
                println!("    Loading from chromosome cache for {}...", cell_type);
                
                let chr_fragments = BamProcessor::extract_fragments_from_chromosome_cache(
                    &sample_dir,
                    &peak_indices_vec,
                    &bg_indices_vec,
                    &peaks,
                    &background,
                )?;
                
                // Store fragments for this type
                for (region_key, frags) in chr_fragments {
                    type_fragments.entry(region_key)
                        .or_insert_with(Vec::new)
                        .extend(frags);
                }
            } else if BamProcessor::has_fragment_cache(&sample_dir) {
                println!("    Loading from per-bin cache for {} (slower)...", cell_type);
                
                // Convert indices to u32 for old cache format
                let peak_bin_vec: Vec<u32> = peak_indices_vec.iter().map(|&x| x as u32).collect();
                let bg_bin_vec: Vec<u32> = bg_indices_vec.iter().map(|&x| x as u32).collect();
                
                // Load peak fragments
                if !peak_bin_vec.is_empty() {
                    let peak_fragments = BamProcessor::extract_fragments_from_cache(
                        &sample_dir,
                        &peak_bin_vec,
                        true
                    )?;
                    
                    for (bin_str, frags) in peak_fragments {
                        let bin_idx: usize = bin_str.parse().unwrap();
                        if bin_idx < peaks.len() {
                            let region = &peaks[bin_idx];
                            let region_key = format!("{}:{}-{}", region.chrom, region.start, region.end);
                            type_fragments.entry(region_key)
                                .or_insert_with(Vec::new)
                                .extend(frags);
                        }
                    }
                }
                
                // Load background fragments
                if !bg_bin_vec.is_empty() {
                    let bg_fragments = BamProcessor::extract_fragments_from_cache(
                        &sample_dir,
                        &bg_bin_vec,
                        false
                    )?;
                    
                    for (bin_str, frags) in bg_fragments {
                        let bin_idx: usize = bin_str.parse().unwrap();
                        if bin_idx < background.len() {
                            let region = &background[bin_idx];
                            let region_key = format!("{}:{}-{}", region.chrom, region.start, region.end);
                            type_fragments.entry(region_key)
                                .or_insert_with(Vec::new)
                                .extend(frags);
                        }
                    }
                }
            } else {
                // This shouldn't happen with staged data
                return Err(anyhow::anyhow!("No fragment cache found for {}. Please re-run staging.", cell_type));
            }
            
            let total_frags: usize = type_fragments.values().map(|v| v.len()).sum();
            println!("    ✓ Loaded {} fragments for {}", total_frags, cell_type);
            
            fragments_by_type.insert(cell_type.clone(), type_fragments);
        }
        
        println!("  ✓ Successfully loaded type-specific fragments (preserving cell type identity!)");
        
        // Calculate total fragments across all types
        let total_fragments: usize = fragments_by_type.values()
            .flat_map(|type_frags| type_frags.values())
            .map(|v| v.len())
            .sum();
        println!("  ✓ Total fragments loaded: {}", total_fragments);
        
        // Run uniform fragment sampling with pre-assigned cell types and per-type fragment pools
        let uniform_sampler = UniformSampler::new(seed);
        
        // Process cells using type-specific fragment pools
        let mut cells = uniform_sampler.sample_fragments_for_cells_with_type_pools(
            cell_regions,
            &fragments_by_type,
            &cell_type_assignments,
            &cells_by_type,
        )?;
        
        // Apply doublet simulation if requested
        let doublet_annotations = if doublet_rate > 0.0 {
            // Create a temporary ScatrsConfig for the doublet simulator (for compatibility)
            // Get unique cell types from assignments
            let cell_types: Vec<String> = stage_config.as_ref().unwrap().samples.clone();
            let temp_config = ScatrsConfig {
                cell_types: cell_types.clone(),
                cell_count: num_cells as u32,
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
            // Collect all fragments from all cells to create ambient pool
            let ambient_simulator = AmbientSimulator::new(seed);
            let mut all_fragments = Vec::new();
            for cell in &cells {
                all_fragments.extend(cell.fragments.clone());
            }
            
            if all_fragments.is_empty() {
                println!("Warning: No fragments available for ambient contamination pool");
            } else {
                println!("  Created ambient pool with {} fragments from all {} cells", 
                    all_fragments.len(), cells.len());
                
                // Get the target fragments per cell from the configuration
                let target_fragments = match &sim_config.fragment_distribution {
                    crate::models::FragmentDistribution::Uniform { fragments_per_cell } => *fragments_per_cell,
                    crate::models::FragmentDistribution::Gaussian { mean_fragments_per_cell, .. } => *mean_fragments_per_cell as u32,
                };
                
                ambient_simulator.add_ambient_contamination(&mut cells, ambient_contamination, &all_fragments, target_fragments)?;
            }
        }
        
        // Write output files
        println!("\n💾 Writing output files...");
        let output_start = std::time::Instant::now();
        
        if compress {
            println!("  Writing compressed fragments file...");
            let fragments_file = output_dir.join("fragments.tsv.gz");
            FragmentWriter::write_fragments_file_gz(&cells, &fragments_file)?;
            println!("    ✓ {}", fragments_file.display());
            
            println!("  Writing compressed BED file...");
            let bed_file = output_dir.join("fragments.bed.gz");
            BedWriter::write_bed_gz(&cells, &bed_file)?;
            println!("    ✓ {}", bed_file.display());
        } else {
            println!("  Writing fragments file...");
            let fragments_file = output_dir.join("fragments.tsv");
            FragmentWriter::write_fragments_file(&cells, &fragments_file)?;
            println!("    ✓ {}", fragments_file.display());
            
            println!("  Writing BED file...");
            let bed_file = output_dir.join("fragments.bed");
            BedWriter::write_bed(&cells, &bed_file)?;
            println!("    ✓ {}", bed_file.display());
        }
        
        // Write metadata and summary
        println!("  Writing metadata files...");
        let metadata_file = output_dir.join("cell_metadata.tsv");
        BedWriter::write_cell_metadata(&cells, &metadata_file)?;
        println!("    ✓ {}", metadata_file.display());
        
        let barcodes_file = output_dir.join("barcodes.tsv");
        FragmentWriter::write_barcodes(&cells, &barcodes_file)?;
        println!("    ✓ {}", barcodes_file.display());
        
        let summary_file = output_dir.join("summary.txt");
        FragmentWriter::write_summary_stats(&cells, &summary_file)?;
        println!("    ✓ {}", summary_file.display());
        
        // Write doublet annotations if any
        if !doublet_annotations.is_empty() {
            println!("  Writing doublet annotations...");
            let doublet_file = output_dir.join("doublet_annotations.tsv");
            println!("    Attempting to write to: {}", doublet_file.display());
            println!("    Directory exists: {}", doublet_file.parent().unwrap().exists());
            FragmentWriter::write_doublet_annotations(&doublet_annotations, &doublet_file)?;
            println!("    File exists after write: {}", doublet_file.exists());
            println!("    ✓ {}", doublet_file.display());
        }
        
        println!("  ✓ All output files written in {:.1}s", output_start.elapsed().as_secs_f64());
        
        println!("\n╔══════════════════════════════════════════════════════╗");
        println!("║          ✅ Simulation complete!                      ║");
        println!("╚══════════════════════════════════════════════════════╝");
        println!("Output directory: {:?}", output_dir);
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
    signal_to_noise: 0.8       # Ratio of peak vs background fragments (0-1)
    doublet_rate: 0.1          # Fraction of cells that are doublets (0-1)
    ambient_contamination: 0.02 # Fraction of ambient DNA contamination (0-1)
    fragment_distribution:
      type: gaussian           # Distribution type: gaussian or uniform
      mean_fragments_per_cell: 10000.0  # Average fragments per cell
      variance: 4000000.0
      min: 1000
      max: 50000
    seed: 42                   # Random seed for reproducibility (null for random)
    compress: false            # Whether to gzip output files

  - name: sim2                 # Second simulation with different parameters
    weight_column: weight2     # Can use different weight columns from CSV
    output_dir: results/sim2/
    cell_count: 500
    signal_to_noise: 0.9
    doublet_rate: 0.05
    ambient_contamination: 0.01
    fragment_distribution:
      type: uniform
      fragments_per_cell: 5000  # Exact fragments per cell for uniform distribution
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
        
        std::fs::write(&output_file, yaml_content)
            .map_err(|e| anyhow::anyhow!("Failed to write config file to {:?}: {}", output_file, e))?;
        
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

/// Assign cell types to cells based on sample weights
fn assign_cell_types_by_weight(
    samples_with_weights: &[(String, f64)],
    num_cells: usize,
    seed: Option<u64>,
) -> Vec<String> {
    use rand::prelude::*;
    
    // Calculate number of cells per type based on weights
    let total_weight: f64 = samples_with_weights.iter().map(|(_, w)| w).sum();
    let mut assignments = Vec::with_capacity(num_cells);
    
    for (sample_name, weight) in samples_with_weights {
        let cells_for_type = ((weight / total_weight) * num_cells as f64).round() as usize;
        for _ in 0..cells_for_type {
            assignments.push(sample_name.clone());
        }
    }
    
    // Adjust for rounding errors
    while assignments.len() < num_cells {
        // Add to the type with highest weight
        let max_sample = samples_with_weights.iter()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(name, _)| name.clone())
            .unwrap_or_else(|| samples_with_weights[0].0.clone());
        assignments.push(max_sample);
    }
    while assignments.len() > num_cells {
        assignments.pop();
    }
    
    // Shuffle for random positioning
    let mut rng = StdRng::seed_from_u64(seed.unwrap_or(42));
    assignments.shuffle(&mut rng);
    
    assignments
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
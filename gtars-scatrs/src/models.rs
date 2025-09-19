use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::path::{Path, PathBuf};
use std::collections::HashMap;

// ============================================================================
// Fragment Model
// ============================================================================

/// Represents a single ATAC-seq fragment
/// 
/// A fragment corresponds to the DNA segment between paired-end reads,
/// representing accessible chromatin captured by the transposase.
/// 
/// # Example
/// ```
/// use gtars::scatrs::ScatrsFragment;
/// 
/// let fragment = ScatrsFragment::new("chr1".to_string(), 1000, 1200)
///     .with_cell_id("cell_001".to_string())
///     .with_quality(30);
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScatrsFragment {
    /// Chromosome name (e.g., "chr1", "chrX")
    pub chrom: String,
    /// Start position (0-based)
    pub start: u64,
    /// End position (exclusive)
    pub end: u64,
    /// Assigned cell barcode
    pub cell_id: String,
    /// DNA strand ('+' or '-')
    pub strand: Option<char>,
    /// Mapping quality score
    pub quality: Option<u8>,
}

impl ScatrsFragment {
    pub fn new(chrom: String, start: u64, end: u64) -> Self {
        Self {
            chrom,
            start,
            end,
            cell_id: String::new(),
            strand: None,
            quality: None,
        }
    }

    pub fn with_cell_id(mut self, cell_id: String) -> Self {
        self.cell_id = cell_id;
        self
    }

    pub fn with_strand(mut self, strand: char) -> Self {
        self.strand = Some(strand);
        self
    }

    pub fn with_quality(mut self, quality: u8) -> Self {
        self.quality = Some(quality);
        self
    }

    pub fn length(&self) -> u64 {
        self.end - self.start
    }

    /// Validate that the fragment has reasonable properties for ATAC-seq data
    /// 
    /// Returns true if the fragment passes all validation checks:
    /// - Fragment length is within expected range (50-2000bp)
    /// - Start position is less than end position
    /// - Chromosome name is not empty
    pub fn validate(&self) -> bool {
        let length = self.length();
        
        // Check basic properties
        if self.chrom.is_empty() {
            return false;
        }
        
        if self.start >= self.end {
            return false;
        }
        
        // Check fragment size is reasonable for ATAC-seq
        // Typical range is 50-2000bp, but we allow a bit more flexibility
        if length < 50 || length > 2000 {
            return false;
        }
        
        true
    }
    
    /// Validate with custom size range
    pub fn validate_with_range(&self, min_length: u64, max_length: u64) -> bool {
        let length = self.length();
        
        if self.chrom.is_empty() || self.start >= self.end {
            return false;
        }
        
        length >= min_length && length <= max_length
    }
}

// ============================================================================
// Region Model
// ============================================================================

/// Genomic interval with overlap and merge operations
/// 
/// Used to represent peaks, background regions, or any genomic interval.
/// Supports efficient overlap detection and merging of adjacent regions.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScatrsRegion {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub name: Option<String>,
    pub score: Option<f64>,
}

impl ScatrsRegion {
    pub fn new(chrom: String, start: u64, end: u64) -> Self {
        Self {
            chrom,
            start,
            end,
            name: None,
            score: None,
        }
    }

    pub fn with_name(mut self, name: String) -> Self {
        self.name = Some(name);
        self
    }

    pub fn with_score(mut self, score: f64) -> Self {
        self.score = Some(score);
        self
    }

    pub fn length(&self) -> u64 {
        self.end - self.start
    }

    pub fn overlaps(&self, other: &ScatrsRegion) -> bool {
        self.chrom == other.chrom && !(self.end <= other.start || self.start >= other.end)
    }

    pub fn merge(&self, other: &ScatrsRegion) -> Option<ScatrsRegion> {
        if self.chrom != other.chrom {
            return None;
        }
        
        Some(ScatrsRegion {
            chrom: self.chrom.clone(),
            start: self.start.min(other.start),
            end: self.end.max(other.end),
            name: self.name.clone().or(other.name.clone()),
            score: match (self.score, other.score) {
                (Some(s1), Some(s2)) => Some(s1.max(s2)),
                (Some(s), None) | (None, Some(s)) => Some(s),
                _ => None,
            },
        })
    }

    pub fn extend(&self, bp: u64) -> ScatrsRegion {
        ScatrsRegion {
            chrom: self.chrom.clone(),
            start: self.start.saturating_sub(bp),
            end: self.end + bp,
            name: self.name.clone(),
            score: self.score,
        }
    }
}

impl PartialEq for ScatrsRegion {
    fn eq(&self, other: &Self) -> bool {
        self.chrom == other.chrom && self.start == other.start && self.end == other.end
    }
}

impl Eq for ScatrsRegion {}

impl PartialOrd for ScatrsRegion {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ScatrsRegion {
    fn cmp(&self, other: &Self) -> Ordering {
        self.chrom.cmp(&other.chrom)
            .then(self.start.cmp(&other.start))
            .then(self.end.cmp(&other.end))
    }
}

// ============================================================================
// Peak Model
// ============================================================================

#[derive(Debug, Clone)]
pub struct Peak {
    pub region: ScatrsRegion,
    pub cell_type: String,
    pub read_count: u32,
    pub signal_value: f64,
    pub p_value: f64,
    pub q_value: f64,
    pub peak_point: Option<u64>,
}

impl Peak {
    pub fn new(
        region: ScatrsRegion,
        cell_type: String,
        signal_value: f64,
        p_value: f64,
        q_value: f64,
    ) -> Self {
        Self {
            region,
            cell_type,
            read_count: 0,
            signal_value,
            p_value,
            q_value,
            peak_point: None,
        }
    }

    pub fn with_read_count(mut self, count: u32) -> Self {
        self.read_count = count;
        self
    }

    pub fn with_peak_point(mut self, point: u64) -> Self {
        self.peak_point = Some(point);
        self
    }
}

// ============================================================================
// Cell Model
// ============================================================================

/// A simulated single cell with fragments and metadata
/// 
/// Contains all fragments assigned to this cell, along with quality metrics
/// and artifact information (doublet status, ambient contamination).
#[derive(Debug)]
pub struct SimulatedCell {
    pub cell_id: String,
    pub cell_type: String,
    pub fragments: Vec<ScatrsFragment>,
    pub metadata: CellMetadata,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CellMetadata {
    pub fragment_count: u32,
    pub peak_fragments: u32,
    pub background_fragments: u32,
    pub tss_enrichment: Option<f64>,
    pub frip_score: Option<f64>,
    pub is_doublet: bool,
    pub doublet_types: Option<Vec<(String, f64)>>,  // Cell types and their proportions
    pub ambient_fragments: u32,
}

impl SimulatedCell {
    pub fn new(cell_id: String, cell_type: String) -> Self {
        Self {
            cell_id,
            cell_type,
            fragments: Vec::new(),
            metadata: CellMetadata::default(),
        }
    }

    pub fn add_fragment(&mut self, fragment: ScatrsFragment) {
        self.fragments.push(fragment);
        self.metadata.fragment_count += 1;
    }

    pub fn add_peak_fragment(&mut self, fragment: ScatrsFragment) {
        self.fragments.push(fragment);
        self.metadata.fragment_count += 1;
        self.metadata.peak_fragments += 1;
    }

    pub fn add_background_fragment(&mut self, fragment: ScatrsFragment) {
        self.fragments.push(fragment);
        self.metadata.fragment_count += 1;
        self.metadata.background_fragments += 1;
    }

    pub fn calculate_frip(&mut self) {
        if self.metadata.fragment_count > 0 {
            self.metadata.frip_score = Some(
                self.metadata.peak_fragments as f64 / self.metadata.fragment_count as f64
            );
        }
    }
}

impl Default for CellMetadata {
    fn default() -> Self {
        Self {
            fragment_count: 0,
            peak_fragments: 0,
            background_fragments: 0,
            tss_enrichment: None,
            frip_score: None,
            is_doublet: false,
            doublet_types: None,
            ambient_fragments: 0,
        }
    }
}

// ============================================================================
// Configuration Model
// ============================================================================

/// Configuration for single-cell simulation
/// 
/// Controls all aspects of the simulation including cell counts,
/// fragment distributions, and artifact rates.
/// 
/// # Example
/// ```toml
/// cell_count = 1000
/// doublet_rate = 0.1
/// ambient_contamination = 0.02
/// signal_to_noise = 0.8
/// ```
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ScatrsConfig {
    pub cell_types: Vec<String>,
    pub cell_count: u32,
    pub signal_to_noise: f64,
    pub extend_peaks: u32,
    pub bin_size: u32,
    pub fragment_distribution: FragmentDistribution,
    pub output_dir: PathBuf,
    pub peak_dir: Option<PathBuf>,
    pub bam_dir: Option<PathBuf>,
    pub blacklist_file: Option<PathBuf>,
    pub genome_file: Option<PathBuf>,
    #[serde(default)]
    pub doublet_rate: f64,  // Fraction of cells that are doublets (0.0-1.0)
    #[serde(default)]
    pub ambient_contamination: f64,  // Fraction of fragments from ambient pool (0.0-1.0)
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(tag = "type")]
#[serde(rename_all = "lowercase")] 
pub enum FragmentDistribution {
    Uniform { 
        fragments_per_cell: u32 
    },
    Gaussian { 
        mean_fragments_per_cell: f64, 
        variance: f64, 
        min: u32, 
        max: u32 
    },
}

impl Default for ScatrsConfig {
    fn default() -> Self {
        Self {
            cell_types: vec!["default".to_string()],
            cell_count: 100,
            signal_to_noise: 0.8,
            extend_peaks: 250,
            bin_size: 5000,
            fragment_distribution: FragmentDistribution::Gaussian {
                mean_fragments_per_cell: 10000.0,
                variance: 2000.0,
                min: 1000,
                max: 50000,
            },
            output_dir: PathBuf::from("output"),
            peak_dir: None,
            bam_dir: None,
            blacklist_file: None,
            genome_file: None,
            doublet_rate: 0.0,  // No doublets by default
            ambient_contamination: 0.0,  // No ambient contamination by default
        }
    }
}

impl ScatrsConfig {
    pub fn from_file(path: &PathBuf) -> anyhow::Result<Self> {
        let content = std::fs::read_to_string(path)?;
        let config: Self = toml::from_str(&content)?;
        Ok(config)
    }

    pub fn to_file(&self, path: &PathBuf) -> anyhow::Result<()> {
        let content = toml::to_string_pretty(self)?;
        std::fs::write(path, content)?;
        Ok(())
    }
}

// ============================================================================
// Manifest Model
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct SampleManifest {
    pub sample_name: String,
    pub bam_file: PathBuf,
    #[serde(default)]
    pub peaks_file: Option<PathBuf>,
    #[serde(default = "default_weight")]
    pub weight: f64,  // Relative weight for cell type proportions
}

fn default_weight() -> f64 {
    1.0
}

// ============================================================================
// Stage Configuration Model
// ============================================================================

/// Minimal configuration file output by the staging step
/// 
/// Contains only essential metadata about staging location and parameters.
/// Sample weights are read directly from the CSV manifest during simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StageConfig {
    /// Version of the staging format
    pub version: String,
    /// Timestamp when staging was performed
    pub timestamp: String,
    /// Sample names that were staged (for validation)
    pub samples: Vec<String>,
    /// Staging parameters used
    pub parameters: StageParameters,
    /// Summary statistics from staging
    pub stats: StageStats,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StageParameters {
    pub extend_peaks: u32,
    pub bin_size: u32,
    pub merge_distance: u32,
    pub blacklist_file: Option<PathBuf>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StageStats {
    pub total_fragments: u64,
    pub total_peaks: Option<u32>,
    pub total_background_regions: Option<u32>,
    pub samples_processed: u32,
}

impl StageConfig {
    /// Helper function to resolve a path relative to a base directory
    /// If the path is already absolute, returns it unchanged
    /// If the path is relative, resolves it relative to the base directory
    fn resolve_config_path(path: &PathBuf, base_dir: &Path) -> PathBuf {
        if path.is_absolute() {
            path.clone()
        } else if path.as_os_str().is_empty() {
            // Handle empty paths by returning them unchanged
            path.clone()
        } else {
            // Resolve relative to config directory
            base_dir.join(path)
        }
    }
    
    pub fn to_file(&self, path: &Path) -> anyhow::Result<()> {
        let content = serde_yaml::to_string(self)?;
        std::fs::write(path, content)
            .map_err(|e| anyhow::anyhow!("Failed to write stage config to {:?}: {}", path, e))?;
        Ok(())
    }
    
    pub fn from_file(path: &Path) -> anyhow::Result<Self> {
        let content = std::fs::read_to_string(path)
            .map_err(|e| anyhow::anyhow!("Failed to read stage config from {:?}: {}", path, e))?;
        let mut config: Self = serde_yaml::from_str(&content)?;
        
        // Get the parent directory of the stage config file
        let config_dir = path.parent()
            .ok_or_else(|| anyhow::anyhow!("Stage config file path has no parent directory: {:?}", path))?;
        
        // Resolve blacklist file path if present
        if let Some(ref blacklist) = config.parameters.blacklist_file {
            config.parameters.blacklist_file = Some(Self::resolve_config_path(blacklist, config_dir));
        }
        
        Ok(config)
    }
}

// ============================================================================
// Unified Configuration Model (NEW)
// ============================================================================

/// Unified configuration for both staging and simulation
/// 
/// A single YAML config that serves as input to both stage and simulate commands.
/// The stage command adds the staging_config path after completion.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UnifiedConfig {
    /// Path to CSV manifest with samples
    pub sample_table: PathBuf,
    
    /// Path to chromosome sizes file (used by both staging and simulation)
    pub chrom_sizes: PathBuf,
    
    /// Optional blacklist file (used by both staging and simulation)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub blacklist: Option<PathBuf>,
    
    /// Staging parameters
    pub staging: StagingParams,
    
    /// Path to staging metadata (added after staging completes)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub staging_config: Option<PathBuf>,
    
    /// Multiple simulation configurations
    pub simulations: Vec<SimulationConfig>,
    
    /// Additional custom fields (for extensibility)
    #[serde(flatten)]
    pub extra: HashMap<String, serde_yaml::Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StagingParams {
    pub output_dir: PathBuf,
    #[serde(default = "default_extend_peaks")]
    pub extend_peaks: u32,
    #[serde(default = "default_bin_size")]
    pub bin_size: u32,
    #[serde(default = "default_merge_distance")]
    pub merge_distance: u32,
    #[serde(default)]
    pub skip_if_exists: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationConfig {
    /// Name of this simulation configuration
    pub name: String,
    
    /// Which weight column to use from the CSV manifest (e.g., "weight", "weight2", etc.)
    /// If None or column doesn't exist, uses equal weights
    #[serde(skip_serializing_if = "Option::is_none")]
    pub weight_column: Option<String>,
    
    pub output_dir: PathBuf,
    pub cell_count: u32,
    #[serde(default = "default_signal_to_noise")]
    pub signal_to_noise: f64,
    #[serde(default)]
    pub doublet_rate: f64,
    #[serde(default)]
    pub ambient_contamination: f64,
    pub fragment_distribution: FragmentDistribution,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub seed: Option<u64>,
    #[serde(default)]
    pub compress: bool,
}

// Default value functions
fn default_extend_peaks() -> u32 { 250 }
fn default_bin_size() -> u32 { 5000 }
fn default_merge_distance() -> u32 { 20 }
fn default_signal_to_noise() -> f64 { 0.8 }

impl UnifiedConfig {
    /// Helper function to resolve a path relative to a base directory
    /// If the path is already absolute, returns it unchanged
    /// If the path is relative, resolves it relative to the base directory
    fn resolve_config_path(path: &PathBuf, base_dir: &Path) -> PathBuf {
        if path.is_absolute() {
            path.clone()
        } else if path.as_os_str().is_empty() {
            // Handle empty paths by returning them unchanged
            path.clone()
        } else {
            // Resolve relative to config directory
            base_dir.join(path)
        }
    }
    
    pub fn from_yaml(path: &Path) -> anyhow::Result<Self> {
        let content = std::fs::read_to_string(path)
            .map_err(|e| anyhow::anyhow!("Failed to read config file from {:?}: {}", path, e))?;
        let mut config: Self = serde_yaml::from_str(&content)
            .map_err(|e| anyhow::anyhow!("Failed to parse YAML config from {:?}: {}", path, e))?;
        
        // Get the parent directory of the config file
        let config_dir = path.parent()
            .ok_or_else(|| anyhow::anyhow!("Config file path has no parent directory"))?;
        
        // Resolve all paths relative to the config file's directory
        config.sample_table = Self::resolve_config_path(&config.sample_table, config_dir);
        config.chrom_sizes = Self::resolve_config_path(&config.chrom_sizes, config_dir);
        
        // Resolve optional blacklist path
        if let Some(ref blacklist) = config.blacklist {
            config.blacklist = Some(Self::resolve_config_path(blacklist, config_dir));
        }
        
        // Resolve staging output directory
        config.staging.output_dir = Self::resolve_config_path(&config.staging.output_dir, config_dir);
        
        // Resolve staging_config path if present (this would have been added by a previous stage run)
        if let Some(ref staging_config) = config.staging_config {
            config.staging_config = Some(Self::resolve_config_path(staging_config, config_dir));
        }
        
        // Resolve paths for each simulation configuration
        for simulation in &mut config.simulations {
            simulation.output_dir = Self::resolve_config_path(&simulation.output_dir, config_dir);
        }
        
        Ok(config)
    }
    
    pub fn to_yaml(&self, path: &Path) -> anyhow::Result<()> {
        let content = serde_yaml::to_string(self)?;
        std::fs::write(path, content)?;
        Ok(())
    }
    
    /// Update the config file with staging_config path after staging
    /// This method only updates the staging_config field without modifying other paths
    pub fn update_staging_config(&mut self, staging_config_path: PathBuf, config_path: &Path) -> anyhow::Result<()> {
        // Read the original YAML file as a mutable Value
        let content = std::fs::read_to_string(config_path)
            .map_err(|e| anyhow::anyhow!("Failed to read config file: {}", e))?;
        let mut yaml_value: serde_yaml::Value = serde_yaml::from_str(&content)
            .map_err(|e| anyhow::anyhow!("Failed to parse YAML: {}", e))?;
        
        // Update only the staging_config field
        if let serde_yaml::Value::Mapping(ref mut map) = yaml_value {
            map.insert(
                serde_yaml::Value::String("staging_config".to_string()),
                serde_yaml::Value::String(staging_config_path.to_string_lossy().to_string())
            );
        }
        
        // Write the updated YAML back to file
        let updated_content = serde_yaml::to_string(&yaml_value)?;
        std::fs::write(config_path, updated_content)?;
        
        // Update our in-memory config as well
        self.staging_config = Some(staging_config_path);
        
        Ok(())
    }
    
    /// Get a specific simulation configuration by name
    pub fn get_simulation(&self, name: &str) -> Option<&SimulationConfig> {
        self.simulations.iter().find(|s| s.name == name)
    }
}

// ============================================================================
// Staged Data Metadata
// ============================================================================

/// Metadata for a staged sample
/// 
/// Contains processing statistics and references to the original BAM file
/// for later fragment extraction during simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleMetadata {
    pub sample_name: String,
    pub bam_file: PathBuf,
    pub total_fragments: u64,
    pub fragments_in_peaks: Option<u64>,
    pub fragments_in_background: Option<u64>,
    pub processing_time_secs: f64,
    pub timestamp: String,
}
use crate::models::{Peak, ScatrsRegion, SimulatedCell, SampleManifest, SampleMetadata};
use std::path::{Path, PathBuf};
use std::fs::File;
use std::io::{Write, BufWriter, BufRead, BufReader};
use anyhow::{Result, Context};
use flate2::write::GzEncoder;
use flate2::Compression;
use csv::Reader;

// ============================================================================
// NarrowPeak Reader
// ============================================================================

pub struct NarrowPeakReader;

impl NarrowPeakReader {
    pub fn read_peaks(path: &Path, cell_type: &str) -> Result<Vec<Peak>> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open narrowPeak file: {:?}", path))?;
        let reader = BufReader::new(file);
        let mut peaks = Vec::new();
        
        for (line_num, line) in reader.lines().enumerate() {
            let line = line.with_context(|| format!("Failed to read line {} from {:?}", line_num + 1, path))?;
            
            // Skip comment lines
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 10 {
                eprintln!("Warning: Line {} has {} fields (expected 10), skipping", line_num + 1, fields.len());
                continue;
            }
            
            let region = ScatrsRegion::new(
                fields[0].to_string(),
                fields[1].parse().with_context(|| format!("Invalid start position on line {}", line_num + 1))?,
                fields[2].parse().with_context(|| format!("Invalid end position on line {}", line_num + 1))?,
            )
            .with_name(fields[3].to_string())
            .with_score(fields[4].parse().unwrap_or(0.0));
            
            let mut peak = Peak::new(
                region,
                cell_type.to_string(),
                fields[6].parse().unwrap_or(0.0),  // signal_value
                fields[7].parse().unwrap_or(1.0),  // p_value
                fields[8].parse().unwrap_or(1.0),  // q_value
            );
            
            // Add peak point if available
            if let Ok(peak_point) = fields[9].parse::<u64>() {
                peak = peak.with_peak_point(peak_point);
            }
            
            peaks.push(peak);
        }
        
        Ok(peaks)
    }
}

// ============================================================================
// Manifest Reader
// ============================================================================

pub struct ManifestReader;

/// Extended sample manifest that can store multiple weight columns
pub struct ExtendedSampleManifest {
    pub sample_name: String,
    pub bam_file: PathBuf,
    pub peaks_file: Option<PathBuf>,
    pub weights: std::collections::HashMap<String, f64>,  // column_name -> weight
}

impl ManifestReader {
    /// Helper function to resolve a path relative to the manifest directory
    fn resolve_manifest_path(path: &PathBuf, manifest_dir: &Path) -> PathBuf {
        if path.is_absolute() {
            path.clone()
        } else if path.as_os_str().is_empty() {
            path.clone()
        } else {
            manifest_dir.join(path)
        }
    }
    
    pub fn read_manifest(manifest_path: &Path) -> Result<Vec<SampleManifest>> {
        let file = File::open(manifest_path)
            .with_context(|| format!("Failed to open manifest file: {:?}", manifest_path))?;
        let reader = BufReader::new(file);
        let mut csv_reader = Reader::from_reader(reader);
        
        // Get the parent directory of the manifest file for resolving relative paths
        let manifest_dir = manifest_path.parent()
            .ok_or_else(|| anyhow::anyhow!("Manifest file path has no parent directory"))?;
        
        let mut samples = Vec::new();
        for result in csv_reader.deserialize() {
            let mut sample: SampleManifest = result
                .with_context(|| "Failed to parse manifest row")?;
            
            // Resolve paths relative to manifest directory
            sample.bam_file = Self::resolve_manifest_path(&sample.bam_file, manifest_dir);
            if let Some(ref peaks_file) = sample.peaks_file {
                sample.peaks_file = Some(Self::resolve_manifest_path(peaks_file, manifest_dir));
            }
            
            // Validate that BAM file exists
            if !sample.bam_file.exists() {
                eprintln!("Warning: BAM file not found for sample {}: {:?}", 
                    sample.sample_name, sample.bam_file);
            }
            
            // Validate peaks file if provided
            if let Some(ref peaks_path) = sample.peaks_file {
                if !peaks_path.exists() {
                    eprintln!("Warning: Peaks file not found for sample {}: {:?}", 
                        sample.sample_name, peaks_path);
                }
            }
            
            samples.push(sample);
        }
        
        if samples.is_empty() {
            anyhow::bail!("No samples found in manifest file");
        }
        
        println!("Loaded {} samples from manifest", samples.len());
        for sample in &samples {
            println!("  - {}: BAM={:?}, Peaks={:?}", 
                sample.sample_name, 
                sample.bam_file.file_name().unwrap_or_default(),
                sample.peaks_file.as_ref().and_then(|p| p.file_name()).unwrap_or_default()
            );
        }
        
        Ok(samples)
    }
    
    /// Read manifest with flexible weight columns
    pub fn read_manifest_with_weights(manifest_path: &Path, weight_column: Option<&str>) -> Result<Vec<(SampleManifest, f64)>> {
        let file = File::open(manifest_path)
            .with_context(|| format!("Failed to open manifest file: {:?}", manifest_path))?;
        let mut reader = csv::Reader::from_reader(BufReader::new(file));
        
        // Get the parent directory of the manifest file for resolving relative paths
        let manifest_dir = manifest_path.parent()
            .ok_or_else(|| anyhow::anyhow!("Manifest file path has no parent directory"))?;
        
        // Get headers to check for weight columns
        let headers = reader.headers()?.clone();
        
        let mut samples_with_weights = Vec::new();
        
        for result in reader.records() {
            let record = result?;
            
            // Parse the basic fields
            let sample_name = record.get(0)
                .ok_or_else(|| anyhow::anyhow!("Missing sample_name"))?
                .to_string();
            let mut bam_file = PathBuf::from(record.get(1)
                .ok_or_else(|| anyhow::anyhow!("Missing bam_file"))?);
            let mut peaks_file = record.get(2)
                .and_then(|s| if s.is_empty() { None } else { Some(PathBuf::from(s)) });
            
            // Resolve paths relative to manifest directory
            bam_file = Self::resolve_manifest_path(&bam_file, manifest_dir);
            if let Some(ref pf) = peaks_file {
                peaks_file = Some(Self::resolve_manifest_path(pf, manifest_dir));
            }
            
            // Get weight from specified column or use default
            let weight = if let Some(col_name) = weight_column {
                // Try to find the specified weight column
                headers.iter()
                    .position(|h| h == col_name)
                    .and_then(|idx| record.get(idx))
                    .and_then(|v| v.parse::<f64>().ok())
                    .unwrap_or(1.0)
            } else {
                // Try default "weight" column, otherwise 1.0
                headers.iter()
                    .position(|h| h == "weight")
                    .and_then(|idx| record.get(idx))
                    .and_then(|v| v.parse::<f64>().ok())
                    .unwrap_or(1.0)
            };
            
            let sample = SampleManifest {
                sample_name,
                bam_file,
                peaks_file,
                weight: weight,  // Store the selected weight
            };
            
            samples_with_weights.push((sample, weight));
        }
        
        if samples_with_weights.is_empty() {
            anyhow::bail!("No samples found in manifest file");
        }
        
        println!("Loaded {} samples with weights from column: {}", 
            samples_with_weights.len(), 
            weight_column.unwrap_or("weight"));
        
        Ok(samples_with_weights)
    }
    
    pub fn create_example_manifest(path: &Path) -> Result<()> {
        let mut writer = csv::Writer::from_path(path)?;
        
        // Write header with multiple weight columns
        writer.write_record(&["sample_name", "bam_file", "peaks_file", "weight", "weight2"])?;
        
        // Write example rows showing different weight distributions
        writer.write_record(&[
            "T_cells",
            "/path/to/tcells.bam",
            "/path/to/tcells_peaks.narrowPeak",
            "2.0",   // Higher weight for first simulation
            "1.0"    // Equal weight for second simulation
        ])?;
        writer.write_record(&[
            "B_cells",
            "/path/to/bcells.bam",
            "",      // Optional peaks file
            "1.0",   // Lower weight for first simulation
            "2.0"    // Higher weight for second simulation
        ])?;
        writer.write_record(&[
            "NK_cells",
            "/path/to/nkcells.bam",
            "/path/to/nkcells_peaks.narrowPeak",
            "0.5",   // Lowest weight for first simulation
            "1.5"    // Medium weight for second simulation
        ])?;
        
        writer.flush()?;
        println!("Created example manifest at: {:?}", path);
        println!("\nNotes:");
        println!("  - Weight columns control cell type proportions in simulations");
        println!("  - Each simulation can use a different weight column via 'weight_column' parameter");
        println!("  - If weight_column is null or missing, uses equal weights (1.0 for all)");
        println!("  - Weights are relative: 2:1:0.5 gives ~57%:29%:14% cell distribution");
        
        Ok(())
    }
}

// ============================================================================
// BED Writer
// ============================================================================

pub struct BedWriter;

impl BedWriter {
    pub fn write_bed(cells: &[SimulatedCell], path: &Path) -> Result<()> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create BED file: {:?}", path))?;
        let mut writer = BufWriter::new(file);
        
        for cell in cells {
            for fragment in &cell.fragments {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t1\t{}",
                    fragment.chrom,
                    fragment.start,
                    fragment.end,
                    cell.cell_id,
                    fragment.strand.unwrap_or('.')
                )?;
            }
        }
        
        writer.flush()?;
        println!("Wrote {} cells to BED file: {:?}", cells.len(), path);
        
        Ok(())
    }
    
    pub fn write_bed_gz(cells: &[SimulatedCell], path: &Path) -> Result<()> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create compressed BED file: {:?}", path))?;
        let encoder = GzEncoder::new(file, Compression::default());
        let mut writer = BufWriter::new(encoder);
        
        for cell in cells {
            for fragment in &cell.fragments {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t1\t{}",
                    fragment.chrom,
                    fragment.start,
                    fragment.end,
                    cell.cell_id,
                    fragment.strand.unwrap_or('.')
                )?;
            }
        }
        
        writer.flush()?;
        println!("Wrote {} cells to compressed BED file: {:?}", cells.len(), path);
        
        Ok(())
    }
    
    pub fn write_cell_metadata(cells: &[SimulatedCell], path: &Path) -> Result<()> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create metadata file: {:?}", path))?;
        let mut writer = BufWriter::new(file);
        
        // Write header
        writeln!(
            writer,
            "cell_id\tcell_type\ttotal_fragments\tpeak_fragments\tbackground_fragments\tfrip_score"
        )?;
        
        // Write cell metadata
        for cell in cells {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{:.4}",
                cell.cell_id,
                cell.cell_type,
                cell.metadata.fragment_count,
                cell.metadata.peak_fragments,
                cell.metadata.background_fragments,
                cell.metadata.frip_score.unwrap_or(0.0)
            )?;
        }
        
        writer.flush()?;
        println!("Wrote cell metadata to: {:?}", path);
        
        Ok(())
    }
}

// ============================================================================
// Fragment Writer
// ============================================================================

pub struct FragmentWriter;

impl FragmentWriter {
    pub fn write_fragments_file(cells: &[SimulatedCell], path: &Path) -> Result<()> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create fragments file: {:?}", path))?;
        let mut writer = BufWriter::new(file);
        
        // Sort fragments by chromosome, start, end, and cell barcode for efficiency
        let mut all_fragments = Vec::new();
        for cell in cells {
            for fragment in &cell.fragments {
                all_fragments.push((
                    &fragment.chrom,
                    fragment.start,
                    fragment.end,
                    &cell.cell_id,
                ));
            }
        }
        
        all_fragments.sort_unstable();
        
        // Write fragments in 10x Genomics format
        for (chrom, start, end, cell_id) in all_fragments {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t1",
                chrom,
                start,
                end,
                cell_id
            )?;
        }
        
        writer.flush()?;
        println!("Wrote fragments file: {:?}", path);
        
        Ok(())
    }
    
    pub fn write_fragments_file_gz(cells: &[SimulatedCell], path: &Path) -> Result<()> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create compressed fragments file: {:?}", path))?;
        let encoder = GzEncoder::new(file, Compression::default());
        let mut writer = BufWriter::new(encoder);
        
        // Sort fragments
        let mut all_fragments = Vec::new();
        for cell in cells {
            for fragment in &cell.fragments {
                all_fragments.push((
                    &fragment.chrom,
                    fragment.start,
                    fragment.end,
                    &cell.cell_id,
                ));
            }
        }
        
        all_fragments.sort_unstable();
        
        // Write sorted fragments
        for (chrom, start, end, cell_id) in all_fragments {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t1",
                chrom,
                start,
                end,
                cell_id
            )?;
        }
        
        writer.flush()?;
        println!("Wrote compressed fragments file: {:?}", path);
        
        Ok(())
    }
    
    pub fn write_barcodes(cells: &[SimulatedCell], path: &Path) -> Result<()> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create barcodes file: {:?}", path))?;
        let mut writer = BufWriter::new(file);
        
        // Write unique cell barcodes
        let mut barcodes: Vec<_> = cells.iter().map(|c| &c.cell_id).collect();
        barcodes.sort_unstable();
        barcodes.dedup();
        
        for barcode in barcodes {
            writeln!(writer, "{}", barcode)?;
        }
        
        writer.flush()?;
        println!("Wrote {} barcodes to: {:?}", cells.len(), path);
        
        Ok(())
    }
    
    pub fn write_doublet_annotations(
        annotations: &[(String, Vec<String>)],
        path: &Path,
    ) -> Result<()> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create doublet annotations file: {:?}", path))?;
        let mut writer = BufWriter::new(file);
        
        // Write header
        writeln!(writer, "cell_id\tcell_type_1\tcell_type_2")?;
        
        // Write annotations
        for (cell_id, types) in annotations {
            if types.len() >= 2 {
                writeln!(writer, "{}\t{}\t{}", cell_id, types[0], types[1])?;
            }
        }
        
        writer.flush()?;
        println!("Wrote {} doublet annotations to: {:?}", annotations.len(), path);
        
        Ok(())
    }
    
    pub fn write_summary_stats(cells: &[SimulatedCell], path: &Path) -> Result<()> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create summary file: {:?}", path))?;
        let mut writer = BufWriter::new(file);
        
        // Calculate summary statistics
        let total_cells = cells.len();
        let total_fragments: u32 = cells.iter().map(|c| c.metadata.fragment_count).sum();
        let avg_fragments = if total_cells > 0 {
            total_fragments as f64 / total_cells as f64
        } else {
            0.0
        };
        
        let total_peak_fragments: u32 = cells.iter().map(|c| c.metadata.peak_fragments).sum();
        let avg_frip: f64 = if total_cells > 0 {
            cells.iter()
                .filter_map(|c| c.metadata.frip_score)
                .sum::<f64>() / total_cells as f64
        } else {
            0.0
        };
        
        // Count cells by type
        let mut cell_type_counts = std::collections::HashMap::new();
        for cell in cells {
            *cell_type_counts.entry(&cell.cell_type).or_insert(0) += 1;
        }
        
        // Write summary
        writeln!(writer, "=== Simulation Summary ===")?;
        writeln!(writer, "Total cells: {}", total_cells)?;
        writeln!(writer, "Total fragments: {}", total_fragments)?;
        writeln!(writer, "Average fragments per cell: {:.2}", avg_fragments)?;
        writeln!(writer, "Total peak fragments: {}", total_peak_fragments)?;
        writeln!(writer, "Average FRiP score: {:.4}", avg_frip)?;
        writeln!(writer)?;
        writeln!(writer, "Cells by type:")?;
        for (cell_type, count) in cell_type_counts {
            writeln!(writer, "  {}: {}", cell_type, count)?;
        }
        
        writer.flush()?;
        println!("Wrote summary statistics to: {:?}", path);
        
        Ok(())
    }
}

// ============================================================================
// Staged Data I/O
// ============================================================================

/// Writer for staged data files
/// 
/// Handles binary serialization of fragment counts and regions for efficient
/// storage and fast loading during simulation.
pub struct StagedDataWriter;

impl StagedDataWriter {
    /// Save fragment counts to binary file
    pub fn save_counts(counts: &[u32], path: &Path) -> Result<()> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create counts file: {:?}", path))?;
        bincode::serialize_into(file, counts)
            .with_context(|| format!("Failed to serialize counts to: {:?}", path))?;
        Ok(())
    }
    
    /// Save regions to binary file
    pub fn save_regions(regions: &[ScatrsRegion], path: &Path) -> Result<()> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create regions file: {:?}", path))?;
        bincode::serialize_into(file, regions)
            .with_context(|| format!("Failed to serialize regions to: {:?}", path))?;
        Ok(())
    }
    
    /// Save sample metadata to YAML file
    pub fn save_sample_metadata(metadata: &SampleMetadata, path: &Path) -> Result<()> {
        let yaml = serde_yaml::to_string(metadata)
            .with_context(|| "Failed to serialize metadata to YAML")?;
        std::fs::write(path, yaml)
            .with_context(|| format!("Failed to write metadata to: {:?}", path))?;
        Ok(())
    }
}

/// Reader for staged data files
/// 
/// Loads binary serialized data created by StagedDataWriter.
pub struct StagedDataReader;

impl StagedDataReader {
    /// Load fragment counts from binary file
    pub fn load_counts(path: &Path) -> Result<Vec<u32>> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open counts file: {:?}", path))?;
        let counts = bincode::deserialize_from(file)
            .with_context(|| format!("Failed to deserialize counts from: {:?}", path))?;
        Ok(counts)
    }
    
    /// Load regions from binary file
    pub fn load_regions(path: &Path) -> Result<Vec<ScatrsRegion>> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open regions file: {:?}", path))?;
        let regions = bincode::deserialize_from(file)
            .with_context(|| format!("Failed to deserialize regions from: {:?}", path))?;
        Ok(regions)
    }
    
    /// Load sample metadata from YAML file
    pub fn load_sample_metadata(path: &Path) -> Result<SampleMetadata> {
        let yaml = std::fs::read_to_string(path)
            .with_context(|| format!("Failed to read metadata from: {:?}", path))?;
        let metadata = serde_yaml::from_str(&yaml)
            .with_context(|| format!("Failed to parse metadata YAML from: {:?}", path))?;
        Ok(metadata)
    }
}
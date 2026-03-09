use std::fs;
use std::path::Path;

use anyhow::{Context, Result};
use clap::ArgMatches;
use serde::Serialize;

use gtars_sc::io::{read_10x, write_10x};
use gtars_sc::rna::qc::compute_rna_qc;
use gtars_sc::rna::pipeline::{run_full_rna_pipeline, run_rna_preprocessing};
use gtars_sc::types::{DownstreamConfig, FeatureType, RnaPipelineConfig};

use super::cli;
use super::output::{write_json, write_table, write_tsv, Format};

pub fn run_sc(matches: &ArgMatches) -> Result<()> {
    match matches.subcommand() {
        Some(("rna", sub)) => run_rna(sub),
        Some(("downstream", sub)) => run_downstream(sub),
        Some(("io", sub)) => run_io(sub),
        _ => unreachable!("sc subcommand not found"),
    }
}

fn run_rna(matches: &ArgMatches) -> Result<()> {
    match matches.subcommand() {
        Some((cli::RNA_PREPROCESS_CMD, sub)) => run_rna_preprocess(sub),
        Some((cli::RNA_QC_CMD, sub)) => run_rna_qc(sub),
        Some((cli::RNA_CONFIG_CMD, sub)) => run_rna_config(sub),
        _ => unreachable!("rna subcommand not found"),
    }
}

fn run_downstream(matches: &ArgMatches) -> Result<()> {
    match matches.subcommand() {
        Some((cli::DOWNSTREAM_ANALYZE_CMD, sub)) => run_downstream_analyze(sub),
        Some((cli::DOWNSTREAM_CLUSTER_CMD, sub)) => run_downstream_cluster(sub),
        _ => unreachable!("downstream subcommand not found"),
    }
}

fn run_io(matches: &ArgMatches) -> Result<()> {
    match matches.subcommand() {
        Some((cli::IO_INSPECT_CMD, sub)) => run_io_inspect(sub),
        _ => unreachable!("io subcommand not found"),
    }
}

// =========================================================================
// RNA handlers
// =========================================================================

fn run_rna_preprocess(matches: &ArgMatches) -> Result<()> {
    let input = matches.get_one::<String>("input").expect("input required");
    let output_dir = matches.get_one::<String>("output").expect("output required");
    let format = Format::from_str(matches.get_one::<String>("format").unwrap());
    let quiet = matches.get_flag("quiet");

    let config = build_rna_config(matches)?;

    if !quiet {
        eprintln!("Reading 10X data from {}", input);
    }

    let fm = read_10x(Path::new(input))
        .with_context(|| format!("reading 10X data from {}", input))?;

    let input_n_cells = fm.n_cells();
    let input_n_features = fm.n_features();

    if !quiet {
        eprintln!(
            "Input: {} cells x {} features",
            input_n_cells, input_n_features
        );
    }

    let result = run_rna_preprocessing(fm, &config)
        .context("running RNA preprocessing pipeline")?;

    let output_path = Path::new(output_dir);

    // Write the filtered/normalized matrix
    write_10x(output_path, &result.matrix)
        .with_context(|| format!("writing matrix to {}", output_dir))?;

    // Write embedding as TSV
    write_embedding_tsv(output_path, &result.embedding.values, &result.matrix.cell_ids)?;

    // Write HVGs
    let hvg_path = output_path.join("hvgs.json");
    let hvg_json = serde_json::to_string_pretty(&result.variable_features)?;
    fs::write(&hvg_path, hvg_json)?;

    // Write metadata
    let metadata = PreprocessMetadata {
        command: "sc rna preprocess".to_string(),
        params: config.clone(),
        input: InputInfo {
            path: input.to_string(),
            n_cells: input_n_cells,
            n_features: input_n_features,
        },
        output: OutputInfo {
            n_cells: result.matrix.n_cells(),
            n_features: result.matrix.n_features(),
            n_hvgs: result.variable_features.len(),
            n_pcs: result.embedding.values.ncols(),
        },
        variance_explained: result.embedding.variance_explained.clone(),
    };
    let meta_json = serde_json::to_string_pretty(&metadata)?;
    fs::write(output_path.join("metadata.json"), meta_json)?;

    if !quiet {
        eprintln!(
            "Output: {} cells x {} features, {} HVGs, {} PCs",
            result.matrix.n_cells(),
            result.matrix.n_features(),
            result.variable_features.len(),
            result.embedding.values.ncols(),
        );
        eprintln!("Written to {}", output_dir);
    }

    // Also print summary to stdout in requested format
    match format {
        Format::Json => write_json(&metadata)?,
        Format::Table => {
            println!("Preprocessing complete:");
            println!(
                "  Input:  {} cells x {} features",
                input_n_cells, input_n_features
            );
            println!(
                "  Output: {} cells x {} features",
                result.matrix.n_cells(),
                result.matrix.n_features()
            );
            println!("  HVGs:   {}", result.variable_features.len());
            println!("  PCs:    {}", result.embedding.values.ncols());
        }
        Format::Tsv => {
            println!("metric\tvalue");
            println!("input_cells\t{}", input_n_cells);
            println!("input_features\t{}", input_n_features);
            println!("output_cells\t{}", result.matrix.n_cells());
            println!("output_features\t{}", result.matrix.n_features());
            println!("n_hvgs\t{}", result.variable_features.len());
            println!("n_pcs\t{}", result.embedding.values.ncols());
        }
    }

    Ok(())
}

fn run_rna_qc(matches: &ArgMatches) -> Result<()> {
    let input = matches.get_one::<String>("input").expect("input required");
    let format = Format::from_str(matches.get_one::<String>("format").unwrap());
    let quiet = matches.get_flag("quiet");

    if !quiet {
        eprintln!("Reading 10X data from {}", input);
    }

    let fm = read_10x(Path::new(input))
        .with_context(|| format!("reading 10X data from {}", input))?;

    let qc = compute_rna_qc(&fm);

    let n_cells = qc.len();
    let median_features = percentile_u32(&qc.iter().map(|m| m.n_features).collect::<Vec<_>>(), 50.0);
    let median_counts = percentile_f64(&qc.iter().map(|m| m.n_counts).collect::<Vec<_>>(), 50.0);
    let median_pct_mt = percentile_f64(&qc.iter().map(|m| m.pct_mt).collect::<Vec<_>>(), 50.0);

    match format {
        Format::Json => {
            let output = QcOutput {
                command: "sc rna qc".to_string(),
                input: input.to_string(),
                summary: QcSummary {
                    n_cells,
                    n_features: fm.n_features(),
                    median_features_per_cell: median_features,
                    median_counts_per_cell: median_counts,
                    median_pct_mt,
                },
                cells: qc
                    .iter()
                    .map(|m| QcCell {
                        cell_id: m.cell_id.clone(),
                        n_features: m.n_features,
                        n_counts: m.n_counts,
                        pct_mt: m.pct_mt,
                    })
                    .collect(),
            };
            write_json(&output)?;
        }
        Format::Tsv => {
            let headers = &["cell_id", "n_features", "n_counts", "pct_mt"];
            let rows: Vec<Vec<String>> = qc
                .iter()
                .map(|m| {
                    vec![
                        m.cell_id.clone(),
                        m.n_features.to_string(),
                        format!("{:.1}", m.n_counts),
                        format!("{:.2}", m.pct_mt),
                    ]
                })
                .collect();
            write_tsv(headers, &rows)?;
        }
        Format::Table => {
            println!(
                "QC Summary: {} cells x {} features",
                n_cells,
                fm.n_features()
            );
            println!("  Median features/cell: {}", median_features);
            println!("  Median counts/cell:   {:.1}", median_counts);
            println!("  Median % MT:          {:.2}%", median_pct_mt);
            println!();

            let headers = &["CELL_ID", "N_FEATURES", "N_COUNTS", "PCT_MT"];
            let rows: Vec<Vec<String>> = qc
                .iter()
                .take(20) // show first 20 for table view
                .map(|m| {
                    vec![
                        m.cell_id.clone(),
                        m.n_features.to_string(),
                        format!("{:.1}", m.n_counts),
                        format!("{:.2}", m.pct_mt),
                    ]
                })
                .collect();
            write_table(headers, &rows)?;
            if n_cells > 20 {
                println!("  ... ({} more cells)", n_cells - 20);
            }
        }
    }

    Ok(())
}

fn run_rna_config(matches: &ArgMatches) -> Result<()> {
    if matches.get_flag("defaults") {
        let preprocess = RnaPipelineConfig::default();
        let downstream = DownstreamConfig::default();
        let output = ConfigDefaults {
            preprocess,
            downstream,
        };
        let json = serde_json::to_string_pretty(&output)?;
        println!("{}", json);
    } else {
        println!("Use --defaults to print default configuration as JSON.");
        println!("Pipe to a file and modify for your dataset:");
        println!("  gtars sc rna config --defaults > pipeline.json");
    }
    Ok(())
}

// =========================================================================
// Downstream handlers
// =========================================================================

fn run_downstream_analyze(matches: &ArgMatches) -> Result<()> {
    let input = matches.get_one::<String>("input").expect("input required");
    let output_dir = matches.get_one::<String>("output").expect("output required");
    let format = Format::from_str(matches.get_one::<String>("format").unwrap());
    let quiet = matches.get_flag("quiet");

    let preprocess_config = RnaPipelineConfig::default(); // will be loaded from metadata
    let downstream_config = build_downstream_config(matches)?;

    if !quiet {
        eprintln!("Reading 10X data from {}", input);
    }

    let fm = read_10x(Path::new(input))
        .with_context(|| format!("reading 10X data from {}", input))?;

    // Check for embedding file — if present, this is a preprocessed directory
    let embedding_path = Path::new(input).join("embedding.tsv");
    if !embedding_path.exists() {
        // Run full pipeline (preprocess + downstream)
        if !quiet {
            eprintln!("No embedding found, running full pipeline");
        }
        let result = run_full_rna_pipeline(fm, &preprocess_config, &downstream_config)
            .context("running full RNA pipeline")?;

        let output_path = Path::new(output_dir);
        write_10x(output_path, &result.preprocessing.matrix)?;
        write_embedding_tsv(
            output_path,
            &result.preprocessing.embedding.values,
            &result.preprocessing.matrix.cell_ids,
        )?;

        // Write clusters
        let cluster_output = ClusterOutput {
            n_clusters: result.clusters.n_clusters,
            modularity: result.clusters.quality,
            assignments: result.clusters.assignments.clone(),
            silhouette_avg: result.silhouette_avg,
        };
        let cluster_json = serde_json::to_string_pretty(&cluster_output)?;
        fs::write(output_path.join("clusters.json"), cluster_json)?;

        // Write markers if present
        if let Some(ref markers) = result.markers {
            let markers_json = serde_json::to_string_pretty(markers)?;
            fs::write(output_path.join("markers.json"), markers_json)?;
        }

        if !quiet {
            eprintln!(
                "Found {} clusters (modularity: {:.4})",
                result.clusters.n_clusters,
                result.clusters.quality.unwrap_or(0.0)
            );
            eprintln!("Written to {}", output_dir);
        }

        match format {
            Format::Json => write_json(&cluster_output)?,
            Format::Table => {
                println!("Clustering complete:");
                println!("  Clusters: {}", result.clusters.n_clusters);
                if let Some(q) = result.clusters.quality {
                    println!("  Modularity: {:.4}", q);
                }
                if let Some(s) = result.silhouette_avg {
                    println!("  Silhouette: {:.4}", s);
                }
                // Print cluster sizes
                let mut sizes: Vec<(u32, usize)> = Vec::new();
                for c in 0..result.clusters.n_clusters {
                    let count = result
                        .clusters
                        .assignments
                        .iter()
                        .filter(|&&a| a == c)
                        .count();
                    sizes.push((c, count));
                }
                sizes.sort_by(|a, b| b.1.cmp(&a.1));
                println!("  Cluster sizes:");
                for (c, n) in &sizes {
                    println!("    {}: {}", c, n);
                }
            }
            Format::Tsv => {
                println!("cell_id\tcluster");
                for (i, &c) in result.clusters.assignments.iter().enumerate() {
                    println!(
                        "{}\t{}",
                        result.preprocessing.cell_metadata[i].cell_id, c
                    );
                }
            }
        }
    } else {
        eprintln!("Preprocessed directory detected. Use 'downstream cluster' for clustering only.");
        std::process::exit(2);
    }

    Ok(())
}

fn run_downstream_cluster(matches: &ArgMatches) -> Result<()> {
    let input = matches.get_one::<String>("input").expect("input required");
    let output_dir = matches.get_one::<String>("output").expect("output required");
    let format = Format::from_str(matches.get_one::<String>("format").unwrap());
    let quiet = matches.get_flag("quiet");

    // Read embedding from preprocessed directory
    let embedding_path = Path::new(input).join("embedding.tsv");
    if !embedding_path.exists() {
        eprintln!("Error: no embedding.tsv found in {}", input);
        eprintln!("Run 'gtars sc rna preprocess' first.");
        std::process::exit(1);
    }

    let (embedding, cell_ids) = read_embedding_tsv(&embedding_path)?;

    let k: usize = matches
        .get_one::<String>("k")
        .unwrap()
        .parse()
        .context("--k must be a positive integer")?;
    let resolution: f64 = matches
        .get_one::<String>("resolution")
        .unwrap()
        .parse()
        .context("--resolution must be a number")?;
    let prune_snn: f64 = matches
        .get_one::<String>("prune-snn")
        .unwrap()
        .parse()
        .context("--prune-snn must be a number")?;

    if !quiet {
        eprintln!(
            "Clustering {} cells (k={}, resolution={}, prune={})",
            cell_ids.len(),
            k,
            resolution,
            prune_snn,
        );
    }

    let knn = gtars_sc::reduce::neighbors::build_knn(&embedding, k)
        .context("building KNN graph")?;
    let snn = gtars_sc::reduce::neighbors::build_snn(&knn, prune_snn);
    let clusters = gtars_sc::cluster::leiden::leiden_clustering(&snn, resolution, 10)
        .context("running Leiden clustering")?;

    let output_path = Path::new(output_dir);
    fs::create_dir_all(output_path)?;

    let cluster_output = ClusterOutput {
        n_clusters: clusters.n_clusters,
        modularity: clusters.quality,
        assignments: clusters.assignments.clone(),
        silhouette_avg: None,
    };
    let cluster_json = serde_json::to_string_pretty(&cluster_output)?;
    fs::write(output_path.join("clusters.json"), &cluster_json)?;

    if !quiet {
        eprintln!(
            "Found {} clusters (modularity: {:.4})",
            clusters.n_clusters,
            clusters.quality.unwrap_or(0.0)
        );
    }

    match format {
        Format::Json => write_json(&cluster_output)?,
        Format::Table => {
            println!("Clusters: {}", clusters.n_clusters);
            if let Some(q) = clusters.quality {
                println!("Modularity: {:.4}", q);
            }
        }
        Format::Tsv => {
            println!("cell_id\tcluster");
            for (i, &c) in clusters.assignments.iter().enumerate() {
                println!("{}\t{}", cell_ids[i], c);
            }
        }
    }

    Ok(())
}

// =========================================================================
// IO handlers
// =========================================================================

fn run_io_inspect(matches: &ArgMatches) -> Result<()> {
    let input = matches.get_one::<String>("input").expect("input required");
    let format = Format::from_str(matches.get_one::<String>("format").unwrap());

    let fm = read_10x(Path::new(input))
        .with_context(|| format!("reading 10X data from {}", input))?;

    let feature_type_str = match &fm.feature_type {
        FeatureType::Gene => "Gene Expression",
        FeatureType::Peak => "Peaks",
        FeatureType::Custom(s) => s.as_str(),
    };

    let sparsity = 1.0 - (fm.matrix.nnz() as f64 / (fm.n_features() as f64 * fm.n_cells() as f64));

    match format {
        Format::Json => {
            let output = InspectOutput {
                input: input.to_string(),
                n_cells: fm.n_cells(),
                n_features: fm.n_features(),
                nnz: fm.matrix.nnz(),
                sparsity,
                feature_type: feature_type_str.to_string(),
            };
            write_json(&output)?;
        }
        Format::Tsv => {
            println!("metric\tvalue");
            println!("n_cells\t{}", fm.n_cells());
            println!("n_features\t{}", fm.n_features());
            println!("nnz\t{}", fm.matrix.nnz());
            println!("sparsity\t{:.6}", sparsity);
            println!("feature_type\t{}", feature_type_str);
        }
        Format::Table => {
            println!("{}", input);
            println!("  Cells:        {}", fm.n_cells());
            println!("  Features:     {}", fm.n_features());
            println!("  Non-zeros:    {}", fm.matrix.nnz());
            println!("  Sparsity:     {:.2}%", sparsity * 100.0);
            println!("  Feature type: {}", feature_type_str);
        }
    }

    Ok(())
}

// =========================================================================
// Helpers
// =========================================================================

fn build_rna_config(matches: &ArgMatches) -> Result<RnaPipelineConfig> {
    let min_features: u32 = matches
        .get_one::<String>("min-features")
        .unwrap()
        .parse()
        .context("--min-features must be a positive integer")?;
    let min_cells: u32 = matches
        .get_one::<String>("min-cells")
        .unwrap()
        .parse()
        .context("--min-cells must be a positive integer")?;
    let max_pct_mt: f64 = matches
        .get_one::<String>("max-pct-mt")
        .unwrap()
        .parse()
        .context("--max-pct-mt must be a number")?;
    let scale_factor: f64 = matches
        .get_one::<String>("scale-factor")
        .unwrap()
        .parse()
        .context("--scale-factor must be a number")?;
    let n_hvgs: usize = matches
        .get_one::<String>("n-hvgs")
        .unwrap()
        .parse()
        .context("--n-hvgs must be a positive integer")?;
    let n_pcs: usize = matches
        .get_one::<String>("n-pcs")
        .unwrap()
        .parse()
        .context("--n-pcs must be a positive integer")?;

    let clip_str = matches.get_one::<String>("clip").unwrap();
    let clip_value = if clip_str == "none" {
        None
    } else {
        Some(
            clip_str
                .parse::<f64>()
                .context("--clip must be a number or 'none'")?,
        )
    };

    Ok(RnaPipelineConfig {
        min_features,
        min_cells,
        max_pct_mt,
        scale_factor,
        n_variable_features: n_hvgs,
        n_pcs,
        clip_value,
    })
}

fn build_downstream_config(matches: &ArgMatches) -> Result<DownstreamConfig> {
    let k: usize = matches
        .get_one::<String>("k")
        .unwrap()
        .parse()
        .context("--k must be a positive integer")?;
    let resolution: f64 = matches
        .get_one::<String>("resolution")
        .unwrap()
        .parse()
        .context("--resolution must be a number")?;
    let prune_snn: f64 = matches
        .get_one::<String>("prune-snn")
        .unwrap()
        .parse()
        .context("--prune-snn must be a number")?;

    Ok(DownstreamConfig {
        k_neighbors: k,
        prune_snn,
        resolution,
        compute_markers: !matches.get_flag("no-markers"),
        compute_silhouette: !matches.get_flag("no-silhouette"),
        ..Default::default()
    })
}

fn write_embedding_tsv(
    dir: &Path,
    embedding: &ndarray::Array2<f64>,
    cell_ids: &[String],
) -> Result<()> {
    use std::io::Write;
    let path = dir.join("embedding.tsv");
    let file = fs::File::create(&path)
        .with_context(|| format!("creating {}", path.display()))?;
    let mut w = std::io::BufWriter::new(file);

    let n_pcs = embedding.ncols();
    // Header
    write!(w, "cell_id")?;
    for pc in 0..n_pcs {
        write!(w, "\tPC_{}", pc + 1)?;
    }
    writeln!(w)?;

    // Data
    for (i, cell_id) in cell_ids.iter().enumerate() {
        write!(w, "{}", cell_id)?;
        for pc in 0..n_pcs {
            write!(w, "\t{:.8}", embedding[[i, pc]])?;
        }
        writeln!(w)?;
    }

    Ok(())
}

fn read_embedding_tsv(path: &Path) -> Result<(ndarray::Array2<f64>, Vec<String>)> {
    use std::io::BufRead;
    let file = fs::File::open(path)
        .with_context(|| format!("opening {}", path.display()))?;
    let reader = std::io::BufReader::new(file);
    let mut lines = reader.lines();

    // Parse header to get number of PCs
    let header = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("empty embedding file"))?
        .context("reading header")?;
    let n_pcs = header.split('\t').count() - 1; // subtract cell_id column

    let mut cell_ids = Vec::new();
    let mut data = Vec::new();

    for line in lines {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < n_pcs + 1 {
            continue;
        }
        cell_ids.push(parts[0].to_string());
        for &val_str in &parts[1..=n_pcs] {
            data.push(val_str.parse::<f64>().context("parsing embedding value")?);
        }
    }

    let n_cells = cell_ids.len();
    let embedding = ndarray::Array2::from_shape_vec((n_cells, n_pcs), data)
        .context("reshaping embedding matrix")?;

    Ok((embedding, cell_ids))
}

fn percentile_u32(values: &[u32], pct: f64) -> u32 {
    if values.is_empty() {
        return 0;
    }
    let mut sorted = values.to_vec();
    sorted.sort();
    let idx = ((pct / 100.0) * (sorted.len() - 1) as f64).round() as usize;
    sorted[idx.min(sorted.len() - 1)]
}

fn percentile_f64(values: &[f64], pct: f64) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let idx = ((pct / 100.0) * (sorted.len() - 1) as f64).round() as usize;
    sorted[idx.min(sorted.len() - 1)]
}

// =========================================================================
// Output structs (serde)
// =========================================================================

#[derive(Serialize)]
struct InspectOutput {
    input: String,
    n_cells: usize,
    n_features: usize,
    nnz: usize,
    sparsity: f64,
    feature_type: String,
}

#[derive(Serialize)]
struct QcOutput {
    command: String,
    input: String,
    summary: QcSummary,
    cells: Vec<QcCell>,
}

#[derive(Serialize)]
struct QcSummary {
    n_cells: usize,
    n_features: usize,
    median_features_per_cell: u32,
    median_counts_per_cell: f64,
    median_pct_mt: f64,
}

#[derive(Serialize)]
struct QcCell {
    cell_id: String,
    n_features: u32,
    n_counts: f64,
    pct_mt: f64,
}

#[derive(Serialize)]
struct PreprocessMetadata {
    command: String,
    params: RnaPipelineConfig,
    input: InputInfo,
    output: OutputInfo,
    variance_explained: Vec<f64>,
}

#[derive(Serialize)]
struct InputInfo {
    path: String,
    n_cells: usize,
    n_features: usize,
}

#[derive(Serialize)]
struct OutputInfo {
    n_cells: usize,
    n_features: usize,
    n_hvgs: usize,
    n_pcs: usize,
}

#[derive(Serialize)]
struct ClusterOutput {
    n_clusters: u32,
    #[serde(skip_serializing_if = "Option::is_none")]
    modularity: Option<f64>,
    assignments: Vec<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    silhouette_avg: Option<f64>,
}

#[derive(Serialize)]
struct ConfigDefaults {
    preprocess: RnaPipelineConfig,
    downstream: DownstreamConfig,
}

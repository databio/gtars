use std::fs;
use std::path::Path;

use anyhow::{Context, Result};
use clap::ArgMatches;
use serde::{Deserialize, Serialize};

use gtars_sc::io::{read_10x, write_10x};
use gtars_sc::rna::qc::compute_rna_qc;
use gtars_sc::rna::pipeline::{run_full_rna_pipeline, run_rna_preprocessing};
use gtars_sc::types::{DownstreamConfig, FeatureType, RnaPipelineConfig};

use super::cli;
use super::error::{report_error, ScError};
use super::output::{write_json, write_json_compact, write_jsonl, write_table, write_tsv, Format};
#[cfg(feature = "sc-parquet")]
use super::output::{write_parquet, ParquetColumn};

pub fn run_sc(matches: &ArgMatches) -> Result<()> {
    // Determine if output is JSON-like (for structured error reporting)
    let json_mode = detect_json_mode(matches);

    let result = match matches.subcommand() {
        Some(("rna", sub)) => run_rna(sub),
        Some(("downstream", sub)) => run_downstream(sub),
        Some(("io", sub)) => run_io(sub),
        _ => unreachable!("sc subcommand not found"),
    };

    if let Err(err) = result {
        report_error(&err, json_mode);
    }

    Ok(())
}

/// Walk the subcommand tree to find a --format flag and check if it's JSON-like.
fn detect_json_mode(matches: &ArgMatches) -> bool {
    fn check_format(m: &ArgMatches) -> Option<bool> {
        m.try_get_one::<String>("format")
            .ok()
            .flatten()
            .map(|f| f == "json" || f == "jsonl")
    }
    if let Some(v) = check_format(matches) {
        return v;
    }
    if let Some((_, sub)) = matches.subcommand() {
        if let Some(v) = check_format(sub) {
            return v;
        }
        if let Some((_, sub2)) = sub.subcommand() {
            if let Some(v) = check_format(sub2) {
                return v;
            }
        }
    }
    false
}

fn run_rna(matches: &ArgMatches) -> Result<()> {
    match matches.subcommand() {
        Some((cli::RNA_PREPROCESS_CMD, sub)) => run_rna_preprocess(sub),
        Some((cli::RNA_QC_CMD, sub)) => run_rna_qc(sub),
        Some((cli::RNA_CONFIG_CMD, sub)) => run_rna_config(sub),
        Some((cli::RNA_FILTER_CMD, sub)) => run_rna_filter(sub),
        Some((cli::RNA_NORMALIZE_CMD, sub)) => run_rna_normalize(sub),
        Some((cli::RNA_HVG_CMD, sub)) => run_rna_hvg(sub),
        Some((cli::RNA_SCALE_CMD, sub)) => run_rna_scale(sub),
        Some((cli::RNA_PCA_CMD, sub)) => run_rna_pca(sub),
        _ => unreachable!("rna subcommand not found"),
    }
}

fn run_downstream(matches: &ArgMatches) -> Result<()> {
    match matches.subcommand() {
        Some((cli::DOWNSTREAM_ANALYZE_CMD, sub)) => run_downstream_analyze(sub),
        Some((cli::DOWNSTREAM_CLUSTER_CMD, sub)) => run_downstream_cluster(sub),
        Some((cli::DOWNSTREAM_MARKERS_CMD, sub)) => run_downstream_markers(sub),
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
    let format = get_format(matches);
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
        Format::Json | Format::JsonCompact | Format::Jsonl => write_json_auto(&format, &metadata)?,
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
        #[cfg(feature = "sc-parquet")]
        Format::Parquet => anyhow::bail!("Parquet output is not supported for preprocess summary; use json or tsv"),
    }

    Ok(())
}

fn run_rna_qc(matches: &ArgMatches) -> Result<()> {
    let input = matches.get_one::<String>("input").expect("input required");
    let format = get_format(matches);
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

    let qc_cells: Vec<QcCell> = qc
        .iter()
        .map(|m| QcCell {
            cell_id: m.cell_id.clone(),
            n_features: m.n_features,
            n_counts: m.n_counts,
            pct_mt: m.pct_mt,
        })
        .collect();

    match format {
        Format::Json | Format::JsonCompact => {
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
                cells: qc_cells,
            };
            write_json_auto(&format, &output)?;
        }
        Format::Jsonl => {
            write_jsonl(&qc_cells)?;
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
        #[cfg(feature = "sc-parquet")]
        Format::Parquet => {
            let path = std::path::Path::new("qc_metrics.parquet");
            let headers = &["cell_id", "n_features", "n_counts", "pct_mt"];
            let columns = vec![
                ParquetColumn::Str(qc.iter().map(|m| m.cell_id.clone()).collect()),
                ParquetColumn::U32(qc.iter().map(|m| m.n_features).collect()),
                ParquetColumn::F64(qc.iter().map(|m| m.n_counts).collect()),
                ParquetColumn::F64(qc.iter().map(|m| m.pct_mt).collect()),
            ];
            write_parquet(path, headers, &columns)?;
            eprintln!("Written to {}", path.display());
        }
    }

    Ok(())
}

fn run_rna_config(matches: &ArgMatches) -> Result<()> {
    if matches.get_flag("schema") {
        let schema = schemars::schema_for!(ConfigDefaults);
        let json = serde_json::to_string_pretty(&schema)?;
        println!("{}", json);
    } else if matches.get_flag("defaults") {
        let preprocess = RnaPipelineConfig::default();
        let downstream = DownstreamConfig::default();
        let output = ConfigDefaults {
            preprocess,
            downstream,
        };
        let yaml = serde_yaml::to_string(&output)?;
        println!("# gtars sc configuration");
        println!("# Generated by: gtars sc rna config --defaults");
        println!("# Save to a file and pass with: --config pipeline.yaml");
        println!();
        print!("{}", yaml);
    } else {
        println!("Use --defaults to print default configuration as YAML.");
        println!("Use --schema to print JSON Schema for config validation.");
        println!();
        println!("  gtars sc rna config --defaults > pipeline.yaml");
        println!("  gtars sc rna config --schema > config_schema.json");
    }
    Ok(())
}

// =========================================================================
// Individual step handlers (Phase 2)
// =========================================================================

fn run_rna_filter(matches: &ArgMatches) -> Result<()> {
    let input = matches.get_one::<String>("input").expect("input required");
    let output_dir = matches.get_one::<String>("output").expect("output required");
    let format = get_format(matches);
    let quiet = matches.get_flag("quiet");

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

    if !quiet {
        eprintln!("Reading 10X data from {}", input);
    }

    let fm = read_10x(Path::new(input))
        .with_context(|| format!("reading 10X data from {}", input))?;

    let input_cells = fm.n_cells();
    let input_features = fm.n_features();

    // Filter genes first, then cells (matching pipeline order)
    use gtars_sc::rna::qc::{filter_genes, filter_rna_cells};

    let fm = filter_genes(&fm, min_cells)?;
    let genes_after = fm.n_features();

    let fm = filter_rna_cells(&fm, min_features, max_pct_mt)?;
    let cells_after = fm.n_cells();

    let output_path = Path::new(output_dir);
    write_10x(output_path, &fm)
        .with_context(|| format!("writing filtered matrix to {}", output_dir))?;

    // Write metadata
    let metadata = FilterMetadata {
        command: "sc rna filter".to_string(),
        input: InputInfo {
            path: input.to_string(),
            n_cells: input_cells,
            n_features: input_features,
        },
        params: FilterParams {
            min_features,
            min_cells,
            max_pct_mt,
        },
        output: FilterOutput {
            n_cells: cells_after,
            n_features: genes_after,
            cells_removed: input_cells - cells_after,
            genes_removed: input_features - genes_after,
        },
    };
    let meta_json = serde_json::to_string_pretty(&metadata)?;
    fs::write(output_path.join("metadata.json"), meta_json)?;

    if !quiet {
        eprintln!(
            "Filtered: {} → {} cells, {} → {} genes",
            input_cells, cells_after, input_features, genes_after
        );
    }

    match format {
        Format::Json | Format::JsonCompact | Format::Jsonl => write_json_auto(&format, &metadata)?,
        Format::Table => {
            println!("Filter complete:");
            println!(
                "  Cells:    {} → {} ({} removed)",
                input_cells,
                cells_after,
                input_cells - cells_after
            );
            println!(
                "  Features: {} → {} ({} removed)",
                input_features,
                genes_after,
                input_features - genes_after
            );
        }
        Format::Tsv => {
            println!("metric\tvalue");
            println!("input_cells\t{}", input_cells);
            println!("output_cells\t{}", cells_after);
            println!("input_features\t{}", input_features);
            println!("output_features\t{}", genes_after);
        }
        #[cfg(feature = "sc-parquet")]
        Format::Parquet => anyhow::bail!("Parquet output is not supported for filter summary; use json or tsv"),
    }

    Ok(())
}

fn run_rna_normalize(matches: &ArgMatches) -> Result<()> {
    let input = matches.get_one::<String>("input").expect("input required");
    let output_dir = matches.get_one::<String>("output").expect("output required");
    let format = get_format(matches);
    let quiet = matches.get_flag("quiet");

    let scale_factor: f64 = matches
        .get_one::<String>("scale-factor")
        .unwrap()
        .parse()
        .context("--scale-factor must be a number")?;

    if !quiet {
        eprintln!("Reading 10X data from {}", input);
    }

    let mut fm = read_10x(Path::new(input))
        .with_context(|| format!("reading 10X data from {}", input))?;

    use gtars_sc::rna::normalize::log_normalize;
    log_normalize(&mut fm, scale_factor)?;

    let output_path = Path::new(output_dir);
    write_10x(output_path, &fm)
        .with_context(|| format!("writing normalized matrix to {}", output_dir))?;

    // Copy hvgs.json forward if present in input dir
    let hvgs_src = Path::new(input).join("hvgs.json");
    if hvgs_src.exists() {
        fs::copy(&hvgs_src, output_path.join("hvgs.json"))?;
    }

    let metadata = NormalizeMetadata {
        command: "sc rna normalize".to_string(),
        scale_factor,
        n_cells: fm.n_cells(),
        n_features: fm.n_features(),
    };
    let meta_json = serde_json::to_string_pretty(&metadata)?;
    fs::write(output_path.join("metadata.json"), meta_json)?;

    if !quiet {
        eprintln!(
            "Normalized {} cells x {} features (scale_factor={})",
            fm.n_cells(),
            fm.n_features(),
            scale_factor,
        );
    }

    match format {
        Format::Json | Format::JsonCompact | Format::Jsonl => write_json_auto(&format, &metadata)?,
        Format::Table => {
            println!("Normalization complete:");
            println!("  Cells:        {}", fm.n_cells());
            println!("  Features:     {}", fm.n_features());
            println!("  Scale factor: {}", scale_factor);
        }
        Format::Tsv => {
            println!("metric\tvalue");
            println!("n_cells\t{}", fm.n_cells());
            println!("n_features\t{}", fm.n_features());
            println!("scale_factor\t{}", scale_factor);
        }
        #[cfg(feature = "sc-parquet")]
        Format::Parquet => anyhow::bail!("Parquet output is not supported for normalize summary; use json or tsv"),
    }

    Ok(())
}

fn run_rna_hvg(matches: &ArgMatches) -> Result<()> {
    let input = matches.get_one::<String>("input").expect("input required");
    let output_dir = matches.get_one::<String>("output");
    let format = get_format(matches);
    let quiet = matches.get_flag("quiet");

    let n_features: usize = matches
        .get_one::<String>("n-features")
        .unwrap()
        .parse()
        .context("--n-features must be a positive integer")?;

    if !quiet {
        eprintln!("Reading 10X data from {}", input);
    }

    let fm = read_10x(Path::new(input))
        .with_context(|| format!("reading 10X data from {}", input))?;

    // Warn if data looks normalized (no integer counts)
    let sample_vals: Vec<f64> = fm
        .matrix
        .iter()
        .take(100)
        .map(|(&v, _)| v)
        .collect();
    let looks_normalized = sample_vals.iter().any(|v| *v > 0.0 && v.fract() != 0.0);
    if looks_normalized && !quiet {
        eprintln!("Warning: input appears to be normalized. HVG selection (VST) expects raw counts.");
        eprintln!("Run 'sc rna hvg' before 'sc rna normalize' for correct results.");
    }

    use gtars_sc::rna::hvg::find_variable_features;
    let hvgs = find_variable_features(&fm, n_features)?;

    if !quiet {
        eprintln!("Selected {} highly variable genes", hvgs.len());
    }

    // If output dir specified, write matrix + hvgs.json
    if let Some(out_dir) = output_dir {
        let output_path = Path::new(out_dir);
        write_10x(output_path, &fm)
            .with_context(|| format!("writing matrix to {}", out_dir))?;

        let hvg_json = serde_json::to_string_pretty(&hvgs)?;
        fs::write(output_path.join("hvgs.json"), &hvg_json)?;

        let metadata = HvgMetadata {
            command: "sc rna hvg".to_string(),
            n_features_input: fm.n_features(),
            n_hvgs: hvgs.len(),
            n_cells: fm.n_cells(),
        };
        let meta_json = serde_json::to_string_pretty(&metadata)?;
        fs::write(output_path.join("metadata.json"), meta_json)?;
    }

    // Output to stdout
    let hvg_output = HvgOutput {
        command: "sc rna hvg".to_string(),
        n_features_input: fm.n_features(),
        n_hvgs: hvgs.len(),
        n_cells: fm.n_cells(),
        hvgs: hvgs.clone(),
    };

    match format {
        Format::Json | Format::JsonCompact => write_json_auto(&format, &hvg_output)?,
        Format::Jsonl => {
            // One gene per line
            for gene in &hvgs {
                println!("{}", serde_json::to_string(gene)?);
            }
        }
        Format::Table => {
            println!(
                "HVG selection: {} → {} variable genes ({} cells)",
                fm.n_features(),
                hvgs.len(),
                fm.n_cells(),
            );
            println!();
            println!("Top 20 HVGs:");
            for (i, gene) in hvgs.iter().take(20).enumerate() {
                println!("  {:>3}. {}", i + 1, gene);
            }
            if hvgs.len() > 20 {
                println!("  ... ({} more)", hvgs.len() - 20);
            }
        }
        Format::Tsv => {
            println!("rank\tgene");
            for (i, gene) in hvgs.iter().enumerate() {
                println!("{}\t{}", i + 1, gene);
            }
        }
        #[cfg(feature = "sc-parquet")]
        Format::Parquet => anyhow::bail!("Parquet output is not supported for HVG list; use json or tsv"),
    }

    Ok(())
}

fn run_rna_scale(matches: &ArgMatches) -> Result<()> {
    let input = matches.get_one::<String>("input").expect("input required");
    let output_dir = matches.get_one::<String>("output").expect("output required");
    let format = get_format(matches);
    let quiet = matches.get_flag("quiet");

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

    // Determine hvgs.json path
    let hvgs_path = if let Some(hvg_path) = matches.get_one::<String>("hvgs") {
        std::path::PathBuf::from(hvg_path)
    } else {
        Path::new(input).join("hvgs.json")
    };

    if !hvgs_path.exists() {
        return Err(ScError::input(format!(
            "No hvgs.json found at {}. Run 'sc rna hvg' first, or pass --hvgs <path>.",
            hvgs_path.display()
        )).into());
    }

    let hvgs_content = fs::read_to_string(&hvgs_path)
        .with_context(|| format!("reading {}", hvgs_path.display()))?;
    let hvgs: Vec<String> =
        serde_json::from_str(&hvgs_content).context("parsing hvgs.json")?;

    if !quiet {
        eprintln!("Reading 10X data from {}", input);
    }

    let fm = read_10x(Path::new(input))
        .with_context(|| format!("reading 10X data from {}", input))?;

    if !quiet {
        eprintln!(
            "Scaling {} HVGs across {} cells (clip={:?})",
            hvgs.len(),
            fm.n_cells(),
            clip_value,
        );
    }

    use gtars_sc::rna::scale::scale_data;
    let scaled = scale_data(&fm, &hvgs, clip_value, None)?;

    // Write scaled dense matrix as gzipped TSV
    let output_path = Path::new(output_dir);
    fs::create_dir_all(output_path)?;

    write_scaled_tsv_gz(output_path, &scaled, &hvgs, &fm.cell_ids)?;

    // Copy hvgs.json forward
    let hvg_json = serde_json::to_string_pretty(&hvgs)?;
    fs::write(output_path.join("hvgs.json"), &hvg_json)?;

    let metadata = ScaleMetadata {
        command: "sc rna scale".to_string(),
        n_features: scaled.nrows(),
        n_cells: scaled.ncols(),
        clip_value,
    };
    let meta_json = serde_json::to_string_pretty(&metadata)?;
    fs::write(output_path.join("metadata.json"), meta_json)?;

    if !quiet {
        eprintln!(
            "Scaled matrix: {} features x {} cells",
            scaled.nrows(),
            scaled.ncols(),
        );
    }

    match format {
        Format::Json | Format::JsonCompact | Format::Jsonl => write_json_auto(&format, &metadata)?,
        Format::Table => {
            println!("Scale complete:");
            println!("  Features: {} (HVGs)", scaled.nrows());
            println!("  Cells:    {}", scaled.ncols());
            if let Some(clip) = clip_value {
                println!("  Clipped:  [-{}, {}]", clip, clip);
            }
        }
        Format::Tsv => {
            println!("metric\tvalue");
            println!("n_features\t{}", scaled.nrows());
            println!("n_cells\t{}", scaled.ncols());
            if let Some(clip) = clip_value {
                println!("clip_value\t{}", clip);
            }
        }
        #[cfg(feature = "sc-parquet")]
        Format::Parquet => anyhow::bail!("Parquet output is not supported for scale summary; use json or tsv"),
    }

    Ok(())
}

fn run_rna_pca(matches: &ArgMatches) -> Result<()> {
    let input = matches.get_one::<String>("input").expect("input required");
    let output_dir = matches.get_one::<String>("output").expect("output required");
    let format = get_format(matches);
    let quiet = matches.get_flag("quiet");

    let n_pcs: usize = matches
        .get_one::<String>("n-pcs")
        .unwrap()
        .parse()
        .context("--n-pcs must be a positive integer")?;

    // Read scaled matrix
    let scaled_path = Path::new(input).join("scaled.tsv.gz");
    if !scaled_path.exists() {
        return Err(ScError::input(format!(
            "No scaled.tsv.gz found in {}. Run 'sc rna scale' first.",
            input,
        )).into());
    }

    if !quiet {
        eprintln!("Reading scaled matrix from {}", scaled_path.display());
    }

    let (scaled, feature_names, cell_ids) = read_scaled_tsv_gz(&scaled_path)?;

    let n_pcs = n_pcs.min(scaled.nrows()).min(scaled.ncols());

    if !quiet {
        eprintln!(
            "Running PCA: {} features x {} cells → {} PCs",
            scaled.nrows(),
            scaled.ncols(),
            n_pcs,
        );
    }

    use gtars_sc::reduce::lanczos::lanczos_svd;
    let svd_result = lanczos_svd(&scaled, n_pcs)?;

    // Cell embeddings = Vt^T * diag(S) (Seurat convention: cells × PCs)
    let n_cells = scaled.ncols();
    let mut embeddings = ndarray::Array2::<f64>::zeros((n_cells, n_pcs));
    for cell in 0..n_cells {
        for pc in 0..n_pcs {
            embeddings[[cell, pc]] = svd_result.vt[[pc, cell]] * svd_result.s[pc];
        }
    }

    let output_path = Path::new(output_dir);
    fs::create_dir_all(output_path)?;

    write_embedding_tsv(output_path, &embeddings, &cell_ids)?;

    // Copy hvgs.json forward if present
    let hvgs_src = Path::new(input).join("hvgs.json");
    if hvgs_src.exists() {
        fs::copy(&hvgs_src, output_path.join("hvgs.json"))?;
    }

    let metadata = PcaMetadata {
        command: "sc rna pca".to_string(),
        n_features: feature_names.len(),
        n_cells,
        n_pcs,
        variance_explained: svd_result.variance_explained.clone(),
    };
    let meta_json = serde_json::to_string_pretty(&metadata)?;
    fs::write(output_path.join("metadata.json"), meta_json)?;

    if !quiet {
        let top_var: f64 = svd_result.variance_explained.iter().take(5).sum();
        eprintln!(
            "PCA complete: {} PCs (top 5 explain {:.1}% variance)",
            n_pcs,
            top_var * 100.0,
        );
    }

    match format {
        Format::Json | Format::JsonCompact | Format::Jsonl => write_json_auto(&format, &metadata)?,
        Format::Table => {
            println!("PCA complete:");
            println!("  Cells:      {}", n_cells);
            println!("  Components: {}", n_pcs);
            println!("  Variance explained (top 10):");
            for (i, ve) in svd_result.variance_explained.iter().take(10).enumerate() {
                println!("    PC_{}: {:.4}%", i + 1, ve * 100.0);
            }
        }
        Format::Tsv => {
            println!("pc\tvariance_explained");
            for (i, ve) in svd_result.variance_explained.iter().enumerate() {
                println!("PC_{}\t{:.6}", i + 1, ve);
            }
        }
        #[cfg(feature = "sc-parquet")]
        Format::Parquet => anyhow::bail!("Parquet output is not supported for PCA summary; use json or tsv"),
    }

    Ok(())
}

// =========================================================================
// Downstream handlers
// =========================================================================

fn run_downstream_analyze(matches: &ArgMatches) -> Result<()> {
    let input = matches.get_one::<String>("input").expect("input required");
    let output_dir = matches.get_one::<String>("output").expect("output required");
    let format = get_format(matches);
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
            Format::Json | Format::JsonCompact | Format::Jsonl => write_json_auto(&format, &cluster_output)?,
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
            #[cfg(feature = "sc-parquet")]
            Format::Parquet => anyhow::bail!("Parquet output is not supported for analyze summary; use json or tsv"),
        }
    } else {
        return Err(ScError::param(
            "Preprocessed directory detected (has embedding.tsv). Use 'downstream cluster' for clustering only."
        ).into());
    }

    Ok(())
}

fn run_downstream_cluster(matches: &ArgMatches) -> Result<()> {
    let input = matches.get_one::<String>("input").expect("input required");
    let output_dir = matches.get_one::<String>("output").expect("output required");
    let format = get_format(matches);
    let quiet = matches.get_flag("quiet");

    // Read embedding from preprocessed directory
    let embedding_path = Path::new(input).join("embedding.tsv");
    if !embedding_path.exists() {
        return Err(ScError::input(format!(
            "No embedding.tsv found in {}. Run 'gtars sc rna preprocess' first.",
            input,
        )).into());
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
        Format::Json | Format::JsonCompact | Format::Jsonl => write_json_auto(&format, &cluster_output)?,
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
        #[cfg(feature = "sc-parquet")]
        Format::Parquet => anyhow::bail!("Parquet output is not supported for cluster summary; use json or tsv"),
    }

    Ok(())
}

fn run_downstream_markers(matches: &ArgMatches) -> Result<()> {
    let input = matches.get_one::<String>("input").expect("input required");
    let format = get_format(matches);
    let quiet = matches.get_flag("quiet");

    let min_pct: f64 = matches
        .get_one::<String>("min-pct")
        .unwrap()
        .parse()
        .context("--min-pct must be a number")?;
    let min_log2fc: f64 = matches
        .get_one::<String>("min-log2fc")
        .unwrap()
        .parse()
        .context("--min-log2fc must be a number")?;

    // Read cluster assignments
    let clusters_path = if let Some(p) = matches.get_one::<String>("clusters") {
        std::path::PathBuf::from(p)
    } else {
        Path::new(input).join("clusters.json")
    };

    if !clusters_path.exists() {
        return Err(ScError::input(format!(
            "No clusters.json found at {}. Run 'sc downstream cluster' first, or pass --clusters <path>.",
            clusters_path.display()
        )).into());
    }

    let clusters_content = fs::read_to_string(&clusters_path)
        .with_context(|| format!("reading {}", clusters_path.display()))?;
    let cluster_data: ClusterOutput =
        serde_json::from_str(&clusters_content).context("parsing clusters.json")?;

    if !quiet {
        eprintln!("Reading 10X data from {}", input);
    }

    let fm = read_10x(Path::new(input))
        .with_context(|| format!("reading 10X data from {}", input))?;

    if fm.n_cells() != cluster_data.assignments.len() {
        return Err(ScError::input(format!(
            "Matrix has {} cells but clusters.json has {} assignments",
            fm.n_cells(),
            cluster_data.assignments.len(),
        )).with_details(serde_json::json!({
            "matrix_cells": fm.n_cells(),
            "cluster_assignments": cluster_data.assignments.len(),
        })).into());
    }

    if !quiet {
        eprintln!(
            "Finding markers for {} clusters ({} cells x {} features)",
            cluster_data.n_clusters,
            fm.n_cells(),
            fm.n_features(),
        );
    }

    use gtars_sc::markers::wilcoxon::find_all_markers;
    use gtars_sc::types::MarkerConfig;

    let marker_config = MarkerConfig {
        min_pct,
        min_log2fc,
        only_positive: true,
    };

    let markers = find_all_markers(&fm, &cluster_data.assignments, &marker_config)?;

    if !quiet {
        eprintln!("Found {} significant markers", markers.len());
    }

    match format {
        Format::Json | Format::JsonCompact => {
            let output = MarkersOutput {
                command: "sc downstream markers".to_string(),
                n_clusters: cluster_data.n_clusters,
                n_markers: markers.len(),
                markers: markers.clone(),
            };
            write_json_auto(&format, &output)?;
        }
        Format::Jsonl => {
            write_jsonl(&markers)?;
        }
        Format::Table => {
            println!(
                "Markers: {} significant genes across {} clusters",
                markers.len(),
                cluster_data.n_clusters,
            );
            println!();
            let headers = &["GENE", "CLUSTER", "AVG_LOG2FC", "P_VAL_ADJ", "PCT_IN", "PCT_OUT"];
            let rows: Vec<Vec<String>> = markers
                .iter()
                .take(50)
                .map(|m| {
                    vec![
                        m.gene.clone(),
                        m.cluster.to_string(),
                        format!("{:.3}", m.avg_log2fc),
                        format!("{:.2e}", m.pval_adj),
                        format!("{:.3}", m.pct_in),
                        format!("{:.3}", m.pct_out),
                    ]
                })
                .collect();
            write_table(headers, &rows)?;
            if markers.len() > 50 {
                println!("  ... ({} more markers)", markers.len() - 50);
            }
        }
        Format::Tsv => {
            let headers = &["gene", "cluster", "avg_log2fc", "pval", "pval_adj", "pct_in", "pct_out"];
            let rows: Vec<Vec<String>> = markers
                .iter()
                .map(|m| {
                    vec![
                        m.gene.clone(),
                        m.cluster.to_string(),
                        format!("{:.6}", m.avg_log2fc),
                        format!("{:.2e}", m.pval),
                        format!("{:.2e}", m.pval_adj),
                        format!("{:.6}", m.pct_in),
                        format!("{:.6}", m.pct_out),
                    ]
                })
                .collect();
            write_tsv(headers, &rows)?;
        }
        #[cfg(feature = "sc-parquet")]
        Format::Parquet => {
            let path = std::path::Path::new("markers.parquet");
            let headers = &["gene", "cluster", "avg_log2fc", "pval", "pval_adj", "pct_in", "pct_out"];
            let columns = vec![
                ParquetColumn::Str(markers.iter().map(|m| m.gene.clone()).collect()),
                ParquetColumn::U32(markers.iter().map(|m| m.cluster).collect()),
                ParquetColumn::F64(markers.iter().map(|m| m.avg_log2fc).collect()),
                ParquetColumn::F64(markers.iter().map(|m| m.pval).collect()),
                ParquetColumn::F64(markers.iter().map(|m| m.pval_adj).collect()),
                ParquetColumn::F64(markers.iter().map(|m| m.pct_in).collect()),
                ParquetColumn::F64(markers.iter().map(|m| m.pct_out).collect()),
            ];
            write_parquet(path, headers, &columns)?;
            eprintln!("Written to {}", path.display());
        }
    }

    Ok(())
}

// =========================================================================
// IO handlers
// =========================================================================

fn run_io_inspect(matches: &ArgMatches) -> Result<()> {
    let input = matches.get_one::<String>("input").expect("input required");
    let format = get_format(matches);

    let fm = read_10x(Path::new(input))
        .with_context(|| format!("reading 10X data from {}", input))?;

    let feature_type_str = match &fm.feature_type {
        FeatureType::Gene => "Gene Expression",
        FeatureType::Peak => "Peaks",
        FeatureType::Custom(s) => s.as_str(),
    };

    let sparsity = 1.0 - (fm.matrix.nnz() as f64 / (fm.n_features() as f64 * fm.n_cells() as f64));

    match format {
        Format::Json | Format::JsonCompact | Format::Jsonl => {
            let output = InspectOutput {
                input: input.to_string(),
                n_cells: fm.n_cells(),
                n_features: fm.n_features(),
                nnz: fm.matrix.nnz(),
                sparsity,
                feature_type: feature_type_str.to_string(),
            };
            write_json_auto(&format, &output)?;
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
        #[cfg(feature = "sc-parquet")]
        Format::Parquet => anyhow::bail!("Parquet output is not supported for inspect; use json or tsv"),
    }

    Ok(())
}

// =========================================================================
// Helpers
// =========================================================================

/// Load config from YAML/JSON file, then override with CLI flags.
/// CLI flags override config values only when explicitly provided (not default).
fn load_config_file(matches: &ArgMatches) -> Result<Option<ConfigDefaults>> {
    let config_path = match matches.get_one::<String>("config") {
        Some(p) => p,
        None => return Ok(None),
    };

    let content = fs::read_to_string(config_path)
        .with_context(|| format!("reading config file {}", config_path))?;

    let config: ConfigDefaults = if config_path.ends_with(".json") {
        serde_json::from_str(&content)
            .with_context(|| format!("parsing JSON config {}", config_path))?
    } else if config_path.ends_with(".toml") {
        toml::from_str(&content)
            .with_context(|| format!("parsing TOML config {}", config_path))?
    } else {
        // YAML for .yaml, .yml, or anything else
        serde_yaml::from_str(&content)
            .with_context(|| format!("parsing YAML config {}", config_path))?
    };

    Ok(Some(config))
}

/// Extract format from CLI args, applying --compact flag.
fn get_format(matches: &ArgMatches) -> Format {
    let format = Format::from_str(matches.get_one::<String>("format").unwrap());
    let compact = matches.try_get_one::<bool>("compact").ok().flatten().copied().unwrap_or(false);
    format.with_compact(compact)
}

/// Write JSON respecting compact mode.
fn write_json_auto<T: Serialize>(format: &Format, value: &T) -> anyhow::Result<()> {
    match format {
        Format::JsonCompact => write_json_compact(value),
        _ => write_json(value),
    }
}

/// Check if a CLI arg was explicitly provided (not just the default value).
fn is_explicit(matches: &ArgMatches, id: &str) -> bool {
    matches.value_source(id) == Some(clap::parser::ValueSource::CommandLine)
}

fn build_rna_config(matches: &ArgMatches) -> Result<RnaPipelineConfig> {
    // Start with defaults, then overlay config file, then explicit CLI flags
    let mut config = if let Some(file_config) = load_config_file(matches)? {
        file_config.preprocess
    } else {
        RnaPipelineConfig::default()
    };

    // Only override from CLI if explicitly provided (not default)
    if is_explicit(matches, "min-features") {
        config.min_features = matches
            .get_one::<String>("min-features")
            .unwrap()
            .parse()
            .context("--min-features must be a positive integer")?;
    }
    if is_explicit(matches, "min-cells") {
        config.min_cells = matches
            .get_one::<String>("min-cells")
            .unwrap()
            .parse()
            .context("--min-cells must be a positive integer")?;
    }
    if is_explicit(matches, "max-pct-mt") {
        config.max_pct_mt = matches
            .get_one::<String>("max-pct-mt")
            .unwrap()
            .parse()
            .context("--max-pct-mt must be a number")?;
    }
    if is_explicit(matches, "scale-factor") {
        config.scale_factor = matches
            .get_one::<String>("scale-factor")
            .unwrap()
            .parse()
            .context("--scale-factor must be a number")?;
    }
    if is_explicit(matches, "n-hvgs") {
        config.n_variable_features = matches
            .get_one::<String>("n-hvgs")
            .unwrap()
            .parse()
            .context("--n-hvgs must be a positive integer")?;
    }
    if is_explicit(matches, "n-pcs") {
        config.n_pcs = matches
            .get_one::<String>("n-pcs")
            .unwrap()
            .parse()
            .context("--n-pcs must be a positive integer")?;
    }
    if is_explicit(matches, "clip") {
        let clip_str = matches.get_one::<String>("clip").unwrap();
        config.clip_value = if clip_str == "none" {
            None
        } else {
            Some(
                clip_str
                    .parse::<f64>()
                    .context("--clip must be a number or 'none'")?,
            )
        };
    }

    Ok(config)
}

fn build_downstream_config(matches: &ArgMatches) -> Result<DownstreamConfig> {
    let mut config = if let Some(file_config) = load_config_file(matches)? {
        file_config.downstream
    } else {
        DownstreamConfig::default()
    };

    if is_explicit(matches, "k") {
        config.k_neighbors = matches
            .get_one::<String>("k")
            .unwrap()
            .parse()
            .context("--k must be a positive integer")?;
    }
    if is_explicit(matches, "resolution") {
        config.resolution = matches
            .get_one::<String>("resolution")
            .unwrap()
            .parse()
            .context("--resolution must be a number")?;
    }
    if is_explicit(matches, "prune-snn") {
        config.prune_snn = matches
            .get_one::<String>("prune-snn")
            .unwrap()
            .parse()
            .context("--prune-snn must be a number")?;
    }
    if is_explicit(matches, "no-markers") {
        config.compute_markers = !matches.get_flag("no-markers");
    }
    if is_explicit(matches, "no-silhouette") {
        config.compute_silhouette = !matches.get_flag("no-silhouette");
    }

    Ok(config)
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

fn write_scaled_tsv_gz(
    dir: &Path,
    scaled: &ndarray::Array2<f64>,
    feature_names: &[String],
    cell_ids: &[String],
) -> Result<()> {
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;

    let path = dir.join("scaled.tsv.gz");
    let file = fs::File::create(&path)
        .with_context(|| format!("creating {}", path.display()))?;
    let gz = GzEncoder::new(file, Compression::default());
    let mut w = std::io::BufWriter::new(gz);

    // Header: gene_name\tcell_1\tcell_2\t...
    write!(w, "gene")?;
    for cell_id in cell_ids {
        write!(w, "\t{}", cell_id)?;
    }
    writeln!(w)?;

    // Each row is a feature (gene), each column is a cell
    // scaled is (n_features × n_cells)
    for (i, gene) in feature_names.iter().enumerate() {
        write!(w, "{}", gene)?;
        for j in 0..scaled.ncols() {
            write!(w, "\t{:.8}", scaled[[i, j]])?;
        }
        writeln!(w)?;
    }

    Ok(())
}

fn read_scaled_tsv_gz(path: &Path) -> Result<(ndarray::Array2<f64>, Vec<String>, Vec<String>)> {
    use flate2::read::GzDecoder;
    use std::io::BufRead;

    let file = fs::File::open(path)
        .with_context(|| format!("opening {}", path.display()))?;
    let gz = GzDecoder::new(file);
    let reader = std::io::BufReader::new(gz);
    let mut lines = reader.lines();

    // Parse header to get cell IDs
    let header = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("empty scaled matrix file"))?
        .context("reading header")?;
    let header_parts: Vec<&str> = header.split('\t').collect();
    let cell_ids: Vec<String> = header_parts[1..].iter().map(|s| s.to_string()).collect();
    let n_cells = cell_ids.len();

    let mut feature_names = Vec::new();
    let mut data = Vec::new();

    for line in lines {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < n_cells + 1 {
            continue;
        }
        feature_names.push(parts[0].to_string());
        for &val_str in &parts[1..=n_cells] {
            data.push(val_str.parse::<f64>().context("parsing scaled value")?);
        }
    }

    let n_features = feature_names.len();
    let matrix = ndarray::Array2::from_shape_vec((n_features, n_cells), data)
        .context("reshaping scaled matrix")?;

    Ok((matrix, feature_names, cell_ids))
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

#[derive(Serialize, Deserialize)]
struct ClusterOutput {
    n_clusters: u32,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    modularity: Option<f64>,
    assignments: Vec<u32>,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    silhouette_avg: Option<f64>,
}

#[derive(Serialize, Deserialize, schemars::JsonSchema)]
struct ConfigDefaults {
    preprocess: RnaPipelineConfig,
    downstream: DownstreamConfig,
}

// Phase 2 output structs

#[derive(Serialize)]
struct FilterMetadata {
    command: String,
    input: InputInfo,
    params: FilterParams,
    output: FilterOutput,
}

#[derive(Serialize)]
struct FilterParams {
    min_features: u32,
    min_cells: u32,
    max_pct_mt: f64,
}

#[derive(Serialize)]
struct FilterOutput {
    n_cells: usize,
    n_features: usize,
    cells_removed: usize,
    genes_removed: usize,
}

#[derive(Serialize)]
struct NormalizeMetadata {
    command: String,
    scale_factor: f64,
    n_cells: usize,
    n_features: usize,
}

#[derive(Serialize)]
struct HvgMetadata {
    command: String,
    n_features_input: usize,
    n_hvgs: usize,
    n_cells: usize,
}

#[derive(Serialize)]
struct HvgOutput {
    command: String,
    n_features_input: usize,
    n_hvgs: usize,
    n_cells: usize,
    hvgs: Vec<String>,
}

#[derive(Serialize)]
struct ScaleMetadata {
    command: String,
    n_features: usize,
    n_cells: usize,
    clip_value: Option<f64>,
}

#[derive(Serialize)]
struct PcaMetadata {
    command: String,
    n_features: usize,
    n_cells: usize,
    n_pcs: usize,
    variance_explained: Vec<f64>,
}

#[derive(Serialize)]
struct MarkersOutput {
    command: String,
    n_clusters: u32,
    n_markers: usize,
    markers: Vec<gtars_sc::types::MarkerResult>,
}

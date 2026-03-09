use clap::{Arg, ArgAction, Command};

pub const SC_CMD: &str = "sc";

// RNA subcommands
const RNA_CMD: &str = "rna";
pub const RNA_PREPROCESS_CMD: &str = "preprocess";
pub const RNA_QC_CMD: &str = "qc";
pub const RNA_CONFIG_CMD: &str = "config";
pub const RNA_FILTER_CMD: &str = "filter";
pub const RNA_NORMALIZE_CMD: &str = "normalize";
pub const RNA_HVG_CMD: &str = "hvg";
pub const RNA_SCALE_CMD: &str = "scale";
pub const RNA_PCA_CMD: &str = "pca";

// Downstream subcommands
const DOWNSTREAM_CMD: &str = "downstream";
pub const DOWNSTREAM_ANALYZE_CMD: &str = "analyze";
pub const DOWNSTREAM_CLUSTER_CMD: &str = "cluster";
pub const DOWNSTREAM_MARKERS_CMD: &str = "markers";

// IO subcommands
const IO_CMD: &str = "io";
pub const IO_INSPECT_CMD: &str = "inspect";

pub fn create_sc_cli() -> Command {
    Command::new(SC_CMD)
        .about("Single-cell analysis pipeline")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(create_rna_cli())
        .subcommand(create_downstream_cli())
        .subcommand(create_io_cli())
}

fn create_rna_cli() -> Command {
    Command::new(RNA_CMD)
        .about("RNA preprocessing commands")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(create_rna_preprocess_cli())
        .subcommand(create_rna_qc_cli())
        .subcommand(create_rna_config_cli())
        .subcommand(create_rna_filter_cli())
        .subcommand(create_rna_normalize_cli())
        .subcommand(create_rna_hvg_cli())
        .subcommand(create_rna_scale_cli())
        .subcommand(create_rna_pca_cli())
}

fn create_downstream_cli() -> Command {
    Command::new(DOWNSTREAM_CMD)
        .about("Post-preprocessing analysis (clustering, markers)")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(create_downstream_analyze_cli())
        .subcommand(create_downstream_cluster_cli())
        .subcommand(create_downstream_markers_cli())
}

fn create_io_cli() -> Command {
    Command::new(IO_CMD)
        .about("Data I/O utilities")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(create_io_inspect_cli())
}

// --- Shared args ---

fn arg_input() -> Arg {
    Arg::new("input")
        .required(true)
        .help("Path to 10X directory (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz)")
}

fn arg_output() -> Arg {
    Arg::new("output")
        .long("output")
        .short('o')
        .help("Output directory (default: stdout for summaries)")
}

fn arg_format() -> Arg {
    #[cfg(not(feature = "sc-parquet"))]
    let values = vec!["table", "json", "jsonl", "tsv"];
    #[cfg(feature = "sc-parquet")]
    let values = vec!["table", "json", "jsonl", "tsv", "parquet"];
    Arg::new("format")
        .long("format")
        .default_value("table")
        .value_parser(values)
        .help("Output format (jsonl: one JSON object per line for streaming)")
}

fn arg_compact() -> Arg {
    Arg::new("compact")
        .long("compact")
        .action(ArgAction::SetTrue)
        .help("Compact JSON output (single line, no pretty-printing)")
}

fn arg_seed() -> Arg {
    Arg::new("seed")
        .long("seed")
        .default_value("42")
        .help("Random seed for reproducibility")
}

fn arg_config() -> Arg {
    Arg::new("config")
        .long("config")
        .short('c')
        .help("Path to YAML/TOML config file (flags override config values)")
}

fn arg_quiet() -> Arg {
    Arg::new("quiet")
        .long("quiet")
        .short('q')
        .action(ArgAction::SetTrue)
        .help("Suppress progress messages on stderr")
}

// --- RNA subcommands ---

fn create_rna_preprocess_cli() -> Command {
    Command::new(RNA_PREPROCESS_CMD)
        .about("Run full RNA preprocessing pipeline (QC -> filter -> normalize -> HVG -> scale -> PCA)")
        .arg(arg_input())
        .arg(arg_output().required(true))
        .arg(arg_format())
        .arg(arg_compact())
        .arg(arg_config())
        .arg(arg_seed())
        .arg(arg_quiet())
        .arg(
            Arg::new("min-features")
                .long("min-features")
                .default_value("200")
                .help("Minimum genes per cell"),
        )
        .arg(
            Arg::new("min-cells")
                .long("min-cells")
                .default_value("3")
                .help("Minimum cells per gene"),
        )
        .arg(
            Arg::new("max-pct-mt")
                .long("max-pct-mt")
                .default_value("5.0")
                .help("Maximum mitochondrial percentage"),
        )
        .arg(
            Arg::new("scale-factor")
                .long("scale-factor")
                .default_value("10000")
                .help("Normalization scale factor"),
        )
        .arg(
            Arg::new("n-hvgs")
                .long("n-hvgs")
                .default_value("2000")
                .help("Number of highly variable genes"),
        )
        .arg(
            Arg::new("n-pcs")
                .long("n-pcs")
                .default_value("50")
                .help("Number of principal components"),
        )
        .arg(
            Arg::new("clip")
                .long("clip")
                .default_value("10.0")
                .help("Scale data clip value (use 'none' to skip clipping)"),
        )
}

fn create_rna_qc_cli() -> Command {
    Command::new(RNA_QC_CMD)
        .about("Compute per-cell QC metrics (n_features, n_counts, pct_mt)")
        .arg(arg_input())
        .arg(arg_format())
        .arg(arg_compact())
        .arg(arg_quiet())
}

fn create_rna_config_cli() -> Command {
    Command::new(RNA_CONFIG_CMD)
        .about("Dump default configuration or JSON Schema")
        .arg(
            Arg::new("defaults")
                .long("defaults")
                .action(ArgAction::SetTrue)
                .help("Print default config as YAML"),
        )
        .arg(
            Arg::new("schema")
                .long("schema")
                .action(ArgAction::SetTrue)
                .help("Print JSON Schema for config validation"),
        )
}

fn create_rna_filter_cli() -> Command {
    Command::new(RNA_FILTER_CMD)
        .about("Filter cells and genes (expects raw counts)")
        .arg(arg_input())
        .arg(arg_output().required(true))
        .arg(arg_format())
        .arg(arg_compact())
        .arg(arg_quiet())
        .arg(
            Arg::new("min-features")
                .long("min-features")
                .default_value("200")
                .help("Minimum genes per cell"),
        )
        .arg(
            Arg::new("min-cells")
                .long("min-cells")
                .default_value("3")
                .help("Minimum cells per gene"),
        )
        .arg(
            Arg::new("max-pct-mt")
                .long("max-pct-mt")
                .default_value("5.0")
                .help("Maximum mitochondrial percentage"),
        )
}

fn create_rna_normalize_cli() -> Command {
    Command::new(RNA_NORMALIZE_CMD)
        .about("Log-normalize count matrix")
        .arg(arg_input())
        .arg(arg_output().required(true))
        .arg(arg_format())
        .arg(arg_compact())
        .arg(arg_quiet())
        .arg(
            Arg::new("scale-factor")
                .long("scale-factor")
                .default_value("10000")
                .help("Normalization scale factor"),
        )
}

fn create_rna_hvg_cli() -> Command {
    Command::new(RNA_HVG_CMD)
        .about("Find highly variable genes (expects raw counts, before normalization)")
        .arg(arg_input())
        .arg(arg_output())
        .arg(arg_format())
        .arg(arg_compact())
        .arg(arg_quiet())
        .arg(
            Arg::new("n-features")
                .long("n-features")
                .default_value("2000")
                .help("Number of highly variable genes to select"),
        )
}

fn create_rna_scale_cli() -> Command {
    Command::new(RNA_SCALE_CMD)
        .about("Scale and center data on HVG subset (produces dense matrix)")
        .arg(arg_input())
        .arg(arg_output().required(true))
        .arg(arg_format())
        .arg(arg_compact())
        .arg(arg_quiet())
        .arg(
            Arg::new("clip")
                .long("clip")
                .default_value("10.0")
                .help("Clip scaled values to [-clip, clip] (use 'none' to skip)"),
        )
        .arg(
            Arg::new("hvgs")
                .long("hvgs")
                .help("Path to hvgs.json (default: <input>/hvgs.json)"),
        )
}

fn create_rna_pca_cli() -> Command {
    Command::new(RNA_PCA_CMD)
        .about("PCA dimensionality reduction on scaled data")
        .arg(arg_input())
        .arg(arg_output().required(true))
        .arg(arg_format())
        .arg(arg_compact())
        .arg(arg_quiet())
        .arg(
            Arg::new("n-pcs")
                .long("n-pcs")
                .default_value("50")
                .help("Number of principal components"),
        )
}

// --- Downstream subcommands ---

fn create_downstream_analyze_cli() -> Command {
    Command::new(DOWNSTREAM_ANALYZE_CMD)
        .about("Run full downstream analysis (KNN -> SNN -> cluster -> markers)")
        .arg(arg_input())
        .arg(arg_output().required(true))
        .arg(arg_format())
        .arg(arg_compact())
        .arg(arg_config())
        .arg(arg_seed())
        .arg(arg_quiet())
        .arg(
            Arg::new("k")
                .long("k")
                .default_value("20")
                .help("Number of nearest neighbors"),
        )
        .arg(
            Arg::new("resolution")
                .long("resolution")
                .default_value("0.8")
                .help("Leiden clustering resolution"),
        )
        .arg(
            Arg::new("prune-snn")
                .long("prune-snn")
                .default_value("0.0667")
                .help("SNN pruning threshold (Jaccard)"),
        )
        .arg(
            Arg::new("no-markers")
                .long("no-markers")
                .action(ArgAction::SetTrue)
                .help("Skip marker gene detection"),
        )
        .arg(
            Arg::new("no-silhouette")
                .long("no-silhouette")
                .action(ArgAction::SetTrue)
                .help("Skip silhouette score computation"),
        )
}

fn create_downstream_cluster_cli() -> Command {
    Command::new(DOWNSTREAM_CLUSTER_CMD)
        .about("Run KNN + SNN + Leiden clustering only (no markers)")
        .arg(arg_input())
        .arg(arg_output().required(true))
        .arg(arg_format())
        .arg(arg_compact())
        .arg(arg_seed())
        .arg(arg_quiet())
        .arg(
            Arg::new("k")
                .long("k")
                .default_value("20")
                .help("Number of nearest neighbors"),
        )
        .arg(
            Arg::new("resolution")
                .long("resolution")
                .default_value("0.8")
                .help("Leiden clustering resolution"),
        )
        .arg(
            Arg::new("prune-snn")
                .long("prune-snn")
                .default_value("0.0667")
                .help("SNN pruning threshold (Jaccard)"),
        )
}

fn create_downstream_markers_cli() -> Command {
    Command::new(DOWNSTREAM_MARKERS_CMD)
        .about("Find marker genes for each cluster (Wilcoxon rank-sum test)")
        .arg(arg_input())
        .arg(arg_format())
        .arg(arg_compact())
        .arg(arg_quiet())
        .arg(
            Arg::new("clusters")
                .long("clusters")
                .help("Path to clusters.json (default: <input>/clusters.json)"),
        )
        .arg(
            Arg::new("min-pct")
                .long("min-pct")
                .default_value("0.1")
                .help("Minimum fraction of cells expressing a gene"),
        )
        .arg(
            Arg::new("min-log2fc")
                .long("min-log2fc")
                .default_value("0.25")
                .help("Minimum log2 fold change"),
        )
        .arg(
            Arg::new("only-positive")
                .long("only-positive")
                .action(ArgAction::SetTrue)
                .default_value("true")
                .help("Only return upregulated markers"),
        )
}

// --- IO subcommands ---

fn create_io_inspect_cli() -> Command {
    Command::new(IO_INSPECT_CMD)
        .about("Show matrix dimensions, sparsity, and feature types")
        .arg(arg_input())
        .arg(arg_format())
        .arg(arg_compact())
}

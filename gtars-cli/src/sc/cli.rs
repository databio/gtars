use clap::{Arg, ArgAction, Command};

pub const SC_CMD: &str = "sc";

// RNA subcommands
const RNA_CMD: &str = "rna";
pub const RNA_PREPROCESS_CMD: &str = "preprocess";
pub const RNA_QC_CMD: &str = "qc";
pub const RNA_CONFIG_CMD: &str = "config";

// Downstream subcommands
const DOWNSTREAM_CMD: &str = "downstream";
pub const DOWNSTREAM_ANALYZE_CMD: &str = "analyze";
pub const DOWNSTREAM_CLUSTER_CMD: &str = "cluster";

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
}

fn create_downstream_cli() -> Command {
    Command::new(DOWNSTREAM_CMD)
        .about("Post-preprocessing analysis (clustering, markers)")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(create_downstream_analyze_cli())
        .subcommand(create_downstream_cluster_cli())
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
    Arg::new("format")
        .long("format")
        .default_value("table")
        .value_parser(["table", "json", "tsv"])
        .help("Output format")
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
        .arg(arg_quiet())
}

fn create_rna_config_cli() -> Command {
    Command::new(RNA_CONFIG_CMD)
        .about("Dump default configuration")
        .arg(
            Arg::new("defaults")
                .long("defaults")
                .action(ArgAction::SetTrue)
                .help("Print default config as YAML"),
        )
}

// --- Downstream subcommands ---

fn create_downstream_analyze_cli() -> Command {
    Command::new(DOWNSTREAM_ANALYZE_CMD)
        .about("Run full downstream analysis (KNN -> SNN -> cluster -> markers)")
        .arg(arg_input())
        .arg(arg_output().required(true))
        .arg(arg_format())
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

// --- IO subcommands ---

fn create_io_inspect_cli() -> Command {
    Command::new(IO_INSPECT_CMD)
        .about("Show matrix dimensions, sparsity, and feature types")
        .arg(arg_input())
        .arg(arg_format())
}

use std::path::Path;
use std::time::Instant;

use anyhow::Result;
use clap::ArgMatches;

use gtars_genomicdist::{GeneModel, SignalMatrix};

/// Derive the default output path: strip `.gz` then append `.bin`.
fn default_output_path(input: &str) -> String {
    let stripped = input.strip_suffix(".gz").unwrap_or(input);
    format!("{}.bin", stripped)
}

pub fn run_prep(matches: &ArgMatches) -> Result<()> {
    let gtf_path = matches.get_one::<String>("gtf");
    let signal_path = matches.get_one::<String>("signal-matrix");
    let output_path = matches.get_one::<String>("output");

    if gtf_path.is_none() && signal_path.is_none() {
        anyhow::bail!("Provide at least one of --gtf or --signal-matrix");
    }

    if let Some(gtf) = gtf_path {
        let out = output_path
            .cloned()
            .unwrap_or_else(|| default_output_path(gtf));

        eprintln!("Parsing GTF: {}", gtf);
        let start = Instant::now();
        let model = GeneModel::from_gtf(gtf.as_str(), true, true)
            .map_err(|e| anyhow::anyhow!("Failed to parse GTF: {}", e))?;
        eprintln!("  parsed in {:.1}s", start.elapsed().as_secs_f64());

        eprintln!("Saving bincode: {}", out);
        let start = Instant::now();
        model
            .save_bin(Path::new(&out))
            .map_err(|e| anyhow::anyhow!("Failed to save bincode: {}", e))?;

        let size = std::fs::metadata(&out)
            .map(|m| m.len())
            .unwrap_or(0);
        eprintln!(
            "  wrote {} ({:.1} MB) in {:.1}s",
            out,
            size as f64 / 1_048_576.0,
            start.elapsed().as_secs_f64()
        );
    }

    if let Some(sm_path) = signal_path {
        let out = output_path
            .cloned()
            .unwrap_or_else(|| default_output_path(sm_path));

        eprintln!("Parsing signal matrix: {}", sm_path);
        let start = Instant::now();
        let sm = SignalMatrix::from_tsv(sm_path)
            .map_err(|e| anyhow::anyhow!("Failed to parse signal matrix: {}", e))?;
        eprintln!("  parsed in {:.1}s", start.elapsed().as_secs_f64());

        eprintln!("Saving bincode: {}", out);
        let start = Instant::now();
        sm.save_bin(Path::new(&out))
            .map_err(|e| anyhow::anyhow!("Failed to save bincode: {}", e))?;

        let size = std::fs::metadata(&out)
            .map(|m| m.len())
            .unwrap_or(0);
        eprintln!(
            "  wrote {} ({:.1} MB) in {:.1}s",
            out,
            size as f64 / 1_048_576.0,
            start.elapsed().as_secs_f64()
        );
    }

    Ok(())
}

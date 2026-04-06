use std::path::Path;
use std::time::Instant;

use anyhow::Result;
use clap::ArgMatches;

use gtars_genomicdist::{GenomicDistAnnotation, SignalMatrix};
use gtars_genomicdist::models::BinaryGenomeAssembly;

/// Derive the default output path: strip `.gz` then append `.bin`.
fn default_output_path(input: &str) -> String {
    let stripped = input.strip_suffix(".gz").unwrap_or(input);
    format!("{}.bin", stripped)
}

pub fn run_prep(matches: &ArgMatches) -> Result<()> {
    let gtf_path = matches.get_one::<String>("gtf");
    let signal_path = matches.get_one::<String>("signal-matrix");
    let fasta_path = matches.get_one::<String>("fasta");
    let output_path = matches.get_one::<String>("output");

    if gtf_path.is_none() && signal_path.is_none() && fasta_path.is_none() {
        anyhow::bail!("Provide at least one of --gtf, --signal-matrix, or --fasta");
    }

    if let Some(gtf) = gtf_path {
        let out = output_path
            .cloned()
            .unwrap_or_else(|| default_output_path(gtf));

        eprintln!("Parsing GTF: {}", gtf);
        let start = Instant::now();
        let ann = GenomicDistAnnotation::from_gtf(gtf)
            .map_err(|e| anyhow::anyhow!("Failed to build GDA: {}", e))?;
        eprintln!(
            "  parsed in {:.1}s ({} genes)",
            start.elapsed().as_secs_f64(),
            ann.gene_model.genes.len(),
        );

        eprintln!("Saving GDA: {}", out);
        let start = Instant::now();
        ann.save_bin(Path::new(&out))
            .map_err(|e| anyhow::anyhow!("Failed to save GDA: {}", e))?;

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

    if let Some(fa) = fasta_path {
        let out = output_path
            .cloned()
            .unwrap_or_else(|| {
                let stripped = fa.strip_suffix(".gz").unwrap_or(fa);
                format!("{}.fab", stripped)
            });

        eprintln!("Converting FASTA to .fab: {}", fa);
        let start = Instant::now();
        BinaryGenomeAssembly::write_from_fasta(Path::new(fa), Path::new(&out))
            .map_err(|e| anyhow::anyhow!("Failed to create .fab: {}", e))?;

        let size = std::fs::metadata(&out)
            .map(|m| m.len())
            .unwrap_or(0);
        eprintln!(
            "  wrote {} ({:.1} GB) in {:.1}s",
            out,
            size as f64 / 1_073_741_824.0,
            start.elapsed().as_secs_f64()
        );
    }

    Ok(())
}

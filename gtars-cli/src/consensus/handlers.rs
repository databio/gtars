use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use anyhow::{Context, Result};
use clap::ArgMatches;

use gtars_core::models::RegionSet;
use gtars_genomicdist::consensus::consensus;

pub fn run_consensus(matches: &ArgMatches) -> Result<()> {
    let bed_paths: Vec<&String> = matches
        .get_many::<String>("beds")
        .expect("--beds is required")
        .collect();

    let min_count: u32 = matches
        .get_one::<String>("min-count")
        .unwrap()
        .parse()
        .context("--min-count must be a positive integer")?;

    let output_path = matches.get_one::<String>("output");

    // Load all BED files
    let sets: Vec<RegionSet> = bed_paths
        .iter()
        .map(|p| {
            RegionSet::try_from(p.as_str())
                .map_err(|e| anyhow::anyhow!("Failed to load BED file {}: {}", p, e))
        })
        .collect::<Result<Vec<_>>>()?;

    eprintln!("Computing consensus across {} BED files...", sets.len());

    let regions = consensus(&sets);

    // Filter by min-count
    let filtered: Vec<_> = regions
        .iter()
        .filter(|r| r.count >= min_count)
        .collect();

    eprintln!(
        "{} consensus regions ({} after --min-count {} filter)",
        regions.len(),
        filtered.len(),
        min_count,
    );

    match output_path {
        Some(p) => {
            let mut file = File::create(Path::new(p))
                .with_context(|| format!("Failed to create output file: {}", p))?;
            for r in &filtered {
                writeln!(file, "{}\t{}\t{}\t{}", r.chr, r.start, r.end, r.count)?;
            }
            eprintln!("Output written to {}", p);
        }
        None => {
            let stdout = io::stdout();
            let mut out = stdout.lock();
            for r in &filtered {
                writeln!(out, "{}\t{}\t{}\t{}", r.chr, r.start, r.end, r.count)?;
            }
        }
    }

    Ok(())
}

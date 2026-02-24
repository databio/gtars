use std::io::{self, Write};
use std::path::Path;

use anyhow::{Context, Result};
use clap::ArgMatches;

use gtars_core::models::RegionSet;
use gtars_core::utils::get_chrom_sizes;
use gtars_genomicdist::IntervalRanges;

pub fn run_ranges(matches: &ArgMatches) -> Result<()> {
    match matches.subcommand() {
        Some(("reduce", m)) => {
            let rs = load_input(m)?;
            let result = rs.reduce();
            write_output(&result, m.get_one::<String>("output"))
        }
        Some(("trim", m)) => {
            let rs = load_input(m)?;
            let cs_path = m
                .get_one::<String>("chrom-sizes")
                .expect("--chrom-sizes is required");
            let chrom_sizes = get_chrom_sizes(cs_path);
            let result = rs.trim(&chrom_sizes);
            write_output(&result, m.get_one::<String>("output"))
        }
        Some(("promoters", m)) => {
            let rs = load_input(m)?;
            let upstream: u32 = m
                .get_one::<String>("upstream")
                .unwrap()
                .parse()
                .context("--upstream must be a positive integer")?;
            let downstream: u32 = m
                .get_one::<String>("downstream")
                .unwrap()
                .parse()
                .context("--downstream must be a positive integer")?;
            let result = rs.promoters(upstream, downstream);
            write_output(&result, m.get_one::<String>("output"))
        }
        Some(("setdiff", m)) => {
            let (a, b) = load_pair(m)?;
            let result = a.setdiff(&b);
            write_output(&result, m.get_one::<String>("output"))
        }
        Some(("pintersect", m)) => {
            let (a, b) = load_pair(m)?;
            let result = a.pintersect(&b);
            write_output(&result, m.get_one::<String>("output"))
        }
        Some(("concat", m)) => {
            let (a, b) = load_pair(m)?;
            let result = a.concat(&b);
            write_output(&result, m.get_one::<String>("output"))
        }
        Some(("union", m)) => {
            let (a, b) = load_pair(m)?;
            let result = a.union(&b);
            write_output(&result, m.get_one::<String>("output"))
        }
        Some(("jaccard", m)) => {
            let (a, b) = load_pair(m)?;
            let j = a.jaccard(&b);
            println!("{}", j);
            Ok(())
        }
        _ => unreachable!("ranges subcommand not found"),
    }
}

fn load_input(matches: &ArgMatches) -> Result<RegionSet> {
    let path = matches
        .get_one::<String>("input")
        .expect("--input is required");
    RegionSet::try_from(path.as_str())
        .map_err(|e| anyhow::anyhow!("Failed to load BED file: {}", e))
}

fn load_pair(matches: &ArgMatches) -> Result<(RegionSet, RegionSet)> {
    let a_path = matches.get_one::<String>("BED_A").expect("-a is required");
    let b_path = matches.get_one::<String>("BED_B").expect("-b is required");
    let a = RegionSet::try_from(a_path.as_str())
        .map_err(|e| anyhow::anyhow!("Failed to load BED file A: {}", e))?;
    let b = RegionSet::try_from(b_path.as_str())
        .map_err(|e| anyhow::anyhow!("Failed to load BED file B: {}", e))?;
    Ok((a, b))
}

fn write_output(rs: &RegionSet, output: Option<&String>) -> Result<()> {
    match output {
        Some(p) => {
            rs.to_bed(Path::new(p))
                .with_context(|| format!("Failed to write output to {}", p))?;
            eprintln!("Output written to {}", p);
        }
        None => {
            let stdout = io::stdout();
            let mut out = stdout.lock();
            for region in &rs.regions {
                writeln!(out, "{}", region.as_string())?;
            }
        }
    }
    Ok(())
}

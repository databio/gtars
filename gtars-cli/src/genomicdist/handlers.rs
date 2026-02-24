use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use anyhow::{Context, Result};
use clap::ArgMatches;
use serde::Serialize;

use gtars_core::models::RegionSet;
use gtars_core::utils::get_chrom_sizes;
use gtars_genomicdist::models::{ChromosomeStatistics, RegionBin, TssIndex};
use gtars_genomicdist::statistics::GenomicIntervalSetStatistics;
use gtars_genomicdist::{
    GeneModel, ExpectedPartitionResult, PartitionResult,
    calc_expected_partitions, calc_partitions, genome_partition_list,
};

#[derive(Serialize)]
struct GenomicDistOutput {
    scalars: Scalars,
    #[serde(skip_serializing_if = "Option::is_none")]
    partitions: Option<PartitionResult>,
    distributions: Distributions,
    #[serde(skip_serializing_if = "Option::is_none")]
    expected_partitions: Option<ExpectedPartitionResult>,
}

#[derive(Serialize)]
struct Scalars {
    number_of_regions: u32,
    mean_region_width: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    median_tss_dist: Option<f64>,
}

#[derive(Serialize)]
struct Distributions {
    widths: Vec<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    tss_distances: Option<Vec<u32>>,
    neighbor_distances: Vec<i64>,
    nearest_neighbors: Vec<u32>,
    region_distribution: Vec<RegionBin>,
    chromosome_stats: HashMap<String, ChromosomeStatistics>,
}

pub fn run_genomicdist(matches: &ArgMatches) -> Result<()> {
    let bed_path = matches
        .get_one::<String>("bed")
        .expect("--bed is required");

    let gtf_path = matches.get_one::<String>("gtf");
    let tss_path = matches.get_one::<String>("tss");
    let chrom_sizes_path = matches.get_one::<String>("chrom-sizes");
    let output_path = matches.get_one::<String>("output");
    let n_bins: u32 = matches
        .get_one::<String>("bins")
        .unwrap()
        .parse()
        .context("--bins must be a positive integer")?;

    // Load BED file
    let rs = RegionSet::try_from(bed_path.as_str())
        .map_err(|e| anyhow::anyhow!("Failed to load BED file: {}", e))?;

    // Optionally load chrom sizes
    let chrom_sizes: Option<HashMap<String, u32>> = chrom_sizes_path.map(|p| get_chrom_sizes(p));

    // --- Unconditional computations ---
    let widths = rs.calc_widths();
    let chromosome_stats = rs.chromosome_statistics();
    let region_dist_map = rs.region_distribution_with_bins(n_bins);
    let neighbor_distances = rs
        .calc_neighbor_distances()
        .map_err(|e| anyhow::anyhow!("Failed to compute neighbor distances: {}", e))?;
    let nearest_neighbors = rs
        .calc_nearest_neighbors()
        .map_err(|e| anyhow::anyhow!("Failed to compute nearest neighbors: {}", e))?;

    // Convert region distribution HashMap to sorted Vec
    let mut region_distribution: Vec<RegionBin> = region_dist_map.into_values().collect();
    region_distribution.sort_by_key(|b| (b.rid, b.chr.clone(), b.start));

    // Scalars
    let number_of_regions = rs.regions.len() as u32;
    let mean_region_width = if widths.is_empty() {
        0.0
    } else {
        widths.iter().map(|&w| w as f64).sum::<f64>() / widths.len() as f64
    };

    // --- Optional: load gene model from GTF ---
    let gene_model: Option<GeneModel> = match gtf_path {
        Some(p) => {
            let model = GeneModel::from_gtf(p.as_str(), true, true)
                .map_err(|e| anyhow::anyhow!("Failed to load GTF: {}", e))?;
            Some(model)
        }
        None => {
            eprintln!("No --gtf provided, skipping partitions.");
            None
        }
    };

    // --- TSS distances ---
    let tss_distances: Option<Vec<u32>> = if let Some(tss_p) = tss_path {
        // Explicit TSS BED file takes priority
        let tss_index = TssIndex::try_from(tss_p.as_str())
            .map_err(|e| anyhow::anyhow!("Failed to load TSS BED file: {}", e))?;
        Some(
            tss_index
                .calc_tss_distances(&rs)
                .map_err(|e| anyhow::anyhow!("Failed to compute TSS distances: {}", e))?,
        )
    } else if let Some(ref model) = gene_model {
        // Use gene starts from GTF as TSS
        let tss_index = TssIndex::try_from(model.genes.inner.clone())
            .map_err(|e| anyhow::anyhow!("Failed to build TSS index from GTF genes: {}", e))?;
        Some(
            tss_index
                .calc_tss_distances(&rs)
                .map_err(|e| anyhow::anyhow!("Failed to compute TSS distances: {}", e))?,
        )
    } else {
        if tss_path.is_none() {
            eprintln!("No --tss or --gtf provided, skipping TSS distances.");
        }
        None
    };

    let median_tss_dist = tss_distances.as_ref().map(|dists| {
        let mut sorted: Vec<f64> = dists.iter().map(|&d| d as f64).collect();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let n = sorted.len();
        if n == 0 {
            0.0
        } else if n % 2 == 0 {
            (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
        } else {
            sorted[n / 2]
        }
    });

    // --- Partitions ---
    let partitions: Option<PartitionResult> = gene_model.as_ref().map(|model| {
        let partition_list = genome_partition_list(model, 100, 2000, chrom_sizes.as_ref());
        calc_partitions(&rs, &partition_list, false)
    });

    let expected_partitions: Option<ExpectedPartitionResult> = match (&gene_model, &chrom_sizes) {
        (Some(model), Some(cs)) => {
            let partition_list = genome_partition_list(model, 100, 2000, Some(cs));
            Some(calc_expected_partitions(&rs, &partition_list, cs, false))
        }
        (Some(_), None) => {
            eprintln!("No --chrom-sizes provided, skipping expected partitions.");
            None
        }
        _ => None,
    };

    // --- Build output ---
    let output = GenomicDistOutput {
        scalars: Scalars {
            number_of_regions,
            mean_region_width,
            median_tss_dist,
        },
        partitions,
        distributions: Distributions {
            widths,
            tss_distances,
            neighbor_distances,
            nearest_neighbors,
            region_distribution,
            chromosome_stats,
        },
        expected_partitions,
    };

    let json = serde_json::to_string_pretty(&output)
        .context("Failed to serialize output to JSON")?;

    match output_path {
        Some(p) => {
            let mut file = File::create(Path::new(p))
                .with_context(|| format!("Failed to create output file: {}", p))?;
            file.write_all(json.as_bytes())?;
            eprintln!("Output written to {}", p);
        }
        None => {
            io::stdout().write_all(json.as_bytes())?;
            println!(); // trailing newline
        }
    }

    Ok(())
}

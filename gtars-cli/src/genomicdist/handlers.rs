use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use anyhow::{Context, Result};
use clap::ArgMatches;
use serde::Serialize;

use gtars_core::models::{Region, RegionSet};
use gtars_core::utils::get_chrom_sizes;
use gtars_genomicdist::models::{ChromosomeStatistics, GenomeAssembly, RegionBin, Strand, TssIndex};
use gtars_genomicdist::statistics::GenomicIntervalSetStatistics;
use gtars_genomicdist::{
    GeneModel, GenomicDistAnnotation, ExpectedPartitionResult, PartitionResult,
    calc_expected_partitions, calc_partitions, genome_partition_list,
    SignalMatrix, calc_summary_signal, ConditionStats,
    calc_gc_content, calc_dinucl_freq, DINUCL_ORDER,
    median_abs_distance,
};

#[derive(Serialize)]
struct GenomicDistOutput {
    scalars: Scalars,
    #[serde(skip_serializing_if = "Option::is_none")]
    partitions: Option<PartitionResult>,
    distributions: Distributions,
    #[serde(skip_serializing_if = "Option::is_none")]
    expected_partitions: Option<ExpectedPartitionResult>,
    #[serde(skip_serializing_if = "Option::is_none")]
    open_signal: Option<OpenSignalOutput>,
    #[serde(skip_serializing_if = "Option::is_none")]
    gc_content: Option<GcContentOutput>,
    #[serde(skip_serializing_if = "Option::is_none")]
    dinucl_freq: Option<DinuclFreqOutput>,
}

#[derive(Serialize)]
struct GcContentOutput {
    /// Mean GC content across all regions (0–1)
    mean: f64,
    /// Per-region GC content values (0–1), one per region in input order
    per_region: Vec<f64>,
}

#[derive(Serialize)]
struct DinuclFreqOutput {
    /// Dinucleotide names in canonical order (matches DINUCL_ORDER)
    dinucleotides: Vec<String>,
    /// `chr_start_end` label per region
    region_labels: Vec<String>,
    /// Per-region matrix: outer is regions, inner is 16 values matching
    /// `dinucleotides` order. Percentages (0–100) by default, or raw counts
    /// if --dinucl-raw-counts flag was passed.
    frequencies: Vec<[f64; 16]>,
    /// Whether `frequencies` are raw counts (true) or percentages (false)
    raw_counts: bool,
}

#[derive(Serialize)]
struct OpenSignalOutput {
    condition_names: Vec<String>,
    matrix_stats: Vec<ConditionStats>,
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
    tss_distances: Option<Vec<i64>>,
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
    let signal_matrix_path = matches.get_one::<String>("signal-matrix");
    let fasta_path = matches.get_one::<String>("fasta");
    let ignore_unk_chroms = matches.get_flag("ignore-unk-chroms");
    let n_bins: u32 = matches
        .get_one::<String>("bins")
        .unwrap()
        .parse()
        .context("--bins must be a positive integer")?;
    let promoter_upstream: u32 = matches
        .get_one::<String>("promoter-upstream")
        .unwrap()
        .parse()
        .context("--promoter-upstream must be a positive integer")?;
    let promoter_downstream: u32 = matches
        .get_one::<String>("promoter-downstream")
        .unwrap()
        .parse()
        .context("--promoter-downstream must be a positive integer")?;

    // Load BED file
    let rs = RegionSet::try_from(bed_path.as_str())
        .map_err(|e| anyhow::anyhow!("Failed to load BED file: {}", e))?;

    // Optionally load chrom sizes from explicit flag
    let explicit_chrom_sizes: Option<HashMap<String, u32>> = chrom_sizes_path.map(|p| get_chrom_sizes(p));

    // --- Unconditional computations ---
    let widths = rs.calc_widths();
    let chromosome_stats = rs.chromosome_statistics();
    let region_dist_map = match explicit_chrom_sizes.as_ref() {
        Some(cs) => rs.region_distribution_with_chrom_sizes(n_bins, cs),
        None => {
            eprintln!(
                "warning: --chrom-sizes not provided; using BED-file-derived bin width."
            );
            eprintln!(
                "         Outputs will NOT be comparable across files or aligned with"
            );
            eprintln!(
                "         reference genome positions. Pass --chrom-sizes <file> for"
            );
            eprintln!("         reference-aligned bins.");
            rs.region_distribution_with_bins(n_bins)
        }
    };
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

    // --- Optional: load gene model from GTF or GDA binary ---
    let gene_model: Option<GeneModel> = match gtf_path {
        Some(p) => {
            if p.ends_with(".bin") {
                let ann = GenomicDistAnnotation::load_bin(p)
                    .map_err(|e| anyhow::anyhow!("Failed to load GDA binary: {}", e))?;
                Some(ann.gene_model)
            } else {
                let model = GeneModel::from_gtf(p.as_str(), true, true)
                    .map_err(|e| anyhow::anyhow!("Failed to load GTF: {}", e))?;
                Some(model)
            }
        }
        None => {
            eprintln!("No --gtf provided, skipping partitions.");
            None
        }
    };

    let chrom_sizes: Option<HashMap<String, u32>> = explicit_chrom_sizes;

    // --- TSS distances (signed: negative = upstream, positive = downstream) ---
    let tss_distances: Option<Vec<i64>> = if let Some(tss_p) = tss_path {
        // Explicit TSS BED file takes priority
        let tss_index = TssIndex::try_from(tss_p.as_str())
            .map_err(|e| anyhow::anyhow!("Failed to load TSS BED file: {}", e))?;
        Some(
            tss_index
                .calc_feature_distances(&rs, gtars_genomicdist::CoordinateMode::Bed)
                .map_err(|e| anyhow::anyhow!("Failed to compute TSS distances: {}", e))?,
        )
    } else if let Some(ref model) = gene_model {
        // Extract actual TSS positions from gene model using strand:
        // Plus/Unstranded → gene start, Minus → gene end
        let tss_regions: Vec<Region> = model.genes.inner.regions.iter()
            .zip(model.genes.strands.iter())
            .map(|(r, strand)| {
                let tss_pos = match strand {
                    Strand::Minus => r.end.saturating_sub(1),
                    _ => r.start,
                };
                Region { chr: r.chr.clone(), start: tss_pos, end: tss_pos + 1, rest: None }
            })
            .collect();
        let tss_rs = RegionSet { regions: tss_regions, header: None, path: None };
        let tss_index = TssIndex::try_from(tss_rs)
            .map_err(|e| anyhow::anyhow!("Failed to build TSS index from GTF genes: {}", e))?;
        Some(
            tss_index
                .calc_feature_distances(&rs, gtars_genomicdist::CoordinateMode::Bed)
                .map_err(|e| anyhow::anyhow!("Failed to compute TSS distances: {}", e))?,
        )
    } else {
        if tss_path.is_none() {
            eprintln!("No --tss or --gtf provided, skipping TSS distances.");
        }
        None
    };

    // Median of absolute distances (for the scalar summary).
    // Filters out i64::MAX sentinels (regions on chromosomes with no TSS features).
    let median_tss_dist = tss_distances.as_ref().and_then(|dists| {
        median_abs_distance(dists)
    });

    // --- Partitions ---
    let partitions: Option<PartitionResult> = gene_model.as_ref().map(|model| {
        let partition_list = genome_partition_list(model, promoter_upstream, promoter_downstream, chrom_sizes.as_ref());
        calc_partitions(&rs, &partition_list, false)
    });

    let expected_partitions: Option<ExpectedPartitionResult> = match (&gene_model, &chrom_sizes) {
        (Some(model), Some(cs)) => {
            let partition_list = genome_partition_list(model, promoter_upstream, promoter_downstream, Some(cs));
            Some(calc_expected_partitions(&rs, &partition_list, cs, false))
        }
        (Some(_), None) => {
            eprintln!("No --chrom-sizes provided, skipping expected partitions.");
            None
        }
        _ => None,
    };

    // --- Optional: open chromatin signal enrichment ---
    let open_signal: Option<OpenSignalOutput> = match signal_matrix_path {
        Some(p) => {
            let sm = if p.ends_with(".bin") {
                SignalMatrix::load_bin(p)
                    .map_err(|e| anyhow::anyhow!("Failed to load signal matrix bincode: {}", e))?
            } else {
                SignalMatrix::from_tsv(p)
                    .map_err(|e| anyhow::anyhow!("Failed to load signal matrix: {}", e))?
            };
            let result = calc_summary_signal(&rs, &sm, gtars_genomicdist::CoordinateMode::Bed)
                .map_err(|e| anyhow::anyhow!("Failed to compute signal summary: {}", e))?;
            Some(OpenSignalOutput {
                condition_names: result.condition_names,
                matrix_stats: result.matrix_stats,
            })
        }
        None => None,
    };

    // --- Optional: GC content + dinucleotide frequencies (require FASTA) ---
    let dinucl_raw_counts = matches.get_flag("dinucl-raw-counts");
    let (gc_content_out, dinucl_freq_out) = match fasta_path {
        Some(p) => {
            let assembly = GenomeAssembly::try_from(p.as_str())
                .map_err(|e| anyhow::anyhow!("Failed to load FASTA: {}", e))?;
            let gc_per_region = calc_gc_content(&rs, &assembly, ignore_unk_chroms)
                .map_err(|e| anyhow::anyhow!("Failed to compute GC content: {}", e))?;
            let gc_mean = if gc_per_region.is_empty() {
                0.0
            } else {
                gc_per_region.iter().sum::<f64>() / gc_per_region.len() as f64
            };
            let gc_out = GcContentOutput {
                mean: gc_mean,
                per_region: gc_per_region,
            };

            let (labels, matrix) = calc_dinucl_freq(&rs, &assembly, dinucl_raw_counts)
                .map_err(|e| anyhow::anyhow!("Failed to compute dinucl freq: {}", e))?;
            let dinucl_out = DinuclFreqOutput {
                dinucleotides: DINUCL_ORDER
                    .iter()
                    .map(|d| d.to_string().unwrap_or_default())
                    .collect(),
                region_labels: labels,
                frequencies: matrix,
                raw_counts: dinucl_raw_counts,
            };
            (Some(gc_out), Some(dinucl_out))
        }
        None => (None, None),
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
        open_signal,
        gc_content: gc_content_out,
        dinucl_freq: dinucl_freq_out,
    };

    let compact = matches.get_flag("compact");
    let json = if compact {
        serde_json::to_string(&output)
    } else {
        serde_json::to_string_pretty(&output)
    }
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

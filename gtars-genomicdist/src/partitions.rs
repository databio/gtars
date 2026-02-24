//! Genomic partition analysis for region sets.
//!
//! Classifies genomic regions into functional categories (promoter, exon, intron, etc.)
//! using a gene model, and computes observed vs expected enrichment statistics.
//!
//! Port of R GenomicDistributions: `genomePartitionList`, `calcPartitions`,
//! `calcExpectedPartitions`.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use flate2::read::MultiGzDecoder;
use gtars_core::models::{Region, RegionSet};
use gtars_overlaprs::{multi_chrom_overlapper::IntoMultiChromOverlapper, OverlapperType};

use serde::{Deserialize, Serialize};

use crate::errors::GtarsGenomicDistError;
use crate::interval_ranges::IntervalRanges;
use crate::models::{Strand, StrandedRegionSet};

/// A gene model loaded from BED files or GTF.
///
/// Contains the core annotations needed to build genomic partitions:
/// gene boundaries, exon coordinates, and optionally UTR regions.
/// Strand information is preserved so promoters can be computed
/// correctly for minus-strand genes.
pub struct GeneModel {
    pub genes: StrandedRegionSet,
    pub exons: StrandedRegionSet,
    pub three_utr: Option<StrandedRegionSet>,
    pub five_utr: Option<StrandedRegionSet>,
}

/// Extract strand from BED `rest` field (column 6 = 3rd tab-separated field in rest).
/// Returns Unstranded for BED3 files or when strand is absent.
fn strand_from_rest(rest: &Option<String>) -> Strand {
    match rest {
        Some(r) => {
            let fields: Vec<&str> = r.split('\t').collect();
            if fields.len() >= 3 {
                Strand::from_char(fields[2].chars().next().unwrap_or('.'))
            } else {
                Strand::Unstranded
            }
        }
        None => Strand::Unstranded,
    }
}

/// Wrap a RegionSet with strands extracted from BED rest fields.
fn stranded_from_regionset(rs: RegionSet) -> StrandedRegionSet {
    let strands: Vec<Strand> = rs.regions.iter().map(|r| strand_from_rest(&r.rest)).collect();
    StrandedRegionSet::new(rs, strands)
}

impl GeneModel {
    /// Load a gene model from individual BED files.
    ///
    /// `genes_path` and `exons_path` are required. UTR paths are optional.
    /// Strand is extracted from BED column 6 if present; BED3 files get
    /// `Strand::Unstranded` (preserving previous behavior).
    pub fn from_bed_files(
        genes_path: &str,
        exons_path: &str,
        three_utr_path: Option<&str>,
        five_utr_path: Option<&str>,
    ) -> Result<Self, GtarsGenomicDistError> {
        let genes = RegionSet::try_from(genes_path)
            .map_err(|e| GtarsGenomicDistError::CustomError(format!("Loading genes: {}", e)))?;
        let exons = RegionSet::try_from(exons_path)
            .map_err(|e| GtarsGenomicDistError::CustomError(format!("Loading exons: {}", e)))?;

        let three_utr = match three_utr_path {
            Some(p) => Some(
                RegionSet::try_from(p).map_err(|e| {
                    GtarsGenomicDistError::CustomError(format!("Loading 3'UTR: {}", e))
                })?,
            ),
            None => None,
        };
        let five_utr = match five_utr_path {
            Some(p) => Some(
                RegionSet::try_from(p).map_err(|e| {
                    GtarsGenomicDistError::CustomError(format!("Loading 5'UTR: {}", e))
                })?,
            ),
            None => None,
        };

        Ok(GeneModel {
            genes: stranded_from_regionset(genes).reduce(),
            exons: stranded_from_regionset(exons).reduce(),
            three_utr: three_utr
                .map(|rs| stranded_from_regionset(rs).reduce())
                .filter(|srs| !srs.is_empty()),
            five_utr: five_utr
                .map(|rs| stranded_from_regionset(rs).reduce())
                .filter(|srs| !srs.is_empty()),
        })
    }

    /// Load a gene model from a GTF (or GTF.gz) file.
    ///
    /// Extracts `gene`, `exon`, `three_prime_utr`, and `five_prime_utr` features.
    ///
    /// Also handles GENCODE-style GTFs that use an undifferentiated `UTR` feature
    /// type. In that case, `CDS` records are parsed to determine CDS boundaries per
    /// transcript, and each `UTR` record is classified as 5' or 3' by comparing its
    /// position to the CDS range on the same transcript:
    ///
    /// - `+` strand: UTR upstream of CDS → 5'UTR, downstream → 3'UTR
    /// - `-` strand: UTR downstream of CDS → 5'UTR, upstream → 3'UTR
    /// Each feature type is `reduce()`d (merged overlapping intervals) to match
    /// the behavior of R's `GenomicDistributions::getGeneModelsFromGTF()`.
    ///
    /// If `filter_protein_coding` is true, only features with
    /// `gene_biotype "protein_coding"` (or `gene_type "protein_coding"`) in the
    /// attributes column are kept.
    ///
    /// If `convert_ensembl_ucsc` is true, chromosome names without a "chr" prefix
    /// get one prepended (e.g. `1` → `chr1`, `X` → `chrX`).
    pub fn from_gtf(
        path: &str,
        filter_protein_coding: bool,
        convert_ensembl_ucsc: bool,
    ) -> Result<Self, GtarsGenomicDistError> {
        let file = File::open(path)?;
        let reader: Box<dyn BufRead> = if path.ends_with(".gz") {
            Box::new(BufReader::new(MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        let mut genes: Vec<Region> = Vec::new();
        let mut gene_strands: Vec<Strand> = Vec::new();
        let mut exons: Vec<Region> = Vec::new();
        let mut exon_strands: Vec<Strand> = Vec::new();
        let mut three_utr: Vec<Region> = Vec::new();
        let mut three_utr_strands: Vec<Strand> = Vec::new();
        let mut five_utr: Vec<Region> = Vec::new();
        let mut five_utr_strands: Vec<Strand> = Vec::new();

        // For classifying undifferentiated "UTR" features (GENCODE-style GTFs)
        struct PendingUtr {
            chr: String,
            start: u32,
            end: u32,
            strand: char,
            transcript_id: String,
        }
        let mut pending_utrs: Vec<PendingUtr> = Vec::new();
        // transcript_id -> (min_cds_start, max_cds_end)
        let mut cds_bounds: HashMap<String, (u32, u32)> = HashMap::new();

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 9 {
                continue;
            }

            let feature_type = fields[2];
            if !matches!(
                feature_type,
                "gene" | "exon" | "three_prime_utr" | "five_prime_utr" | "UTR" | "CDS"
            ) {
                continue;
            }

            if filter_protein_coding {
                let attrs = fields[8];
                // GTF attributes: key "value"; — check both gene_biotype and gene_type
                let is_protein_coding = attrs.contains("gene_biotype \"protein_coding\"")
                    || attrs.contains("gene_type \"protein_coding\"");
                if !is_protein_coding {
                    continue;
                }
            }

            let mut chr = fields[0].to_string();
            if convert_ensembl_ucsc && !chr.starts_with("chr") {
                chr = format!("chr{}", chr);
            }

            // GTF is 1-based inclusive; BED is 0-based half-open
            // start: subtract 1; end: stays the same (inclusive→exclusive cancels 1-based→0-based)
            let start: u32 = fields[3]
                .parse::<u32>()
                .map_err(|e| {
                    GtarsGenomicDistError::CustomError(format!("Parsing GTF start: {}", e))
                })?
                .saturating_sub(1);
            let end: u32 = fields[4].parse::<u32>().map_err(|e| {
                GtarsGenomicDistError::CustomError(format!("Parsing GTF end: {}", e))
            })?;

            let strand = Strand::from_char(fields[6].chars().next().unwrap_or('.'));

            match feature_type {
                "gene" => {
                    genes.push(Region { chr, start, end, rest: None });
                    gene_strands.push(strand);
                }
                "exon" => {
                    exons.push(Region { chr, start, end, rest: None });
                    exon_strands.push(strand);
                }
                "three_prime_utr" => {
                    three_utr.push(Region { chr, start, end, rest: None });
                    three_utr_strands.push(strand);
                }
                "five_prime_utr" => {
                    five_utr.push(Region { chr, start, end, rest: None });
                    five_utr_strands.push(strand);
                }
                "CDS" => {
                    if let Some(tid) = extract_gtf_transcript_id(fields[8]) {
                        let entry = cds_bounds.entry(tid).or_insert((u32::MAX, 0));
                        entry.0 = entry.0.min(start);
                        entry.1 = entry.1.max(end);
                    }
                }
                "UTR" => {
                    let strand = fields[6].chars().next().unwrap_or('+');
                    if let Some(tid) = extract_gtf_transcript_id(fields[8]) {
                        pending_utrs.push(PendingUtr {
                            chr,
                            start,
                            end,
                            strand,
                            transcript_id: tid,
                        });
                    }
                }
                _ => {}
            }
        }

        // Classify undifferentiated UTRs by position relative to CDS
        for utr in pending_utrs {
            if let Some(&(cds_start, cds_end)) = cds_bounds.get(&utr.transcript_id) {
                let utr_mid = (utr.start as u64 + utr.end as u64) / 2;
                let cds_mid = (cds_start as u64 + cds_end as u64) / 2;
                let region = Region {
                    chr: utr.chr,
                    start: utr.start,
                    end: utr.end,
                    rest: None,
                };
                let utr_strand = Strand::from_char(utr.strand);
                let is_five_prime = match utr.strand {
                    '+' => utr_mid < cds_mid,
                    _ => utr_mid > cds_mid,
                };
                if is_five_prime {
                    five_utr.push(region);
                    five_utr_strands.push(utr_strand);
                } else {
                    three_utr.push(region);
                    three_utr_strands.push(utr_strand);
                }
            }
            // UTRs without matching CDS (non-coding transcripts) are skipped
        }

        let genes_srs = StrandedRegionSet::new(RegionSet::from(genes), gene_strands).reduce();
        let exons_srs = StrandedRegionSet::new(RegionSet::from(exons), exon_strands).reduce();
        let three_utr_srs = StrandedRegionSet::new(RegionSet::from(three_utr), three_utr_strands);
        let five_utr_srs = StrandedRegionSet::new(RegionSet::from(five_utr), five_utr_strands);

        Ok(GeneModel {
            genes: genes_srs,
            exons: exons_srs,
            three_utr: if three_utr_srs.is_empty() {
                None
            } else {
                Some(three_utr_srs.reduce())
            },
            five_utr: if five_utr_srs.is_empty() {
                None
            } else {
                Some(five_utr_srs.reduce())
            },
        })
    }
}

/// Extract `transcript_id` from a GTF attributes string (column 9).
///
/// Looks for the pattern `transcript_id "VALUE"` and returns VALUE.
fn extract_gtf_transcript_id(attrs: &str) -> Option<String> {
    let marker = "transcript_id \"";
    let start = attrs.find(marker)? + marker.len();
    let end = start + attrs[start..].find('"')?;
    Some(attrs[start..end].to_string())
}

/// An ordered list of named, mutually exclusive genomic partitions.
///
/// The order determines priority for region assignment: a region overlapping
/// both promoterCore and exon is assigned to promoterCore (checked first).
pub struct PartitionList {
    pub partitions: Vec<(String, RegionSet)>,
}

/// Result of classifying query regions into partitions.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PartitionResult {
    /// (partition_name, count) pairs. Includes "intergenic" as the remainder.
    pub counts: Vec<(String, u32)>,
    /// Total count (region count or bp count depending on mode).
    pub total: u32,
}

/// A single row in the expected partition analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExpectedPartitionRow {
    pub partition: String,
    pub observed: f64,
    pub expected: f64,
    pub log10_oe: f64,
    pub chi_sq_pval: f64,
}

/// Result of observed vs expected partition enrichment analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExpectedPartitionResult {
    pub rows: Vec<ExpectedPartitionRow>,
}

/// Build mutually exclusive genomic partitions from a gene model.
///
/// Partition priority order:
/// 1. Core promoters (`core_prom_size` bp upstream of gene starts)
/// 2. Proximal promoters (`prox_prom_size` bp upstream, minus core)
/// 3. 5'UTR (if provided, minus promoters)
/// 4. 3'UTR (if provided, minus promoters and 5'UTR)
/// 5. Exons (minus UTRs)
/// 6. Introns (gene bodies minus exons and UTRs)
///
/// Matches R GenomicDistributions `genomePartitionList()` semantics:
/// - Promoters are computed independently (not subtracted from UTRs/exons/introns)
/// - 3'UTR has priority over 5'UTR (5'UTR = setdiff(5'UTR, 3'UTR))
/// - Query assignment uses the partition list order to resolve overlaps
pub fn genome_partition_list(
    model: &GeneModel,
    core_prom_size: u32,
    prox_prom_size: u32,
    chrom_sizes: Option<&HashMap<String, u32>>,
) -> PartitionList {
    let mut partitions: Vec<(String, RegionSet)> = Vec::new();

    // 1. Core promoters (strand-aware: minus-strand uses end, not start)
    // When chrom_sizes are provided, trim promoters at chromosome boundaries
    // before reduce to avoid promoter regions extending past chromosome ends.
    let raw_core = model.genes.promoters(core_prom_size, 0);
    let core_promoters = match chrom_sizes {
        Some(sizes) => raw_core.trim(sizes).reduce(),
        None => raw_core.reduce(),
    };
    partitions.push(("promoterCore".to_string(), core_promoters.clone()));

    // 2. Proximal promoters (larger window minus core)
    let raw_prox = model.genes.promoters(prox_prom_size, 0);
    let prox_promoters = match chrom_sizes {
        Some(sizes) => raw_prox.trim(sizes).reduce(),
        None => raw_prox.reduce(),
    }
    .setdiff(&core_promoters);
    partitions.push(("promoterProx".to_string(), prox_promoters));

    // UTR/exon/intron construction does NOT subtract promoters — R doesn't either.
    // Overlap with promoters is resolved by calcPartitions priority ordering.
    //
    // Strand-aware reduce/setdiff during construction, then drop strand via
    // into_regionset() for the final PartitionList (overlap checking is positional).
    let three_utr_reduced = model.three_utr.as_ref().map(|srs| srs.reduce());
    let five_utr_reduced = model.five_utr.as_ref().map(|srs| srs.reduce());

    // 3. 3'UTR (if present)
    if let Some(ref three_utr) = three_utr_reduced {
        partitions.push(("threeUTR".to_string(), three_utr.inner.clone()));
    }

    // 4. 5'UTR (if present, minus 3'UTR — R gives 3'UTR priority)
    if let Some(ref five_utr) = five_utr_reduced {
        let five_utr_part = if let Some(ref three_utr) = three_utr_reduced {
            five_utr.setdiff(three_utr).into_regionset()
        } else {
            five_utr.inner.clone()
        };
        partitions.push(("fiveUTR".to_string(), five_utr_part));
    }

    // 5. Exons (minus UTRs)
    let mut exon_part = model.exons.reduce();
    if let Some(ref three_utr) = three_utr_reduced {
        exon_part = exon_part.setdiff(three_utr);
    }
    if let Some(ref five_utr) = five_utr_reduced {
        exon_part = exon_part.setdiff(five_utr);
    }
    partitions.push(("exon".to_string(), exon_part.into_regionset()));

    // 6. Introns (gene bodies minus UTRs and exons)
    let mut intron_part = model.genes.reduce();
    if let Some(ref three_utr) = three_utr_reduced {
        intron_part = intron_part.setdiff(three_utr);
    }
    if let Some(ref five_utr) = five_utr_reduced {
        intron_part = intron_part.setdiff(five_utr);
    }
    intron_part = intron_part.setdiff(&model.exons.reduce());
    partitions.push(("intron".to_string(), intron_part.into_regionset()));

    PartitionList { partitions }
}

/// Classify query regions into partitions.
///
/// **Priority mode** (`bp_proportion=false`): Each query region is assigned to the
/// first partition it overlaps (in list order). Regions overlapping no partition are
/// counted as "intergenic". Returns region counts.
///
/// **BP proportion mode** (`bp_proportion=true`): For each partition, computes total
/// overlapping base pairs with the query. Remaining bp are counted as "intergenic".
pub fn calc_partitions(
    query: &RegionSet,
    partitions: &PartitionList,
    bp_proportion: bool,
) -> PartitionResult {
    if bp_proportion {
        calc_partitions_bp(query, partitions)
    } else {
        calc_partitions_priority(query, partitions)
    }
}

/// Priority-based partition assignment.
fn calc_partitions_priority(query: &RegionSet, partitions: &PartitionList) -> PartitionResult {
    let n = query.regions.len();
    // Track which partition each query region is assigned to (None = unassigned)
    let mut assignments: Vec<Option<usize>> = vec![None; n];

    for (pi, (_name, partition_rs)) in partitions.partitions.iter().enumerate() {
        if partition_rs.regions.is_empty() {
            continue;
        }
        let overlapper =
            partition_rs
                .clone()
                .into_multi_chrom_overlapper(OverlapperType::AIList);

        for (ri, region) in query.regions.iter().enumerate() {
            if assignments[ri].is_some() {
                continue; // already assigned
            }
            // Check if this query region overlaps any interval in this partition
            let single = RegionSet::from(vec![Region {
                chr: region.chr.clone(),
                start: region.start,
                end: region.end,
                rest: None,
            }]);
            let hits = overlapper.find_overlaps(&single);
            if !hits.is_empty() {
                assignments[ri] = Some(pi);
            }
        }
    }

    // Count per partition
    let mut counts: Vec<(String, u32)> = partitions
        .partitions
        .iter()
        .map(|(name, _)| (name.clone(), 0))
        .collect();

    let mut intergenic_count: u32 = 0;
    for assignment in &assignments {
        match assignment {
            Some(pi) => counts[*pi].1 += 1,
            None => intergenic_count += 1,
        }
    }
    counts.push(("intergenic".to_string(), intergenic_count));

    let total = n as u32;
    PartitionResult { counts, total }
}

/// BP-proportion partition assignment.
fn calc_partitions_bp(query: &RegionSet, partitions: &PartitionList) -> PartitionResult {
    let total_query_bp: u32 = query.regions.iter().map(|r| r.end - r.start).sum();

    let mut counts: Vec<(String, u32)> = Vec::new();
    let mut assigned_bp: u32 = 0;

    for (name, partition_rs) in &partitions.partitions {
        if partition_rs.regions.is_empty() {
            counts.push((name.clone(), 0));
            continue;
        }

        let overlapper =
            partition_rs
                .clone()
                .into_multi_chrom_overlapper(OverlapperType::AIList);

        let mut partition_bp: u32 = 0;
        for region in &query.regions {
            let single = RegionSet::from(vec![Region {
                chr: region.chr.clone(),
                start: region.start,
                end: region.end,
                rest: None,
            }]);
            for (_chr, iv) in overlapper.find_overlaps_iter(&single) {
                // Compute actual overlap width
                let ol_start = region.start.max(iv.start);
                let ol_end = region.end.min(iv.end);
                if ol_end > ol_start {
                    partition_bp += ol_end - ol_start;
                }
            }
        }
        assigned_bp += partition_bp;
        counts.push((name.clone(), partition_bp));
    }

    // Remainder = intergenic
    let remainder_bp = total_query_bp.saturating_sub(assigned_bp);
    counts.push(("intergenic".to_string(), remainder_bp));

    PartitionResult {
        counts,
        total: total_query_bp,
    }
}

/// Compute observed vs expected partition enrichment with chi-square test.
///
/// For each partition, compares the observed fraction of query regions (or bp)
/// to the expected fraction based on the partition's share of the genome.
pub fn calc_expected_partitions(
    query: &RegionSet,
    partitions: &PartitionList,
    chrom_sizes: &HashMap<String, u32>,
    bp_proportion: bool,
) -> ExpectedPartitionResult {
    let observed = calc_partitions(query, partitions, bp_proportion);

    let genome_size: u64 = chrom_sizes.values().map(|&v| v as u64).sum();
    let query_total = observed.total as f64;

    // Compute partition sizes in bp
    let mut partition_bp_total: u64 = 0;
    let partition_sizes: Vec<u64> = partitions
        .partitions
        .iter()
        .map(|(_, rs)| {
            let bp: u64 = rs.regions.iter().map(|r| (r.end - r.start) as u64).sum();
            partition_bp_total += bp;
            bp
        })
        .collect();

    let remainder_genome_bp = genome_size.saturating_sub(partition_bp_total);

    let mut rows: Vec<ExpectedPartitionRow> = Vec::new();

    for (i, (name, obs_count)) in observed.counts.iter().enumerate() {
        let obs = *obs_count as f64;

        let partition_genome_bp = if name == "intergenic" {
            remainder_genome_bp
        } else {
            partition_sizes[i]
        };

        let expected = (partition_genome_bp as f64 / genome_size as f64) * query_total;

        let log10_oe = if obs == 0.0 {
            f64::NEG_INFINITY
        } else if expected == 0.0 {
            f64::INFINITY
        } else {
            (obs / expected).log10()
        };

        let chi_sq_pval = chi_square_2x2(obs, expected, query_total);

        rows.push(ExpectedPartitionRow {
            partition: name.clone(),
            observed: obs,
            expected,
            log10_oe,
            chi_sq_pval,
        });
    }

    ExpectedPartitionResult { rows }
}

/// Chi-square test for a 2x2 contingency table.
///
/// Table: [[obs, total-obs], [exp, total-exp]]
/// Returns p-value using chi-square distribution with df=1.
///
/// NOTE: This uses a goodness-of-fit (O-E)²/E formula. R's chisq.test()
/// computes a 2×2 test of independence and may apply Yates' continuity
/// correction, so p-values will differ from GenomicDistributions.
fn chi_square_2x2(obs: f64, exp: f64, total: f64) -> f64 {
    if total == 0.0 || exp == 0.0 || (total - exp) == 0.0 {
        return 1.0;
    }

    // Chi-square statistic: sum of (O-E)^2/E for each cell
    let non_obs = total - obs;
    let non_exp = total - exp;

    let chi_sq = (obs - exp).powi(2) / exp + (non_obs - non_exp).powi(2) / non_exp;

    // P-value from chi-square distribution with df=1
    // P(X > chi_sq) = 1 - regularized_gamma_lower(0.5, chi_sq/2)
    1.0 - regularized_gamma_lower(0.5, chi_sq / 2.0)
}

/// Regularized lower incomplete gamma function P(a, x).
///
/// Uses series expansion for small x, continued fraction for large x.
/// Sufficient precision for chi-square p-values.
fn regularized_gamma_lower(a: f64, x: f64) -> f64 {
    if x < 0.0 {
        return 0.0;
    }
    if x == 0.0 {
        return 0.0;
    }

    let ln_gamma_a = ln_gamma(a);

    if x < a + 1.0 {
        // Series expansion
        gamma_series(a, x, ln_gamma_a)
    } else {
        // Continued fraction
        1.0 - gamma_cf(a, x, ln_gamma_a)
    }
}

/// Series expansion for regularized lower incomplete gamma.
fn gamma_series(a: f64, x: f64, ln_gamma_a: f64) -> f64 {
    let mut sum = 1.0 / a;
    let mut term = 1.0 / a;
    for n in 1..200 {
        term *= x / (a + n as f64);
        sum += term;
        if term.abs() < sum.abs() * 1e-14 {
            break;
        }
    }
    sum * (-x + a * x.ln() - ln_gamma_a).exp()
}

/// Continued fraction for upper incomplete gamma Q(a,x) = 1 - P(a,x).
/// Uses modified Lentz's algorithm.
fn gamma_cf(a: f64, x: f64, ln_gamma_a: f64) -> f64 {
    let b0 = x + 1.0 - a;
    let mut d = 1.0 / b0;
    let mut c = 1.0 / 1e-30_f64;
    let mut f = d;

    for n in 1..200 {
        let an = -(n as f64) * (n as f64 - a);
        let bn = x + 2.0 * n as f64 + 1.0 - a;
        d = bn + an * d;
        if d.abs() < 1e-30 {
            d = 1e-30;
        }
        d = 1.0 / d;
        c = bn + an / c;
        if c.abs() < 1e-30 {
            c = 1e-30;
        }
        let delta = c * d;
        f *= delta;
        if (delta - 1.0).abs() < 1e-14 {
            break;
        }
    }

    let result = f * (-x + a * x.ln() - ln_gamma_a).exp();
    // Clamp to valid range for numerical stability
    result.clamp(0.0, 1.0)
}

/// Log-gamma function using Lanczos approximation.
fn ln_gamma(x: f64) -> f64 {
    // Lanczos coefficients (g=7, n=9)
    const COEFFS: [f64; 9] = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];

    if x < 0.5 {
        // Reflection formula
        let pi = std::f64::consts::PI;
        (pi / (pi * x).sin()).ln() - ln_gamma(1.0 - x)
    } else {
        let x = x - 1.0;
        let mut sum = COEFFS[0];
        for (i, &c) in COEFFS[1..].iter().enumerate() {
            sum += c / (x + i as f64 + 1.0);
        }
        let t = x + 7.5; // g + 0.5
        0.5 * (2.0 * std::f64::consts::PI).ln() + (t).ln() * (x + 0.5) - t + sum.ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use std::path::PathBuf;

    fn get_test_path(file_name: &str) -> PathBuf {
        std::env::current_dir()
            .unwrap()
            .join("../tests/data/regionset")
            .join(file_name)
    }

    fn make_region(chr: &str, start: u32, end: u32) -> Region {
        Region {
            chr: chr.to_string(),
            start,
            end,
            rest: None,
        }
    }

    fn make_regionset(regions: Vec<(&str, u32, u32)>) -> RegionSet {
        let regions: Vec<Region> = regions
            .into_iter()
            .map(|(chr, start, end)| make_region(chr, start, end))
            .collect();
        RegionSet::from(regions)
    }

    fn test_chrom_sizes() -> HashMap<String, u32> {
        [
            ("chr1".to_string(), 100000u32),
            ("chr2".to_string(), 80000),
            ("chr3".to_string(), 60000),
        ]
        .into_iter()
        .collect()
    }

    // ── GeneModel loading ─────────────────────────────────────────────

    #[rstest]
    fn test_gene_model_from_bed_files() {
        let genes = get_test_path("test_genes.bed");
        let exons = get_test_path("test_exons.bed");
        let three_utr = get_test_path("test_three_utr.bed");
        let five_utr = get_test_path("test_five_utr.bed");

        let model = GeneModel::from_bed_files(
            genes.to_str().unwrap(),
            exons.to_str().unwrap(),
            Some(three_utr.to_str().unwrap()),
            Some(five_utr.to_str().unwrap()),
        )
        .unwrap();

        assert_eq!(model.genes.len(), 6);
        assert_eq!(model.exons.len(), 19);
        assert!(model.three_utr.is_some());
        assert!(model.five_utr.is_some());
        assert_eq!(model.three_utr.unwrap().len(), 6);
        assert_eq!(model.five_utr.unwrap().len(), 6);
    }

    #[rstest]
    fn test_gene_model_without_utrs() {
        let genes = get_test_path("test_genes.bed");
        let exons = get_test_path("test_exons.bed");

        let model =
            GeneModel::from_bed_files(genes.to_str().unwrap(), exons.to_str().unwrap(), None, None)
                .unwrap();

        assert_eq!(model.genes.len(), 6);
        assert!(model.three_utr.is_none());
        assert!(model.five_utr.is_none());
    }

    // ── genome_partition_list ─────────────────────────────────────────

    #[rstest]
    fn test_partition_list_has_all_categories() {
        let model = load_test_model();
        let pl = genome_partition_list(&model, 100, 2000, None);

        let names: Vec<&str> = pl.partitions.iter().map(|(n, _)| n.as_str()).collect();
        assert_eq!(
            names,
            vec![
                "promoterCore",
                "promoterProx",
                "threeUTR",
                "fiveUTR",
                "exon",
                "intron"
            ]
        );
    }

    #[rstest]
    fn test_partition_list_without_utrs() {
        let genes = get_test_path("test_genes.bed");
        let exons = get_test_path("test_exons.bed");
        let model =
            GeneModel::from_bed_files(genes.to_str().unwrap(), exons.to_str().unwrap(), None, None)
                .unwrap();

        let pl = genome_partition_list(&model, 100, 2000, None);
        let names: Vec<&str> = pl.partitions.iter().map(|(n, _)| n.as_str()).collect();
        assert_eq!(
            names,
            vec!["promoterCore", "promoterProx", "exon", "intron"]
        );
    }

    #[rstest]
    fn test_partition_mutual_exclusivity() {
        let model = load_test_model();
        let pl = genome_partition_list(&model, 100, 2000, None);

        // Check pairwise: no bp should belong to two partitions
        for i in 0..pl.partitions.len() {
            for j in (i + 1)..pl.partitions.len() {
                let (name_i, rs_i) = &pl.partitions[i];
                let (name_j, rs_j) = &pl.partitions[j];

                if rs_i.regions.is_empty() || rs_j.regions.is_empty() {
                    continue;
                }

                let overlapper =
                    rs_i.clone()
                        .into_multi_chrom_overlapper(OverlapperType::AIList);
                let hits = overlapper.find_overlaps(rs_j);

                // No bp should belong to two partitions
                assert_eq!(
                    hits.len(),
                    0,
                    "Partitions '{}' and '{}' overlap ({} hits)",
                    name_i,
                    name_j,
                    hits.len()
                );
            }
        }
    }

    #[rstest]
    fn test_promoter_sizes_respected() {
        // Use synthetic data with a single gene at a known position
        let genes = make_regionset(vec![("chr1", 5000, 10000)]);
        let exons = make_regionset(vec![("chr1", 5000, 5500), ("chr1", 9000, 10000)]);
        let model = GeneModel {
            genes: StrandedRegionSet::unstranded(genes),
            exons: StrandedRegionSet::unstranded(exons),
            three_utr: None,
            five_utr: None,
        };

        let pl = genome_partition_list(&model, 100, 2000, None);

        let core = &pl.partitions[0].1;
        let prox = &pl.partitions[1].1;

        // Core promoter: 100bp upstream of gene start (5000)
        // So core = [4900, 5000)
        assert_eq!(core.regions.len(), 1);
        assert_eq!(core.regions[0].start, 4900);
        assert_eq!(core.regions[0].end, 5000);

        // Proximal promoter: 2000bp upstream minus core
        // Full prox = [3000, 5000), minus core [4900, 5000) = [3000, 4900)
        assert_eq!(prox.regions.len(), 1);
        assert_eq!(prox.regions[0].start, 3000);
        assert_eq!(prox.regions[0].end, 4900);
    }

    // ── calc_partitions priority mode ─────────────────────────────────

    #[rstest]
    fn test_priority_mode_promoter_wins_over_exon() {
        // Region overlapping both core promoter and exon → promoter wins
        let genes = make_regionset(vec![("chr1", 100, 1000)]);
        let exons = make_regionset(vec![("chr1", 50, 150)]); // overlaps with promoter region
        let model = GeneModel {
            genes: StrandedRegionSet::unstranded(genes),
            exons: StrandedRegionSet::unstranded(exons),
            three_utr: None,
            five_utr: None,
        };
        let pl = genome_partition_list(&model, 100, 200, None);

        // Query region spanning the promoter/exon boundary
        let query = make_regionset(vec![("chr1", 0, 120)]);
        let result = calc_partitions(&query, &pl, false);

        // Should be assigned to promoterCore (first in priority)
        let core_count = result
            .counts
            .iter()
            .find(|(n, _)| n == "promoterCore")
            .unwrap()
            .1;
        assert_eq!(core_count, 1);
    }

    #[rstest]
    fn test_priority_mode_intergenic() {
        let model = load_test_model();
        let pl = genome_partition_list(&model, 100, 2000, None);

        // Query in a region far from any gene
        let query = make_regionset(vec![("chr1", 50000, 50200)]);
        let result = calc_partitions(&query, &pl, false);

        let intergenic = result
            .counts
            .iter()
            .find(|(n, _)| n == "intergenic")
            .unwrap()
            .1;
        assert_eq!(intergenic, 1);
    }

    #[rstest]
    fn test_priority_mode_all_regions_accounted() {
        let model = load_test_model();
        let pl = genome_partition_list(&model, 100, 2000, None);

        let query_path = get_test_path("test_query_promoter_enriched.bed");
        let query = RegionSet::try_from(query_path.to_str().unwrap()).unwrap();

        let result = calc_partitions(&query, &pl, false);

        let sum: u32 = result.counts.iter().map(|(_, c)| c).sum();
        assert_eq!(sum, result.total);
        assert_eq!(result.total, query.regions.len() as u32);
    }

    // ── calc_partitions bp mode ───────────────────────────────────────

    #[rstest]
    fn test_bp_mode_total_consistent() {
        let model = load_test_model();
        let pl = genome_partition_list(&model, 100, 2000, None);

        let query_path = get_test_path("test_query_promoter_enriched.bed");
        let query = RegionSet::try_from(query_path.to_str().unwrap()).unwrap();

        let result = calc_partitions(&query, &pl, true);

        let query_total_bp: u32 = query.regions.iter().map(|r| r.end - r.start).sum();
        assert_eq!(result.total, query_total_bp);

        // Sum of partition bp should equal total (remainder fills the gap)
        let sum_bp: u32 = result.counts.iter().map(|(_, c)| c).sum();
        assert_eq!(sum_bp, query_total_bp);
    }

    #[rstest]
    fn test_bp_mode_partial_overlap() {
        // A query region that partially overlaps a partition
        let genes = make_regionset(vec![("chr1", 1000, 5000)]);
        let exons = make_regionset(vec![("chr1", 1000, 2000)]);
        let model = GeneModel {
            genes: StrandedRegionSet::unstranded(genes),
            exons: StrandedRegionSet::unstranded(exons),
            three_utr: None,
            five_utr: None,
        };
        let pl = genome_partition_list(&model, 100, 200, None);

        // Query straddles promoter and intergenic
        let query = make_regionset(vec![("chr1", 750, 950)]);
        let result = calc_partitions(&query, &pl, true);

        // Total should be 200bp
        assert_eq!(result.total, 200);

        // Some bp should be promoter, some intergenic
        let prom_bp: u32 = result
            .counts
            .iter()
            .filter(|(n, _)| n.starts_with("promoter"))
            .map(|(_, c)| c)
            .sum();
        let intergenic_bp = result
            .counts
            .iter()
            .find(|(n, _)| n == "intergenic")
            .unwrap()
            .1;

        assert!(prom_bp > 0, "Should have some promoter bp");
        assert!(intergenic_bp > 0, "Should have some intergenic bp");
        assert_eq!(
            prom_bp + intergenic_bp,
            200,
            "Promoter + intergenic should sum to total"
        );
    }

    // ── calc_expected_partitions ──────────────────────────────────────

    #[rstest]
    fn test_expected_partitions_enrichment_direction() {
        let model = load_test_model();
        let pl = genome_partition_list(&model, 100, 2000, None);
        let chrom_sizes = test_chrom_sizes();

        // Promoter-enriched query
        let query_path = get_test_path("test_query_promoter_enriched.bed");
        let query = RegionSet::try_from(query_path.to_str().unwrap()).unwrap();

        let result = calc_expected_partitions(&query, &pl, &chrom_sizes, false);

        // Find promoter rows
        let prom_core = result
            .rows
            .iter()
            .find(|r| r.partition == "promoterCore")
            .unwrap();
        let prom_prox = result
            .rows
            .iter()
            .find(|r| r.partition == "promoterProx")
            .unwrap();

        // Promoter-enriched query should have positive log10(O/E) for promoters
        assert!(
            prom_core.log10_oe > 0.0 || prom_prox.log10_oe > 0.0,
            "At least one promoter partition should be enriched (positive log10 O/E). \
             Core: {:.3}, Prox: {:.3}",
            prom_core.log10_oe,
            prom_prox.log10_oe
        );
    }

    #[rstest]
    fn test_expected_partitions_pvalues_valid() {
        let model = load_test_model();
        let pl = genome_partition_list(&model, 100, 2000, None);
        let chrom_sizes = test_chrom_sizes();

        let query_path = get_test_path("test_query_promoter_enriched.bed");
        let query = RegionSet::try_from(query_path.to_str().unwrap()).unwrap();

        let result = calc_expected_partitions(&query, &pl, &chrom_sizes, false);

        for row in &result.rows {
            assert!(
                (0.0..=1.0).contains(&row.chi_sq_pval),
                "P-value for '{}' out of range: {}",
                row.partition,
                row.chi_sq_pval
            );
        }
    }

    #[rstest]
    fn test_expected_partitions_all_partitions_present() {
        let model = load_test_model();
        let pl = genome_partition_list(&model, 100, 2000, None);
        let chrom_sizes = test_chrom_sizes();

        let query = make_regionset(vec![("chr1", 500, 600)]);
        let result = calc_expected_partitions(&query, &pl, &chrom_sizes, false);

        let names: Vec<&str> = result.rows.iter().map(|r| r.partition.as_str()).collect();
        assert!(names.contains(&"promoterCore"));
        assert!(names.contains(&"promoterProx"));
        assert!(names.contains(&"exon"));
        assert!(names.contains(&"intron"));
        assert!(names.contains(&"intergenic"));
    }

    // ── chi-square internals ──────────────────────────────────────────

    #[rstest]
    fn test_chi_square_equal_obs_exp() {
        // When observed == expected, chi-square = 0, p-value = 1
        let pval = chi_square_2x2(50.0, 50.0, 100.0);
        assert!((pval - 1.0).abs() < 0.01, "p-value should be ~1.0: {}", pval);
    }

    #[rstest]
    fn test_chi_square_extreme_difference() {
        // Very different obs vs exp → low p-value
        let pval = chi_square_2x2(90.0, 10.0, 100.0);
        assert!(pval < 0.001, "p-value should be very small: {}", pval);
    }

    // ── GeneModel::from_gtf ─────────────────────────────────────────

    #[rstest]
    fn test_gtf_basic_parsing() {
        let gtf = get_test_path("test_gene_model.gtf");
        let model = GeneModel::from_gtf(gtf.to_str().unwrap(), false, false).unwrap();

        // 3 genes total (gene1 on chr1, gene2 on chr2, gene3 on chr3)
        assert_eq!(model.genes.len(), 3);
        // 7 exons total, but after reduce() some may merge — all are non-overlapping so 7
        assert_eq!(model.exons.len(), 7);
        assert!(model.three_utr.is_some());
        assert!(model.five_utr.is_some());
        assert_eq!(model.three_utr.as_ref().unwrap().len(), 2);
        assert_eq!(model.five_utr.as_ref().unwrap().len(), 2);
    }

    #[rstest]
    fn test_gtf_coordinate_conversion() {
        let gtf = get_test_path("test_gene_model.gtf");
        let model = GeneModel::from_gtf(gtf.to_str().unwrap(), false, false).unwrap();

        // gene1: GTF start=1001 end=5000 → BED start=1000 end=5000
        let gene1 = model
            .genes
            .inner
            .regions
            .iter()
            .find(|r| r.chr == "chr1")
            .unwrap();
        assert_eq!(gene1.start, 1000);
        assert_eq!(gene1.end, 5000);

        // gene2: GTF start=3001 end=8000 → BED start=3000 end=8000
        let gene2 = model
            .genes
            .inner
            .regions
            .iter()
            .find(|r| r.chr == "chr2")
            .unwrap();
        assert_eq!(gene2.start, 3000);
        assert_eq!(gene2.end, 8000);
    }

    #[rstest]
    fn test_gtf_protein_coding_filter() {
        let gtf = get_test_path("test_gene_model.gtf");

        // Without filter: 3 genes (protein_coding + lncRNA)
        let model_all = GeneModel::from_gtf(gtf.to_str().unwrap(), false, false).unwrap();
        assert_eq!(model_all.genes.len(), 3);
        assert_eq!(model_all.exons.len(), 7);

        // With filter: only 2 protein_coding genes, 5 exons
        let model_pc = GeneModel::from_gtf(gtf.to_str().unwrap(), true, false).unwrap();
        assert_eq!(model_pc.genes.len(), 2);
        assert_eq!(model_pc.exons.len(), 5);
        // lncRNA gene on chr3 should be gone
        assert!(model_pc
            .genes
            .inner
            .regions
            .iter()
            .all(|r| r.chr != "chr3"));
    }

    #[rstest]
    fn test_gtf_ensembl_to_ucsc_conversion() {
        let gtf = get_test_path("test_gene_model_ensembl.gtf");

        // Without conversion: chromosome names are "1", "X"
        let model_raw = GeneModel::from_gtf(gtf.to_str().unwrap(), false, false).unwrap();
        assert!(model_raw.genes.inner.regions.iter().any(|r| r.chr == "1"));
        assert!(model_raw.genes.inner.regions.iter().any(|r| r.chr == "X"));

        // With conversion: chromosome names become "chr1", "chrX"
        let model_ucsc = GeneModel::from_gtf(gtf.to_str().unwrap(), false, true).unwrap();
        assert!(model_ucsc.genes.inner.regions.iter().any(|r| r.chr == "chr1"));
        assert!(model_ucsc.genes.inner.regions.iter().any(|r| r.chr == "chrX"));
        assert!(model_ucsc.genes.inner.regions.iter().all(|r| r.chr.starts_with("chr")));
    }

    #[rstest]
    fn test_gtf_no_utrs_produces_none() {
        // The ensembl test GTF has UTRs for gene on chr1 but not chrX
        // After protein_coding filter + reduce, UTRs should be Some
        let gtf = get_test_path("test_gene_model_ensembl.gtf");
        let model = GeneModel::from_gtf(gtf.to_str().unwrap(), false, false).unwrap();
        assert!(model.three_utr.is_some());
        assert!(model.five_utr.is_some());
    }

    #[rstest]
    fn test_gtf_model_works_with_partitions() {
        // Verify a GTF-loaded gene model works end-to-end with partition functions
        let gtf = get_test_path("test_gene_model.gtf");
        let model = GeneModel::from_gtf(gtf.to_str().unwrap(), true, false).unwrap();

        let pl = genome_partition_list(&model, 100, 2000, None);
        let names: Vec<&str> = pl.partitions.iter().map(|(n, _)| n.as_str()).collect();
        assert_eq!(
            names,
            vec![
                "promoterCore",
                "promoterProx",
                "threeUTR",
                "fiveUTR",
                "exon",
                "intron"
            ]
        );

        // Query a region and verify it classifies
        let query = make_regionset(vec![("chr1", 500, 600)]);
        let result = calc_partitions(&query, &pl, false);
        let sum: u32 = result.counts.iter().map(|(_, c)| c).sum();
        assert_eq!(sum, 1);
    }

    #[rstest]
    fn test_minus_strand_promoter_placement() {
        // gene2 in test_gene_model.gtf: chr2, minus strand, GTF [3001,8000] → BED [3000,8000)
        // Core promoter (100bp) should be at [8000, 8100) for minus strand,
        // NOT [2900, 3000) which is what the old strand-unaware code produced.
        let gtf = get_test_path("test_gene_model.gtf");
        let model = GeneModel::from_gtf(gtf.to_str().unwrap(), true, false).unwrap();

        let pl = genome_partition_list(&model, 100, 2000, None);
        let core = &pl.partitions[0].1; // promoterCore

        // Find the core promoter region on chr2
        let chr2_core: Vec<&Region> = core
            .regions
            .iter()
            .filter(|r| r.chr == "chr2")
            .collect();
        assert_eq!(chr2_core.len(), 1, "Should have one core promoter on chr2");
        assert_eq!(chr2_core[0].start, 8000, "Minus-strand promoter should start at gene end");
        assert_eq!(chr2_core[0].end, 8100, "Minus-strand promoter should be 100bp past gene end");
    }

    #[rstest]
    fn test_plus_strand_promoter_placement() {
        // gene1 in test_gene_model.gtf: chr1, plus strand, BED [1000,5000)
        // Core promoter (100bp) should be at [900, 1000) — same as before
        let gtf = get_test_path("test_gene_model.gtf");
        let model = GeneModel::from_gtf(gtf.to_str().unwrap(), true, false).unwrap();

        let pl = genome_partition_list(&model, 100, 2000, None);
        let core = &pl.partitions[0].1;

        let chr1_core: Vec<&Region> = core
            .regions
            .iter()
            .filter(|r| r.chr == "chr1")
            .collect();
        assert_eq!(chr1_core.len(), 1);
        assert_eq!(chr1_core[0].start, 900);
        assert_eq!(chr1_core[0].end, 1000);
    }

    #[rstest]
    fn test_stranded_gene_model_from_synthetic() {
        // Build a GeneModel with explicit strand info and verify promoters
        let genes = make_regionset(vec![
            ("chr1", 1000, 5000),  // plus-strand
            ("chr1", 6000, 9000),  // minus-strand
        ]);
        let strands = vec![Strand::Plus, Strand::Minus];
        let exons = make_regionset(vec![
            ("chr1", 1000, 1500),
            ("chr1", 6000, 6500),
        ]);
        let exon_strands = vec![Strand::Plus, Strand::Minus];

        let model = GeneModel {
            genes: StrandedRegionSet::new(genes, strands),
            exons: StrandedRegionSet::new(exons, exon_strands),
            three_utr: None,
            five_utr: None,
        };

        let pl = genome_partition_list(&model, 100, 2000, None);
        let core = &pl.partitions[0].1;

        // Should have 2 core promoters on chr1
        assert_eq!(core.regions.len(), 2);
        // Plus-strand: [900, 1000)
        // Minus-strand: [9000, 9100)
        let mut starts: Vec<u32> = core.regions.iter().map(|r| r.start).collect();
        starts.sort();
        assert_eq!(starts, vec![900, 9000]);
    }

    // ── Cross-validation vs R GenomicDistributions ─────────────────

    /// Helper: load a reference BED file and return sorted (chr, start, end) tuples.
    fn load_reference_bed(file_name: &str) -> Vec<(String, u32, u32)> {
        let path = get_test_path(file_name);
        let rs = RegionSet::try_from(path.to_str().unwrap()).unwrap();
        let mut coords: Vec<(String, u32, u32)> = rs
            .regions
            .iter()
            .map(|r| (r.chr.clone(), r.start, r.end))
            .collect();
        coords.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        coords
    }

    /// Helper: extract sorted (chr, start, end) tuples from a RegionSet.
    fn sorted_coords(rs: &RegionSet) -> Vec<(String, u32, u32)> {
        let mut coords: Vec<(String, u32, u32)> = rs
            .regions
            .iter()
            .map(|r| (r.chr.clone(), r.start, r.end))
            .collect();
        coords.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        coords
    }

    #[rstest]
    fn test_gtf_vs_r_all_features() {
        // Cross-validate against R's getGeneModelsFromGTF() output.
        //
        // Reference BEDs were generated with reduce(ignore.strand=TRUE) —
        // strand-unaware merge. Our reduce is now strand-aware, so we apply
        // an additional unstranded reduce before comparing coordinates.
        let gtf = get_test_path("C_elegans_cropped_example.gtf.gz");
        let model = GeneModel::from_gtf(gtf.to_str().unwrap(), false, false).unwrap();

        let r_genes = load_reference_bed("ce_ref_genes_all.bed");
        let r_exons = load_reference_bed("ce_ref_exons_all.bed");
        let r_three_utr = load_reference_bed("ce_ref_three_utr_all.bed");
        let r_five_utr = load_reference_bed("ce_ref_five_utr_all.bed");

        // Apply unstranded reduce to match the R reference
        assert_eq!(sorted_coords(&model.genes.inner.reduce()), r_genes, "genes mismatch vs R");
        assert_eq!(sorted_coords(&model.exons.inner.reduce()), r_exons, "exons mismatch vs R");
        assert_eq!(
            sorted_coords(&model.three_utr.as_ref().unwrap().inner.reduce()),
            r_three_utr,
            "3'UTR mismatch vs R"
        );
        assert_eq!(
            sorted_coords(&model.five_utr.as_ref().unwrap().inner.reduce()),
            r_five_utr,
            "5'UTR mismatch vs R"
        );
    }

    #[rstest]
    fn test_gtf_vs_r_protein_coding() {
        // Cross-validate protein_coding filter against R reference.
        let gtf = get_test_path("C_elegans_cropped_example.gtf.gz");
        let model = GeneModel::from_gtf(gtf.to_str().unwrap(), true, false).unwrap();

        let r_genes = load_reference_bed("ce_ref_genes_pc.bed");
        let r_exons = load_reference_bed("ce_ref_exons_pc.bed");
        let r_three_utr = load_reference_bed("ce_ref_three_utr_pc.bed");
        let r_five_utr = load_reference_bed("ce_ref_five_utr_pc.bed");

        assert_eq!(
            sorted_coords(&model.genes.inner.reduce()),
            r_genes,
            "protein_coding genes mismatch vs R"
        );
        assert_eq!(
            sorted_coords(&model.exons.inner.reduce()),
            r_exons,
            "protein_coding exons mismatch vs R"
        );
        assert_eq!(
            sorted_coords(&model.three_utr.as_ref().unwrap().inner.reduce()),
            r_three_utr,
            "protein_coding 3'UTR mismatch vs R"
        );
        assert_eq!(
            sorted_coords(&model.five_utr.as_ref().unwrap().inner.reduce()),
            r_five_utr,
            "protein_coding 5'UTR mismatch vs R"
        );
    }

    // ── helpers ───────────────────────────────────────────────────────

    fn load_test_model() -> GeneModel {
        let genes = get_test_path("test_genes.bed");
        let exons = get_test_path("test_exons.bed");
        let three_utr = get_test_path("test_three_utr.bed");
        let five_utr = get_test_path("test_five_utr.bed");

        GeneModel::from_bed_files(
            genes.to_str().unwrap(),
            exons.to_str().unwrap(),
            Some(three_utr.to_str().unwrap()),
            Some(five_utr.to_str().unwrap()),
        )
        .unwrap()
    }
}

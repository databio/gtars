# gtars-genomicdist

Rust port of the R [GenomicDistributions](https://code.databio.org/GenomicDistributions/) package — plus a handful of extra calculations and binary serialization formats used by BEDbase. The crate provides:

- **Summary statistics** and per-chromosome distributions (`GenomicIntervalSetStatistics` trait).
- **Multi-set algebra** over a `RegionSetList` (`RegionSetListOps`, `pairwise_jaccard`). Single-set GRanges-style operations (`reduce`, `setdiff`, `promoters`, ...) now live on `RegionSet` in [gtars-core](core.md).
- **Genomic partitioning** using a gene model (promoter / UTR / exon / intron / intergenic).
- **Signal matrix overlap and summary** (`calc_summary_signal`, TSV + packed-binary I/O).
- **Consensus region calling** across multiple input sets.
- **Nearest-TSS / nearest-feature distance** calculations (`TssIndex`).
- **GC content** and **dinucleotide frequency** against a reference genome.
- **BED file classification** (optional `bedclassifier` feature, polars-backed).
- The **GDA** binary gene-model format and the **.fab** binary FASTA format for fast mmap-backed sequence lookup.

All operations take a `&RegionSet` (from `gtars-core`) as input and either extend it via trait implementations or run as free functions. If you're new to the crate, start with the `GenomicIntervalSetStatistics` trait (plus the interval set operations on `RegionSet` in gtars-core); together those cover the GRanges/GenomicDistributions surface area most users need.

!!! tip "Coming from R GenomicDistributions?"
    Most function names follow the pattern `calc_*` where R used `calc*`, and behavior is matched where practical. A few deliberate divergences are documented inline — the most notable is that midpoint calculations default to BED (floor) conventions rather than GRanges banker's rounding, and spacing/nearest-neighbor calculations *exclude* single-region chromosomes rather than emitting sentinel values.

## Installation

```toml
[dependencies]
gtars-genomicdist = "0.8"
```

### Feature flags

| flag | effect |
|---|---|
| *(default)* | All traits, statistics, partitions, signal, consensus, and TSS calculations. |
| `bedclassifier` | Enables `classify_bed()` and pulls in `polars`. Transitively enables `gtars-core/dataframe`. |

Enable the classifier feature:

```toml
[dependencies]
gtars-genomicdist = { version = "0.8", features = ["bedclassifier"] }
```

## Core types

All of these are re-exported at the crate root, so you typically import them directly:

```rust
use gtars_genomicdist::{
    // Statistics
    GenomicIntervalSetStatistics,
    // Multi-set algebra
    RegionSetListOps, pairwise_jaccard,
    // Strand-aware wrappers (SortedRegionSet is re-exported from gtars-core)
    SortedRegionSet, Strand, StrandedRegionSet,
    // Partitions
    GeneModel, PartitionList, PartitionResult,
    ExpectedPartitionResult, ExpectedPartitionRow,
    calc_expected_partitions, calc_partitions, genome_partition_list,
    // Signal
    SignalMatrix, SignalSummaryResult, ConditionStats,
    calc_summary_signal,
    // Consensus
    ConsensusRegion, consensus,
    // Sequence stats
    calc_gc_content, calc_dinucl_freq, DINUCL_ORDER,
    // Utilities
    chrom_karyotype_key, median_abs_distance,
    // Gene model serialization
    GenomicDistAnnotation,
    // BED classifier (feature-gated)
    #[cfg(feature = "bedclassifier")] classify_bed,
    // Re-export from gtars-core
    CoordinateMode,
};
```

The submodules (`statistics`, `region_set_list_ops`, `stranded_region_set`, `partitions`, `signal`, `consensus`, `models`, `asset`, `bed_classifier`, `errors`, `utils`) are all `pub` as well, so both `gtars_genomicdist::calc_partitions` and `gtars_genomicdist::partitions::calc_partitions` resolve to the same symbol.

## Statistics

The `GenomicIntervalSetStatistics` trait extends `RegionSet` with per-region and per-pair quantities — direct ports of the R GenomicDistributions functions (`calcWidth`, `calcNeighborDist`, `calcNearestNeighbors`, `calcChromBins`). These return vectors or per-chromosome tables, one number per peak or per gap. Use them when you want to *see* each value — plot a histogram, feed them to a downstream test, render a per-chromosome heatmap.

Bring the trait into scope and the methods appear on any `RegionSet`.

### Per-region and per-pair quantities

Direct ports of R GenomicDistributions. Each returns a vector or per-chromosome table.

```rust
use gtars_core::models::RegionSet;
use gtars_genomicdist::GenomicIntervalSetStatistics;

let peaks = RegionSet::try_from("peaks.bed")?;

// Per-chromosome summary: count, min/max/mean/median width, chr bounds
for (chr, stats) in peaks.chromosome_statistics() {
    println!(
        "{chr}: {} regions, median width {}",
        stats.number_of_regions,
        stats.median_region_length
    );
}

// Region-width distribution (one value per region)
let widths: Vec<u32> = peaks.calc_widths();

// Signed inter-region gaps on each chromosome (only positive gaps are returned)
let gaps: Vec<i64> = peaks.calc_neighbor_distances()?;

// Nearest-neighbor distance per region (multi-region chromosomes only)
let nn: Vec<u32> = peaks.calc_nearest_neighbors()?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

| method | question it answers | returns |
|---|---|---|
| `chromosome_statistics()` | "What do the per-chromosome counts and width distributions look like?" | `HashMap<String, ChromosomeStatistics>` |
| `region_distribution_with_bins(n_bins)` | "How many peaks fall in each of `n_bins` bins sized by the longest chromosome?" | `HashMap<String, RegionBin>` |
| `region_distribution_with_chrom_sizes(n_bins, chrom_sizes)` | Same, but bins are sized per-chromosome by actual chromosome length — matches R `getGenomeBins` | `HashMap<String, RegionBin>` |
| `calc_widths()` | "How wide is each peak?" — port of R `calcWidth()` | `Vec<u32>` |
| `calc_neighbor_distances()` | "What's the bp gap between every pair of consecutive peaks on each chromosome?" | `Result<Vec<i64>>` |
| `calc_nearest_neighbors()` | "How far is each peak from its closest neighbor?" — port of R `calcNearestNeighbors()` | `Result<Vec<u32>>` |

### Output structs

```rust
pub struct ChromosomeStatistics {
    pub chromosome: String,
    pub number_of_regions: u32,
    pub start_nucleotide_position: u32,  // leftmost start
    pub end_nucleotide_position: u32,    // rightmost end
    pub minimum_region_length: u32,
    pub maximum_region_length: u32,
    pub mean_region_length: f64,
    pub median_region_length: f64,
}

pub struct RegionBin { pub chr: String, pub start: u32, pub end: u32, pub n: u32, pub rid: u32 }
```

!!! warning "Output length for spacing/nearest calculations"
    `calc_neighbor_distances` and `calc_nearest_neighbors` **skip single-region chromosomes** (matching R GenomicDistributions behavior). The output length is therefore **not 1:1 with the input region count** — it's the number of multi-region chromosomes' regions. No sentinel values are emitted. If you need 1:1 alignment, filter your input to multi-region chromosomes first.

!!! warning "`n_bins` is a target, not a total"
    In `region_distribution_with_chrom_sizes`, `n_bins` is the target bin count for the **longest** chromosome in `chrom_sizes`. Bin width is derived from that, and every chromosome is tiled at the same bp width — so shorter chromosomes get fewer bins and the total bin count is `sum(min(n_bins, ceil(chrom_size / bin_width)))`, which can exceed `n_bins` substantially. No chromosome gets more than `n_bins` bins: the last bin absorbs any leftover tail. Regions on chromosomes missing from `chrom_sizes`, or whose midpoint lies past the stated chromosome end, are skipped. To target a specific bin width in bp, set `n_bins = max_chrom_len / desired_bp`.

## Interval set algebra

GRanges/IRanges-style set algebra now lives in **gtars-core**: most operations are inherent methods on `RegionSet`, and the two-set operations `setdiff`, `intersect`, `jaccard`, `coverage`, and `overlap_coefficient` come from the `gtars_core::models::IntervalSetOps` trait. All operations return new `RegionSet`s (immutable pattern) and are strand-unaware by default — use `StrandedRegionSet` if you need strand-aware promoters, reduce, or setdiff.

```rust
use gtars_core::models::{IntervalSetOps, RegionSet};

let a = RegionSet::try_from("peaks_a.bed")?;
let b = RegionSet::try_from("peaks_b.bed")?;

let merged   = a.reduce();                         // merge overlapping intervals
let common   = a.intersect(&b);                    // range-level intersection
let a_only   = a.setdiff(&b);                      // remove b from a
let combined = a.union(&b);                        // union + reduce
let jac      = a.jaccard(&b);                      // bp-Jaccard similarity
let closest  = a.closest(&b);                      // Vec<(a_idx, b_idx, dist)>
let clustered = a.cluster(1000);                   // per-region cluster id
# Ok::<(), Box<dyn std::error::Error>>(())
```

### Method reference

| method | description |
|---|---|
| `trim(chrom_sizes)` | clamp regions to `[0, chrom_size)`; drop regions on unknown chromosomes |
| `promoters(upstream, downstream)` | `[start - upstream, start + downstream)` per region |
| `reduce()` | merge overlapping/adjacent intervals per chromosome |
| `setdiff(other)` | remove `other` from self (`IntervalSetOps`) |
| `pintersect(other)` | *pairwise* (by index) intersection |
| `intersect(other)` | range-level intersection (`IntervalSetOps`) |
| `concat(other)` | concatenate without merging |
| `union(other)` | `concat(other).reduce()` |
| `jaccard(other)` | bp-level Jaccard `|A ∩ B| / |A ∪ B|` (`IntervalSetOps`) |
| `coverage(other)` | fraction of `self` bp covered by `other` (`IntervalSetOps`) |
| `overlap_coefficient(other)` | `|A ∩ B| / min(|A|, |B|)` (`IntervalSetOps`) |
| `shift(offset)` | translate by signed bp offset (saturating at 0) |
| `flank(width, use_start, both)` | upstream/downstream/both-side flanks |
| `resize(width, fix)` | fixed width anchored at `"start"`, `"end"`, or `"center"` |
| `narrow(start, end, width)` | GRanges-style narrow with 1-based relative coords |
| `disjoin()` | tile into non-overlapping pieces at every boundary |
| `gaps(chrom_sizes)` | peak-free regions, including leading/trailing/whole-chrom gaps |
| `closest(other)` | `Vec<(self_idx, other_idx, signed_dist)>` |
| `cluster(max_gap)` | `Vec<u32>` cluster ids in original order |

All-vs-all intersection fragments (`intersect_all`) are provided by the overlap indexes in [gtars-overlaprs](overlaprs.md).

!!! note "`rest` fields are dropped"
    Operations that merge or synthesize new intervals (reduce, setdiff, promoters, etc.) produce regions with `rest: None`. There is no unambiguous way to carry the original metadata through a merge, so the contract is: use the interval set operations for coordinate-only work.

### `pairwise_jaccard`

A standalone helper for computing the full N×N Jaccard matrix over a slice of region sets, optimized to pre-reduce each set once and walk pairs linearly:

```rust
use gtars_genomicdist::pairwise_jaccard;

let sets = vec![rs1, rs2, rs3];
let jac = pairwise_jaccard(&sets); // Vec<f64> of length 9 (row-major)

for i in 0..3 {
    for j in 0..3 {
        print!("{:.3} ", jac[i * 3 + j]);
    }
    println!();
}
```

### RegionSetListOps

`RegionSetListOps` implements the same set-algebra operations on a `RegionSetList` by index — useful for bindings (wasm/Python/R) that want to operate on pairs without copying whole `RegionSet`s across an FFI boundary.

```rust
use gtars_core::models::RegionSetList;
use gtars_genomicdist::RegionSetListOps;

let rsl = RegionSetList::try_from(vec!["a.bed", "b.bed", "c.bed"])?;

let jac_01 = rsl.jaccard_at(0, 1);                   // Option<f64>
let union_02 = rsl.union_at(0, 2);                   // Option<RegionSet>
let drop_one = rsl.union_except(1);                  // union of all but index 1
let (full_union, loo) = rsl.bulk_union_except()      // all leave-one-out unions in O(n)
    .ok_or("empty list")?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

Methods: `pintersect_at`, `pintersect_count`, `jaccard_at`, `union_at`, `setdiff_at`, `region_count`, `union_except`, `bulk_union_except`, `union_all`, `intersect_all`.

## Partitions

Ports of `genomePartitionList()`, `calcPartitions()`, and `calcExpectedPartitions()` from R GenomicDistributions. You load a gene model from BED files or GTF, build an ordered `PartitionList` of promoter / UTR / exon / intron categories, and classify your query regions into mutually exclusive buckets (`intergenic` is added as the implicit remainder).

```rust
use gtars_core::models::RegionSet;
use gtars_genomicdist::{
    GeneModel, calc_partitions, calc_expected_partitions, genome_partition_list,
};
use std::collections::HashMap;

// 1. Load the gene model (choose one)
let model = GeneModel::from_bed_files(
    "genes.bed",
    "exons.bed",
    Some("three_utr.bed"),  // optional
    Some("five_utr.bed"),   // optional
)?;

// Or, from a GTF:
// let model = GeneModel::from_gtf("gencode.v47.gtf.gz", true, true)?;

// 2. Build partition list (core/prox promoter sizes in bp; chrom_sizes optional)
let partitions = genome_partition_list(&model, 100, 2000, None);

// 3. Classify query regions
let query = RegionSet::try_from("peaks.bed")?;

let result = calc_partitions(&query, &partitions, /* bp_proportion = */ false);
for (name, count) in &result.counts {
    let pct = *count as f64 / result.total as f64 * 100.0;
    println!("{:<15} {:>6}  ({:.1}%)", name, count, pct);
}

// 4. Observed vs expected, with chi-square p-value
let chrom_sizes: HashMap<String, u32> = [("chr1".into(), 248_956_422)].into_iter().collect();
let expected = calc_expected_partitions(&query, &partitions, &chrom_sizes, false);
for row in &expected.rows {
    println!(
        "{:<15} obs={:>6}  exp={:>10.1}  log10(O/E)={:+.2}  p={:.2e}",
        row.partition, row.observed, row.expected, row.log10_oe, row.chi_sq_pval
    );
}
# Ok::<(), Box<dyn std::error::Error>>(())
```

### Priority order

`genome_partition_list` produces partitions in this order, which is also the priority for classification:

1. **`promoterCore`** — `core_prom_size` bp upstream of each gene start (strand-aware).
2. **`promoterProx`** — `prox_prom_size` bp upstream, minus core.
3. **`threeUTR`** — if UTR files are provided.
4. **`fiveUTR`** — minus 3'UTR (R gives 3'UTR priority).
5. **`exon`** — minus UTRs.
6. **`intron`** — gene bodies minus UTRs and exons.
7. **`intergenic`** — implicit remainder, emitted by `calc_partitions`.

`calc_partitions` has two modes:

- **Priority mode (`bp_proportion = false`)** — each query region is assigned to the first partition it overlaps; remainders are counted as `intergenic`. Mutually exclusive.
- **BP proportion mode (`bp_proportion = true`)** — for each partition, computes total overlapping base pairs of query. Not mutually exclusive: a query region overlapping multiple partitions contributes bp to each (matching R behavior).

!!! note "Chi-square p-values differ from R"
    `calc_expected_partitions` uses a goodness-of-fit `(O-E)²/E` formula. R's `chisq.test()` computes a 2×2 test of independence with optional Yates correction, so p-values will not match R GenomicDistributions output byte-for-byte. The `log10_oe` column is directly comparable.

### `GeneModel::from_gtf`

GTF loading handles GENCODE-style files that use an undifferentiated `UTR` feature type by parsing `CDS` records to infer which UTRs are 5' vs 3'. Two flags:

- `filter_protein_coding` — keep only features with `gene_biotype "protein_coding"` or `gene_type "protein_coding"` in the attributes column.
- `convert_ensembl_ucsc` — prepend `chr` to chromosome names that don't already have it (so Ensembl `1` becomes UCSC `chr1`).

## Signal matrix overlap

`SignalMatrix` holds a matrix of signal values across genomic regions × conditions (e.g. a peak × cell-type matrix of ChIP-seq intensities). `calc_summary_signal` overlaps a query region set against the matrix, aggregates by MAX per query region, and returns Tukey boxplot statistics per condition.

```rust
use gtars_core::models::{RegionSet, CoordinateMode};
use gtars_genomicdist::{SignalMatrix, calc_summary_signal};

// Load from a TSV where col 0 is "chr_start_end", remaining cols are condition values
let sm = SignalMatrix::from_tsv("signal_matrix.tsv")?;

// Or load from the packed binary format (produced by sm.save_bin())
let sm = SignalMatrix::load_bin("signal_matrix.sigm")?;
// ...or from an in-memory byte slice (for wasm):
// let sm = SignalMatrix::load_bin_from_bytes(&bytes)?;

let query = RegionSet::try_from("peaks.bed")?;
let summary = calc_summary_signal(&query, &sm, CoordinateMode::Bed)?;

for stats in &summary.matrix_stats {
    println!(
        "{}: median={:.2}  [{:.2}–{:.2}]",
        stats.condition, stats.median, stats.lower_hinge, stats.upper_hinge
    );
}
# Ok::<(), Box<dyn std::error::Error>>(())
```

`SignalSummaryResult` contains per-query-region max signal vectors (`signal_matrix`), per-condition Tukey stats (`matrix_stats`), and `condition_names`. The packed binary format (`.sigm`, magic `SIGM`) is produced by `sm.save_bin()` — substantially faster to load than TSV for large matrices.

## Consensus regions

Given N region sets, `consensus` produces the reduced union annotated with how many input sets overlap each union region. Useful for filtering by replicate support threshold:

```rust
use gtars_core::models::RegionSet;
use gtars_genomicdist::{consensus, ConsensusRegion};

let reps = [
    RegionSet::try_from("rep1.bed")?,
    RegionSet::try_from("rep2.bed")?,
    RegionSet::try_from("rep3.bed")?,
];

let cons: Vec<ConsensusRegion> = consensus(&reps);

// Keep regions present in ≥ 2/3 replicates
let robust: Vec<_> = cons.into_iter().filter(|c| c.count >= 2).collect();
println!("{} robust regions", robust.len());
# Ok::<(), Box<dyn std::error::Error>>(())
```

`ConsensusRegion` has `chr`, `start`, `end`, and `count` — the number of input sets overlapping that union region.

## TSS / nearest-feature distances

`TssIndex` builds a per-chromosome sorted vector of TSS midpoints from a `RegionSet`, enabling O(R · log M) distance queries instead of the naive O(R · M):

```rust
use gtars_core::models::{RegionSet, CoordinateMode};
use gtars_genomicdist::models::TssIndex;

let tss = TssIndex::try_from("gencode_tss.bed")?;

let peaks = RegionSet::try_from("peaks.bed")?;

// Unsigned nearest-feature distance per query region
let dists: Vec<u32> = tss.calc_tss_distances(&peaks, CoordinateMode::Bed)?;

// Signed distances (positive = feature downstream, negative = upstream)
let signed: Vec<i64> = tss.calc_feature_distances(&peaks, CoordinateMode::Bed)?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

Regions on chromosomes missing from the TSS index are padded with `u32::MAX` / `i64::MAX` sentinels (one per query region). Use `median_abs_distance(&signed)` from the `utils` module to summarize while filtering those sentinels.

## GC content and dinucleotide frequencies

Both functions take anything implementing `SequenceAccess` — currently `GenomeAssembly` (in-memory HashMap, built from a FASTA) or `BinaryGenomeAssembly` (mmap-backed `.fab` binary).

```rust
use gtars_core::models::RegionSet;
use gtars_genomicdist::{calc_gc_content, calc_dinucl_freq, DINUCL_ORDER};
use gtars_genomicdist::models::{GenomeAssembly, BinaryGenomeAssembly};

// In-memory loader — ~2s to build for hg38 but then gives zero-copy slices
let genome = GenomeAssembly::try_from("hg38.fa")?;

// Or the mmap .fab loader — instant construction, zero-copy per region
let genome = BinaryGenomeAssembly::from_file("hg38.fab".as_ref())?;

let peaks = RegionSet::try_from("peaks.bed")?;

// GC content per region, 0.0–1.0
let gc: Vec<f64> = calc_gc_content(&peaks, &genome, /* ignore_unk_chroms = */ true)?;

// Dinucleotide frequencies: 16 columns per region in DINUCL_ORDER
// raw_counts=false → percentages (matches R default)
let (labels, matrix) = calc_dinucl_freq(&peaks, &genome, false, true)?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

Create a `.fab` file from a regular FASTA once up-front:

```rust
use gtars_genomicdist::models::BinaryGenomeAssembly;

BinaryGenomeAssembly::write_from_fasta(
    "hg38.fa".as_ref(),
    "hg38.fab".as_ref(),
)?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

Or via the CLI: `gtars prep --fasta hg38.fa` (see the [gtars-cli](cli.md) page).

## BED classification

Under the `bedclassifier` feature, `classify_bed` inspects the column layout of a BED file and assigns one of several UCSC / ENCODE format classifications (BED3, BED6+, narrowPeak, broadPeak, gappedPeak, RNA elements, etc.).

```rust
#[cfg(feature = "bedclassifier")]
{
    use gtars_core::models::RegionSet;
    use gtars_genomicdist::bed_classifier::classify_bed;

    let rs = RegionSet::try_from("unknown_format.bed").unwrap();
    let c = classify_bed(&rs).unwrap();
    println!(
        "format={}  bed_compliance={}  compliant_cols={}",
        c.data_format, c.bed_compliance, c.compliant_columns
    );
}
```

## GDA binary format

`GenomicDistAnnotation` is a serializable wrapper around `GeneModel` in the GDA (Genomic Dist Annotation) binary format. Chrom sizes are **not** stored in the GDA file — they come from a separate source (refgenie, `.chrom.sizes` file, or an API).

```rust
use gtars_genomicdist::GenomicDistAnnotation;

// Build from GTF once, save as GDA
let gda = GenomicDistAnnotation::from_gtf("gencode.v47.gtf.gz")?;
gda.save_bin("gencode.v47.gda")?;

// Load from disk
let gda = GenomicDistAnnotation::load_bin("gencode.v47.gda")?;

// Or from memory (wasm / API)
// let gda = GenomicDistAnnotation::load_bin_from_bytes(&bytes)?;

let partitions = gtars_genomicdist::genome_partition_list(&gda.gene_model, 100, 2000, None);
# Ok::<(), Box<dyn std::error::Error>>(())
```

## Utilities

From `gtars_genomicdist::utils`:

- **`partition_genome_into_bins(chrom_sizes, n_bins)`** — tile all chromosomes with fixed-width bins (bin width = longest chrom / n_bins, floored, min 1 bp). Returns a `RegionSet`.
- **`median_abs_distance(&[i64])`** — median absolute value, filtering `i64::MAX` sentinels produced by `TssIndex::calc_feature_distances`.
- **`chrom_karyotype_key(chr)`** — sort key producing the standard karyotype order: `chr1, chr2, …, chr22, chrX, chrY, chrM, chrUn_*`. Works with or without the `chr` prefix.

## Strand-aware wrappers

Two wrappers extend `RegionSet` with stronger invariants:

- **`SortedRegionSet`** — a newtype guaranteeing `(chr, start)` sort order (defined in gtars-core, re-exported here). Constructed via `SortedRegionSet::new(rs)`, which sorts in place (move, no clone). Downstream code that requires sorted input can accept `&SortedRegionSet` to avoid re-sorting on every call.
- **`StrandedRegionSet`** — pairs a `RegionSet` with a parallel `Vec<Strand>`. Strand-aware `promoters_stranded`, `reduce`, `setdiff`, and `trim` are methods on this type and are used internally by `genome_partition_list` to produce correct partitions for minus-strand genes.

```rust
use gtars_core::models::RegionSet;
use gtars_genomicdist::{SortedRegionSet, StrandedRegionSet, Strand};

let rs = RegionSet::try_from("peaks.bed")?;
let sorted = SortedRegionSet::new(rs); // in-place sort

let stranded = StrandedRegionSet::new(
    sorted.0.clone(),
    vec![Strand::Plus; sorted.0.regions.len()],
);
# Ok::<(), Box<dyn std::error::Error>>(())
```

`Strand::from_char('+' | '-' | _)` converts a single character to `Plus` / `Minus` / `Unstranded`.

## Errors

All operations that can fail return `Result<T, GtarsGenomicDistError>`. Variants:

- `CustomError(String)` — general I/O and parsing wrapper.
- `GCContentError(chr, start, end, msg)` — sequence lookup failed during GC calculation.
- `TSSContentError(msg)` — no TSS features found for the region set.
- `SignalMatrixError(msg)` — signal-matrix parse failure.
- `Io(std::io::Error)` — transparent I/O error.

The `bedclassifier` feature has its own `BedClassifierError` enum for format-classification failures.

## Where to go next

- **[gtars-core](core.md)** — `RegionSet`, `RegionSetList`, and `CoordinateMode`, which this crate consumes.
- **[gtars-lola](lola.md)** — LOLA enrichment is built on the `RegionSet` interval operations (for universe construction) and the IGD index.
- **[gtars-overlaprs](overlaprs.md)** — the overlap-detection engine used internally by `calc_partitions` and `consensus`.
- **[gtars CLI](cli.md)** — `gtars genomicdist` subcommands for running these analyses from the command line.

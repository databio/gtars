# gtars-r

R bindings for gtars — high-performance genomic interval analysis in R, backed by Rust via `extendr`. The package is `gtars`.

## Installation

```r
install.packages("remotes")
remotes::install_github("databio/gtars", subdir = "gtars-r")
```

Or install from a local clone of the repo:

```console
R CMD INSTALL gtars-r
```

Building the package requires a Rust toolchain (`cargo`, `rustc`). See the [gtars-r README](https://github.com/databio/gtars/blob/dev/gtars-r/README.md) for detailed setup notes.

## Loading

```r
library(gtars)
```

## Submodule pages

| page | covers |
|---|---|
| [**RegionSet & RegionSetList**](r/regionset.md) | S4 classes for genomic regions, interval set algebra, overlap queries, conversion to/from `GRanges`, pairwise Jaccard. |
| [**GenomicDistributions**](r/genomicdist.md) | Drop-in replacements for R GenomicDistributions — `calcWidth`, `calcGCContent`, `calcDinuclFreq`, `genomePartitionList`, `calcPartitions`, `calcExpectedPartitions`, `calcSummarySignal`, `calcFeatureDist`, etc. |
| [**LOLA enrichment**](r/lola.md) | Drop-in replacement for R LOLA — `loadRegionDB`, `runLOLA`, `checkUniverseAppropriateness`, `redefineUserSets`, `buildRestrictedUniverse`. |
| [**IGD indexing**](r/igd.md) | Low-level `igd_create` and `igd_search` for building and querying IGD databases directly. |
| [**RefgetStore**](r/refgetstore.md) | GA4GH refget protocol client. |

## Polymorphic inputs

A defining feature of the R gtars API: nearly every function that takes a query accepts **any** of these input types interchangeably:

- a `RegionSet` S4 object,
- a `GRanges` object (from `GenomicRanges`),
- a file path to a BED / BED.gz / narrowPeak file,
- a `data.frame` with `chr` / `start` / `end` columns (and optional `strand`).

This means you can drop gtars into an existing Bioconductor workflow without explicit conversion:

```r
library(GenomicRanges)
library(gtars)

gr <- GRanges(...)

widths <- calcWidth(gr)                        # GRanges works directly
merged <- reduce("peaks.bed")                   # file path dispatches to gtars
jac <- jaccard(gr, "other_peaks.bed")           # mix GRanges with file path
```

## Quick start

```r
library(gtars)

# 1. Load regions
peaks <- RegionSet("peaks.bed")
length(peaks)
show(peaks)

# 2. Summary statistics
calcWidth(peaks) |> summary()
chromosomeStatistics(peaks)

# 3. Interval set algebra
merged <- reduce(peaks)
proms  <- promoters(peaks, upstream = 2000L, downstream = 200L)

# 4. Convert back to GRanges for downstream Bioconductor work
#    (GenomicRanges only needed here — it's a bridge, not a core dep)
library(GenomicRanges)
gr <- as_granges(peaks)
```

## LOLA enrichment quick start

```r
library(gtars)

regionDB <- loadRegionDB("LOLACore/hg38")

results <- runLOLA(
  userSets      = list(RegionSet("peaks_a.bed"), RegionSet("peaks_b.bed")),
  userUniverse  = "universe.bed",
  regionDB      = regionDB,
  redefineUserSets = TRUE,
)

head(results[order(-pValueLog)], 20)
```

## refget quick start

```r
library(gtars)

readable <- "ACGTACGT"
gtars::sha512t24u_digest(readable)
gtars::md5_digest(readable)

store <- refget_store("raw")
```

## IGD quick start

```r
library(gtars)

# Build an IGD index from a directory of BED files
igd_create(output_path = "./igd_out", filelist = "./my_beds/", db_name = "my_peaks")

# Query the index
hits <- igd_search(
  database_path = "./igd_out/my_peaks.igd",
  query_path    = "query.bed",
)
```

## Coordinate-system notes

The R package uses **0-based half-open** coordinates internally (BED convention) to match the Rust core. When you pass a `GRanges` object to a gtars function, coordinates are converted from 1-based closed to 0-based half-open automatically. Converting back via `as_granges(rs)` does the reverse.

If you're building a `RegionSet` from a `data.frame` directly, the `start` and `end` columns are expected to already be in 0-based half-open form.

## See also

- **[gtars-core](core.md)** and **[gtars-genomicdist](genomicdist.md)** — Rust reference pages with the algorithmic detail and caveats that apply identically in R.
- **[R GenomicDistributions](https://code.databio.org/GenomicDistributions/)** and **[R LOLA](https://www.bioconductor.org/packages/release/bioc/html/LOLA.html)** — the original packages gtars R is a port of.

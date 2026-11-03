# GenomicDistributions in R (gtars)

The `gtars` R package ships **drop-in replacements** for the core functions in the [R GenomicDistributions](https://code.databio.org/GenomicDistributions/) package, backed by the Rust [`gtars-genomicdist`](../genomicdist.md) crate. The function names, signatures, and return shapes match R GenomicDistributions closely enough that most analyses should port over with one-line changes.

Every function accepts polymorphic query input — `RegionSet`, `GRanges`, file path, or `data.frame` — so you can mix gtars with existing Bioconductor workflows without explicit conversion.

For algorithmic detail and caveats that apply identically in R, see [gtars-genomicdist](../genomicdist.md).

## Drop-in replacements for GenomicDistributions

| gtars R function | R GenomicDistributions equivalent | notes |
|---|---|---|
| `calcWidth(query)` | `calcWidth()` | delegates to `widths()` S4 method |
| `calcNeighborDist(query)` | `calcNeighborDist()` | signed inter-region gaps per chromosome |
| `calcNearestNeighbors(query)` | `calcNearestNeighbors()` | per-region min neighbor distance |
| `regionDistribution(query, nBins, chromSizes)` | `calcChromBins()` | bin counts; pass `chromSizes` for reference-aligned bins |
| `calcGCContent(query, ref, ignoreUnkChroms)` | `calcGCContent()` | requires a `GenomeAssembly` pointer, not BSgenome |
| `calcDinuclFreq(query, ref, rawCounts, ignoreUnknownChroms)` | `calcDinuclFreq()` | returns a data.frame with region + 16 dinucleotide columns |
| `genomePartitionList(...)` | `genomePartitionList()` | build partition list from gene model components |
| `partitionListFromGTF(path, ...)` | *(not in original)* | convenience — loads gene model from GTF in one step |
| `calcPartitions(query, partitionList, bpProportion)` | `calcPartitions()` | priority-based (or bp-proportion) partition classification |
| `calcExpectedPartitions(query, partitionList, genomeSize)` | `calcExpectedPartitions()` | observed vs. expected with chi-square p-values |
| `calcSummarySignal(query, signalMatrix)` | `calcSummarySignal()` | overlap query regions with a signal matrix |
| `calcFeatureDist(query, features)` | `calcFeatureDist()` | signed distance to nearest feature |
| `calcTSSDist(query, features)` | *(not in original)* | absolute distance variant |
| `loadGenomeAssembly(fasta_path)` | — | replaces BSgenome loading (FASTA directly) |

Lower-level snake_case loaders are also exported for pre-built `gtars prep` assets: `load_binary_genome_assembly()` (`.fab` genome), `load_gda_bin()` with `gda_partition_list()`, `gda_tss_index()`, and `gda_gene_model()`, and `load_signal_matrix_bin()` / `load_signal_matrix_tsv()`.

## Porting an analysis

Most GenomicDistributions scripts port by adding `library(gtars)` after `library(GenomicDistributions)` — gtars masks the functions with faster implementations:

```r
library(GenomicDistributions)  # optional, for plotting helpers
library(gtars)                  # masks calcWidth, calcGCContent, etc.
library(GenomicRanges)

query <- GRanges(...)           # or a file path, or a data.frame
widths <- calcWidth(query)      # now backed by Rust

# Use existing plotGenomicDist* helpers on the results
plotChromBins(regionDistribution(query, chromSizes = hg38))
```

## Basic statistics

```r
library(gtars)

# Works with any query input: RegionSet, GRanges, file path, or data.frame
query <- "peaks.bed"

w  <- calcWidth(query)                  # numeric vector of widths
nd <- calcNeighborDist(query)           # signed gaps per chromosome
nn <- calcNearestNeighbors(query)       # min neighbor distance per region
```

!!! warning "`calcNeighborDist` / `calcNearestNeighbors` output length"
    Both functions **skip chromosomes with only one region** (matching R GenomicDistributions). The returned vector is shorter than `length(query)` whenever any chromosome has a single peak, and is **not aligned 1:1** with input regions.

## Region distribution

```r
# Without chromSizes: bin size derived from the query itself
# (not comparable across files)
rd <- regionDistribution(query, nBins = 250L)

# With chromSizes: reference-aligned bins (comparable across files,
# aligned with reference genome positions)
chromSizes <- c(chr1 = 248956422L, chr2 = 242193529L)
rd <- regionDistribution(query, nBins = 250L, chromSizes = chromSizes)
```

Returns a data.table compatible with R GenomicDistributions' `plotChromBins`.

!!! warning "`nBins` is a target, not a total"
    When `chromSizes` is provided, `nBins` is the target bin count for the **longest** chromosome in `chromSizes`. Bin width is derived as `max(chromSizes) %/% nBins` (floored, minimum 1 bp), and every chromosome is tiled at the same bp width, so shorter chromosomes get proportionally fewer bins. No chromosome gets more than `nBins` bins: a midpoint in the leftover tail is folded into the last bin, which runs to the chromosome end. The genome-wide bin count can still substantially exceed `nBins` when `chromSizes` has many entries. Only bins that contain at least one region midpoint are returned, and regions on chromosomes missing from `chromSizes` (or whose midpoint lies past the chromosome end) are skipped. To target a specific bin width in bp instead, pass `nBins = max_chrom_len %/% desired_bp`.

## GC content and dinucleotides

Unlike the original R GenomicDistributions (which uses BSgenome), gtars loads a reference genome directly from a FASTA file:

```r
genome <- loadGenomeAssembly("hg38.fa")

gc <- calcGCContent(query, genome, ignoreUnkChroms = TRUE)
# numeric vector of GC content per region, 0.0 to 1.0

dinucl <- calcDinuclFreq(query, genome, rawCounts = FALSE)
# data.frame with columns: region + 16 dinucleotide columns
# rawCounts = FALSE → percentages (each row sums to 100)
# rawCounts = TRUE  → integer counts

# Pooled global counts across all regions
pooled <- colSums(calcDinuclFreq(query, genome, rawCounts = TRUE)[, -1])
```

## Partitions

Build a `PartitionList` of promoter / UTR / exon / intron / intergenic categories from gene model components, then classify your query. You can either supply each component explicitly or load them in one call from a GTF.

### From gene model components

```r
# From GRanges, file paths, or data.frames — any combination works
pl <- genomePartitionList(
  genesGR      = "genes.bed",
  exonsGR      = "exons.bed",
  threeUTRGR   = "three_utr.bed",
  fiveUTRGR    = "five_utr.bed",
  corePromSize = 100L,
  proxPromSize = 2000L,
  chromSizes   = c(chr1 = 248956422L, chr2 = 242193529L),  # optional
)
```

Strand information from the inputs (GRanges strand column, or the `strand` column of a data.frame) is used for strand-aware promoter and reduce operations — critical for getting correct promoters on minus-strand genes.

### From a GTF

```r
# Single-call convenience — loads gene model and builds partitions
pl <- partitionListFromGTF(
  "gencode.v47.annotation.gtf.gz",
  filterProteinCoding = TRUE,
  convertEnsemblUCSC  = FALSE,
  corePromSize        = 100L,
  proxPromSize        = 2000L,
  chromSizes          = chromSizes,
)
```

### Classifying query regions

```r
# Priority mode — each query region assigned to first matching partition
counts <- calcPartitions(query, pl, bpProportion = FALSE)
# data.frame with partition, Freq columns + attr(counts, "total")

total <- attr(counts, "total")
counts$pct <- counts$Freq / total * 100
print(counts)

# Observed vs expected with chi-square p-values
chromSizes <- c(chr1 = 248956422L, chr2 = 242193529L)
expected <- calcExpectedPartitions(query, pl, chromSizes, bpProportion = FALSE)
# data.frame with partition, observed, expected, log10OE, pvalue columns
```

!!! warning "p-values differ slightly from R GenomicDistributions"
    `calcExpectedPartitions` uses a goodness-of-fit `(O-E)²/E` formula. R `chisq.test()` computes a 2×2 test of independence with optional Yates correction, so p-values will not match R GenomicDistributions output byte-for-byte. The `log10OE` column is directly comparable.

## Signal matrix overlap

```r
# signalMatrix is a data.frame / data.table where column 1 is the region ID
# in "chr_start_end" format and the remaining columns are condition names
# with numeric signal values.
signalMatrix <- data.table::fread("signal_matrix.tsv")

result <- calcSummarySignal(query, signalMatrix)
# result is a list compatible with GenomicDistributions' plotSummarySignal:
#   $signalSummaryMatrix — data.table with queryPeak + one column per condition
#   $matrixStats — data.frame with 5 rows (Tukey stats) and one column per condition

# Row names of matrixStats are the Tukey boxplot quantities:
#   lowerWhisker, lowerHinge, median, upperHinge, upperWhisker
print(result$matrixStats)
```

Plug directly into R GenomicDistributions' `plotSummarySignal()` without further conversion.

## Distances to TSS / nearest features

```r
tss <- "gencode.v47.tss.bed"

# Signed distance (positive = feature downstream, negative = upstream)
d_signed <- calcFeatureDist(query, tss)

# Absolute distance
d_abs <- calcTSSDist(query, tss)
```

Both accept RegionSet, GRanges, file path, or data.frame for either argument.

## End-to-end example

```r
library(gtars)

# 1. Load inputs once
peaks <- RegionSet("peaks.bed")
genome <- loadGenomeAssembly("hg38.fa")
chromSizes <- c(chr1 = 248956422L, chr2 = 242193529L)
pl <- partitionListFromGTF("gencode.v47.gtf.gz", corePromSize = 100L,
                            proxPromSize = 2000L, chromSizes = chromSizes)

# 2. Widths and GC content
cat(sprintf("%d peaks, median width %.0f\n",
            length(peaks), median(calcWidth(peaks))))

gc <- calcGCContent(peaks, genome, ignoreUnkChroms = TRUE)
cat(sprintf("Mean GC content: %.3f\n", mean(gc)))

# 3. Partition enrichment
expected <- calcExpectedPartitions(peaks, pl, chromSizes, bpProportion = FALSE)
print(expected)

# 4. Distance to nearest TSS
tssDist <- calcTSSDist(peaks, "gencode.v47.tss.bed")
cat(sprintf("Median distance to nearest TSS: %d bp\n", median(tssDist)))

# 5. Convert to GRanges for downstream Bioconductor work
#    (only step that needs GenomicRanges — it's a bridge, not a core dep)
library(GenomicRanges)
gr <- as_granges(peaks)
```

## See also

- **[RegionSet & RegionSetList (R)](regionset.md)** — S4 class reference and method list.
- **[R LOLA interface](lola.md)** — enrichment testing built on top of these types.
- **[gtars-genomicdist](../genomicdist.md)** — Rust reference with all the algorithmic detail.
- **[R GenomicDistributions](https://code.databio.org/GenomicDistributions/)** — the original package gtars is a port of. Plotting helpers from GenomicDistributions still work on gtars return values.

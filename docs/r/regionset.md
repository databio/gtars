# RegionSet & RegionSetList (R)

The `gtars` R package exposes S4 classes `RegionSet` and `RegionSetList` — the R-idiomatic equivalents of Bioconductor's `GRanges` and `GRangesList`, backed by a high-performance Rust implementation via `externalptr`.

Unlike most gtars bindings, the R API is **polymorphic on input type**: nearly every function that takes a query will accept a `RegionSet`, a `GRanges` object, a file path, or a `data.frame` interchangeably. This lets you mix existing Bioconductor workflows with gtars without converting explicitly.

For the underlying semantics see [gtars-core](../core.md) (core types) and [gtars-genomicdist](../genomicdist.md) (statistics, set algebra).

## Loading the package

```r
library(gtars)
```

## `RegionSet`

An S4 class representing a set of genomic regions. Internally wraps a Rust `RegionSet` pointer plus an R-side strand vector.

### Construction

`RegionSet()` accepts any of these inputs:

```r
# 1. File path (BED, bed.gz, narrowPeak, etc.)
rs <- RegionSet("peaks.bed")

# 2. GRanges object (requires GenomicRanges)
library(GenomicRanges)
gr <- GRanges(seqnames = c("chr1", "chr1"),
              ranges = IRanges(start = c(100, 500), end = c(200, 600)),
              strand = c("+", "-"))
rs <- RegionSet(gr)   # coords converted from 1-based closed to 0-based half-open

# 3. data.frame with chr/start/end (and optional strand)
df <- data.frame(chr = c("chr1", "chr1"),
                 start = c(100, 500),
                 end = c(200, 600),
                 strand = c("+", "-"))
rs <- RegionSet(df)

# 4. An existing RegionSet (returned as-is)
rs2 <- RegionSet(rs)
```

`as_regionset()` is an alias for `RegionSet()`.

!!! warning "0-based half-open, BED convention"
    R gtars uses **0-based half-open** coordinates internally, matching BED convention (not GRanges' 1-based closed convention). Constructor/accessor functions convert automatically when bridging to `GRanges`; if you're building a RegionSet from a data.frame directly, the `start` and `end` columns are expected in 0-based half-open form.

### Basic S4 methods

```r
length(rs)              # number of regions
show(rs)                # pretty-print summary (first 5 regions)
rs                      # same — calls show()

as.data.frame(rs)       # convert to a data.frame with chr/start/end/strand
rs[1:10]                # subset by index (numeric, logical, or character)
rs[as.data.frame(rs)$strand == "+"]   # filter by strand
```

### Converting to `GRanges`

```r
gr <- as_granges(rs)
# Coordinates converted from 0-based half-open to 1-based closed.
# Strand information is preserved.
```

Requires the `GenomicRanges` package to be installed.

### Statistics methods

All of these accept `RegionSet`, `GRanges`, file path, or data.frame as input:

```r
widths(rs)                     # numeric vector of region widths (end - start)
neighborDistances(rs)          # signed gaps between consecutive regions per chromosome
nearestNeighbors(rs)           # per-region min neighbor distance
chromosomeStatistics(rs)       # data.frame, one row per chromosome

distribution(rs, nBins = 250)                     # bin counts across the genome
distribution(rs, nBins = 250, chromSizes = hg38)  # with reference chrom sizes

clusterRegions(rs, maxGap = 5000L)  # cluster id per region
```

!!! warning "Output length for `neighborDistances` / `nearestNeighbors`"
    Both skip chromosomes with only one region (matching R GenomicDistributions). The returned vector is **not aligned 1:1** with the input — it's shorter than `length(rs)` whenever any chromosome has a single peak.

!!! warning "`nBins` is a target, not a total"
    In `distribution(chromSizes = ...)`, `nBins` is the target bin count for the **longest** chromosome in `chromSizes`. Bin width is derived as `max(chromSizes) %/% nBins` (floored, minimum 1 bp), and every chromosome is tiled at the same bp width, so shorter chromosomes get proportionally fewer bins. No chromosome gets more than `nBins` bins: a midpoint in the leftover tail is folded into the last bin, which runs to the chromosome end. The genome-wide bin count can still substantially exceed `nBins` when `chromSizes` has many entries. Only bins that contain at least one region midpoint are returned, and regions on chromosomes missing from `chromSizes` (or whose midpoint lies past the chromosome end) are skipped. To target a specific bin width in bp instead, pass `nBins = max_chrom_len %/% desired_bp`.

### Interval set algebra

R gtars overrides the standard Bioconductor generics (`union`, `intersect`, `setdiff`, `reduce`, `promoters`, `shift`, `flank`, `resize`, `narrow`, `disjoin`, `gaps`, `findOverlaps`, `countOverlaps`) plus adds gtars-specific ones (`trim`, `pintersect`, `concat`, `jaccard`).

The unary methods (`reduce`, `trim`, `promoters`, `shift`, `flank`, `resize`, `narrow`, `disjoin`, `gaps`) dispatch to gtars for a `RegionSet`, file path, or data.frame. A `GRanges` input keeps its own IRanges method and returns a `GRanges`; wrap it in `RegionSet()` first to use gtars.

```r
merged   <- reduce(rs)
trimmed  <- trim(rs, chromSizes = hg38)
proms    <- promoters(rs, upstream = 2000L, downstream = 200L)
shifted  <- shift(rs, shift = 100L)
resized  <- resize(rs, width = 500L, fix = "center")
disjoint <- disjoin(rs)
gapped   <- gaps(rs, chrom_sizes = hg38)

# Binary operations
u  <- union(rs1, rs2)
i  <- intersect(rs1, rs2)
d  <- setdiff(rs1, rs2)
pi <- pintersect(rs1, rs2)   # pairwise (by index)
c  <- concat(rs1, rs2)

j  <- jaccard(rs1, rs2)      # scalar Jaccard similarity

# Overlap queries
hits  <- findOverlaps(rs1, rs2)
counts <- countOverlaps(rs1, rs2)
```

Methods dispatch on `("RegionSet", "RegionSet")`, `("RegionSet", "ANY")`, and `("ANY", "RegionSet")` so you can mix in `GRanges`, file paths, or data.frames on either side:

```r
# Works — second argument is a file path
j <- jaccard(rs, "other_peaks.bed")

# Works — first argument is a GRanges
merged <- union(gr, rs)
```

### Consensus regions

```r
cons <- consensus(list(rep1_rs, rep2_rs, rep3_rs))
# Returns a data.frame with chr, start, end, count columns.
# 'count' is the number of input sets overlapping each union region.

# Keep regions present in ≥ 2/3 replicates
robust <- cons[cons$count >= 2, ]
```

## `RegionSetList`

An S4 class for collections of `RegionSet`s — the gtars equivalent of `GRangesList`. Provides the single most efficient path for operating on many region sets at once, since pointers are passed between R and Rust without copying.

### Construction

```r
# Variadic: any mix of RegionSet / file path / GRanges / data.frame
rsl <- RegionSetList(
  RegionSet("rep1.bed"),
  RegionSet("rep2.bed"),
  RegionSet("rep3.bed"),
)

# Equivalent: single list argument
rsl <- RegionSetList(list(rs1, rs2, rs3))

# Or directly from file paths — each is auto-wrapped via RegionSet()
rsl <- RegionSetList("rep1.bed", "rep2.bed", "rep3.bed")

# Empty list
empty <- RegionSetList()
```

### S4 methods

```r
length(rsl)        # number of region sets
show(rsl)          # pretty-print with sizes
rsl[[1]]           # extract a single RegionSet by index (1-based)
names(rsl)         # character vector of names, or NULL
```

### Operations

```r
# Flatten into a single RegionSet (no merge/dedup)
flat <- concat(rsl)
# Apply reduce() on the result if you want the union
merged_all <- reduce(concat(rsl))

# Full N x N Jaccard similarity matrix
jac <- pairwise_jaccard(rsl)
# Returns a symmetric numeric matrix with 1.0 on the diagonal.
```

## Coordinate-system notes

| system | representation | gtars handling |
|---|---|---|
| BED / gtars | 0-based half-open | **internal** — always used by `RegionSet` |
| GRanges / IRanges | 1-based closed | converted in both directions by `RegionSet(gr)` / `as_granges(rs)` |
| data.frame input | 0-based half-open | passed through as-is — make sure your `start`/`end` columns are BED-style |

If you're working entirely in GRanges-land, use `as_granges()` to convert back before handing off to Bioconductor workflows. If you're reading BED files and writing results back to BED, stay in `RegionSet` to avoid double conversion.

## See also

- **[R GenomicDistributions wrappers](genomicdist.md)** — `calcWidth`, `calcGCContent`, `calcPartitions`, etc. — the drop-in API for users migrating from R GenomicDistributions.
- **[R LOLA interface](lola.md)** — `loadRegionDB`, `runLOLA`, `checkUniverseAppropriateness`, etc.
- **[R IGD interface](igd.md)** — `igd_create` / `igd_search` for low-level IGD access.
- **[gtars-core](../core.md)** — the underlying Rust types.
- **[Core models tour](../regionSet.md)** — Python + Rust walkthrough of the same concepts.

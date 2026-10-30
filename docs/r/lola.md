# LOLA enrichment in R (gtars)

The `gtars` R package ships a **drop-in replacement** for the [LOLA](https://www.bioconductor.org/packages/release/bioc/html/LOLA.html) R/Bioconductor package, powered by the Rust [`gtars-lola`](../lola.md) crate. Function names and signatures match LOLA closely enough that existing scripts port over with one-line changes.

For the underlying statistical reference see [gtars-lola](../lola.md) — Fisher's exact test on a hypergeometric survival function, CMLE odds ratio via Brent's method, per-user-set Benjamini-Hochberg FDR correction, and 1-based TSV output compatible with R LOLA's `writeCombinedEnrichment`.

## Loading the database

`loadRegionDB()` reads a standard LOLA-format folder (one directory per collection, each with `collection.txt`, `index.txt`, and a `regions/` subdirectory of BED files):

```r
library(gtars)

# Load the full database
regionDB <- loadRegionDB("LOLACore/hg38")

# Load only specific collections
regionDB <- loadRegionDB("LOLACore/hg38", collections = c("encode_tfbs", "roadmap_epigenomics"))

# Per-collection file cap (useful for smoke tests)
regionDB <- loadRegionDB("LOLACore/hg38", limit = 10L)
```

### From a list of BED files

If you don't have the full LOLA folder layout, `loadRegionDBFromBeds()` builds a minimal region database from an arbitrary list of BED paths:

```r
regionDB <- loadRegionDBFromBeds(
  bedFiles  = c("file1.bed", "file2.bed", "file3.bed"),
  filenames = c("Condition A", "Condition B", "Condition C"),  # optional display names
)
```

### Accessing annotations

```r
# Per-file annotation table (matches R LOLA's regionDB$regionAnno)
regionAnno <- regionDBAnno(regionDB)
# data.table with columns: filename, cellType, description, tissue,
# dataSource, antibody, treatment, collection

# Collection-level annotation
collectionAnno <- regionDBCollectionAnno(regionDB)
# data.table with columns: collectionname, collector, date, source, description

# List region set filenames
listRegionSets(regionDB)
listRegionSets(regionDB, collections = "encode_tfbs")  # filtered
```

### Extracting region sets

```r
# Get all region sets as a RegionSetList
all_sets <- getRegionSets(regionDB)

# Or a subset by 1-based indices
subset <- getRegionSets(regionDB, indices = c(1, 5, 12))

# RegionSetList is the same S4 class described in the RegionSet page
length(subset)
names(subset)
rs1 <- subset[[1]]  # individual RegionSet
```

## Running enrichment

```r
runLOLA(
  userSets,           # RegionSet, or list of RegionSets (or GRanges / file paths — polymorphic)
  userUniverse,       # RegionSet or GRanges representing the background
  regionDB,           # from loadRegionDB
  minOverlap = 1,     # minimum bp overlap to count as overlapping
  cores = 1,          # reserved for future parallelism — currently ignored
  redefineUserSets = FALSE,  # automatically rewrite user sets against universe
  direction = "enrichment",  # "enrichment" or "depletion"
)
```

**Returns a `data.table`** with one row per `(user_set, db_set)` pair, matching the R LOLA schema: `userSet`, `dbSet`, `collection`, `pValueLog`, `oddsRatio`, `support`, `rnkPV`, `rnkOR`, `rnkSup`, `maxRnk`, `meanRnk`, `b`, `c`, `d`, `description`, `cellType`, `tissue`, `antibody`, `treatment`, `dataSource`, `filename`, `qValue`, `size`.

- `userSet` and `dbSet` are **1-based** in the output data.table (R convention), even though they are 0-based internally.
- `pValueLog` is `-log10(p)` from Fisher's exact test, capped at ~322.
- `oddsRatio` is the CMLE odds ratio — matches `fisher.test()$estimate`, not the simple `(a·d)/(b·c)` point estimate.
- `qValue` is Benjamini-Hochberg adjusted, computed **per user set independently**.
- Rows are sorted by descending `pValueLog`, then ascending `meanRnk` — identical to R LOLA output order.

### Basic example

```r
library(gtars)

regionDB <- loadRegionDB("LOLACore/hg38")

# Polymorphic inputs — RegionSet, GRanges, file path all accepted
userSets <- list(
  RegionSet("condition_a.bed"),
  RegionSet("condition_b.bed"),
)
universe <- "universe.bed"  # file path works directly

results <- runLOLA(userSets, universe, regionDB,
                    minOverlap = 1,
                    direction = "enrichment")

# Top 20 hits for user set 1
head(results[userSet == 1], 20)
```

### Two-condition comparison

```r
# Differential enrichment — use the union of both conditions as the universe
restricted <- buildRestrictedUniverse(list(peaks_a, peaks_b))

results <- runLOLA(list(peaks_a, peaks_b), restricted, regionDB)

# Top hits unique to condition A
head(results[userSet == 1 & qValue < 0.05][order(-pValueLog)], 20)
```

## Universe diagnostics

LOLA results are only meaningful if your universe **contains** your user sets: every user region should overlap at least one universe region.

### `checkUniverseAppropriateness`

```r
report <- checkUniverseAppropriateness(userSets, userUniverse)
# data.frame with columns:
#   userSet (1-based), totalRegions, regionsInUniverse, coverage, manyToMany
# Warnings are emitted via warning() for low coverage or many-to-many mappings.
```

Warnings fire when coverage is below 50% (severe) or 90% (moderate), or when any user region overlaps more than one universe region.

### `redefineUserSets`

Rewrite each user set in terms of universe regions — eliminates many-to-many mapping artifacts.

```r
redefined <- redefineUserSets(userSets, userUniverse)
# Returns a list of RegionSet objects.

# Use directly
results <- runLOLA(redefined, userUniverse, regionDB)

# Or pass redefineUserSets = TRUE to runLOLA to do this automatically
results <- runLOLA(userSets, userUniverse, regionDB, redefineUserSets = TRUE)
```

### `buildRestrictedUniverse`

For differential enrichment: build a universe that is exactly the union of all user sets, disjoined into non-overlapping pieces (R LOLA's `disjoin(unlist(userSets))`):

```r
restricted <- buildRestrictedUniverse(list(peaks_a, peaks_b, peaks_c))
# Returns a RegionSet
```

Unlike the other functions on this page, `buildRestrictedUniverse` only accepts `RegionSet` objects (not file paths or `GRanges`), so wrap inputs with `RegionSet()` first.

## Porting from R LOLA

Most LOLA scripts port by changing one line — replacing `library(LOLA)` with `library(gtars)`. The core API is compatible:

| R LOLA | gtars R | notes |
|---|---|---|
| `loadRegionDB()` | `loadRegionDB()` | same signature |
| `runLOLA()` | `runLOLA()` | same signature; returns `data.table` |
| `checkUniverseAppropriateness()` | `checkUniverseAppropriateness()` | same signature |
| `redefineUserSets()` | `redefineUserSets()` | same signature; returns list of `RegionSet` |
| `writeCombinedEnrichment()` | *(use `data.table::fwrite`)* | output table is already in the right format |
| `extractEnrichmentOverlaps()` | *(not implemented)* | file an issue if needed |
| `getRegionFile()` | `getRegionSets(regionDB, index)` | returns a `RegionSetList`, not `GRanges`; take `[[1]]` and call `as_granges()` to convert |

The notable differences:

- **Return type is `data.table`, not `data.frame`.** R LOLA returns a plain data.frame; gtars returns a `data.table` by default. Convert with `as.data.frame()` if you prefer.
- **Faster p-value and odds ratio.** The p-value uses the survival function (no cancellation at small tails), and the odds ratio is the CMLE (not the point estimate). Expect tail significance to be reported more precisely.

## End-to-end example

```r
library(gtars)
library(data.table)

# 1. Load the database once
regionDB <- loadRegionDB("LOLACore/hg38")
cat(sprintf("Loaded %d region sets\n", length(listRegionSets(regionDB))))

# 2. Load user sets and universe — polymorphic inputs
peaks_a <- "peaks_condition_a.bed"
peaks_b <- "peaks_condition_b.bed"
universe <- "universe.bed"

# 3. Universe sanity check
diag <- checkUniverseAppropriateness(list(peaks_a, peaks_b), universe)
print(diag)

# 4. Run enrichment with auto-redefinition
results <- runLOLA(
  list(peaks_a, peaks_b), universe, regionDB,
  minOverlap = 1,
  redefineUserSets = TRUE,
  direction = "enrichment",
)

# 5. Top hits per user set
for (us in unique(results$userSet)) {
  cat(sprintf("\n=== User set %d — top 10 ===\n", us))
  top <- head(
    results[userSet == us][order(-pValueLog)],
    10,
  )
  print(top[, .(filename, cellType, pValueLog, oddsRatio, support, qValue)])
}

# 6. Write to TSV matching R LOLA's writeCombinedEnrichment format
fwrite(results, "lola_results.tsv", sep = "\t")
```

## See also

- **[RegionSet & RegionSetList (R)](regionset.md)** — the S4 classes used for user sets and extracted database sets.
- **[R GenomicDistributions wrappers](genomicdist.md)** — statistics and partition analysis.
- **[R IGD interface](igd.md)** — IGD is the overlap index that backs a `RegionDB`.
- **[gtars-lola](../lola.md)** — Rust reference with the full statistical detail.
- **[R LOLA package](https://www.bioconductor.org/packages/release/bioc/html/LOLA.html)** — the original Bioconductor package this is a port of.

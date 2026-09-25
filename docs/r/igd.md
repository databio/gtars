# IGD in R (gtars)

Low-level R bindings for IGD ([Indexed Genomic Data](../igd.md)) — a fast binary index over collections of BED files, used internally by gtars LOLA and by the BEDbase retrieval pipeline. The R API exposes two functions: `igd_create` (build an index from BED files on disk) and `igd_search` (query an existing index with a BED file).

If you're running LOLA enrichment, you usually don't need to touch IGD directly — `loadRegionDB()` in [`r/lola.md`](lola.md) builds an IGD index internally from the region database. Use this page if you want to index a custom collection of BED files and query it directly.

## Creating an IGD database

```r
library(gtars)

igd_create(
  output_path = "path/to/output/",
  filelist    = "path/to/bed/files/",
  db_name     = "igd_database",  # default
)
```

**Arguments:**

- `output_path` — directory where the IGD database files will be written. Created if it does not exist.
- `filelist` — one of:
    - a path to a `.txt` file listing BED file paths (one per line; the `.txt` extension is how it is recognized),
    - a path to a directory containing BED files (all `.bed` and `.gz` files in it are indexed),
    - `"-"` or `"stdin"` to read paths from standard input.
- `db_name` — prefix for the output filenames (default `"igd_database"`). The function produces two files: `<db_name>.igd` (the index) and `<db_name>.tsv` (file metadata).

Returns `NULL` invisibly on success. Errors are raised via `stop()` on invalid input.

### Example

```r
# From a directory of BED files
igd_create(
  output_path = "./igd_out",
  filelist    = "./my_bed_files/",
  db_name     = "my_peaks",
)

# Produces:
# ./igd_out/my_peaks.igd
# ./igd_out/my_peaks.tsv

# From an explicit file list
writeLines(
  c("./rep1.bed", "./rep2.bed", "./rep3.bed"),
  "bed_list.txt",
)
igd_create(
  output_path = "./igd_out",
  filelist    = "bed_list.txt",
  db_name     = "replicates",
)
```

## Searching an IGD database

```r
hits <- igd_search(
  database_path = "path/to/my_peaks.igd",
  query_path    = "path/to/query.bed",
)
```

**Arguments:**

- `database_path` — path to an existing `.igd` file (produced by `igd_create` or any other IGD-compatible tool).
- `query_path` — path to a BED file containing the query regions.

**Returns** a `data.frame` with one row per database file that has at least one hit, and columns `filename`, `numRegions` (regions in that file), and `hits` (number of overlaps with the query).

### Example

```r
hits <- igd_search(
  database_path = "./igd_out/my_peaks.igd",
  query_path    = "query.bed",
)

head(hits)

# Typical analysis: rank database files by overlap count
hits[order(-hits$hits), ]
```

## Relationship to other gtars R functions

- **`loadRegionDB()` in `lola.R`** — builds an IGD index internally for a LOLA region database. The IGD is wrapped inside a `RegionDB` pointer and used by `runLOLA()`. Users running enrichment typically don't need to call `igd_create` / `igd_search` separately.
- **`countOverlaps()` / `findOverlaps()` on `RegionSet`** — ad-hoc overlap queries between two `RegionSet` objects use the AIList index from [`gtars-overlaprs`](../overlaprs.md), not IGD. Use those methods when you have both the query and subject in memory.
- **IGD is preferred when** you have a large static collection of BED files that you'll query many times against different user sets, and you want the index persisted on disk.

!!! note "Naming convention"
    Unlike the rest of the R gtars API, which uses camelCase (`runLOLA`, `calcPartitions`, etc.), the IGD functions use snake_case (`igd_create`, `igd_search`). This is a legacy from the direct CLI mapping and may be aligned in a future release.

## See also

- **[gtars-igd](../igd.md)** — the underlying Rust crate, with the binary format reference and performance characteristics.
- **[R LOLA interface](lola.md)** — `loadRegionDB()` and `runLOLA()`, which use IGD internally.
- **[gtars CLI](../cli.md)** — `gtars igd create` and `gtars igd search` commands, equivalent to these R functions.

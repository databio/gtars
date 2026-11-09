# BLAST over a RefgetStore in R

The gtars R package builds [BLAST](https://blast.ncbi.nlm.nih.gov/) databases from a [RefgetStore](refgetstore.md), searches a store, and tells you which genomes a sequence could come from. The `blast_*` functions have the same names and options as [`gtars blast`](../blast.md) on the command line (long flag names, with `-` written as `_`). Tables come back as data.frames and the rest as named lists.

## Make a store

```r
library(gtars)

set.seed(7)
chroms <- vapply(1:3, function(i) paste(sample(c("A", "C", "G", "T"), 2000, TRUE), collapse = ""), "")

work <- file.path(tempdir(), "blast_demo"); dir.create(work)
store <- refget_store_on_disk(file.path(work, "store"))
for (g in list(a = chroms[1:2], b = chroms[2:3])) {
  fa <- tempfile(tmpdir = work, fileext = ".fa")
  writeLines(as.vector(rbind(paste0(">s", seq_along(g)), g)), fa)
  add_fasta(store, fa)
}
writeLines(c(">q1", substr(chroms[1], 501, 900)), file.path(work, "query.fa"))
```

## Build a database

```r
s <- blast_build(file.path(work, "store"), file.path(work, "db"), all = TRUE)
s$preflight$subjects
blast_verify(file.path(work, "db"), file.path(work, "store"))$ok
```

Choose collections with `all = TRUE`, or `collection = "ucsc:hg38"`, `namespace`, `taxon` and the other selection arguments (a value or a vector). `dry_run = TRUE` returns what a build would do without writing.

### Progress and Ctrl-C

By default `progress = TRUE` prints each step with `message()`, so `suppressMessages()` silences it. Pass a function to get each event as a named list (`event`, its fields, and `message`):

```r
blast_build(file.path(work, "store"), file.path(work, "db2"), all = TRUE,
            progress = function(ev) if (ev$event == "sequences") cat(ev$done, "of", ev$total, "\n"))
```

Ctrl-C cancels the build, removes partial output, and signals the condition `gtars_rgblast_cancelled`.

## Search

`blast_search()` runs rgblastn, a Rust port of blastn, in process on your selection of the store. It needs no database and no BLAST+ install, and covers tasks `megablast`, `blastn` and `blastn-short` with `-outfmt 6`.

```r
r <- blast_search(file.path(work, "store"), query = file.path(work, "query.fa"), all = TRUE,
                  outfmt = "6")
cat(r$output)
```

With `out = NULL` the hits come back in `r$output`. Give `out = "hits.tsv"` to write a file instead (the run metadata goes next to it in `hits.tsv.meta.json`). More blastn options go in `args = c("-evalue", "1e-20")`.

For any other blastn option, run BLAST+ on a database you built: `blast_search(query = ..., engine = "blastn", db = file.path(work, "db"))`. blastn comes from `blastn =`, else `$GTARS_RGBLAST_BLASTN`, else the `PATH`.

## Errors are classed conditions

Every error inherits from `gtars_rgblast_error` and from one of `gtars_rgblast_not_found`, `_invalid_input`, `_conflict`, `_io`, `_network`, `_data` or `_cancelled`. A missing BLAST+ also has `gtars_rgblast_engine_not_found`, and a blastn run that fails has `gtars_rgblast_engine_failed` (with `returncode`). Catch them by class:

```r
tryCatch(
  blast_search(query = file.path(work, "query.fa"), engine = "blastn",
               db = file.path(work, "db"), blastn = "/nonexistent/blastn"),
  gtars_rgblast_engine_not_found = function(e) message("no BLAST+: ", conditionMessage(e))
)
```

## Where did this sequence come from?

`blast_origin()` searches one representative assembly per species, then the rest of each species that matched, and returns the best hit per genome:

```r
r <- blast_origin(file.path(work, "store"), query = file.path(work, "query.fa"), all = TRUE)
r$summary          # one row per query: status, taxon, top_genome, pident, qcov, ...
r$report           # the full report as a nested list
r$meta             # how the run went
```

To start from hits you already have, use `blast_report("hits.tsv")`, which reads `hits.tsv.meta.json` from `blast_search()`, or `blast_report("hits.tsv", db = "db/hg38")` for hits from a database. `blast_reps(store, all = TRUE)$groups` shows the representative of each species.

## Views and N counts

`blast_view()` writes one searchable view per assembly, and `blast_nstats()` counts N bases per collection (a data.frame).

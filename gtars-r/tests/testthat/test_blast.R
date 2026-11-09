library(testthat)
library(gtars)

fasta_dir <- file.path(test_path(), '..', '..', '..', 'tests', 'data', 'fasta')
base_fa <- file.path(fasta_dir, 'base.fa')
names_fa <- file.path(fasta_dir, 'different_names.fa')

skip_if_no_fasta <- function() {
  skip_if_not(file.exists(base_fa) && file.exists(names_fa), "Test FASTA files not found")
}

# A store with base.fa, different_names.fa (same sequences, other names) and a collection of
# its own. Returns the store path and the three collection digests.
make_store <- function() {
  dir <- tempfile("blast_test_")
  dir.create(dir)
  uniq <- file.path(dir, "uniq.fa")
  writeLines(c(">u1", "ACGTACGTACGTTTGACCA", ">u2", "CCCGGGATATATTT"), uniq)
  path <- file.path(dir, "store")
  store <- refget_store_on_disk(path)
  digests <- c()
  for (src in c(base_fa, names_fa, uniq)) {
    # copies: importing writes an .rgsi next to the FASTA
    fa <- file.path(dir, basename(src))
    if (!file.exists(fa)) file.copy(src, fa)
    add_fasta(store, fa)
    digests <- c(digests, digest_fasta(fa)@digest)
  }
  names(digests) <- c("base", "other", "uniq")
  list(path = path, dir = dir, digests = digests)
}

with_env <- function(code) {
  old <- Sys.getenv(c("GTARS_RGBLAST_CACHE", "SOURCE_DATE_EPOCH"), unset = NA)
  Sys.setenv(GTARS_RGBLAST_CACHE = tempfile("nstats_cache_"), SOURCE_DATE_EPOCH = "0")
  on.exit({
    for (k in names(old)) if (is.na(old[[k]])) Sys.unsetenv(k) else do.call(Sys.setenv, setNames(list(old[[k]]), k))
  })
  force(code)
}

test_that("blast_build writes a database and verify passes", {
  skip_if_no_fasta()
  with_env({
    s <- make_store()
    out <- file.path(s$dir, "db")
    r <- blast_build(s$path, out, all = TRUE, threads = 2L, progress = FALSE)
    expect_type(r, "list")
    expect_gt(r$preflight$subjects, 0)
    expect_true(all(file.exists(unlist(r$files))))
    v <- blast_verify(out, s$path, threads = 2L, progress = FALSE)
    expect_true(v$ok)
    expect_equal(v$ids_failed, 0)
  })
})

test_that("blast_build dry_run writes nothing", {
  skip_if_no_fasta()
  with_env({
    s <- make_store()
    out <- file.path(s$dir, "db")
    r <- blast_build(s$path, out, all = TRUE, dry_run = TRUE, progress = FALSE)
    expect_gt(r$subjects, 0)
    expect_length(Sys.glob(paste0(out, "*")), 0L)
  })
})

test_that("a progress function receives opened first and finished last", {
  skip_if_no_fasta()
  with_env({
    s <- make_store()
    events <- character()
    blast_build(s$path, file.path(s$dir, "db"), all = TRUE, threads = 2L,
                progress = function(ev) {
                  expect_true(is.character(ev$message))
                  events <<- c(events, ev$event)
                })
    expect_equal(events[1], "opened")
    expect_equal(events[length(events)], "finished")
    expect_true("planned" %in% events)
  })
})

test_that("progress = TRUE prints messages", {
  skip_if_no_fasta()
  with_env({
    s <- make_store()
    expect_message(blast_build(s$path, file.path(s$dir, "db"), all = TRUE, threads = 2L))
    expect_silent(suppressMessages(
      blast_build(s$path, file.path(s$dir, "db2"), all = TRUE, threads = 2L)))
  })
})

test_that("suppressMessages() silences progress = TRUE; withCallingHandlers sees it", {
  skip_if_no_fasta()
  with_env({
    s <- make_store()
    seen <- 0L
    withCallingHandlers(
      blast_build(s$path, file.path(s$dir, "db"), all = TRUE, threads = 2L),
      message = function(m) { seen <<- seen + 1L; invokeRestart("muffleMessage") })
    expect_gt(seen, 0L)
  })
})

test_that("an error in the progress function cancels and removes partial output", {
  skip_if_no_fasta()
  with_env({
    s <- make_store()
    out <- file.path(s$dir, "db")
    expect_error(
      blast_build(s$path, out, all = TRUE,
                  progress = function(ev) if (ev$event == "planned") stop("stop here")),
      "stop here")
    expect_length(Sys.glob(paste0(out, "*")), 0L)
  })
})

test_that("an interrupt cancels the build and removes partial output", {
  skip_if_no_fasta()
  with_env({
    s <- make_store()
    out <- file.path(s$dir, "db")
    cond <- tryCatch(
      blast_build(s$path, out, all = TRUE,
                  # what R signals on Ctrl-C (a real SIGINT would stop testthat itself)
                  progress = function(ev) if (ev$event == "planned")
                    signalCondition(structure(class = c("interrupt", "condition"), list()))),
      gtars_rgblast_cancelled = function(e) e)
    expect_s3_class(cond, "gtars_rgblast_cancelled")
    expect_length(Sys.glob(paste0(out, "*")), 0L)
  })
})

test_that("errors are classed conditions", {
  skip_if_no_fasta()
  with_env({
    s <- make_store()
    out <- file.path(s$dir, "db")
    expect_error(blast_build(s$path, out, all = TRUE, collection = "x", progress = FALSE),
                 class = "gtars_rgblast_invalid_input")
    caught <- tryCatch(
      blast_build(s$path, out, all = TRUE, collection = "x", progress = FALSE),
      gtars_rgblast_invalid_input = function(e) conditionMessage(e))
    expect_match(caught, "cannot be used with")
    expect_error(blast_build(s$path, out, collection = "no-such-collection", progress = FALSE),
                 class = "gtars_rgblast_not_found")
    blast_build(s$path, out, all = TRUE, threads = 2L, progress = FALSE)
    expect_error(blast_build(s$path, out, all = TRUE, threads = 2L, progress = FALSE),
                 class = "gtars_rgblast_conflict")
    expect_error(blast_build(s$path, out, all = TRUE, force = TRUE, threads = 2L, progress = FALSE),
                 NA)
    expect_error(blast_build(s$path, out, all = TRUE, fetch_jobs = 4, progress = FALSE),
                 "only applies to remote stores")
    expect_error(blast_build(s$path, out, all = TRUE, progress = 3),
                 class = "gtars_rgblast_error")
  })
})

test_that("blast_view writes views; no mode is invalid input", {
  skip_if_no_fasta()
  with_env({
    s <- make_store()
    out <- file.path(s$dir, "db")
    blast_build(s$path, out, all = TRUE, threads = 2L, progress = FALSE)
    v <- blast_view(out, per_collection = TRUE, progress = FALSE)
    expect_s3_class(v$views, "data.frame")
    expect_equal(nrow(v$views), 3L)
    expect_error(blast_view(out, progress = FALSE), class = "gtars_rgblast_invalid_input")
  })
})

test_that("blast_nstats returns one row per collection", {
  skip_if_no_fasta()
  with_env({
    s <- make_store()
    df <- blast_nstats(s$path, all = TRUE, threads = 2L, progress = FALSE)
    expect_s3_class(df, "data.frame")
    expect_equal(nrow(df), 3L)
    expect_true(all(df$n_fraction >= 0 & df$n_fraction <= 1))
    expect_true(all(c("digest", "bases", "n", "n_fraction") %in% names(df)))
  })
})

# ---------------------------------------------------------------------------------------------
# search, report, origin, reps, cache

origin_dir <- file.path(test_path(), '..', '..', '..', 'gtars-rgblast', 'tests', 'data', 'fixtures', 'origin')

skip_if_no_origin <- function() {
  skip_if_not(file.exists(file.path(origin_dir, "small_jungle_hits.tsv")), "origin fixtures not found")
}

# The blastn for --blastn: PATH, GTARS_RGBLAST_BIN, or a test-only shim script that runs docker
# image ncbi/blast:2.17.0 with the temp and working directories mounted at their own paths.
blastn_path <- function() {
  if (nzchar(Sys.which("blastn"))) return(unname(Sys.which("blastn")))
  d <- Sys.getenv("GTARS_RGBLAST_BIN")
  if (nzchar(d) && file.exists(file.path(d, "blastn"))) return(file.path(d, "blastn"))
  image <- Sys.getenv("GTARS_RGBLAST_DOCKER", "ncbi/blast:2.17.0")
  if (nzchar(Sys.which("docker")) &&
      system2("docker", c("image", "inspect", image), stdout = FALSE, stderr = FALSE) == 0) {
    dirs <- unique(normalizePath(c(tempdir(), dirname(tempdir()), getwd())))
    ids <- system2("id", c("-u"), stdout = TRUE)
    gid <- system2("id", c("-g"), stdout = TRUE)
    shim <- tempfile("blastn_")
    writeLines(c("#!/bin/sh", paste0(
      "exec docker run --rm -i -u ", ids, ":", gid, " ",
      paste0("-v '", dirs, ":", dirs, "'", collapse = " "),
      " -w \"$PWD\" ", image, " blastn \"$@\"")), shim)
    Sys.chmod(shim, "755")
    return(shim)
  }
  NULL
}

# Two collections of random 2 kb sequences (they share one) and a query cut from the first.
make_dna_store <- function() {
  dir <- tempfile("blast_dna_")
  dir.create(dir)
  set.seed(7)
  seqs <- vapply(1:3, function(i) paste(sample(c("A", "C", "G", "T"), 2000, replace = TRUE), collapse = ""), "")
  store <- refget_store_on_disk(file.path(dir, "store"))
  for (i in 1:2) {
    recs <- if (i == 1) seqs[1:2] else seqs[2:3]
    fa <- file.path(dir, paste0("c", i, ".fa"))
    writeLines(as.vector(rbind(paste0(">s", i, "_", seq_along(recs)), recs)), fa)
    add_fasta(store, fa)
  }
  q <- file.path(dir, "q.fa")
  writeLines(c(">q1", substr(seqs[1], 501, 900)), q)
  list(path = file.path(dir, "store"), query = q, dir = dir)
}

test_that("blast_report on the fixture returns the TSV's columns and rows", {
  skip_if_no_origin()
  r <- blast_report(file.path(origin_dir, "small_jungle_hits.tsv"),
                    db = file.path(origin_dir, "small_jungle"),
                    query = file.path(origin_dir, "small_jungle_queries.fa"))
  want <- read.delim(file.path(origin_dir, "expected_summary.tsv"), colClasses = "character",
                     check.names = FALSE, comment.char = "")
  names(want) <- sub("^#", "", names(want))
  expect_s3_class(r$summary, "data.frame")
  expect_equal(names(r$summary), names(want))
  expect_equal(nrow(r$summary), nrow(want))
  expect_equal(as.character(r$summary$top_genome), want$top_genome)
  expect_equal(as.character(r$summary$taxid), want$taxid)
  expect_type(r$report, "list")
  out <- tempfile(fileext = ".json")
  blast_report(file.path(origin_dir, "small_jungle_hits.tsv"), db = file.path(origin_dir, "small_jungle"),
               output = out, format = "json")
  expect_match(readLines(out, 1L), "^\\{")
})

test_that("blast_report signals classed conditions", {
  skip_if_no_origin()
  hits <- file.path(origin_dir, "small_jungle_hits.tsv")
  expect_error(blast_report(hits), class = "gtars_rgblast_invalid_input")
  expect_error(blast_report(hits, db = file.path(origin_dir, "small_jungle"), format = "pdf"),
               class = "gtars_rgblast_invalid_input")
})

test_that("blast_report_fields lists the appended columns", {
  expect_equal(blast_report_fields(), c("qlen", "score", "sallseqid", "salltitles", "staxids"))
})

test_that("blast_reps on a store", {
  skip_if_no_fasta()
  with_env({
    s <- make_store()
    r <- blast_reps(s$path, all = TRUE)
    expect_s3_class(r$groups, "data.frame")
    expect_gt(nrow(r$groups), 0L)
    expect_error(blast_reps(s$path), class = "gtars_rgblast_invalid_input")
  })
})

test_that("blast_search without an engine signals engine_not_found", {
  skip_if_no_fasta()
  with_env({
    s <- make_store()
    err <- tryCatch(
      blast_search(query = base_fa, engine = "blastn", db = tempfile("db_"),
                   blastn = "/nonexistent/blastn"),
      gtars_rgblast_engine_not_found = function(e) e)
    expect_s3_class(err, "gtars_rgblast_engine_not_found")
    expect_s3_class(err, "gtars_rgblast_not_found")
    expect_error(blast_search(s$path, query = base_fa, all = TRUE, engine = "nope"),
                 class = "gtars_rgblast_invalid_input")
    expect_error(blast_search(query = base_fa, engine = "blastn"),
                 class = "gtars_rgblast_invalid_input")
  })
})

test_that("blast_search and blast_origin with blastn on a database", {
  bp <- blastn_path()
  skip_if(is.null(bp), "blastn not found")
  with_env({
    d <- make_dna_store()
    db <- file.path(tempfile("db_"), "db")
    dir.create(dirname(db))
    blast_build(d$path, db, all = TRUE, progress = FALSE)
    run <- function(...) blast_search(query = d$query, outfmt = "6", engine = "blastn", db = db,
                                      blastn = bp, progress = FALSE, ...)
    r1 <- run()
    expect_equal(r1$exit_code, 0)
    expect_true(nzchar(r1$output))
    expect_equal(normalizePath(r1$db), normalizePath(db))
    expect_error(run(args = c("-task", "nonsense")), class = "gtars_rgblast_engine_failed")
    o <- blast_origin(query = d$query, engine = "blastn", db = db, blastn = bp, progress = FALSE)
    expect_equal(as.character(o$summary$query), "q1")
    expect_true(all(c("summary", "report", "meta") %in% names(o)))
    expect_error(blast_origin(query = d$query, engine = "blastn", db = db, staged = TRUE),
                 class = "gtars_rgblast_invalid_input")
  })
})

test_that("blast_search and blast_origin run with the default engine, rgblastn", {
  with_env({
    d <- make_dna_store()
    r <- blast_search(d$path, query = d$query, all = TRUE, outfmt = "6", progress = FALSE)
    expect_equal(r$exit_code, 0)
    expect_match(r$output, "^q1\t")
    o <- blast_origin(d$path, query = d$query, all = TRUE, progress = FALSE)
    expect_equal(as.character(o$summary$query), "q1")
    expect_error(blast_search(d$path, query = d$query, all = TRUE, blastn = "x"),
                 class = "gtars_rgblast_invalid_input")
  })
})

# BLAST over a RefgetStore: blast_build(), blast_verify(), blast_view(), blast_nstats(),
# blast_search(), blast_report(), blast_origin() and blast_reps(). Each runs the `gtars blast`
# verb of the same name. Options passed in `...` are that verb's long flags
# with `-` written as `_` (`gtars blast build --help` lists them): TRUE passes a flag, NULL and
# FALSE leave it out, a vector repeats it, and any other value is passed as text.
#
# Errors are classed conditions: every one inherits from `gtars_rgblast_error` and from one of
# `gtars_rgblast_not_found`, `_invalid_input`, `_conflict`, `_io`, `_network`, `_data` or
# `_cancelled`. A missing engine also has `gtars_rgblast_engine_not_found`, and a failing search
# engine `gtars_rgblast_engine_failed`.

# Raise the condition a Rust wrapper described, or return its value.
.blast_call <- function(res) {
  if (isTRUE(res$ok)) return(res$value)
  stop(errorCondition(res$message, class = c(res$class, "gtars_rgblast_error"),
                      call = sys.call(-1)))
}

# Run one verb. The library call runs on a worker thread; this loop, on R's main thread, runs
# `progress` for each event, so `suppressMessages()` and `tryCatch()` see what it raises. An
# error in `progress` cancels the job (partial output is removed) and propagates unchanged;
# Ctrl-C cancels it and signals `gtars_rgblast_cancelled`.
.blast_run <- function(argv, progress) {
  fn <- if (is.function(progress)) {
    progress
  } else if (is.null(progress) || isFALSE(progress)) {
    NULL
  } else if (isTRUE(progress)) {
    function(ev) message(ev$message)
  } else {
    stop(errorCondition("progress must be TRUE, FALSE, NULL or a function",
                        class = c("gtars_rgblast_invalid_input", "gtars_rgblast_error")))
  }
  job <- .blast_call(blast_start_rs(argv))
  finished <- FALSE
  on.exit(if (!finished) blast_cancel_rs(job), add = TRUE)
  tryCatch(
    repeat {
      r <- blast_poll_rs(job, 100)
      if (identical(r$kind, "event")) {
        if (!is.null(fn)) fn(r$event)
      } else if (identical(r$kind, "done")) {
        finished <- TRUE
        return(.blast_call(r$result))
      }
    },
    interrupt = function(i) {
      blast_cancel_rs(job)
      finished <<- TRUE
      stop(errorCondition("interrupted", class = c("gtars_rgblast_cancelled", "gtars_rgblast_error")))
    }
  )
}

# One option value as text: paths have `~` expanded, numbers are written without exponents.
.blast_text <- function(x) {
  if (is.numeric(x)) return(format(x, scientific = FALSE, trim = TRUE, digits = 15))
  x <- as.character(x)
  ifelse(grepl("^https?://", x), x, path.expand(x))
}

# Named options as flags.
.blast_flags <- function(opts) {
  if (length(opts) && (is.null(names(opts)) || any(names(opts) == ""))) {
    stop(errorCondition("options must be named (for example all = TRUE)",
                        class = c("gtars_rgblast_invalid_input", "gtars_rgblast_error")))
  }
  argv <- character()
  for (key in names(opts)) {
    v <- opts[[key]]
    if (is.null(v) || isFALSE(v)) next
    flag <- paste0("--", gsub("_", "-", key))
    if (isTRUE(v)) {
      argv <- c(argv, flag)
    } else {
      argv <- c(argv, paste0(flag, "=", .blast_text(v)))
    }
  }
  argv
}

# Engine arguments (-query, -out, -outfmt, then the rest), after "--".
.blast_store <- function(store) {
  if (is.null(store)) character() else c("--store", .blast_text(store))
}

.blast_engine_args <- function(query = NULL, out = NULL, outfmt = NULL, args = character()) {
  a <- character()
  if (!is.null(query)) a <- c(a, "-query", .blast_text(query))
  if (!is.null(out)) a <- c(a, "-out", .blast_text(out))
  if (!is.null(outfmt)) a <- c(a, "-outfmt", as.character(outfmt))
  a <- c(a, as.character(args))
  if (length(a)) c("--", a) else character()
}

# A list of row lists (one per record) as a data.frame; zero rows give a zero-row data.frame.
.blast_df <- function(rows) {
  if (length(rows) == 0L) return(data.frame())
  cols <- names(rows[[1L]])
  out <- lapply(cols, function(k) {
    v <- lapply(rows, function(r) { x <- r[[k]]; if (is.null(x)) NA else x })
    if (all(vapply(v, length, 1L) == 1L)) unlist(v, use.names = FALSE) else I(v)
  })
  names(out) <- cols
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

# A report result: the summary rows as a data.frame next to the report document.
.blast_report_result <- function(res) {
  df <- .blast_df(res$summary)
  if (ncol(df)) df <- df[, blast_summary_columns_rs(), drop = FALSE]
  res$summary <- df
  res
}

#' BLAST databases from a RefgetStore
#'
#' \code{blast_build()} writes a makeblastdb-compatible version 5 nucleotide database straight
#' from a store; \code{blast_verify()} checks one against its store; \code{blast_view()} writes
#' subset views (\code{.nal} + \code{.bsl}); \code{blast_nstats()} counts N bases per
#' collection.
#'
#' @param store Local store directory or \code{http(s)} store URL.
#' @param out Database base path (\code{makeblastdb -out}).
#' @param db Database base path.
#' @param ... Options: the verb's long flags with \code{-} written as \code{_}, for example
#'   \code{all = TRUE}, \code{collection = c("a", "b")}, \code{threads = 4}. See
#'   \code{gtars blast build --help}.
#' @param progress \code{TRUE} (messages, silenced by \code{suppressMessages()}), \code{FALSE},
#'   or a function called with each event as a named list (\code{event}, its fields, and
#'   \code{message}). An error in the function, or Ctrl-C, cancels the work and removes
#'   partial output.
#' @return The verb's summary as a named list; \code{blast_nstats()} a data.frame, one row
#'   per collection.
#' @section Conditions:
#'   \code{gtars_rgblast_invalid_input} (bad or conflicting options),
#'   \code{gtars_rgblast_not_found} (unknown store or collection),
#'   \code{gtars_rgblast_conflict} (the database exists; use \code{force = TRUE}),
#'   \code{gtars_rgblast_io}, \code{gtars_rgblast_network}, \code{gtars_rgblast_cancelled}.
#' @examples
#' \dontrun{
#' s <- blast_build("store", "db/hg38", collection = "ucsc:hg38")
#' blast_verify("db/hg38", "store")$ok
#' }
#' @export
blast_build <- function(store, out, ..., progress = TRUE) {
  .blast_run(c("build", "--store", .blast_text(store), "--out", .blast_text(out),
               .blast_flags(list(...))), progress)
}

#' @rdname blast_build
#' @export
blast_verify <- function(db, store, ..., progress = FALSE) {
  .blast_run(c("verify", "--db", .blast_text(db), "--store", .blast_text(store),
               .blast_flags(list(...))), progress)
}

#' @rdname blast_build
#' @export
blast_view <- function(db, ..., progress = FALSE) {
  res <- .blast_run(c("view", "--db", .blast_text(db), .blast_flags(list(...))), progress)
  res$views <- .blast_df(res$views)
  res
}

#' @rdname blast_build
#' @export
blast_nstats <- function(store, ..., progress = FALSE) {
  res <- .blast_run(c("nstats", "--store", .blast_text(store), .blast_flags(list(...))), progress)
  df <- .blast_df(res$collections)
  attr(df, "sequences") <- res$sequences
  attr(df, "excluded") <- res$excluded
  df
}

#' The columns \code{blast_search()} adds to the search output
#'
#' The \code{-outfmt 6} fields that \code{blast_search()} appends so that
#' \code{blast_report()} can group the hits.
#'
#' @return A character vector.
#' @examples
#' blast_report_fields()
#' @export
blast_report_fields <- function() blast_report_fields_rs()

#' Search a RefgetStore and report where hits come from
#'
#' \code{blast_search()} searches the selected collections of a store with rgblastn (or, with
#' \code{engine = "blastn", db = ...} and no store, a database you built with stock blastn);
#' \code{blast_report()} groups
#' tabular hits by genome and species; \code{blast_origin()} searches and reports in one step;
#' \code{blast_reps()} lists each species' representative.
#'
#' @inheritParams blast_build
#' @param query Query FASTA file.
#' @param out Engine output file; without it the output is returned as \code{output}.
#' @param outfmt Engine \code{-outfmt}.
#' @param args More engine arguments, as a character vector.
#' @param hits Tabular hits file(s).
#' @return \code{blast_search()}: the run metadata, with \code{output}.
#'   \code{blast_report()} and \code{blast_origin()}: \code{summary} (a data.frame, one row per
#'   query), \code{report}, and for \code{blast_origin()} \code{meta}.
#'   \code{blast_reps()}: \code{groups} and \code{summary}.
#' @export
blast_search <- function(store = NULL, query = NULL, out = NULL, outfmt = NULL,
                         args = character(), ..., progress = FALSE) {
  res <- .blast_run(c("search", .blast_store(store), .blast_flags(list(...)),
                      .blast_engine_args(query, out, outfmt, args)), progress)
  if (!identical(res$exit_code, 0)) {
    tail <- utils::tail(as.character(res$engine_stderr), 5L)
    stop(errorCondition(
      paste0(res$engine$name, " exited with status ", res$exit_code,
             if (length(tail)) paste0(":\n", paste(tail, collapse = "\n")) else ""),
      returncode = res$exit_code,
      class = c("gtars_rgblast_engine_failed", "gtars_rgblast_error")))
  }
  res$engine_stderr <- NULL
  res
}

#' @rdname blast_search
#' @export
blast_report <- function(hits, ..., progress = FALSE) {
  res <- .blast_run(c("report", .blast_text(hits), .blast_flags(list(...))), progress)
  .blast_report_result(res)
}

#' @rdname blast_search
#' @export
blast_origin <- function(store = NULL, query = NULL, args = character(), ..., progress = FALSE) {
  lead <- c("origin", .blast_store(store))
  if (!is.null(query)) lead <- c(lead, "--query", .blast_text(query))
  res <- .blast_run(c(lead, .blast_flags(list(...)), .blast_engine_args(args = args)), progress)
  .blast_report_result(res)
}

#' @rdname blast_search
#' @export
blast_reps <- function(store, ..., progress = FALSE) {
  res <- .blast_run(c("reps", "--store", .blast_text(store), .blast_flags(list(...))), progress)
  res$groups <- .blast_df(res$groups)
  res
}

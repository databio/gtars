#' @useDynLib gtars, .registration = TRUE
NULL

#' A fast, Rust-based findOverlaps analog.
#'
#' @param query   A GRanges object
#' @param subject A Tokenizer object (your Rust-based subject)
#' @param ...     Further arguments for consistency with findOverlaps signature
#'
#' @return A Hits object containing the query/subject index matches.
#'
findOverlapsFast <- function(query, subject, ...) {
    if (!inherits(query, "GRanges")) {
        stop("`query` must be a GRanges object.")
    }
    if (!inherits(subject, "Tokenizer")) {
        stop("`subject` must be a Tokenizer object.")
    }

    chrs   <- as.character(seqnames(query))
    starts <- start(query)
    ends   <- end(query)

    # returns a list of 2 integer vectors: [queryHits, subjectHits]
    overlap_indices <- subject$tokenize(
        chrs   = chrs,
        starts = starts,
        ends   = ends
    )
    q_hits <- overlap_indices[[1]]
    s_hits <- overlap_indices[[2]]

    n_universe <- subject$universe_length()

    # Construct a standard Hits object with S4Vectors
    # (nLnode = length of the query, nRnode = total intervals in subject’s universe)
    hits <- S4Vectors::Hits(
        from       = q_hits,
        to         = s_hits,
        nLnode     = length(query),
        nRnode     = n_universe,
        sort.by.query = TRUE
    )

    hits
}
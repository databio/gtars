library(GenomicRanges)
library(rtracklayer)

# First create our GRanges objects
set_A <- GRanges(
  seqnames = "chr1",
  ranges = IRanges(
    start = c(1, 4, 8, 12, 15, 20, 25),
    end = c(3, 6, 10, 14, 17, 22, 27)
  )
)

set_AA <- GRanges(
  seqnames = "chr1",
  ranges = IRanges(
    start = c(2, 4, 8),
    end = c(3, 6, 10)
  )
)


set_B <- GRangesList(
  group1 = GRanges(
    seqnames = "chr1",
    ranges = IRanges(
      start = c(2, 7, 12, 16, 21),
      end = c(4, 9, 15, 18, 23)
    )
  ),
  group2 = GRanges(
    seqnames = "chr1",
    ranges = IRanges(
      start = c(5, 11, 16, 19, 24),
      end = c(7, 13, 18, 21, 26)
    )
  ),
  group3 = GRanges(
    seqnames = "chr1",
    ranges = IRanges(
      start = c(3, 8, 13, 17, 22),
      end = c(5, 10, 15, 19, 24)
    )
  )
)


export(set_A, '/Users/sam/Documents/Work/gtars/bindings/r/tests/set_A.bed', format="BED")
export(set_AA, '/Users/sam/Documents/Work/gtars/bindings/r/tests/set_AA.bed', format="BED" )

# rextendr::document()

gtars_create <- gtars::r_igd_create('/Users/sam/Documents/Work/episcope/.test/igd/', '/Users/sam/Documents/Work/episcope/.test/test_paths.txt')
gtars_count <- gtars::r_igd_search(database_path = '/Users/sam/Documents/Work/episcope/.test/igd/igd_database.igd', query_path = '/Users/sam/Documents/Work/episcope/.test/set_A.bed')

userSets_beds <- c('/Users/sam/Documents/Work/episcope/.test/set_A.bed', '/Users/sam/Documents/Work/episcope/.test/set_AA.bed')
db_path <- '/Users/sam/Documents/Work/episcope/.test/igd/igd_database.igd'


## test lapply
r_igd_search_rev <- function(query_path = query_path, database_path = database_path) {
  gtars::r_igd_search(database_path = database_path, query_path = query_path)
}

lapply(userSets_beds, r_igd_search_rev, db_path)
lapply(geneSetDatabaseOverlaps, function(x) as.numeric(as.character(x[," number of hits"])))
       
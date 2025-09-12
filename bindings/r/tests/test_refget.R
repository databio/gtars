library(rextendr)

setwd("/Users/sam/Documents/Work/gtars/bindings/r")
rextendr::clean()
rextendr::document()

readable <- 'ATCG'
sha512t24u_digest <- gtars::r_sha512t24u_digest(readable)
md5_digest <- gtars::r_md5_digest(readable)

fasta_path <- "../../gtars/tests/data/fasta/base.fa"
result <- digest_fasta_raw(fasta_path)

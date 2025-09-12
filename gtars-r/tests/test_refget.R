library(rextendr)

# setwd('../')
setwd('/Users/sam/Documents/Work/gtars/gtars-r')
rextendr::clean()
rextendr::document()

readable <- 'ATCG'
gtars::sha512t24u_digest(readable)
gtars::md5_digest(readable)


fasta_path <- '../../gtars/tests/data/fasta/base.fa'
result_raw <- gtars::digest_fasta_raw(fasta_path)
result <- gtars::digest_fasta(fasta_path)

result
result@sequences
result@sequences[[1]]@metadata
result@lvl1



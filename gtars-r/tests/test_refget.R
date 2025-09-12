library(rextendr)

# setwd('../')
setwd('/Users/sam/Documents/Work/gtars/gtars-r')
rextendr::clean()
rextendr::document()

readable <- 'ATCG'
sha512t24u_digest <- gtars::r_sha512t24u_digest(readable)
md5_digest <- gtars::r_md5_digest(readable)


fasta_path <- '../../gtars/tests/data/fasta/base.fa'
# result <- digest_fasta_raw(fasta_path)
result <- digest_fasta(fasta_path)

result
result@sequences
result@sequences[[1]]@metadata
result@lvl1



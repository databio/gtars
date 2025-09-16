library(rextendr)

# setwd('../')
setwd('/Users/sam/Documents/Work/gtars/gtars-r')
rextendr::clean()
rextendr::document()

readable <- 'ATCG'
gtars::sha512t24u_digest(readable)
gtars::md5_digest(readable)


fasta_path <- '../../gtars/tests/data/fasta/base.fa'
result <- gtars::digest_fasta(fasta_path)

result
result@sequences
result@sequences[[1]]@metadata
result@sequences[[1]]@data
result@lvl1

store <- global_refget_store('raw')
store

import_fasta(store, fasta_path)
seq <- get_sequence_by_id(store, 'iYtREV555dUFKg2_agSJW6suquUyPpMw')
get_substring(store, 'iYtREV555dUFKg2_agSJW6suquUyPpMw', 0, 4)
seq2 <- get_sequence_by_collection_and_name(store, result@digest, result@sequences[[2]]@metadata@name)

length(store)



# Create temporary directory and files
temp_dir <- tempdir()
temp_bed_path <- file.path(temp_dir, "test.bed")
temp_output_fa_path <- file.path(temp_dir, "output.fa")


write_store_to_directory(store, temp_dir, "sequences/%s2/%s.seq")
store_load <- load_from_directory(temp_dir)


fasta_content <- paste(
  ">chr1",
  "ATGCATGCATGC",
  ">chr2", 
  "GGGGAAAA",
  sep = "\n"
)

temp_fasta_path <- file.path(temp_dir, "test.fa")
writeLines(fasta_content, temp_fasta_path)

# Create store and import FASTA
store <- global_refget_store("encoded")
import_fasta(store, temp_fasta_path)
result <- digest_fasta(temp_fasta_path)
collection_digest <- result@digest

chr1_sha <- gtars::sha512t24u_digest("ATGCATGCATGC")
chr1_md5 <- gtars::md5_digest("ATGCATGCATGC")
chr2_sha <- gtars::sha512t24u_digest("GGGGAAAA")
chr2_md5 <- gtars::md5_digest("GGGGAAAA")

# Create BED file content
bed_content <- paste(
  "chr1\t0\t5",
  "chr1\t8\t12", 
  "chr2\t0\t4",
  "chr_nonexistent\t10\t20",
  "chr1\t-5\t100",
  sep = "\n"
)


# Write BED file
writeLines(bed_content, temp_bed_path)


# Call the BED file function
get_seqs_bed_file(store, collection_digest, temp_bed_path, temp_output_fa_path)

# Read and check the output
output_content <- readLines(temp_output_fa_path)
cat("Output FASTA content:\n")
cat(paste(output_content, collapse = "\n"))

# Expected output should be:
expected_content <- paste(
  paste0(">chr1 12 dna2bit ", chr1_sha, " ", chr1_md5),
  "ATGCAATGC",
  paste0(">chr2 8 dna2bit ", chr2_sha, " ", chr2_md5), 
  "GGGG",
  sep = "\n"
)

cat("\nExpected content:\n")
cat(expected_content)


vec_result <- get_seqs_bed_file_to_vec(store, collection_digest, temp_bed_path)
vec_result[[1]]@sequence
vec_result[[1]]@chrom_name
vec_result[[1]]@start
vec_result[[1]]@end

vec_result[[2]]@sequence
vec_result[[2]]@chrom_name
vec_result[[2]]@start
vec_result[[2]]@end

vec_result[[3]]@sequence
vec_result[[3]]@chrom_name
vec_result[[3]]@start
vec_result[[3]]@end

# Create expected results as R objects (not S4 classes for simplicity)
expected_vec <- list(
  list(sequence = "ATGCA", chrom_name = "chr1", start = 0L, end = 5L),
  list(sequence = "ATGC", chrom_name = "chr1", start = 8L, end = 12L),
  list(sequence = "GGGG", chrom_name = "chr2", start = 0L, end = 4L)
)


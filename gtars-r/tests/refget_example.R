library(rextendr)

### either use rextendr or build/instsall the package
setwd('../')
# setwd('/Users/sam/Documents/Work/gtars/gtars-r')
rextendr::clean()
rextendr::document()

### string to digest
readable <- 'ATCG'
gtars::sha512t24u_digest(readable)
gtars::md5_digest(readable)

### digest fasta file
fasta_path <- '../../gtars/tests/data/fasta/base.fa'
result <- gtars::digest_fasta(fasta_path)

result
result@sequences
result@sequences[[1]]@metadata
result@sequences[[1]]@data
result@lvl1

### refget store
store <- global_refget_store('raw')
store

### import fasta file into store and retrieve sequences from store
import_fasta(store, fasta_path)
seq <- get_sequence_by_id(store, result@sequences[[1]]@metadata@sha512t24u)
seq@metadata
seq2 <- get_sequence_by_collection_and_name(store, result@digest, result@sequences[[2]]@metadata@name)
seq2@metadata
get_substring(store, result@sequences[[1]]@metadata@sha512t24u, 0, 4)

### save and load store
temp_dir <- tempdir()
temp_fasta_path <- file.path(temp_dir, 'test.fa')
temp_bed_path <- file.path(temp_dir, 'test.bed')
temp_output_fa_path <- file.path(temp_dir, 'output.fa')

write_store_to_directory(store, temp_dir, 'sequences/%s2/%s.seq')
store_load <- load_from_directory(temp_dir)

### extract BED file sequences from store
fasta_content <- paste(
  '>chr1',
  'ATGCATGCATGC',
  '>chr2', 
  'GGGGAAAA',
  sep = '\n'
)
writeLines(fasta_content, temp_fasta_path)

store <- global_refget_store('encoded')
import_fasta(store, temp_fasta_path)
result <- digest_fasta(temp_fasta_path)

chr1_sha <- gtars::sha512t24u_digest('ATGCATGCATGC')
chr1_md5 <- gtars::md5_digest('ATGCATGCATGC')
chr2_sha <- gtars::sha512t24u_digest('GGGGAAAA')
chr2_md5 <- gtars::md5_digest('GGGGAAAA')

bed_content <- paste(
  'chr1\t0\t5',
  'chr1\t8\t12', 
  'chr2\t0\t4',
  'chr_nonexistent\t10\t20',
  'chr1\t-5\t100',
  sep = '\n'
)
writeLines(bed_content, temp_bed_path)

### extract BED file sequences from store as FASTA
get_seqs_bed_file(store, result@digest, temp_bed_path, temp_output_fa_path)

output_content <- readLines(temp_output_fa_path)
cat('Output FASTA content:\n')
cat(paste(output_content, collapse = '\n'))

expected_content <- paste(
  paste0('>chr1 12 dna2bit ', chr1_sha, ' ', chr1_md5),
  'ATGCAATGC',
  paste0('>chr2 8 dna2bit ', chr2_sha, ' ', chr2_md5), 
  'GGGG',
  sep = '\n'
)
cat('\nExpected content:\n')
cat(expected_content)

### extract BED file sequences from store into memory
vec_result <- get_seqs_bed_file_to_vec(store, result@digest, temp_bed_path)
output_content <- data.frame(rbind(
  list(sequence = vec_result[[1]]@sequence, chrom_name = vec_result[[1]]@chrom_name, start = vec_result[[1]]@start, end = vec_result[[1]]@end),
  list(sequence = vec_result[[2]]@sequence, chrom_name = vec_result[[2]]@chrom_name, start = vec_result[[2]]@start, end = vec_result[[2]]@end),
  list(sequence = vec_result[[3]]@sequence, chrom_name = vec_result[[3]]@chrom_name, start = vec_result[[3]]@start, end = vec_result[[3]]@end)
))
print(output_content)

expected_vec <- data.frame(rbind(
  list(sequence = 'ATGCA', chrom_name = 'chr1', start = 0L, end = 5L),
  list(sequence = 'ATGC', chrom_name = 'chr1', start = 8L, end = 12L),
  list(sequence = 'GGGG', chrom_name = 'chr2', start = 0L, end = 4L)
))
print(expected_vec)

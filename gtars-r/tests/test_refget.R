library(testthat)
library(gtars)

test_that("digest functions work correctly", {
  readable <- 'ATCG'
  
  sha_result <- gtars::sha512t24u_digest(readable)
  expect_type(sha_result, "character")
  expect_true(nchar(sha_result) > 0)
  
  md5_result <- gtars::md5_digest(readable)
  expect_type(md5_result, "character")
  expect_true(nchar(md5_result) > 0)
})

test_that("fasta file digestion works", {
  fasta_path <- '../../gtars/tests/data/fasta/base.fa'
  skip_if_not(file.exists(fasta_path), "Test FASTA file not found")
  
  result <- gtars::digest_fasta(fasta_path)
  
  expect_s4_class(result, "FastaDigest")
  expect_true(length(result@sequences) > 0)
  expect_s4_class(result@sequences[[1]], "Sequence")
  expect_s4_class(result@sequences[[1]]@metadata, "SequenceMetadata")
  expect_type(result@sequences[[1]]@data, "raw")
  expect_type(result@lvl1, "character")
})

test_that("refget store initialization works", {
  store <- refget_store('raw')
  expect_s4_class(store, "RefgetStore")
})

test_that("fasta import and sequence retrieval works", {
  fasta_path <- '../../gtars/tests/data/fasta/base.fa'
  skip_if_not(file.exists(fasta_path), "Test FASTA file not found")
  
  store <- refget_store('raw')
  result <- gtars::digest_fasta(fasta_path)
  import_fasta(store, fasta_path)
  
  seq <- get_sequence(store, result@sequences[[1]]@metadata@sha512t24u)
  expect_s4_class(seq, "Sequence")
  expect_s4_class(seq@metadata, "SequenceMetadata")

  seq2 <- get_sequence_by_name(store, result@digest, result@sequences[[2]]@metadata@name)
  expect_s4_class(seq2, "Sequence")
  expect_s4_class(seq2@metadata, "SequenceMetadata")
  
  substring <- get_substring(store, result@sequences[[1]]@metadata@sha512t24u, 0, 4)
  expect_type(substring, "character")
  expect_true(nchar(substring) == 4)
})

test_that("store save and load works", {
  fasta_path <- '../../gtars/tests/data/fasta/base.fa'
  skip_if_not(file.exists(fasta_path), "Test FASTA file not found")
  
  temp_dir <- tempdir()
  store <- refget_store('raw')
  import_fasta(store, fasta_path)
  
  write_store_to_directory(store, temp_dir, 'sequences/%s2/%s.seq')
  expect_true(dir.exists(temp_dir))
  
  store_load <- open_local(temp_dir)
  expect_s4_class(store_load, "RefgetStore")
})

test_that("BED file sequence extraction to FASTA works", {
  temp_dir <- tempdir()
  temp_fasta_path <- file.path(temp_dir, 'test.fa')
  temp_bed_path <- file.path(temp_dir, 'test.bed')
  temp_output_fa_path <- file.path(temp_dir, 'output.fa')
  
  fasta_content <- paste(
    '>chr1',
    'ATGCATGCATGC',
    '>chr2', 
    'GGGGAAAA',
    sep = '\n'
  )
  writeLines(fasta_content, temp_fasta_path)
  
  store <- refget_store('encoded')
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
  
  get_seqs_bed_file(store, result@digest, temp_bed_path, temp_output_fa_path)
  
  expect_true(file.exists(temp_output_fa_path))
  output_content <- readLines(temp_output_fa_path)
  expect_true(length(output_content) > 0)
  
  expected_content <- paste(
    paste0('>chr1 12 dna2bit ', chr1_sha, ' ', chr1_md5),
    'ATGCAATGC',
    paste0('>chr2 8 dna2bit ', chr2_sha, ' ', chr2_md5), 
    'GGGG',
    sep = '\n'
  )
  
  expect_equal(paste(output_content, collapse = '\n'), expected_content)
  
  unlink(c(temp_fasta_path, temp_bed_path, temp_output_fa_path))
})

test_that("BED file sequence extraction to vector works", {
  temp_dir <- tempdir()
  temp_fasta_path <- file.path(temp_dir, 'test_vec.fa')
  temp_bed_path <- file.path(temp_dir, 'test_vec.bed')
  
  fasta_content <- paste(
    '>chr1',
    'ATGCATGCATGC',
    '>chr2', 
    'GGGGAAAA',
    sep = '\n'
  )
  writeLines(fasta_content, temp_fasta_path)
  
  store <- refget_store('encoded')
  import_fasta(store, temp_fasta_path)
  result <- digest_fasta(temp_fasta_path)
  
  bed_content <- paste(
    'chr1\t0\t5',
    'chr1\t8\t12', 
    'chr2\t0\t4',
    'chr_nonexistent\t10\t20',
    'chr1\t-5\t100',
    sep = '\n'
  )
  writeLines(bed_content, temp_bed_path)
  
  vec_result <- get_seqs_bed_file_to_vec(store, result@digest, temp_bed_path)
  
  expect_equal(length(vec_result), 3)
  
  output_content <- data.frame(
    sequence = sapply(vec_result, function(x) x@sequence),
    chrom_name = sapply(vec_result, function(x) x@chrom_name),
    start = sapply(vec_result, function(x) x@start),
    end = sapply(vec_result, function(x) x@end),
    stringsAsFactors = FALSE
  )
  
  expected_vec <- data.frame(
    sequence = c('ATGCA', 'ATGC', 'GGGG'),
    chrom_name = c('chr1', 'chr1', 'chr2'),
    start = c(0L, 8L, 0L),
    end = c(5L, 12L, 4L),
    stringsAsFactors = FALSE
  )
  
  expect_equal(output_content, expected_vec)
  
  unlink(c(temp_fasta_path, temp_bed_path))
})

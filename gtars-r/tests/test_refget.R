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

  expect_s4_class(result, "SequenceCollection")
  expect_true(length(result@sequences) > 0)
  expect_s4_class(result@sequences[[1]], "SequenceRecord")
  expect_s4_class(result@sequences[[1]]@metadata, "SequenceMetadata")
  expect_type(result@sequences[[1]]@data, "character")
  expect_s4_class(result@lvl1, "SeqColDigestLvl1")
})

test_that("refget store initialization works", {
  store <- refget_store('raw')
  expect_s4_class(store, "RefgetStore")

  # Test default mode (encoded)
  store2 <- refget_store()
  expect_s4_class(store2, "RefgetStore")
})

test_that("fasta import and sequence retrieval works", {
  fasta_path <- '../../gtars/tests/data/fasta/base.fa'
  skip_if_not(file.exists(fasta_path), "Test FASTA file not found")

  store <- refget_store('raw')
  result <- gtars::digest_fasta(fasta_path)
  import_fasta(store, fasta_path)

  seq <- get_sequence(store, result@sequences[[1]]@metadata@sha512t24u)
  expect_s4_class(seq, "SequenceRecord")
  expect_s4_class(seq@metadata, "SequenceMetadata")

  seq2 <- get_sequence_by_name(store, result@digest, result@sequences[[2]]@metadata@name)
  expect_s4_class(seq2, "SequenceRecord")
  expect_s4_class(seq2@metadata, "SequenceMetadata")

  substring <- get_substring(store, result@sequences[[1]]@metadata@sha512t24u, 0, 4)
  expect_type(substring, "character")
  expect_true(nchar(substring) == 4)
})

test_that("store save and load works", {
  fasta_path <- '../../gtars/tests/data/fasta/base.fa'
  skip_if_not(file.exists(fasta_path), "Test FASTA file not found")

  temp_dir <- tempdir()
  store_dir <- file.path(temp_dir, "refget_test_store")
  unlink(store_dir, recursive = TRUE)

  store <- refget_store('raw')
  import_fasta(store, fasta_path)

  write_store_to_directory(store, store_dir, 'sequences/%s2/%s.seq')
  expect_true(dir.exists(store_dir))

  store_load <- refget_store_open_local(store_dir)
  expect_s4_class(store_load, "RefgetStore")

  unlink(store_dir, recursive = TRUE)
})

test_that("store configuration functions work", {
  store <- refget_store('raw')

  # Test storage mode
  expect_equal(storage_mode(store), "raw")
  enable_encoding(store)
  expect_equal(storage_mode(store), "encoded")
  disable_encoding(store)
  expect_equal(storage_mode(store), "raw")

  # Test quiet mode
  expect_false(is_quiet(store))
  set_quiet(store, TRUE)
  expect_true(is_quiet(store))
  set_quiet(store, FALSE)
  expect_false(is_quiet(store))

  # Test persistence
  expect_false(is_persisting(store))
})

test_that("collection listing and stats work", {
  fasta_path <- '../../gtars/tests/data/fasta/base.fa'
  skip_if_not(file.exists(fasta_path), "Test FASTA file not found")

  store <- refget_store()
  import_fasta(store, fasta_path)

  # Test list_collections
  collections <- list_collections(store)
  expect_true(length(collections) > 0)
  expect_s4_class(collections[[1]], "SequenceCollectionMetadata")

  # Test list_sequences
  sequences <- list_sequences(store)
  expect_true(length(sequences) > 0)
  expect_s4_class(sequences[[1]], "SequenceMetadata")

  # Test stats
  store_stats <- stats(store)
  expect_type(store_stats, "list")
  expect_true("n_sequences" %in% names(store_stats))
  expect_true("n_collections" %in% names(store_stats))
})

test_that("seqcol level1 and level2 work", {
  fasta_path <- '../../gtars/tests/data/fasta/base.fa'
  skip_if_not(file.exists(fasta_path), "Test FASTA file not found")

  store <- refget_store()
  import_fasta(store, fasta_path)

  result <- digest_fasta(fasta_path)
  collection_digest <- result@digest

  # Test get_level1
  lvl1 <- get_level1(store, collection_digest)
  expect_type(lvl1, "list")
  expect_true("names" %in% names(lvl1))
  expect_true("lengths" %in% names(lvl1))
  expect_true("sequences" %in% names(lvl1))

  # Test get_level2
  lvl2 <- get_level2(store, collection_digest)
  expect_type(lvl2, "list")
  expect_true("names" %in% names(lvl2))
  expect_true(is.character(lvl2$names))
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

test_that("sequence alias operations work", {
  store <- refget_store()

  # Create test sequence
  temp_dir <- tempdir()
  temp_fasta_path <- file.path(temp_dir, 'alias_test.fa')
  writeLines(c('>seq1', 'ATGCATGC'), temp_fasta_path)
  import_fasta(store, temp_fasta_path)
  result <- digest_fasta(temp_fasta_path)
  seq_digest <- result@sequences[[1]]@metadata@sha512t24u

  # Add alias
  add_sequence_alias(store, "test_ns", "my_alias", seq_digest)

  # List namespaces
  namespaces <- list_sequence_alias_namespaces(store)
  expect_true("test_ns" %in% namespaces)

  # List aliases in namespace
  aliases <- list_sequence_aliases(store, "test_ns")
  expect_true("my_alias" %in% aliases)

  # Get by alias
  seq_by_alias <- get_sequence_metadata_by_alias(store, "test_ns", "my_alias")
  expect_s4_class(seq_by_alias, "SequenceMetadata")

  # Get aliases for sequence
  aliases_for_seq <- get_aliases_for_sequence(store, seq_digest)
  expect_true(length(aliases_for_seq) > 0)

  # Remove alias
  removed <- remove_sequence_alias(store, "test_ns", "my_alias")
  expect_true(removed)

  unlink(temp_fasta_path)
})

test_that("collection alias operations work", {
  store <- refget_store()

  # Create test collection
  temp_dir <- tempdir()
  temp_fasta_path <- file.path(temp_dir, 'coll_alias_test.fa')
  writeLines(c('>seq1', 'ATGCATGC'), temp_fasta_path)
  import_fasta(store, temp_fasta_path)
  result <- digest_fasta(temp_fasta_path)
  coll_digest <- result@digest

  # Add alias
  add_collection_alias(store, "test_ns", "my_collection", coll_digest)

  # List namespaces
  namespaces <- list_collection_alias_namespaces(store)
  expect_true("test_ns" %in% namespaces)

  # Get by alias
  coll_by_alias <- get_collection_metadata_by_alias(store, "test_ns", "my_collection")
  expect_s4_class(coll_by_alias, "SequenceCollectionMetadata")

  # Remove alias
  removed <- remove_collection_alias(store, "test_ns", "my_collection")
  expect_true(removed)

  unlink(temp_fasta_path)
})

test_that("compare function works", {
  store <- refget_store()

  # Create two test collections
  temp_dir <- tempdir()
  temp_fasta1 <- file.path(temp_dir, 'compare_test1.fa')
  temp_fasta2 <- file.path(temp_dir, 'compare_test2.fa')
  writeLines(c('>seq1', 'ATGCATGC'), temp_fasta1)
  writeLines(c('>seq1', 'ATGCATGC', '>seq2', 'GGGG'), temp_fasta2)

  import_fasta(store, temp_fasta1)
  import_fasta(store, temp_fasta2)

  result1 <- digest_fasta(temp_fasta1)
  result2 <- digest_fasta(temp_fasta2)

  comparison <- compare(store, result1@digest, result2@digest)
  expect_type(comparison, "list")
  expect_true("digests" %in% names(comparison))
  expect_true("attributes" %in% names(comparison))

  unlink(c(temp_fasta1, temp_fasta2))
})

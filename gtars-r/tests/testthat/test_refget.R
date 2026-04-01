library(testthat)
library(gtars)

fasta_path <- file.path(
  test_path(), '..', '..', '..', 'tests', 'data', 'fasta', 'base.fa'
)

skip_if_no_fasta <- function() {
  skip_if_not(file.exists(fasta_path), "Test FASTA file not found")
}

# =========================================================================
# Digest functions
# =========================================================================

test_that("sha512t24u_digest and md5_digest work", {
  sha_result <- gtars::sha512t24u_digest("ATCG")
  expect_type(sha_result, "character")
  expect_true(nchar(sha_result) > 0)

  md5_result <- gtars::md5_digest("ATCG")
  expect_type(md5_result, "character")
  expect_true(nchar(md5_result) > 0)
})

test_that("digest_fasta returns SequenceCollection", {
  skip_if_no_fasta()
  result <- digest_fasta(fasta_path)

  expect_s4_class(result, "SequenceCollection")
  expect_true(length(result@sequences) > 0)
  expect_s4_class(result@sequences[[1]], "SequenceRecord")
  expect_s4_class(result@sequences[[1]]@metadata, "SequenceMetadata")
  expect_type(result@sequences[[1]]@data, "character")
  expect_s4_class(result@lvl1, "SeqColDigestLvl1")
})

# =========================================================================
# Store constructors
# =========================================================================

test_that("refget_store creates in-memory store", {
  store <- refget_store("raw")
  expect_s4_class(store, "RefgetStore")

  store2 <- refget_store()
  expect_s4_class(store2, "RefgetStore")
})

test_that("refget_store_on_disk creates disk-backed store", {
  disk_path <- file.path(tempdir(), "test_disk_store")
  on.exit(unlink(disk_path, recursive = TRUE))

  store <- refget_store_on_disk(disk_path)
  expect_s4_class(store, "RefgetStore")
})

# =========================================================================
# Store configuration
# =========================================================================

test_that("encoding mode toggling works", {
  store <- refget_store("raw")
  expect_equal(storage_mode(store), "raw")

  enable_encoding(store)
  expect_equal(storage_mode(store), "encoded")

  disable_encoding(store)
  expect_equal(storage_mode(store), "raw")
})

test_that("quiet mode toggling works", {
  store <- refget_store()
  expect_false(is_quiet(store))

  set_quiet(store, TRUE)
  expect_true(is_quiet(store))

  set_quiet(store, FALSE)
  expect_false(is_quiet(store))
})

test_that("persistence toggling works", {
  store <- refget_store()
  expect_false(is_persisting(store))

  persist_path <- file.path(tempdir(), "test_persist")
  on.exit(unlink(persist_path, recursive = TRUE))

  enable_persistence(store, persist_path)
  expect_true(is_persisting(store))
  expect_true(grepl("test_persist", cache_path(store)))

  disable_persistence(store)
  expect_false(is_persisting(store))
})

# =========================================================================
# Import and sequence retrieval
# =========================================================================

test_that("import_fasta and get_sequence work", {
  skip_if_no_fasta()
  store <- refget_store("raw")
  result <- digest_fasta(fasta_path)
  import_fasta(store, fasta_path)

  seq <- get_sequence(store, result@sequences[[1]]@metadata@sha512t24u)
  expect_s4_class(seq, "SequenceRecord")
  expect_s4_class(seq@metadata, "SequenceMetadata")
  expect_true(nchar(seq@data) > 0)
})

test_that("get_sequence_by_name works", {
  skip_if_no_fasta()
  store <- refget_store("raw")
  result <- digest_fasta(fasta_path)
  import_fasta(store, fasta_path)

  seq <- get_sequence_by_name(store, result@digest, result@sequences[[2]]@metadata@name)
  expect_s4_class(seq, "SequenceRecord")
  expect_equal(seq@metadata@sha512t24u, result@sequences[[2]]@metadata@sha512t24u)
})

test_that("get_sequence_metadata works", {
  skip_if_no_fasta()
  store <- refget_store()
  result <- digest_fasta(fasta_path)
  import_fasta(store, fasta_path)

  meta <- get_sequence_metadata(store, result@sequences[[1]]@metadata@sha512t24u)
  expect_s4_class(meta, "SequenceMetadata")
  expect_equal(meta@name, result@sequences[[1]]@metadata@name)
  expect_type(meta@description, "character")
})

test_that("get_substring works", {
  skip_if_no_fasta()
  store <- refget_store("raw")
  result <- digest_fasta(fasta_path)
  import_fasta(store, fasta_path)

  substring <- get_substring(store, result@sequences[[1]]@metadata@sha512t24u, 0, 4)
  expect_type(substring, "character")
  expect_equal(nchar(substring), 4)
})

test_that("backwards compat aliases work", {
  skip_if_no_fasta()
  store <- refget_store("raw")
  result <- digest_fasta(fasta_path)
  import_fasta(store, fasta_path)
  digest <- result@sequences[[1]]@metadata@sha512t24u

  rec2 <- get_sequence_by_id(store, digest)
  expect_s4_class(rec2, "SequenceRecord")

  rec3 <- get_sequence_by_collection_and_name(
    store, result@digest, result@sequences[[1]]@metadata@name
  )
  expect_s4_class(rec3, "SequenceRecord")
})

# =========================================================================
# Collection operations
# =========================================================================

test_that("list_collections and list_sequences work", {
  skip_if_no_fasta()
  store <- refget_store()
  import_fasta(store, fasta_path)

  collections <- list_collections(store)
  expect_true(length(collections) > 0)
  expect_s4_class(collections[[1]], "SequenceCollectionMetadata")

  sequences <- list_sequences(store)
  expect_true(length(sequences) > 0)
  expect_s4_class(sequences[[1]], "SequenceMetadata")
})

test_that("get_collection works", {
  skip_if_no_fasta()
  store <- refget_store()
  result <- digest_fasta(fasta_path)
  import_fasta(store, fasta_path)

  coll <- get_collection(store, result@digest)
  expect_s4_class(coll, "SequenceCollection")
  expect_equal(coll@digest, result@digest)
})

test_that("get_collection_metadata works", {
  skip_if_no_fasta()
  store <- refget_store()
  result <- digest_fasta(fasta_path)
  import_fasta(store, fasta_path)

  meta <- get_collection_metadata(store, result@digest)
  expect_s4_class(meta, "SequenceCollectionMetadata")
  expect_equal(meta@digest, result@digest)
  expect_true(meta@n_sequences > 0)
})

test_that("is_collection_loaded works", {
  skip_if_no_fasta()
  store <- refget_store()
  result <- digest_fasta(fasta_path)
  import_fasta(store, fasta_path)

  expect_true(is_collection_loaded(store, result@digest))
})

test_that("iter_collections and iter_sequences work", {
  skip_if_no_fasta()
  store <- refget_store()
  import_fasta(store, fasta_path)

  colls <- iter_collections(store)
  expect_type(colls, "list")
  expect_s4_class(colls[[1]], "SequenceCollection")

  seqs <- iter_sequences(store)
  expect_type(seqs, "list")
  expect_s4_class(seqs[[1]], "SequenceRecord")
})

test_that("stats works", {
  skip_if_no_fasta()
  store <- refget_store()
  import_fasta(store, fasta_path)

  st <- stats(store)
  expect_type(st, "list")
  expect_true("n_sequences" %in% names(st))
  expect_true("n_collections" %in% names(st))
  expect_true("storage_mode" %in% names(st))
  expect_true(st$n_sequences > 0)
})

# =========================================================================
# Seqcol spec operations
# =========================================================================

test_that("get_level1 and get_level2 work", {
  skip_if_no_fasta()
  store <- refget_store()
  import_fasta(store, fasta_path)
  result <- digest_fasta(fasta_path)

  lvl1 <- get_level1(store, result@digest)
  expect_type(lvl1, "list")
  expect_true(all(c("names", "lengths", "sequences") %in% names(lvl1)))

  lvl2 <- get_level2(store, result@digest)
  expect_type(lvl2, "list")
  expect_true("names" %in% names(lvl2))
  expect_true(is.character(lvl2$names))
})

test_that("ancillary digests toggling works", {
  store <- refget_store()
  # Default may be TRUE or FALSE depending on build; just test toggling
  enable_ancillary_digests(store)
  expect_true(has_ancillary_digests(store))

  disable_ancillary_digests(store)
  expect_false(has_ancillary_digests(store))

  enable_ancillary_digests(store)
  expect_true(has_ancillary_digests(store))
})

test_that("attribute index and lookup work", {
  skip_if_no_fasta()
  store <- refget_store()
  import_fasta(store, fasta_path)
  result <- digest_fasta(fasta_path)

  # Test toggling
  expect_false(has_attribute_index(store))
  enable_attribute_index(store)
  expect_true(has_attribute_index(store))
  disable_attribute_index(store)
  expect_false(has_attribute_index(store))

  # find_collections_by_attribute works with brute-force (index disabled)
  lvl1 <- get_level1(store, result@digest)
  found <- find_collections_by_attribute(store, "names", lvl1$names)
  expect_true(length(found) > 0)

  attr_val <- get_attribute(store, "names", lvl1$names)
  expect_false(is.null(attr_val))
})

# =========================================================================
# Compare collections
# =========================================================================

test_that("compare works", {
  skip_if_no_fasta()
  store <- refget_store()

  temp_fasta1 <- file.path(tempdir(), "compare_test1.fa")
  temp_fasta2 <- file.path(tempdir(), "compare_test2.fa")
  on.exit(unlink(c(temp_fasta1, temp_fasta2)))

  writeLines(c(">seq1", "ATGCATGC"), temp_fasta1)
  writeLines(c(">seq1", "ATGCATGC", ">seq2", "GGGG"), temp_fasta2)

  import_fasta(store, temp_fasta1)
  import_fasta(store, temp_fasta2)

  result1 <- digest_fasta(temp_fasta1)
  result2 <- digest_fasta(temp_fasta2)

  comparison <- compare(store, result1@digest, result2@digest)
  expect_type(comparison, "list")
  expect_true("digests" %in% names(comparison))
  expect_true("attributes" %in% names(comparison))
})

# =========================================================================
# Write/Export operations
# =========================================================================

test_that("write_store_to_directory and reload work", {
  skip_if_no_fasta()
  store <- refget_store("raw")
  import_fasta(store, fasta_path)

  store_dir <- file.path(tempdir(), "refget_test_store")
  on.exit(unlink(store_dir, recursive = TRUE))
  unlink(store_dir, recursive = TRUE)

  write_store_to_directory(store, store_dir, "sequences/%s2/%s.seq")
  expect_true(dir.exists(store_dir))

  store_load <- refget_store_open_local(store_dir)
  expect_s4_class(store_load, "RefgetStore")

  st_orig <- stats(store)
  st_load <- stats(store_load)
  expect_equal(st_load$n_sequences, st_orig$n_sequences)
})

test_that("load_from_directory compat alias works", {
  skip_if_no_fasta()
  store <- refget_store("raw")
  import_fasta(store, fasta_path)

  store_dir <- file.path(tempdir(), "refget_compat_store")
  on.exit(unlink(store_dir, recursive = TRUE))
  unlink(store_dir, recursive = TRUE)

  write_store_to_directory(store, store_dir, "sequences/%s2/%s.seq")
  store_compat <- load_from_directory(store_dir)
  expect_s4_class(store_compat, "RefgetStore")
})

test_that("persistence round-trip with write_store works", {
  skip_if_no_fasta()
  persist_dir <- file.path(tempdir(), "refget_persist_rt")
  on.exit(unlink(persist_dir, recursive = TRUE))
  unlink(persist_dir, recursive = TRUE)

  store <- refget_store()
  set_quiet(store, TRUE)
  enable_persistence(store, persist_dir)
  import_fasta(store, fasta_path)
  write_store(store)

  store2 <- refget_store_open_local(persist_dir)
  st2 <- stats(store2)
  expect_true(st2$n_sequences > 0)
})

test_that("export_fasta creates valid FASTA", {
  skip_if_no_fasta()
  store <- refget_store()
  set_quiet(store, TRUE)
  result <- digest_fasta(fasta_path)
  import_fasta(store, fasta_path)

  export_path <- file.path(tempdir(), "refget_export.fa")
  on.exit(unlink(export_path))

  export_fasta(store, result@digest, export_path)
  expect_true(file.exists(export_path))
  lines <- readLines(export_path)
  expect_true(any(grepl("^>", lines)))
})

test_that("export_fasta with sequence_names subset works", {
  skip_if_no_fasta()
  store <- refget_store()
  set_quiet(store, TRUE)
  result <- digest_fasta(fasta_path)
  import_fasta(store, fasta_path)

  export_path <- file.path(tempdir(), "refget_export_subset.fa")
  on.exit(unlink(export_path))

  seq_names <- sapply(result@sequences[1:2], function(s) s@metadata@name)
  export_fasta(store, result@digest, export_path,
               sequence_names = seq_names, line_width = 60L)
  expect_true(file.exists(export_path))
})

test_that("export_fasta_by_digests works", {
  skip_if_no_fasta()
  store <- refget_store()
  set_quiet(store, TRUE)
  result <- digest_fasta(fasta_path)
  import_fasta(store, fasta_path)

  export_path <- file.path(tempdir(), "refget_export_digests.fa")
  on.exit(unlink(export_path))

  digests <- sapply(result@sequences[1:2], function(s) s@metadata@sha512t24u)
  export_fasta_by_digests(store, digests, export_path, line_width = 70L)
  expect_true(file.exists(export_path))
})

# =========================================================================
# BED file sequence extraction
# =========================================================================

test_that("BED file sequence extraction to FASTA works", {
  skip_if_no_fasta()
  store <- refget_store("encoded")
  import_fasta(store, fasta_path)
  result <- digest_fasta(fasta_path)

  temp_bed <- file.path(tempdir(), "refget_test.bed")
  temp_out <- file.path(tempdir(), "refget_output.fa")
  on.exit(unlink(c(temp_bed, temp_out)))

  writeLines(c("chr1\t0\t4", "chr2\t0\t4"), temp_bed)

  get_seqs_bed_file(store, result@digest, temp_bed, temp_out)
  expect_true(file.exists(temp_out))
  output_content <- readLines(temp_out)
  expect_true(length(output_content) > 0)
})

test_that("BED file sequence extraction to vector works", {
  skip_if_no_fasta()
  store <- refget_store("encoded")
  import_fasta(store, fasta_path)
  result <- digest_fasta(fasta_path)

  temp_bed <- file.path(tempdir(), "refget_test_vec.bed")
  on.exit(unlink(temp_bed))

  writeLines(c("chr1\t0\t4", "chr2\t0\t4"), temp_bed)

  vec <- get_seqs_bed_file_to_vec(store, result@digest, temp_bed)
  expect_type(vec, "list")
  expect_equal(length(vec), 2)
  expect_s4_class(vec[[1]], "RetrievedSequence")
  expect_true(nchar(vec[[1]]@sequence) > 0)
})

# =========================================================================
# Sequence aliases
# =========================================================================

test_that("sequence alias operations work", {
  skip_if_no_fasta()
  store <- refget_store()
  import_fasta(store, fasta_path)
  result <- digest_fasta(fasta_path)
  seq_digest <- result@sequences[[1]]@metadata@sha512t24u

  add_sequence_alias(store, "test_ns", "my_alias", seq_digest)

  namespaces <- list_sequence_alias_namespaces(store)
  expect_true("test_ns" %in% namespaces)

  aliases <- list_sequence_aliases(store, "test_ns")
  expect_true("my_alias" %in% aliases)

  seq_by_alias <- get_sequence_metadata_by_alias(store, "test_ns", "my_alias")
  expect_s4_class(seq_by_alias, "SequenceMetadata")

  aliases_for_seq <- get_aliases_for_sequence(store, seq_digest)
  expect_true(length(aliases_for_seq) > 0)

  removed <- remove_sequence_alias(store, "test_ns", "my_alias")
  expect_true(removed)
})

# =========================================================================
# Collection aliases
# =========================================================================

test_that("collection alias operations work", {
  skip_if_no_fasta()
  store <- refget_store()
  import_fasta(store, fasta_path)
  result <- digest_fasta(fasta_path)
  coll_digest <- result@digest

  add_collection_alias(store, "test_ns", "my_collection", coll_digest)

  namespaces <- list_collection_alias_namespaces(store)
  expect_true("test_ns" %in% namespaces)

  aliases <- list_collection_aliases(store, "test_ns")
  expect_true("my_collection" %in% aliases)

  coll_by_alias <- get_collection_metadata_by_alias(store, "test_ns", "my_collection")
  expect_s4_class(coll_by_alias, "SequenceCollectionMetadata")

  aliases_for_coll <- get_aliases_for_collection(store, coll_digest)
  expect_true(length(aliases_for_coll) > 0)

  removed <- remove_collection_alias(store, "test_ns", "my_collection")
  expect_true(removed)
})

# =========================================================================
# FHR metadata
# =========================================================================

test_that("FHR metadata operations work", {
  skip_if_no_fasta()
  skip_if_not_installed("jsonlite")

  store <- refget_store()
  import_fasta(store, fasta_path)
  result <- digest_fasta(fasta_path)
  coll_digest <- result@digest

  metadata <- list(organism = "test", assembly = "base")
  set_fhr_metadata(store, coll_digest, metadata)

  retrieved <- get_fhr_metadata(store, coll_digest)
  expect_false(is.null(retrieved))
  expect_equal(retrieved$organism, "test")
  expect_equal(retrieved$assembly, "base")

  fhr_list <- list_fhr_metadata(store)
  expect_true(length(fhr_list) > 0)

  remove_fhr_metadata(store, coll_digest)
  after_remove <- get_fhr_metadata(store, coll_digest)
  expect_null(after_remove)
})

test_that("FHR metadata with raw JSON string works", {
  skip_if_no_fasta()
  skip_if_not_installed("jsonlite")

  store <- refget_store()
  import_fasta(store, fasta_path)
  result <- digest_fasta(fasta_path)

  set_fhr_metadata(store, result@digest, '{"source": "test"}')
  raw_meta <- get_fhr_metadata(store, result@digest)
  expect_equal(raw_meta$source, "test")
  remove_fhr_metadata(store, result@digest)
})

# =========================================================================
# Show methods (just verify no errors)
# =========================================================================

test_that("show methods do not error", {
  skip_if_no_fasta()
  store <- refget_store()
  import_fasta(store, fasta_path)
  result <- digest_fasta(fasta_path)

  expect_no_error(capture.output(show(store)))
  expect_no_error(capture.output(show(result)))
  expect_no_error(capture.output(show(result@lvl1)))
  expect_no_error(capture.output(show(result@sequences[[1]])))
  expect_no_error(capture.output(show(result@sequences[[1]]@metadata)))

  coll_meta <- get_collection_metadata(store, result@digest)
  expect_no_error(capture.output(show(coll_meta)))
})

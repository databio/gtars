**Example: Full `RefgetStore` Usage**

```r
library(gtars)

temp_dir <- tempdir()

# 1. Prepare a dummy FASTA file
temp_fasta_path <- file.path(temp_dir, "source.fa")
fasta_content <- paste(
  ">chr1",
  "ATGCATGCATGCAGTCGTAGC",
  ">chr2",
  "GGGGAAAA",
  sep = "\n"
)
writeLines(fasta_content, temp_fasta_path)

# 2. Digest the FASTA to get collection info and digest
collection <- digest_fasta(temp_fasta_path)
collection_digest <- collection@digest
cat(sprintf("Source FASTA digested. Collection digest: %s\n", collection_digest))

# 3. Create an in-memory RefgetStore in Encoded mode
store <- refget_store("encoded")
print(store)

# 4. Add the FASTA to the store
added <- add_fasta(store, temp_fasta_path)
cat(sprintf("FASTA added to the store as collection %s\n", added$digest))

# 5. Get a sequence by its digest (from the first sequence in the collection)
seq_digest_chr1 <- collection[1]@metadata@sha512t24u
record_chr1 <- get_sequence(store, seq_digest_chr1)
if (!is.null(record_chr1)) {
  cat(sprintf("Retrieved sequence by digest: %s, length %s\n",
              record_chr1@metadata@name, record_chr1@metadata@length))
  cat(sprintf("  Sequence (full): %s\n",
              get_substring(store, seq_digest_chr1, 0, record_chr1@metadata@length)))
}

# Or look it up by collection and name
record_chr2 <- get_sequence_by_name(store, collection_digest, "chr2")

# 6. Get a substring (0-based, half-open)
sub_seq <- get_substring(store, seq_digest_chr1, 5, 15)
cat(sprintf("Substring from chr1[5:15]: %s\n", sub_seq))

# 7. Prepare a BED file for region retrieval
temp_bed_path <- file.path(temp_dir, "test.bed")
bed_content <- paste(
  "chr1\t0\t10",
  "chr2\t2\t6",
  "chr_nonexistent\t0\t5",
  sep = "\n"
)
writeLines(bed_content, temp_bed_path)

# 8. Retrieve sequences from the BED file into memory.
# Regions on unknown sequences (chr_nonexistent) are skipped with a warning.
retrieved_list <- get_seqs_bed_file_to_vec(store, collection_digest, temp_bed_path)
cat("Retrieved sequences from BED file (as list):\n")
for (rs in retrieved_list) {
  print(rs)
}

# The same, as a data.frame with columns sequence, chrom_name, start, end
retrieved_df <- get_seqs_bed_file_to_df(store, collection_digest, temp_bed_path)

# 9. Write BED regions to a new FASTA file.
# Unlike step 8, this stops with an error on unknown sequences,
# so give it a BED file with valid regions only.
valid_bed_path <- file.path(temp_dir, "valid.bed")
writeLines(paste("chr1\t0\t10", "chr2\t2\t6", sep = "\n"), valid_bed_path)
temp_output_fa_path <- file.path(temp_dir, "output.fa")
get_seqs_bed_file(store, collection_digest, valid_bed_path, temp_output_fa_path)
cat(sprintf("Retrieved sequences from BED file written to: %s\n", temp_output_fa_path))

# 10. Export the whole collection (or selected sequences) to FASTA
export_fasta(store, collection_digest, file.path(temp_dir, "all.fa"))
export_fasta(store, collection_digest, file.path(temp_dir, "chr1.fa"),
             sequence_names = c("chr1"), line_width = 60)

# 11. Write store to a new directory.
# In the template, %s2 is the first 2 characters of the digest and %s the full digest.
temp_saved_store_path <- file.path(temp_dir, "my_refget_store")
write_store_to_directory(store, temp_saved_store_path, "sequences/%s2/%s.seq")
cat(sprintf("Store saved to: %s\n", temp_saved_store_path))

# 12. Open the store from the directory
store_load <- refget_store_open_local(temp_saved_store_path)
cat(sprintf("Store successfully loaded from: %s\n", temp_saved_store_path))
print(store_load)
```

A saved store can also be served over HTTP and opened remotely:

```r
remote <- refget_store_open_remote(refget_default_cache_dir(), "https://example.com/store",
                                   fetch_jobs = 32L)

# Download many sequences into the shared cache at once, without loading them
digests <- vapply(list_sequences(remote), function(m) m@sha512t24u, character(1))
counters <- cache_sequences(remote, digests, progress = TRUE)

fetch_stats(remote)                       # requests, retries, bytes, cache hits
export_fasta(remote, "all.fa", all = TRUE, write_cache = TRUE)
refget_cache_info(refget_default_cache_dir())
refget_clear_cache(refget_default_cache_dir())   # deletes everything in the cache
```

`cache_path` is a cache root shared by every remote store: downloaded
sequences live in `sequences/` (named by digest, so remotes share them) and each
remote's metadata in `remotes/<key>/`. `refget_default_cache_dir()` is
`REFGET_CACHE_DIR`, else `$XDG_CACHE_HOME/refget`, else `~/.cache/refget`.
Remote stores need a server with HTTP range support: sequence files over 8 MiB
download as parallel 8 MiB ranges, and partial reads (`get_substring`) read only
the bytes they need. Running the same download again fetches only the sequences
still missing. See the [Python reference](../python/refget-store.md) for the
cache layout and download behavior.

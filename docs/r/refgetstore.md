**Example: Full `GlobalRefgetStore` Usage**

``` R
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

# 3. Initialize GlobalRefgetStore in Encoded mode
store <- refget_store("encoded")
cat(sprintf("Initialized store: %s\n", store))

# 4. Import FASTA into the store
import_fasta(store, temp_fasta_path)
cat(sprintf("FASTA imported into the store.\n"))

# 5. Get a sequence by its ID (using the digest from the first sequence in collection)
seq_digest_chr1 = collection[1]@metadata@sha512t24u
record_chr1 = get_sequence_by_id(store, seq_digest_chr1)
if (!is.null(record_chr1)) {
  cat(sprintf("Retrieved sequence by ID: %s, length %s\n", 
              record_chr1@metadata@name, record_chr1@metadata@length))
  cat(sprintf("  Sequence (full): %s\n", 
              get_substring(store, seq_digest_chr1, 0, record_chr1@metadata@length)))
}

# 6. Get a substring
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

# 8. Retrieve sequences from BED file to a list
retrieved_list <- get_seqs_bed_file_to_vec(store, collection_digest, temp_bed_path)
cat("Retrieved sequences from BED file (as list):\n")
for (rs in retrieved_list) {
  print(rs)
}

# 9. Retrieve sequences from BED file and write to new FASTA
temp_output_fa_path <- file.path(temp_dir, "output.fa")
get_seqs_bed_file(store, collection_digest, temp_bed_path, temp_output_fa_path)
cat(sprintf("Retrieved sequences from BED file written to: %s\n", temp_output_fa_path))

# 10. Write store to a new directory
temp_saved_store_path <- file.path(temp_dir, "my_refget_store")
write_store_to_directory(store, temp_saved_store_path, "{digest_prefix}/{digest}.gz")
cat(sprintf("Store saved to: %s\n", temp_saved_store_path))

# 11. Load store from the directory
store_load <- load_local(temp_saved_store_path)
cat(sprintf("Store successfully loaded from: %s\n", temp_saved_store_path))
```

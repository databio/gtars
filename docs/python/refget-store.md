# RefgetStore Python Reference

This is a Python-specific reference guide that provides quick examples for using the `RefgetStore` class from the `gtars.refget` module. For detailed information about the underlying RefgetStore file format specification, see [RefgetStore format specification](https://docs.refgenie.org/refget/reference/refgetstore-format/).

## Creating and Populating a Store

```python
from gtars.refget import RefgetStore

# Create a new store (defaults to Encoded mode, space-efficient)
store = RefgetStore.in_memory()
print(f"Initialized store: {store}")

# Add sequences from a FASTA file
store.add_sequence_collection_from_fasta("genome.fa")

# Inspect what's in the store
sequence_records = store.iter_sequences()
sequence_metadata = store.list_sequences()

# list_collections() returns a page of results plus pagination info
page = store.list_collections()
for coll in page["results"]:
    print(f"{coll.digest}: {coll.n_sequences} sequences")

# Access individual sequences
first_seq = sequence_records[0]
print(f"First sequence: {first_seq.metadata.name}")

# Decode sequence data to string
if first_seq.sequence:
    decoded = first_seq.decode()
    print(f"Sequence: {decoded}")
```

## Saving and Loading Local Stores

```python
import os

# Save the store to disk
store_path = "my_refget_store"
store.write_store_to_dir(store_path, "sequences/%s2/%s.seq")

# Open a local store
loaded_store = RefgetStore.open_local(store_path)
```

## Loading Remote Stores

You can open stores from remote URLs (HTTP/HTTPS). The small index files are
downloaded into a local cache directory; sequence data stays on the server
until you ask for it:

```python
# Open a remote store; index files are cached in cache_dir
cache_dir = "local_cache"
remote_url = "https://refget-server.example.com/hg38"
remote_store = RefgetStore.open_remote(cache_dir, remote_url)

# Collections and sequence metadata are loaded lazily
coll = remote_store.list_collections()["results"][0]
remote_store.load_collection(coll.digest)
record = remote_store.get_sequence_by_name(coll.digest, "chr1")
digest = record.metadata.sha512t24u

# Get a substring. This reads only the bytes for the region with an HTTP
# range request; nothing is written to the sequence cache.
substring = remote_store.get_substring(digest, 0, 1000)
print(f"First 1000 bases: {substring[:50]}...")

# For many reads on one sequence, download and cache it once
remote_store.load_sequence(digest)

# Iterate over the sequences the store knows about (metadata only)
for seq_meta in remote_store:
    print(f"{seq_meta.name}: {seq_meta.length} bp")
```

## Working with Collections

```python
# Get the first collection in the store
collection = store.list_collections()["results"][0]

# Get a sequence by collection and name
record = store.get_sequence_by_name(
    collection.digest,
    "chr1"
)

# Export entire collection to FASTA
store.export_fasta(
    collection.digest,
    "output.fa",
    sequence_names=None,  # None = all sequences
    line_width=80
)

# Export specific sequences from a collection
store.export_fasta(
    collection.digest,
    "chr1_and_chr2.fa",
    sequence_names=["chr1", "chr2"],
    line_width=80
)
```

## Extracting Regions from BED Files

```python
# Get sequences for regions defined in a BED file
retrieved_seqs = store.substrings_from_regions(
    collection.digest,
    "regions.bed"
)

for seq in retrieved_seqs:
    print(f"{seq.chrom_name}:{seq.start}-{seq.end} = {seq.sequence}")

# Export BED regions to a FASTA file
store.export_fasta_from_regions(
    collection.digest,
    "regions.bed",
    "output_regions.fa"
)
```

## Local HTTP Server Example

For testing remote loading locally, you can serve a store directory. The server
must support HTTP range requests, because `get_substring()` asks for just the
bytes it needs. Python's built-in `python -m http.server` does **not** support
them, and `get_substring()` fails against it with "Remote server did not honor
Range header". Use a range-capable static server instead, for example:

```bash
# In the directory containing your refget store
npx http-server -p 8200
```

Then connect to it:

```python
remote_store = RefgetStore.open_remote(
    "local_cache",
    "http://localhost:8200/my_refget_store"
)

# Use it like any other store
coll = remote_store.list_collections()["results"][0]
remote_store.load_collection(coll.digest)
seq_digest = remote_store.get_sequence_by_name(coll.digest, "chr1").metadata.sha512t24u
substring = remote_store.get_substring(seq_digest, 0, 100)
```

The [RefgetStore tutorial](refgetstore.ipynb) includes a small pure-Python
range-capable server if you prefer not to use Node.

## More Information

For a comprehensive tutorial with detailed examples, see [refgetstore.ipynb](refgetstore.ipynb).

For the full API documentation, visit the [gtars repository](https://github.com/databio/gtars).
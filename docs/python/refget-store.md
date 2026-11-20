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

## Packed layout (alpha)

A standard on-disk store writes one `.seq` file per sequence. A store with
hundreds of thousands of sequences then means hundreds of thousands of files,
which is slow on network filesystems. The packed layout writes sequences into
a few large pack files under `packs/` instead, and finds each sequence through
an index file, `sequences.pack.idx`.

```python
# Create a packed store
store = RefgetStore.on_disk("my_packed_store", packed=True)
store.add_sequence_collection_from_fasta("genome.fa")
print(store.layout)  # "packed"

# Removing collections leaves dead bytes in the packs; compact() reclaims them
store.compact()  # {"packs_before": ..., "packs_after": ..., "bytes_reclaimed": ..., "packs_pruned": 0}
```

`store.compact(prune_unreferenced=True)` also deletes pack files that no index
entry points to, such as packs left by a writer that crashed. Use it only when
no other process is writing to the store: a pack another writer has sealed but
not yet committed looks the same, and would be deleted too.

A store's layout is fixed when it is created. Opening an existing standard
store with `packed=True` raises an error; `open_local` opens either kind.

A standard store refuses an import that would leave it with more than
500,000 `.seq` files. To go past that, use the packed layout, or raise the
limit with `standalone_file_limit=` on `add_sequence_collection_from_fasta`
or `add_sequence_collections_from_fastas`. `force=True` does not skip the
limit.

Reads from a packed store can go through a block cache, which is off by
default. `store.set_pack_cache(128)` turns it on with 128 MB; `0` turns it
off. It helps clustered reads (for example a sorted BED) on a cold network
file system, and hurts scattered reads on stores larger than the cache.

The packed layout is alpha: its on-disk format may change before it is
declared stable.

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

# Export an entire collection to FASTA. The default label is the sequence
# digest (SQ....); label="{name}" writes chromosome names instead.
store.export_fasta("output.fa", collections=[collection.digest], label="{name}")

# Export specific sequences from a collection
store.export_fasta(
    "chr1_and_chr2.fa",
    collections=[collection.digest],
    names=["chr1", "chr2"],
    label="{name}",
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

# Export BED regions to a FASTA file. Headers look like >chr1:0-10 by
# default; pass label="{seq}" to use sequence digests instead.
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
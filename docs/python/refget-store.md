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

You can open stores from remote URLs (HTTP/HTTPS). The small metadata files
download at once, in parallel; sequence data stays on the server until you ask
for it, and then comes down over many connections at once:

```python
from gtars.refget import RefgetStore, default_cache_dir

remote_url = "https://refget-server.example.com/hg38"
remote_store = RefgetStore.open_remote(default_cache_dir(), remote_url)

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

# Download many sequences into the cache at once, without loading them
digests = [m.sha512t24u for m in remote_store.list_sequences()]
counters = remote_store.cache_sequences(digests, progress=True)
print(counters["bytes_downloaded"], counters["cache_hits"])

# Iterate over the sequences the store knows about (metadata only)
for seq_meta in remote_store:
    print(f"{seq_meta.name}: {seq_meta.length} bp")
```

### Download options

`open_remote(cache_path, remote_url, *, fetch_jobs=32)`. `fetch_jobs` is how
many HTTP requests run at once. S3 serves about 1-2 MB/s per connection, so
speed comes from many connections.

Every download works the same way:

- Sequence files over 8 MiB download as parallel 8 MiB ranges, so one big
  chromosome comes down over many connections. The server must support HTTP
  range requests; a server that answers a ranged request with the whole file
  (HTTP 200) is refused.
- A failed request is retried up to 5 times. Waits double each time (from
  0.5 s, at most 8 s, with jitter) and follow `Retry-After`.
- Every sequence that came over the network is checked against its digest
  before it is used or cached.
- HTTP 404, 410 and 403 count as "not found" (S3 answers 403 for a missing
  file when listing is denied).

`store.fetch_stats()` returns what the store has done so far: requests,
retries, bytes downloaded, cache hits and local reads.
`cache_sequences(digests, *, progress=False)` downloads many sequences into the
cache at once; `progress=True` prints a line to stderr every 5 s.

### The shared cache

`cache_path` is a cache root shared by every remote store. The default,
`default_cache_dir()`, is `REFGET_CACHE_DIR`, else `$XDG_CACHE_HOME/refget`,
else `~/.cache/refget`:

```
<cache root>/
  CACHEDIR.TAG                         tells backup tools to skip the folder
  sequences/<tag>/<d0d1>/<digest>.seq  downloaded sequences, shared by every remote
  remotes/<key>/                       one folder of metadata per remote URL
  tmp/                                 files being written, renamed into sequences/ when done
```

Sequence files are named by their digest and the way they are packed (`<tag>`
is the alphabet, or `raw` or `zstd`), so two remotes that both hold hg38 share
one copy, and a cached sequence never goes out of date. `load_sequence`,
`load_all_sequences` and `cache_sequences` write sequences to the cache.
Partial reads, streaming and (unless you pass `write_cache=True`)
`export_fasta` do not.

Metadata in `remotes/<key>/` stays valid while the matching digest in the
remote's `rgstore.json` is unchanged; `rgstore.json` itself is fetched on every
open. `cache_info(root)` summarizes the cache (sequence files and bytes, and the
bytes of each remote's metadata), and `clear_cache(root)` deletes everything in
it and returns the bytes freed.

A sequence reaches `sequences/` only after it is complete and checked, so an
interrupted download leaves nothing half-written there. Running the same
download again fetches only the sequences still missing.

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

For testing remote loading locally, you can serve a store directory with any
static server that supports HTTP range requests: sequence files over 8 MiB
download in ranges, and partial reads (`get_substring()`, `stream_sequence()`)
read only the bytes they need. Python's built-in `python -m http.server`
ignores ranges, so use a range-capable server, for example:

```bash
# In the directory containing your refget store
npx http-server -p 8200
```

Then connect to it:

```python
remote_store = RefgetStore.open_remote(
    "local_cache",  # a cache root
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
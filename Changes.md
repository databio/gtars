# Changelog

## Unreleased

- CI (untested): install ncbi-blast+ for BLAST tests; add a `blast` test job and feature checks.

### Breaking

- `gtars-refget`: `ComparisonDigests.a` is now `Option<String>` (was
  `String`), matching `b`. A side compared from a level-2 body has no known
  digest, so either digest can be `None` (JSON `null`).
- `gtars-refget`: remote stores share one download cache. `open_remote(cache_root,
  url, &FetchConfig)` takes a cache root (default `default_cache_root()`:
  `REFGET_CACHE_DIR`, else `$XDG_CACHE_HOME/refget`, else `~/.cache/refget`)
  instead of a per-remote folder, and a `FetchConfig`. Sequences live in
  `sequences/<tag>/<d0d1>/<digest>.seq`, shared by every remote; metadata in
  `remotes/<first 16 hex of sha256(url)>/`. Old per-origin caches are not read.
  Python and R keep their positional `(cache_path, remote_url)` arguments; the
  folder is now a cache root.
- `gtars-refget`: a cached metadata file is refetched whenever the manifest
  field that covers it is missing or changed (`collections/*.rgsi` follow both
  `collections_digest` and `sequences_digest`). `sequences.rgsi` and
  `collections.rgci` are checked against the manifest's sha256.
- `gtars-refget`: `ExportSummary` gains `fetch` (and is no longer `Eq`);
  `ExportOptions` gains `write_cache`. The `http` feature
  now implies `filesystem`.
- `gtars-refget`: FASTA export is one method with two independent choices.
  `export_fasta(selection, output, opts)` replaces `export_fasta` and
  `export_fasta_by_digests`. `ExportSelection` picks the records (one or more
  collections with optional names, a list of sequence digests, or `All`), and
  `ExportOptions.label` (a `LabelTemplate`) picks how they are labeled.
  The default label is `{seq}` (`SQ.<digest>`), so callers that want
  chromosome names must pass `{name}`. FASTA descriptions are no longer written.
  The same change applies to the Python (`export_fasta(output_path,
  collections=..., seq_digests=..., all=..., label=...)`), R (`export_fasta`)
  and CLI (`gtars refget export`) interfaces.
- `gtars-refget`: `export_fasta_from_regions` takes a label template and
  appends `:start-end` to each label. It now returns an `ExportSummary`.

### Added

- `gtars-refget`: new `compare_level2_json` compares two level-2 sequence
  collection bodies given as JSON, with no store or digest needed. Exposed as
  `compare_seqcols` in Python (`gtars.refget`) and `compareSeqcols` in WASM.
- `gtars-rgblast`: `BlastCommand`, one clap enum for every verb (`build`, `verify`,
  `view`, `nstats`, `search`, `report`, `origin`, `reps`) with
  `run_code(prog)`, and `cli::run_json(argv, progress)`, which runs a verb from an
  argument list and returns its result as JSON. The `rgmakeblastdb` and `rgblast`
  binaries mount `BlastCommand` and the bindings call `run_json`, so `gtars blast`,
  the binaries, Python and R share flags, conversions and exit codes (a closed
  pipe now exits 0).
- `gtars blast <verb>` in the `gtars` CLI (feature `blast`, on by default).
- Python `gtars.blast`: `build`, `verify`, `view`, `nstats`, `search`, `report`,
  `origin` and `reps`, with the verbs' long flags as keywords, results as dicts, a
  `progress` callback (called on the calling thread, Ctrl-C cancels),
  `REPORT_FIELDS`, and exceptions tagged with a `kind` (`FileNotFoundError`,
  `ValueError` and other built-ins, plus `EngineNotFoundError` and `EngineError`).
- R `blast_build()`, `blast_verify()`, `blast_view()`, `blast_nstats()`,
  `blast_search()`, `blast_report()`, `blast_origin()`,
  `blast_reps()` and `blast_report_fields()`, with errors as classed conditions
  (`gtars_rgblast_*`, including `_engine_not_found` and `_engine_failed`).
- `search` and `origin` run rgblastn in process on the store by default (no
  BLAST+ install); `--engine blastn --db BASE` runs BLAST+ `blastn` on a database
  you built.
- `gtars-rgblast`: `report::rank::SUMMARY_COLUMNS`, the TSV summary's column order,
  so bindings return rows in the same order.
- Docs: `blast.md`, the Python and R BLAST tutorials and the BLAST API reference.

- `gtars-refget`: parallel remote fetch. Remote stores download over a pool of
  keep-alive connections (`FetchConfig::jobs`, default 32), download sequence
  files over 8 MiB as parallel 8 MiB ranges (the server must support HTTP Range
  requests), retry transient failures with jittered backoff and `Retry-After`,
  treat 403 as "not found or access denied" (S3 answers 403 for a missing key),
  and hash every downloaded sequence against its digest before use. Metadata
  (collection files, alias tables, FHR sidecars) downloads in parallel batches.
- `gtars-refget`: `prefetch(digests, PrefetchOptions)` streams the stored bytes
  of many sequences in parallel, in input or any order; `cache_sequences` fills
  the cache. A sequence is published to the cache only when complete, so a rerun
  downloads only the sequences still missing. `load_collections`, `fetch_stats`,
  `cache_info`, `clear_cache(root)` (deletes everything in the cache) and
  `normalize_remote_url` are new; `RemoteMissing` lets callers tell a missing
  remote file from other errors.
- `gtars refget export --remote URL` with `--cache-dir`, `--resume` and
  `--fetch-jobs`. Exports read ahead in parallel and print a `network:` summary
  for remote stores.
- Python and R: `open_remote(..., fetch_jobs)`, `cache_sequences(digests,
  progress)`, `fetch_stats`, `load_collections`, `export_fasta(...,
  write_cache)`, and cache helpers (`default_cache_dir`, `cache_info`,
  `clear_cache(root)`; `refget_*` in R). The Python binding releases the GIL
  during downloads.
- `gtars-refget`: label templates with the fields `{seq}`, `{seq:NS}`,
  `{col}`, `{col:NS}` and `{name}`, for example `{col:ucsc}/{name}`. A missing
  alias falls back to the digest. A label repeated for the same sequence is
  written once; a label repeated for different sequences is an error raised
  before any output.
- `gtars-refget`: export every sequence in a store (`ExportSelection::All`,
  `--all`, `all=True`), and filter by molecule (`dna` or `protein`).
- Exports write to stdout when the output path is `-`. `gtars refget export`
  defaults to stdout, takes `--collection` more than once, and exits quietly
  when the reader closes the pipe. Python raises `BrokenPipeError` in that case.
- New crate `gtars-rgblast` with three binaries: `rgmakeblastdb` (`build`, `verify`,
  `view`, `nstats`) writes BLAST databases straight from a store;
  `rgblastn` searches a store in process (tasks `megablast`, `blastn` and
  `blastn-short`, `-outfmt 6`); `rgblast` (`search`, `report`, `origin`, `reps`)
  says which species and genomes a query comes from. Remote stores are pulled
  into the shared cache first (`--cache-dir`, `--resume`, `--fetch-jobs`).
- `gtars-refget`: `ReadonlyRefgetStore::sequence_relative_path` is now public
  (the store-relative path of a sequence's `.seq` record).

### Fixed

- `gtars-refget`: the streaming FASTA hasher (`FastaStreamHasher`) now reads
  multi-member gzip, including bgzip (BGZF) and concatenated `.gz` files. It
  used to decode only the first member.
- `gtars-refget`: `pull_fhr(KeepOurs)` loads sidecars that are already cached
  into memory. `pull_aliases` and `pull_fhr` count only a missing file as
  `not_found`; network errors propagate.
- `gtars-refget`: exports stream one sequence at a time and decode in chunks.
  Peak memory is one sequence's stored bytes plus a 1 MB chunk, instead of
  every exported sequence held in memory and never released. Exports leave
  the store's loaded state unchanged and never write a remote store's `.seq`
  files into the local cache.

## 0.10.0 - 2026-09-05

A `gtars-refget` store-safety and store-operations release. Component
versions: `gtars-refget` 0.10.0, `gtars-vrs` 0.7.1 (re-pinned to refget
0.10), `gtars-overlaprs` 0.6.1, `gtars` 0.10.0, `gtars-cli` 0.10.0,
`gtars-python` 0.10.0, `gtars-r` 0.10.0, `gtars-node` 0.8.0,
`gtars-wasm` 0.9.2.

### Breaking

- `gtars-refget`: renamed the two RAM-residency fields on `StoreStats` to say
  what they actually measure: `n_sequences_loaded` -> `n_sequences_in_memory`
  and `n_collections_loaded` -> `n_collections_in_memory`. The old names read
  as per-run ingest counters, which they never were: `n_sequences_in_memory`
  is structurally always 0 on a disk-backed store, and
  `n_collections_in_memory` counts whatever happens to be resident in RAM when
  `stats()` is called (it resets on process start and also counts collections
  merely touched by a read). The old keys are removed outright — no aliases —
  across the Rust API, the Python/Node/R bindings, and the `refget store stats`
  JSON. Callers reading the old keys now get a `KeyError`/`undefined` rather
  than a silently wrong number. In JS the fields are `nSequencesInMemory` /
  `nCollectionsInMemory`.
- `gtars-refget`: `add_sequence_collections_from_fastas` now returns an
  `ImportReport` instead of `Vec<(SequenceCollectionMetadata, bool)>`. The
  per-file results moved to `report.collections`. The singular
  `add_sequence_collection_from_fasta` is unchanged.
- `gtars-refget`: `get_sequence_by_name` now returns an owned `SequenceRecord`
  instead of `&SequenceRecord`. Rust callers holding a borrow must adjust;
  bytes are `Arc`-shared, so the change is not a copy.
- `gtars-refget`: `total_disk_size()` renamed to `logical_sequence_bytes()`
  (cheap, no I/O, sequence content only). `actual_disk_usage()` remains the
  exact on-disk footprint.
- `gtars-refget`: the `Tombstones` type is gone; it existed only to support
  the merge-at-commit scheme that this release replaces (see Fixed).
- `gtars-overlaprs`: removed the never-constructed `MultiChromOverlapperError`
  enum (and the crate's `thiserror` dependency). No code path ever produced
  it.

### Added

- `gtars-refget`: `ImportReport` provides genuine per-run ingest counters —
  `n_sequences_written`, `n_sequences_deduped`, and `n_collections_new` — which
  is what build reports actually want. `written + deduped` equals the number of
  sequence records seen across all processed files. Exposed in Python as the
  `ImportReport` class, and printed by `gtars refget build` as an
  "Ingested this run" summary line.
- `gtars-cli`: new `gtars refget export` subcommand
  (`-s <store> -o out.fa [-c <digest|ns:alias>] [--names ...] [-w N]`).
  `-c` accepts a collection alias as well as a digest and may be omitted when
  the store holds exactly one collection; `.gz` output is gzip-compressed;
  only the requested collection (and names) is loaded, never the whole store.
- `gtars-refget`: unwrapped FASTA export. `line_width` of `Some(0)` (CLI
  `-w 0`) writes each sequence on one line for k-mer tools such as GGCAT and
  SSHash; `None` still wraps at 80. This also removes a latent `chunks(0)`
  panic.
- `gtars-refget`: `FastaImportOptions::collection_alias(namespace, alias)`
  registers a collection alias (e.g. `ucsc:hg38`) at import time. Exposed as
  `gtars refget build --collection-alias ns:alias` and in the Python import
  API. Opt-in; nothing is inferred from filenames.
- `gtars-refget`: `StoreMetadata.logical_sequence_bytes` is computed at
  index-write time and cached in `rgstore.json`, so a store's size can be read
  from the manifest without loading the sequence index. Surfaced via
  `store_metadata()` and `stats()` in Rust, Python, Node, and R.
- `gtars-refget`: Encoded-mode stores now accept sequences that are already
  packed to the alphabet's bit width at insert time (the pattern panget uses),
  in addition to ASCII.

### Fixed

- `gtars-refget`: **concurrent writers no longer lose data.** `RefgetStore`
  took no lock and rewrote its index files from the snapshot loaded at open,
  so two writers were last-writer-wins; on 2026-07-23 four concurrent genome
  inits silently dropped a collection from both indexes. Writes now take an
  exclusive `<store>/.rgstore.lock` (an `O_EXCL` lockfile with heartbeat and
  steal-by-rename, not `flock`, which degrades silently on Lustre), and every
  index, manifest, and alias file is written via write-temp / fsync / rename /
  fsync-parent. No on-disk format change; existing stores and S3 copies stay
  readable.
- `gtars-refget`: **commit writes only this handle's changes.** A commit used
  to write back the handle's entire in-memory state, so a handle opened before
  another process's `remove_collection` resurrected the removed rows after the
  `.seq` files were already unlinked, leaving dangling references. Commit now
  takes the lock, re-reads the index from disk, applies only this handle's own
  additions and removals, and writes. `remove_collection` holds the lock
  across the orphan scan, index write, and unlinks, and writes the index
  before unlinking so a crash leaves recoverable garbage rather than dangling
  rows. Regression tests cover cross-process removal, concurrent adds, and
  torn-read safety.
- `gtars-refget`: orphan-sequence GC is linear instead of quadratic.
  `remove_collection(.., remove_orphan_sequences=true)` and
  `remove_orphan_seq_files` rescanned the whole `md5_lookup` once per orphan;
  on a 15M-sequence store, removing a 736k-orphan collection ran 78 minutes
  and deleted nothing. Now one scan per call. Shard-directory `rmdir` is also
  hoisted out of the per-file loop, halving syscalls on network filesystems.
- `gtars-refget`: a conflicting `--collection-alias` is now rejected before
  the collection is committed. Previously the collection was inserted and its
  `.rgsi` written, then the import failed, leaving an orphaned per-collection
  index file.
- `gtars-refget`: `export_fasta` and `get_collection` now use the collection's
  own sequence names and descriptions for FASTA headers instead of the
  first-imported label of a shared sequence (#270).
- `gtars-refget`: `get_sequence_by_name` now returns the collection's own
  name and description for the sequence, not the first-imported label of a
  sequence shared across collections (same class as #270).
- `gtars-refget`: read-path syscalls trimmed. Per-store open-fd cache is
  sized from the soft `RLIMIT_NOFILE` (bounded), and the per-read `exists()`
  pre-check is replaced by inspecting the open result, removing one `statx`
  per read.
- `gtars-refget`: the FHR metadata sidecar is not rewritten when its content
  is unchanged.
- `gtars-wasm`: `encodedByteRange` return type in the JS typings corrected.

### Internal

- Workspace dev profile uses `debug = "line-tables-only"`: panics and test
  backtraces keep file/line, `target/` shrinks substantially.

## 0.9.0 - 2026-06-12

The gtars 0.9.0 release. Highlights:

- New `gtars-vrs` crate: GA4GH VRS allele identifiers from VCF/HGVS, with an
  HGVS parser/AST, allele normalization, and transcript-anchored mapping.
- New Node.js bindings (`gtars-node`) and R bindings (`gtars-r`).
- Major `gtars-refget` overhaul: on-disk sequence store plus a new binary
  transcript store and coordinate mapper.
- `gtars-genomicdist`: binary FASTA (`.fab`) genome format with zero-copy
  mmap access; new stranded region-set operations.
- New BAM QC tooling in `gtars-uniwig`; overlap-engine rewrite in
  `gtars-overlaprs`; expanded WASM and Python bindings.
- A batch of correctness-audit fixes across VRS, refget, bamqc, and the
  Python HGVS AST.

See the pull request and per-crate changelogs for full details.

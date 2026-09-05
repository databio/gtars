# Changelog

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
- `gtars-refget`: import-time aliases (from `--collection-alias` and from
  FASTA-header namespaces) are published in the same commit as the indexes
  instead of one TSV rewrite per alias, so a reader can no longer see an
  alias whose collection is not yet in `collections.rgci`, and a failed import
  leaves no dangling alias. `FastaImportOptions::force` now also overrides an
  alias another writer published under a different digest, as documented.
- `gtars-refget`: stale-lock breaking is serialized through an `O_EXCL`
  `.rgstore.lock.break` marker and re-verifies the exact lock instance before
  renaming it away; two contenders judging the same lock stale could
  previously both acquire.
- `gtars-refget`: `.seq` payloads are published by temp-file + rename rather
  than an in-place `File::create`, so a second process writing the same
  digest can no longer truncate a live file or expose a partial one.
- `gtars-refget`: cleanup after a failed import unlinks only the `.seq` files
  that import staged. On a reopened store it previously classified every
  pre-existing sequence as an orphan and deleted them.
- `gtars-refget`: `remove_collection` also removes aliases for the collection
  that another process published after this handle opened, and a commit
  refuses to publish a collection whose sequences were removed from the index
  by a concurrent orphan GC while the import ran.
- `gtars-refget`: `fhr_digest` ignores `.rgstore.tmp.*` files in `fhr/`.
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

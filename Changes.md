# Changelog

## Unreleased

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

### Added

- `gtars-refget`: `ImportReport` provides genuine per-run ingest counters —
  `n_sequences_written`, `n_sequences_deduped`, and `n_collections_new` — which
  is what build reports actually want. `written + deduped` equals the number of
  sequence records seen across all processed files. Exposed in Python as the
  `ImportReport` class, and printed by `gtars refget build` as an
  "Ingested this run" summary line.

### Fixed

- `gtars-refget`: `export_fasta` and `get_collection` now use the collection's
  own sequence names and descriptions for FASTA headers instead of the
  first-imported label of a shared sequence (#270).

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

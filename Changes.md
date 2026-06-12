# Changelog

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

# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### RefgetStore (`gtars-refget`)
- Fix residues that did not round-trip through `StorageMode::Encoded`:
  - `D` and `H` in the IUPAC alphabet decoded as `H` and `V`. Fixed in the decoder, so existing stores read correctly without re-importing.
  - `U` in the IUPAC alphabet was stored as `T`. `U` now has its own code.
  - `U` (selenocysteine) and `B` in proteins were stored as `A`. The protein alphabet now holds `U`, `B`, `Z`, `O` and `J`.
- Choose the smallest alphabet that holds every character in a sequence, instead of assuming each alphabet contains the one before it
- Encoding a character that the chosen alphabet cannot hold is now an error instead of silently storing it as another residue (`encode_sequence`, `SequenceEncoder::finalize` and `set_encoding_mode` now return `Result`)
- The 3-bit alphabet stores `U` in place of `X`, so RNA fits in 3 bits. Sequences with `X` now use the protein alphabet.
- Add `verify` to check that every stored sequence still hashes to its digest (Rust, CLI `gtars refget verify`, Python `store.verify()`, R `verify_store()`)
- Faster alphabet detection and encoding

Stores built with 0.11.0 or earlier that contain `U` or `B` in proteins, or `U` in nucleotide sequences, must be re-imported. Run `gtars refget verify <store>` to list affected sequences.

## [0.11.0] -- 2026-09-22

### RefgetStore (`gtars-refget`)
- Add `StorageMode::Zstd` for compressed novel/unplaced sequence storage
- Lazy-load sequence index on first read (`list_sequences()` now returns `Result`)
- Add content-based sequence name matching
- Add resident overlay for RAM-cached sequence reads
- Detect gzip by magic bytes instead of file extension
- Fix remote listing, export, and WASM build
- Fix encoded sequence ingestion

### Python bindings (`gtars` PyPI package)
- Add RefgetStore docstrings visible in Python (fixes #272)
- Polish refget Python docstrings

## [0.10.0] -- 2026-09-05
- refget: writers now take an exclusive store lock and commit only their own changes; concurrent imports no longer lose collections, and a stale handle no longer resurrects removed rows
- refget: all index, manifest, alias, and sequence files are written atomically (write-temp then rename)
- refget: fix cleanup after a failed import deleting pre-existing sequences on a reopened store
- refget: orphan-sequence GC is linear instead of quadratic; large removals finish in seconds instead of hours
- refget: `export_fasta` and `get_sequence_by_name` use the collection's own names and descriptions [#270](https://github.com/databio/gtars/issues/270)
- refget: add `ImportReport` with per-run counters (`n_sequences_written`, `n_sequences_deduped`, `n_collections_new`); `add_sequence_collections_from_fastas` now returns it
- refget: add import-time collection aliases (`--collection-alias ns:alias`, `FastaImportOptions::collection_alias`, Python `collection_alias=`)
- refget: add `logical_sequence_bytes` to the store manifest; rename `total_disk_size()` to `logical_sequence_bytes()`
- refget: rename `StoreStats` fields `n_sequences_loaded`/`n_collections_loaded` to `n_sequences_in_memory`/`n_collections_in_memory` across Rust, Python, R, and Node
- refget: encoded stores accept pre-packed sequences at insert time
- cli: add `gtars refget export` (`-c <digest|ns:alias>`, `-w 0` for unwrapped FASTA, `.gz` output)
- overlaprs: remove unused `MultiChromOverlapperError`
- wasm: fix `encodedByteRange` return type in typings

## [0.9.1] -- 2026-06-17
Point release of refget 0.9.1, python 0.9.2, R 0.9.2, node 0.7.1, wasm 0.9.1 [#268](https://github.com/databio/gtars/pull/268).
- refget: unify remote on-demand sequence retrieval into three flows: `get_substring` (byte-range read), `stream_sequence` (streaming reader), `load_sequence` (download and cache)
- refget: bulk substring reads on a remote-only sequence download and cache the `.seq` once instead of one request per range
- python, R, node: expose `stream_sequence`, `load_sequence`, and `load_all_sequences`; node no longer auto-downloads whole chromosomes
- wasm: add `RemoteRefgetStore` JS layer (HTTP byte-range + OPFS caching) via `@databio/gtars/remote`

## [0.9.0] -- 2026-06-15
- add `gtars-vrs` crate: GA4GH VRS allele identifiers from VCF/HGVS, with an HGVS parser, allele normalization, and transcript-anchored mapping
- add Node.js bindings (`gtars-node`)
- refget: on-disk sequence store overhaul; add binary transcript store and coordinate mapper
- genomicdist: add binary FASTA (`.fab`) format with zero-copy mmap access; add stranded region-set operations
- uniwig: add BAM QC tooling
- overlaprs: overlap-engine rewrite
- expanded WASM and Python bindings
- correctness fixes across vrs, refget, bamqc, and the Python HGVS AST

## [0.8.0] -- 2026-03-19
- genomicdist: port of genomic distributions [#231](https://github.com/databio/gtars/pull/231), bug fixes [#244](https://github.com/databio/gtars/pull/244)
- uniwig: add streaming implementation alongside the batch parallel one [#236](https://github.com/databio/gtars/pull/236)
- refget: RefgetStore updates [#237](https://github.com/databio/gtars/pull/237)
- add `gtars-lola` crate; IGD supports an in-memory database [#242](https://github.com/databio/gtars/pull/242)
- expose more functionality via Python [#241](https://github.com/databio/gtars/pull/241)

## [0.7.0] -- 2026-02-24
- update all crates to Rust edition 2024
- refget: add sequence and collection aliases, ancillary digests, seqcol comparison, and FHR metadata; sequences can exist independently of a collection
- refget: add streaming FASTA digest computation in WASM
- uniwig: fix off-by-one errors; add variable-step format
- refactor R bindings; RefgetStore uses `ExternalPtr` so R garbage collection works
- `RegionSet` returns an error instead of panicking on bad input

## [0.6.0] -- 2026-01-26
- refget: add FAI (FASTA index) computation and metadata storage
- refget: add remote-access RefgetStore with lazy loading and caching; `SequenceRecord` is now an enum (Stub vs Full)
- refget: add `on_disk()`, `in_memory()`, `load_remote()` constructors and export methods

## [0.5.3] -- 2026-01-21
- remove reqwest dependency from gtars-core
- merge the gdrs package into gtars
- genomicdist: add bed classifier

## [0.5.2] -- 2025-11-04
- switch overlap operations to `MultiChromOverlapper`
- genomicdist: more tests, fix region distribution plot, fix function naming in wasm bindings
- load/write functions accept `AsRef<Path>`

## [0.5.1] -- 2025-10-08
- introduce web assembly bindings
- introduce R bindings

## [0.5.0] -- 2025-09-25
- converted the codebase into a workspace with multiple crates

## [0.4.0] -- 2025-09-10
- added `get_max_end_per_chr` and `nucleotides_length` methods to `RegionSet`

## [0.3.0] -- 2025-07-30
- move digests functionality to refget
- add RefgetStore to refget and its associated python bindings
- integrate support for `bits` + backend types for tokenizers (AIList or BITS)
- reworked the tokenization CLI to support the new `bits` and `backend` options

## [0.2.5] -- 2025-04-06
- Rework tokenizer API to be more consistent with the HuggingFace tokenizers API.
- Updates to `RegionSet` to improve performance and usability.
- Added file_digest function to RegionSet struct
- Fixed reqwest error in R bindings
- Fixed [#107](https://github.com/databio/gtars/issues/107)

## [0.2.4] -- 2025-03-05
- Attempt to fix failing python bindings in CI linux [#104](https://github.com/databio/gtars/issues/104)

## [0.2.3] -- 2025-03-05
- Improved RegionSet, by adding a multiple new methods: `to_bed`, `to_bed_gz`, `to_bigbed`, `identifier()`, and others.
- Fixed allowed `fasta_digest` to accept `Path` or `bytes` [#93](https://github.com/databio/gtars/issues/93)

## [0.2.2] -- 2025-02-18
- fix [#90](https://github.com/databio/gtars/issues/90)
- fix [#89](https://github.com/databio/gtars/issues/89)


## [0.2.1] -- 2025-02-11 
- allow comments at the beginning of fragment files
- bump bigtools to 0.5.5, fixing [#74](https://github.com/databio/gtars/issues/74) and [#77](https://github.com/databio/gtars/issues/77)

## [0.2.0] -- 2025-01-13 
- add position shift workflow for bam to bw (was previously added for bam to bed)
- add scaling argument for bam to bw workflow [#53](https://github.com/databio/gtars/issues/53)
- fix accumulation issue for bam workflow [#56](https://github.com/databio/gtars/issues/56)
- fix wiggle file (core) beginning at 0 [#43](https://github.com/databio/gtars/issues/43)
- fix npy file (end) using start instead of end [#61](https://github.com/databio/gtars/issues/61)
- force zoom to 1 for bed/narrowPeak to bw [#34](https://github.com/databio/gtars/issues/34)
- fix IGD overlap issue [#45](https://github.com/databio/gtars/issues/45)
- add ga4gh refget digest functionality [#58](https://github.com/databio/gtars/issues/58)
- fix wig and npy inconsistency [#64](https://github.com/databio/gtars/issues/64)
- fix narrowPeak to bw zoom  [#34](https://github.com/databio/gtars/issues/34)
- fix bed to bw fileheader consistency issue  [#52](https://github.com/databio/gtars/issues/52)
- change npy metadata file structure [#65](https://github.com/databio/gtars/issues/65)

## [0.1.2]
- add position shift workflow for `bam` to `bw` (was previously added for `bam` to `bed`)
- add scaling argument for `bam` to `bw` workflow [#53](https://github.com/databio/gtars/issues/53)
- fix accumulation issue for `bam` workflow [#56](https://github.com/databio/gtars/issues/56)
- fix wiggle file (core) beginning at 0 [#43](https://github.com/databio/gtars/issues/43)
- fix npy file (end) using start instead of end [#61](https://github.com/databio/gtars/issues/61)
- force zoom to 1 for bed/narrowPeak to bw [#34](https://github.com/databio/gtars/issues/34)
- fix IGD overlap issue [#45](https://github.com/databio/gtars/issues/45)
- add ga4gh refget digest functionality [#58](https://github.com/databio/gtars/pull/58)

## [0.1.1] -- 2024-12-03
- hot fix for broken python bindings; remove IGD from the python bindings for now

## [0.1.0] -- 2024-12-03
- Rust implementation of `uniwig` that expands on the C++ version
  - Uniwig now accepts a single sorted  `.bed` file, `.narrowPeak` file, or `.bam` file.
  - Outputs now include  `.wig`, `.npy`, `.bedGraph`, and `.bw`
  - Accumulations can now be counted via `.narrowPeak` scoring
- Rust implementation of `igd` ported from the C version (experimental).
- Region scoring matrix calculation for region clustering
- Fragment file splitter for pseudobulking

## [0.0.15] -- 2024-07-29
-  added meta tokenization tools and a new `MetaTokenizer` struct that can be used to tokenize regions using the meta-token strategy.
-  added some annotations to the `pyo3` `#[pyclass]` and `#[pymethods]` attributes to make the python bindings more readable.

## [0.0.14] -- 2024-06-11
- renamed repository to `gtars` to better reflect the project's goals.

## [0.0.13] -- 2024-06-03
- implemented a fragment file tokenizer that will generate `.gtok` files directly from `fragments.tsv.gz` files.
- fix an off-by-one error in the `region-to-id` maps in the `Universe` structs. This was leading to critical bugs in our models.

## [0.0.12] -- 2024-05-28
- optimize creation of `PyRegionSet` to reduce expensive cloning of `Universe` structs.

## [0.0.11] -- 2024-05-22
- redesigned API for the tokenizers to better emulate the huggingface tokenizers API.
- implemented new traits for tokenizers to allow for more flexibility when creating new tokenizers.
- bumped the version `pyo3` to `0.21.0`
- added `rust-numpy` dependency to the python bindings for exporting tokenized regions as numpy arrays.
- overall stability improvements to the tokenizers and the python bindings.

## [0.0.10] -- 2024-01-24
- update file format specifications

## [0.0.9] -- 2024-01-22
- start working on the concept of a `.gtok` file-format to store tokenized regions
- added basic readers and writers for this format

## [0.0.8] -- 2024-01-17
- add a new `ids_as_strs` getter to the `TokenizedRegionSet` struct so that we can get the ids as strings quickly, this is meant mostly for interface with geniml.

## [0.0.7] -- 2023-11-30
- move things around based on rust club feedback

## [0.0.6] -- 2024-02-20
- update python bindings to support the module/submodule structure (https://github.com/PyO3/pyo3/issues/759#issuecomment-1828431711)
- change name of some submodules
- remove `consts` submodule, just add to base
- expose a `__version__` attribute in the python bindings

## [0.0.5] -- 2024-02-19
- add many "core utils"
- move `gtokenizers` into this package inside `gtars::tokenizers`
- create `tokenize` cli
- add tests for core utils and tokenizers
- RegionSet is now backed by a polars DataFrame
- new python bindings for core utils and tokenizers

## [0.0.4] -- 2023-11-06
- add type annotations to the python bindings

## [0.0.3] -- 2023-11-06
- work on python bindings initialization

## [0.0.2] -- 2023-09-20
- prepare for first release

## [0.0.1] -- 2023-08-15
- initial setup of repository
- two main wrappers: 1) wrapper binary crate, and 2) wrapper library crate
- `gtars` can be used as a library crate. or as a command line tool
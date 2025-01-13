# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] 
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

## [0.1.1]
- hot fix for broken python bindings; remove IGD from the python bindings for now

## [0.1.0]
- Rust implementation of `uniwig` that expands on the C++ version
  - Uniwig now accepts a single sorted  `.bed` file, `.narrowPeak` file, or `.bam` file.
  - Outputs now include  `.wig`, `.npy`, `.bedGraph`, and `.bw`
  - Accumulations can now be counted via `.narrowPeak` scoring
- Rust implementation of `igd` ported from the C version (experimental).
- Region scoring matrix calculation for region clustering
- Fragment file splitter for pseudobulking

## [0.0.15]
-  added meta tokenization tools and a new `MetaTokenizer` struct that can be used to tokenize regions using the meta-token strategy.
-  added some annotations to the `pyo3` `#[pyclass]` and `#[pymethods]` attributes to make the python bindings more readable.

## [0.0.14]
- renamed repository to `gtars` to better reflect the project's goals.

## [0.0.13]
- implemented a fragment file tokenizer that will generate `.gtok` files directly from `fragments.tsv.gz` files.
- fix an off-by-one error in the `region-to-id` maps in the `Universe` structs. This was leading to critical bugs in our models.

## [0.0.12]
- optimize creation of `PyRegionSet` to reduce expensive cloning of `Universe` structs.

## [0.0.11]
- redesigned API for the tokenizers to better emulate the huggingface tokenizers API.
- implemented new traits for tokenizers to allow for more flexibility when creating new tokenizers.
- bumped the version `pyo3` to `0.21.0`
- added `rust-numpy` dependency to the python bindings for exporting tokenized regions as numpy arrays.
- overall stability improvements to the tokenizers and the python bindings.

## [0.0.10]
- update file format specifications

## [0.0.9]
- start working on the concept of a `.gtok` file-format to store tokenized regions
- added basic readers and writers for this format

## [0.0.8]
- add a new `ids_as_strs` getter to the `TokenizedRegionSet` struct so that we can get the ids as strings quickly, this is meant mostly for interface with geniml.

## [0.0.7]
- move things around based on rust club feedback

## [0.0.6]
- update python bindings to support the module/submodule structure (https://github.com/PyO3/pyo3/issues/759#issuecomment-1828431711)
- change name of some submodules
- remove `consts` submodule, just add to base
- expose a `__version__` attribute in the python bindings

## [0.0.5]
- add many "core utils"
- move `gtokenizers` into this package inside `gtars::tokenizers`
- create `tokenize` cli
- add tests for core utils and tokenizers
- RegionSet is now backed by a polars DataFrame
- new python bindings for core utils and tokenizers

## [0.0.4]
- add type annotations to the python bindings

## [0.0.3]
- work on python bindings initialization

## [0.0.2]
- prepare for first release

## [0.0.1]
- initial setup of repository
- two main wrappers: 1) wrapper binary crate, and 2) wrapper library crate
- `gtars` can be used as a library crate. or as a command line tool
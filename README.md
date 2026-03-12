[![codecov](https://codecov.io/gh/databio/gtars/branch/master/graph/badge.svg)](https://codecov.io/gh/databio/gtars)
[![crates.io](https://img.shields.io/crates/v/gtars?&logo=rust)](https://crates.io/crates/gtars)

<h1 align="center">
<img src="https://raw.githubusercontent.com/databio/gtars/master/gtars_logo.png" alt="gtars logo" height="100px">
</h1>


`gtars` is a rust project that provides a set of tools for working with genomic interval data. It includes modules for genomic distribution analysis (`genomicdist`), locus overlap enrichment analysis (`lola`), integrated genome database overlap queries (`igd`), sequence collection management (`refget`), and more. Its primary goal is to provide processors for our python package, [`geniml`](https:github.com/databio/geniml), a library for machine learning on genomic intervals. However, it can be used as a standalone library for working with genomic intervals as well. For more information, see the [public-facing documentation](https://docs.bedbase.org/gtars/) (under construction).


`gtars` provides these things:

1. A set of rust crates.
2. A command-line interface, written in rust.
3. A Python package that provides Python bindings to the rust crates.
4. An R package that provides R bindings to the rust crates.

## Repository organization (for developers)

This repository is a work in progress, and still in early development. This repo is organized like as a [workspace](https://doc.rust-lang.org/cargo/reference/workspaces.html). More specifically:

1. Each piece of core functionality is implemented as a separate rust crate and is mostly independent.
2. Common functionality (structs, traits, helpers) are stored in a `gtars-core` crate.
3. Python bindings are stored in `gtars-py`. They pull in the necessary rust crates and provide a Pythonic interface.
4. A command-line interface is implemented in the `gtars-cli` crate.

## Installation

To install `gtars`, you must first [install the rust toolchain](https://www.rust-lang.org/tools/install).

### Command-line interface

You may build the cli binary locally by navigating to `gtars-cli` and using `cargo build --release`. This will create a binary in `target/release/gtars` at the top level of the workspace. You can then add this to your path, or run it directly.

Alternatively, you can run `cargo install --path gtars-cli` from the top level of the workspace. This will install the binary to your cargo bin directory (usually `~/.cargo/bin`).

We feature-gate binary dependencies maximize compatibility and minimize install size. You can specify features during installation like so:

```
cargo install --path gtars-cli gtars-cli --features "uniwig tokenizers"
```

Finally, you can download precompiled binaries from the [releases page](https://github.com/databio/gtars/releases).

### Python bindings

You can install the Python bindings via pip. First, ensure you have a recent version of pip installed. Then run:

```bash
pip install gtars
```

Then, you can use it in Python like so:

```python
from gtars import __version__
print(__version__)
```

## Usage

`gtars` provides several useful tools. There are 3 ways to use `gtars`.

### 1. From Python

Using bindings, you can call some `gtars` functions from within Python.

### 2. From the CLI

To see the available tools you can use from the CLI run `gtars --help`. To see the help for a specific tool, run `gtars <tool> --help`.

Available subcommands:

| Subcommand | Description |
|---|---|
| `genomicdist` | Compute genomic distribution statistics for a BED file |
| `prep` | Pre-serialize GTF gene models or signal matrices to binary for fast loading |
| `ranges` | Interval set algebra operations on BED files (reduce, trim, promoters, setdiff, pintersect, concat, union, jaccard) |
| `consensus` | Compute consensus regions across multiple BED files |

#### Preparing reference files

Pre-compile reference files to binary for fast repeated loading. This is optional but recommended when running `genomicdist` repeatedly against the same references.

```bash
# Pre-compile a GTF gene model
gtars prep --gtf gencode.v47.annotation.gtf.gz

# Pre-compile an open signal matrix
gtars prep --signal-matrix openSignalMatrix_hg38.txt
```

Output defaults to the input path with `.bin` appended (stripping `.gz` first). Use `-o` to specify a custom output path.

#### Computing genomic distributions

```bash
gtars genomicdist \
  --bed query.bed \
  --gtf gencode.v47.annotation.gtf.bin \
  --tss tss.bed \
  --chrom-sizes hg38.chrom.sizes \
  --signal-matrix openSignalMatrix_hg38.txt.bin \
  --output result.json
```

All flags except `--bed` are optional. Omit any flag to skip that analysis:

| Flag | Required | Description |
|---|---|---|
| `--bed` | yes | Input BED file |
| `--gtf` | no | GTF/GTF.gz or pre-compiled `.bin` — enables partitions and TSS distances |
| `--tss` | no | TSS BED file — overrides GTF-derived TSS for distance calculation |
| `--chrom-sizes` | no | Chrom sizes file — enables expected partitions |
| `--signal-matrix` | no | Signal matrix TSV or pre-compiled `.bin` — enables open chromatin enrichment |
| `--bins` | no | Number of bins for region distribution (default: 250) |
| `--promoter-upstream` | no | Upstream distance from TSS for promoter regions (default: 200) |
| `--promoter-downstream` | no | Downstream distance from TSS for promoter regions (default: 2000) |
| `--output` | no | Output JSON path (default: stdout) |
| `--compact` | no | Compact JSON output (default: pretty-printed) |

### 3. As a rust library

You can link `gtars` as a library in your rust project. To do so, add the following to your `Cargo.toml` file:

```toml
[dependencies]
gtars = { git = "https://github.com/databio/gtars/gtars" }
```

We wall off crates using features, so you will need to enable the features you want. For example, to use the overlap tool:

```toml
[dependencies]
gtars = { git = "https://github.com/databio/gtars/gtars", features = ["overlaprs"] }
```

Then, in your rust code, you can use it like so:

```rust
use gtars::overlaprs::{ ... };
```
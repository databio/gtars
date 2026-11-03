<p align="center">
<h1><img align="center" src="img/gtars_logo.svg" class="img-header" height="100"></h1>
</p>


<p align="center">
<a href="https://pypi.org/project/gtars"><img src="https://img.shields.io/pypi/v/gtars" alt=""></a>
<a href="https://crates.io/crates/gtars"><img src="https://img.shields.io/crates/v/gtars?&logo=rust" alt="crates.io"></a>
<a href="https://github.com/databio/gtars"><img src="https://img.shields.io/badge/source-github-354a75?logo=github"></a>
</p>



## Introduction

`gtars` is a high-performance toolkit for *genomic tools and algorithms in Rust*. Built with Rust for speed and reliability, gtars provides core utilities for machine learning on genomic intervals for the [geniml](https://github.com/databio/geniml) Python package. It also provides lots of utility as a standalone library for alternative downstream use cases.

## Installation

### Rust Library

Gtars uses a feature-flag system to allow you to include only the modules you need. Add to your `Cargo.toml`:

```toml
[dependencies]
# Install specific features
gtars = { version = "0.5", features = ["overlaprs", "tokenizers"] }

# Or install from GitHub
gtars = { git = "https://github.com/databio/gtars", features = ["overlaprs", "tokenizers"] }
```

Modules:

- `core` - Core functionality and data structures
- `tokenizers` - Genomic region tokenizers
- `io` - I/O utilities
- `refget` - Reference sequence access
- `overlaprs` - Overlap operations
- `uniwig` - Coverage computation
- `igd` - Interval search
- `bbcache` - BED file caching
- `scoring` - Fragment scoring
- `fragsplit` - Fragment splitting


Example combinations:

```toml
# For machine learning tasks
gtars = { version = "0.5", features = ["tokenizers", "core"] }

# For genomic analysis
gtars = { version = "0.5", features = ["overlaprs", "uniwig", "scoring"] }

# For data access
gtars = { version = "0.5", features = ["refget", "bbcache", "io"] }
```

### Python Package

```bash
pip install gtars
```

See further documentation under Python bindings.

### Command-Line Interface

Install from source:
```bash
git clone https://github.com/databio/gtars
cd gtars
cargo install --path gtars-cli
```

Or download precompiled binaries from the [releases page](https://github.com/databio/gtars/releases).


## Development

Run tests with `cargo test` from the workspace root. Please see [CONTRIBUTING.md](https://github.com/databio/gtars/blob/master/CONTRIBUTING.md) for development guidelines.

## Module organization

`gtars` is organized into modules. The modules section gives an [overview of each module](modules.md).

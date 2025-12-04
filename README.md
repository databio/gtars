[![codecov](https://codecov.io/gh/databio/gtars/branch/master/graph/badge.svg)](https://codecov.io/gh/databio/gtars)
[![crates.io](https://img.shields.io/crates/v/gtars?&logo=rust)](https://crates.io/crates/gtars)

<h1 align="center">
<img src="docs/gtars_logo_new_with_words.png" alt="gtars logo" height="100px">
</h1>


`gtars` is a rust project that provides a set of tools for working with genomic interval data. Its primary goal is to provide processors for our python package, [`geniml`](https:github.com/databio/geniml), a library for machine learning on genomic intervals. However, it can be used as a standalone library for working with genomic intervals as well. For more information, see the [public-facing documentation](https://docs.bedbase.org/gtars/) (under construction).


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

`gtars` is available in multiple forms to support different environments. Choose the installation method that best fits your needs:

### Quick Install

#### Python Package (Recommended for most users)
```bash
pip install gtars
```

#### Command Line Interface (CLI)
```bash
cargo install gtars-cli --all-features
```

### Detailed Installation Options

#### Prerequisites
- For CLI installation: [Rust toolchain](https://www.rust-lang.org/tools/install)
- For Python: pip (comes with Python)

#### 1. Python Package
The easiest way to use gtars is through the Python package:

```bash
pip install gtars
```

Example usage:
```python
from gtars.tokenizers import Tokenizer
tokenizer = Tokenizer.from_pretrained("databio/atacformer-base-hg38")
```

#### 2. Command Line Interface (CLI)
Install the CLI tool using cargo:

```bash
# Basic installation
cargo install gtars-cli

# With specific features
cargo install gtars-cli --features "uniwig tokenizers"

# With all features
cargo install gtars-cli --all-features
```

Verify installation:
```bash
gtars --help
```

#### 3. Rust Library
Add to your `Cargo.toml`:

```toml
[dependencies]
gtars = { version = "0.5.0", features = ["uniwig", "tokenizers"] }
```

#### 4. JavaScript/WebAssembly
```bash
npm install @databio/gtars-js
```

Example usage:
```javascript
import { Overlapper } from '@databio/gtars-js';
const overlapper = new Overlapper(universe, 'ailist');
```

#### 5. Precompiled Binaries
Download ready-to-use binaries from the [releases page](https://github.com/databio/gtars/releases).

## Usage

`gtars` provides several useful tools accessible through different interfaces:

### 1. Python Package

After installation via pip, you can use gtars functions directly in Python:

```python
from gtars import __version__
from gtars.tokenizers import Tokenizer

print(__version__)
tokenizer = Tokenizer.from_pretrained("databio/atacformer-base-hg38")
```

### 2. Command Line Interface

Once installed, use the CLI to access various tools:

```bash
# View available tools
gtars --help

# Get help for a specific tool
gtars <tool> --help
```

### 3. Rust Library

Include gtars in your Rust project by adding to `Cargo.toml`:

```toml
[dependencies]
gtars = { version = "0.5.0", features = ["overlaprs"] }
```

Then use in your Rust code:

```rust
use gtars::overlaprs::{ ... };
```

### 4. JavaScript/WebAssembly

After npm installation, use in JavaScript/TypeScript:

```javascript
import { Overlapper } from '@databio/gtars-js';
// Your code here
```
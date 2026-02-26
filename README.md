[![codecov](https://codecov.io/gh/databio/gtars/branch/master/graph/badge.svg)](https://codecov.io/gh/databio/gtars)
[![crates.io](https://img.shields.io/crates/v/gtars?&logo=rust)](https://crates.io/crates/gtars)

<h1 align="center">
<img src="https://raw.githubusercontent.com/databio/gtars/master/gtars_logo.png" alt="gtars logo" height="100px">
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

To install `gtars`, you must first [install the rust toolchain](https://www.rust-lang.org/tools/install).

### Command-line interface

You may build the cli binary locally by navigating to `gtars-cli` and using `cargo build --release`. This will create a binary in `target/release/gtars` at the top level of the workspace. You can then add this to your path, or run it directly.

#### Python Package (Recommended for most users)
```bash
pip install gtars
```

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

Example usage:
```python
from gtars.tokenizers import Tokenizer
tokenizer = Tokenizer.from_pretrained("databio/atacformer-base-hg38")
```

### Dev Python bindings




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
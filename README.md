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

To install `gtars`, you must have the rust toolchain installed. You can install it by [following the instructions](https://www.rust-lang.org/tools/install).

You may build the cli binary locally by navigating to `gtars-cli` and using `cargo build --release`. This will create a binary in `target/release/gtars` at the top level of the workspace. You can then add this to your path, or run it directly.

## Usage

`gtars` provides several useful tools. There are 3 ways to use `gtars`. 

### 1. From Python

Using bindings, you can call some `gtars` functions from within Python.

### 2. From the CLI

To see the available tools you can use from the CLI run `gtars --help`. To see the help for a specific tool, run `gtars <tool> --help`.

### 3. As a rust library

You can link `gtars` as a library in your rust project. To do so, add the following to your `Cargo.toml` file:

```toml
[dependencies]
gtars = { git = "https://github.com/databio/gtars" }
```

## Testing

To run the tests, run `cargo test --workspace`.

### Refget tests

The default tests for this module are designed to run quickly on tiny fasta files.
To run the test on a full-scale fasta file, you can look at `test_loading_large_fasta_file`.
This is large test, which is ignored by default, so it doesn't run in the typical `cargo test`. 
To run just this large test on a fasta file, try something like this:

```
FASTA_PATH=tests/data/subset.fa.gz cargo test tests::test_loading_large_fasta_file -- --nocapture --ignored
FASTA_PATH=`refgenie seek test/fasta` cargo test tests::test_loading_large_fasta_file -- --nocapture --ignored
```

## Contributing

### New internal library crate tools

If you'd like to add a new tool, you can do so by creating a new crate at the `gtars` workspace root. This can be done via `cargo new --lib <tool_name>`.

### New public library crate tools

If you want this to be available to users of `gtars`, you can add it to the `gtars` library crate as well. To do so, add your new crate as a dependency in `Cargo.toml`:

```toml
[dependencies]
<tool_name> = { path = "../<tool_name>", optional = true }
```

Then create a new feature to `Cargo.toml`:

```toml
[features]
<tool_name> = ["dep:<tool_name>"]
```

And finally, re-export this tool in `gtars/src/lib.rs`:

```rust
#[cfg(feature = "<tool_name>")]
#[doc(inline)]
pub use <tool_name> as <tool_name>;

```

### New binary crate tools

Finally, if you want to have command-line functionality, you can add it to the `gtars-cli` binary crate. This requires ____ steps:

1. Create a new module inside `gtars-cli/src/main.rs`:

```rust
mod new_tool;
```

2. Use `clap` to define your command-line interface inside `gtars-cli/src/new_tool/cli.rs`:

```rust
use clap::{App, Arg};

pub fn make_new_tool_cli() -> App<'static> {
    App::new("new_tool")
        .about("Does something new")
        .arg(
            Arg::new("input")
                .about("Input file")
                .required(true)
                .index(1),
        )
}
```

3. Write your logic in a wrapper function. This will live inside the `handlers` module of `gtars-cli/src/new_tool`:

```rust
// top of file:
use tool_name::{ ... }

// inside the module:
pub fn new_tool_wrapper() -> Result<(), Box<dyn Error>> {
    // your logic here
}
```

4. Bring this into the `gtars-cli/src/main.rs` file:

```rust
mod new_tool;
```

Please make sure you update the changelog and bump the version number in `Cargo.toml` when you add a new tool.

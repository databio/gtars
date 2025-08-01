[![codecov](https://codecov.io/gh/databio/gtars/branch/master/graph/badge.svg)](https://codecov.io/gh/databio/gtars)
[![crates.io](https://img.shields.io/crates/v/gtars?&logo=rust)](https://crates.io/crates/gtars)

<h1 align="center">
<img src="gtars/docs/gtars_logo_new_with_words.png" alt="gtars logo" height="100px">
</h1>


`gtars` is a rust crate that provides a set of tools for working with genomic interval data. Its primary goal is to provide processors for our python package, [`geniml`](https:github.com/databio/geniml), a library for machine learning on genomic intervals. However, it can be used as a standalone library for working with genomic intervals as well. For more information, see the [public-facing documentation](https://docs.bedbase.org/gtars/) (under construction).


`gtars` provides these things:

1. A rust library crate.
2. A command-line interface, written in rust.
3. A Python package that provides Python bindings to the rust library.
4. An R package that provides R bindings to the rust library

## Repository organization (for developers)

This repository is a work in progress, and still in early development. This repo is organized like so:

1. The main gtars rust package (in subfolder `/gtars`), which contains two crates:
    A. A rust library crate (`/gtars/lib.rs`) that provides functions, traits, and structs for working with genomic interval data.
    B. A rust binary crate (in `/gtars/main.rs`), a small, wrapper command-line interface for the library crate.
2. Python bindings (in `/bindings/python`), which consists of a rust package with a library crate (no binary crate) and Python package.
3. R bindings (in `/bindinds/r`), which consists of an R package.

## Installation

To install `gtars`, you must have the rust toolchain installed. You can install it by [following the instructions](https://www.rust-lang.org/tools/install).

You may build the binary locally using `cargo build --release`. This will create a binary in `target/release/gtars`. You can then add this to your path, or run it directly.

## Usage

`gtars` provides several useful tools. There are 3 ways to use `gtars`. 

### 1. From R/Python

Using bindings, you can call some `gtars` functions from within R or Python.

### 2. From the CLI

To see the available tools you can use from the CLI run `gtars --help`. To see the help for a specific tool, run `gtars <tool> --help`.

### 3. As a rust library

You can link `gtars` as a library in your rust project. To do so, add the following to your `Cargo.toml` file:

```toml
[dependencies]
gtars = { git = "https://github.com/databio/gtars" }
```

## Testing

To run the tests, run `cargo test`.

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

If you'd like to add a new tool, you can do so by creating a new module within the src folder.

### New public library crate tools

If you want this to be available to users of `gtars`, you can add it to the `gtars` library crate as well. To do so, add the following to `src/lib.rs`:
```rust
pub mod <tool_name>;
```

### New binary crate tools

Finally, if you want to have command-line functionality, you can add it to the `gtars` binary crate. This requires two steps:

1. Create a new `cli` using `clap` inside the `interfaces` module of `src/cli.rs`:

```rust
pub fn make_new_tool_cli() -> Command {

}
```

2. Write your logic in a wrapper function. This will live inside the `functions` module of `src/cli.rs`:

```rust
// top of file:
use tool_name::{ ... }

// inside the module:
pub fn new_tool_wrapper() -> Result<(), Box<dyn Error>> {
    // your logic here
}
```

Please make sure you update the changelog and bump the version number in `Cargo.toml` when you add a new tool.

### VSCode users

If you are using VSCode, make sure you link to the `Cargo.toml` inside the `.vscode` folder, so that `rust-analyzer` can link it all together:
```json
{
    "rust-analyzer.linkedProjects": [
        "./vocab/Cargo.toml",
        "./Cargo.toml"
        "./new-tool/Cargo.toml"
    ]
}
```

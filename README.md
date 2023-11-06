# genimtools
`genimtools` is a rust crate that provides a set of tools for working with genomic interval data. Its primary goal is to provide processors for our python package, [`geniml`](https:github.com/databio/geniml), a libary for machine learning on genomic intervals. However, it can be used as a standalone library for working with genomic intervals as well.

## Installation
To install `genimtools`, you must have the rust toolchain installed. You can install it by following the instructions [here](https://www.rust-lang.org/tools/install).

## Usage
`genimtools` is very early in development, and as such, it does not have a lot of functionality yet. However, it does have a few useful tools. To see the available tools, run `genimtools --help`. To see the help for a specific tool, run `genimtools <tool> --help`.

Alternatively, you can link `genimtools` as a library in your rust project. To do so, add the following to your `Cargo.toml` file:
```toml
[dependencies]
genimtools = { git = "https://github.com/databio/genimtools" }
```

## Testing
To run the tests, run `cargo test`.


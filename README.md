<h1 align="center">
<img src="genimtools/docs/logo.svg" alt="genimtools logo" height="100px">
</h1>

`genimtools` is a rust crate that provides a set of tools for working with genomic interval data. Its primary goal is to provide processors for our python package, [`geniml`](https:github.com/databio/geniml), a libary for machine learning on genomic intervals. However, it can be used as a standalone library for working with genomic intervals as well.

`genimtools` is three things:

1. A rust crate that provides functions, traits, and structs for working with genomic interval data.
2. A command-line tool that acts as a command-line interface for the rust crate.
3. A python package that provides python bindings for the rust crate, so that it can be used within python.

This repository is a work in progress, and still in early development. As such, it is not yet available on `crates.io`. However, you can still use it by building it locally.

## Installation
To install `genimtools`, you must have the rust toolchain installed. You can install it by following the instructions [here](https://www.rust-lang.org/tools/install).

You may build the binary locally using `cargo build --release`. This will create a binary in `target/release/genimtools`. You can then add this to your path, or run it directly.

## Usage
`genimtools` is very early in development, and as such, it does not have a lot of functionality yet. However, it does have a few useful tools. To see the available tools, run `genimtools --help`. To see the help for a specific tool, run `genimtools <tool> --help`.

Alternatively, you can link `genimtools` as a library in your rust project. To do so, add the following to your `Cargo.toml` file:
```toml
[dependencies]
genimtools = { git = "https://github.com/databio/genimtools" }
```

## Testing
To run the tests, run `cargo test`.

## Contributing
### New internal library crate tools
If you'd like to add a new tool, you can do so by creating a new module within the src folder.

### New public library crate tools
If you want this to be available to users of `genimtools`, you can add it to the `genimtools` library crate as well. To do so, add the following to `src/lib.rs`:
```rust
pub mod <tool_name>;
```

### New binary crate tools
Finally, if you want to have command-line functionality, you can add it to the `genimtools` binary crate. This requires two steps:
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
# gtars

This is an R package that wraps the `gtars` Rust crate so you can call gtars code from R.

## Brief instructions

To install the development version, you'll have to build it locally. First, make sure src/rust/Cargo.toml includes the following dependencies:

``` 
gtars-core = { git = "https://github.com/databio/gtars", branch = "r-workspace" }
io = { git = "https://github.com/databio/gtars", branch = "r-workspace" }
igd = { git = "https://github.com/databio/gtars", branch = "r-workspace" }
refget = { git = "https://github.com/databio/gtars", branch = "r-workspace" }
tokenizers = { git = "https://github.com/databio/gtars", branch = "r-workspace", features = ["huggingface"] }
```

and comments these dependencies out:
```
# gtars-core = { path = "../../../gtars-core" }
# io = { path = "../../../io" }
# igd = { path = "../../../igd" }
# refget = { path = "../../../refget" }
# tokenizers = { path = "../../../tokenizers", features = ["huggingface"] }
```

If you are updating these gtars dependencies and would like to test them with R bindings, you will need to push the changes to the appropriate branch first, or stick to the local dependency paths and use `rextendr::clean()` and `rextendr::document()` rather than building the package source. You can build R bindings like this:

```console
cd bindings
R CMD build r
```

Install the package that was just built:

```console
R CMD INSTALL gtars_0.5.0.tar.gz
```

Some demo code is located at tests/refget_example.R.


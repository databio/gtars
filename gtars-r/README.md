# gtars

This is an R package that wraps the `gtars` Rust crate so you can call gtars code from R.

## Brief instructions

To install the development version, use the `remotes` package:

``` R
install.packages("remotes")
remotes::install_github("databio/gtars", ref = "dev", subdir = "gtars-r")
```

You can also build it locally. First, make sure src/rust/Cargo.toml includes the following dependencies:

``` 
gtars-core = { git = "https://github.com/databio/gtars", branch = "dev" }
io = { git = "https://github.com/databio/gtars", branch = "dev" }
igd = { git = "https://github.com/databio/gtars", branch = "dev" }
refget = { git = "https://github.com/databio/gtars", branch = "dev" }
tokenizers = { git = "https://github.com/databio/gtars", branch = "dev", features = ["huggingface"] }
```

and comments these dependencies out:
```
# gtars-core = { path = "../../../gtars-core" }
# io = { path = "../../../io" }
# igd = { path = "../../../igd" }
# refget = { path = "../../../refget" }
# tokenizers = { path = "../../../tokenizers", features = ["huggingface"] }
```

If you are updating these gtars dependencies and would like to test them with R bindings, you will need to push the changes to the appropriate branch first, or stick to the local dependency paths and use `rextendr::clean()` and `rextendr::document()` rather than building the package source. You can install the R package from source:
```console
R CMD INSTALL gtars-r
```

Some demo code is located at tests/refget_example.R.


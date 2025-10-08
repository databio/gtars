# gtars

This is an R package that wraps the `gtars` Rust crate so you can call gtars code from R.

## Brief instructions

To install the development version from GitHub, use the `remotes` package:

``` R
install.packages("remotes")
remotes::install_github("databio/gtars", ref = "dev", subdir = "gtars-r")
```

You can install the R package locally from the repo directory:

``` console
R CMD INSTALL gtars-r
```

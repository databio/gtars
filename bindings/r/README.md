# gtars

This is an R package that wraps the `gtars` Rust crate so you can call gtars code from R.

## Brief instructions

To install the development version, you'll have to build it locally. Build R bindings like this:

```console
cd bindings
R CMD build r
```

Then install the package that was just built:

```console
R CMD INSTALL gtars_0.0.1.tar.gz
```
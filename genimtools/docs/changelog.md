# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.0.7]
- move things around based on rust club feedback

## [0.0.6]
- update python bindings to support the module/submodule structure (https://github.com/PyO3/pyo3/issues/759#issuecomment-1828431711)
- change name of some submodules
- remove `consts` submodule, just add to base
- expose a `__version__` attribute in the python bindings

## [0.0.5]
- add many "core utils"
- move `gtokenizers` into this package inside `genimtools::tokenizers`
- create `tokenize` cli
- add tests for core utils and tokenizers
- RegionSet is now backed by a polars DataFrame
- new python bindings for core utils and tokenizers

## [0.0.4]
- add type annotations to the python bindings

## [0.0.3]
- work on python bindings initialization

## [0.0.2]
- prepare for first release

## [0.0.1]
- initial setup of repository
- two main wrappers: 1) wrapper binary crate, and 2) wrapper library crate
- `genimtools` can be used as a library crate. or as a command line tool
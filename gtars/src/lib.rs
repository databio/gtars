//! # gtars: *<small>Performance-critical tools to manipulate, analyze, and process genomic interval data. </small>*
//!
//! `gtars` is a rust crate that provides a set of tools for working with genomic interval data. Its primary goal is to provide
//! processors for our python package, [`geniml`](https:github.com/databio/geniml), a library for machine learning on genomic intervals.
//! However, it can be used as a standalone library for working with genomic intervals as well.
//!
pub mod ailist;
pub mod common;
pub mod io;
pub mod tokenizers;
pub mod uniwig;

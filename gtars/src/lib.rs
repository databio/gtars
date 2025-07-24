//! # gtars: *<small>Performance-critical tools to manipulate, analyze, and process genomic interval data. </small>*
//!
//! `gtars` is a rust crate that provides a set of tools for working with genomic interval data. Its primary goal is to provide
//! processors for our python package, [`geniml`](https:github.com/databio/geniml), a library for machine learning on genomic intervals.
//! However, it can be used as a standalone library for working with genomic intervals as well.
//!
//! There are several modules in this crate. The most comprehensive is the [tokenizers] modules which houses genomic region tokenizers
//! for use as pre-processors to machine learning pipelines.
//!
pub mod bbcache;
pub mod common;
pub mod fragsplit;
pub mod igd;
pub mod io;
pub mod refget;
pub mod scoring;
pub mod overlap;
pub mod tokenizers;
pub mod uniwig;

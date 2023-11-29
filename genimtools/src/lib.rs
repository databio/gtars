//! # Genimtools: *<small>Performance-critical tools to manipulate, analyze, and process genomic interval data. </small>*
//!
//! `genimtools` is a rust crate that provides a set of tools for working with genomic interval data. Its primary goal is to provide 
//! processors for our python package, [`geniml`](https:github.com/databio/geniml), a libary for machine learning on genomic intervals. 
//! However, it can be used as a standalone library for working with genomic intervals as well.
//!
pub mod common;
pub mod tokenizers;
pub mod uniwig;
pub mod vocab;

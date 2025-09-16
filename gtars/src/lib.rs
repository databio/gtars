//! # gtars: *<small>Performance-critical tools to manipulate, analyze, and process genomic interval data. </small>*
//!
//! `gtars` is a rust crate that provides a set of tools for working with genomic interval data. Its primary goal is to provide
//! processors for our python package, [`geniml`](https:github.com/databio/geniml), a library for machine learning on genomic intervals.
//! However, it can be used as a standalone library for working with genomic intervals as well.
//!
//! There are several modules in this crate. The most comprehensive is the [tokenizers] modules which houses genomic region tokenizers
//! for use as pre-processors to machine learning pipelines.
//!
#[cfg(feature = "core")]
#[doc(inline)]
pub use gtars_core as core;

#[cfg(feature = "tokenizers")]
#[doc(inline)]
pub use tokenizers;

#[cfg(feature = "io")]
#[doc(inline)]
pub use io;

#[cfg(feature = "refget")]
#[doc(inline)]
pub use refget;

#[cfg(feature = "overlaprs")]
#[doc(inline)]
pub use overlaprs;

#[cfg(feature = "uniwig")]
#[doc(inline)]
pub use uniwig;

#[cfg(feature = "igd")]
#[doc(inline)]
pub use igd;

#[cfg(feature = "bbcache")]
#[doc(inline)]
pub use bbcache;

#[cfg(feature = "scoring")]
#[doc(inline)]
pub use scoring;

#[cfg(feature = "fragsplit")]
#[doc(inline)]
pub use fragsplit;

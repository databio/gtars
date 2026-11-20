//! The rgblastn engine: a port of NCBI BLAST+ 2.17.0's nucleotide search core
//! (`algo/blast/core`) that searches a `SubjectPlan` through `subject::StoreSubjects` (the
//! plan over the store's mapped files). It never imports `blastdb`. Only `search` (the entry
//! point and its types), `options` and `encoding` are public.
//!
//! Each module header names the C source file(s) and functions it ports.
//!
//! The ported functions keep the C's loop shapes, argument lists and literal constants so a
//! reader can diff them against the source line by line; the lints below would make them
//! drift from it.
#![allow(
    clippy::needless_range_loop,
    clippy::too_many_arguments,
    clippy::collapsible_if,
    clippy::collapsible_else_if,
    clippy::manual_range_contains,
    clippy::if_same_then_else,
    clippy::approx_constant,
    clippy::mut_range_bound,
    clippy::comparison_chain,
    clippy::manual_memcpy,
    clippy::explicit_counter_loop,
    clippy::neg_cmp_op_on_partial_ord,
    clippy::nonminimal_bool,
    clippy::identity_op,
    clippy::erasing_op,
    clippy::precedence,
    clippy::bool_to_int_with_if,
    clippy::needless_late_init,
    clippy::manual_swap,
    clippy::int_plus_one,
    clippy::assign_op_pattern,
    clippy::unnecessary_cast,
    clippy::type_complexity
)]

pub(crate) mod dust;
pub mod encoding;
pub(crate) mod gapalign;
pub(crate) mod greedy;
pub(crate) mod hsp;
pub(crate) mod itree;
pub(crate) mod lookup;
pub mod options;
pub(crate) mod prelim;
pub(crate) mod query;
pub(crate) mod scan;
pub mod search;
pub(crate) mod stat;
pub(crate) mod subject;
pub(crate) mod traceback;
pub(crate) mod ungapped;

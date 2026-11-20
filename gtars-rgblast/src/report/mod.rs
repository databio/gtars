//! `rgblast report`: blastn hits turned into "best hit per genome, grouped by taxon, one line
//! per query".
//!
//! - [`parse`]: tabular 6/7 reader;
//! - [`fanout`]: `sallseqid` (or a subject digest) to (collection, name) pairs;
//! - [`rank`]: best hit per genome, ties, the call and its status;
//! - [`reduce`]: dropping the HSPs the report cannot use, as they stream in;
//! - [`taxonomy`]: the optional NCBI taxdump (lineage, LCA, rank folding);
//! - [`output`]: TSV, JSON and the detail table.

pub mod fanout;
pub mod output;
pub mod parse;
pub mod rank;
pub mod reduce;
pub mod taxonomy;

pub use fanout::DbMeta;
pub use parse::{Hit, HitsInput};
pub use rank::{
    Call, GenomeHit, QueryInput, QueryResult, RankOptions, Report, ReportMeta, StageHits,
    SummaryRow, TaxId, TaxonGroup, build_report,
};
pub use taxonomy::Taxonomy;

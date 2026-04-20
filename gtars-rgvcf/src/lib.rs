//! Reader for the rgvcf binary variant format (refget-anchored, sites-only).
//!
//! See `rgvcf/SPEC.md` for the on-disk format. This crate exposes a
//! read-only mmap-based parser with per-chromosome iteration and
//! sidecar lookup for non-SNV / multi-allelic sites.
//!
//! Typical use:
//! ```no_run
//! use gtars_rgvcf::{RgvcfReader, Record};
//! let reader = RgvcfReader::open("foo.rgvcf").unwrap();
//! assert_eq!(reader.collection_digest().len(), 32);
//! for chrom in reader.chromosomes() {
//!     for rec in reader.iter_chrom(&chrom.name).unwrap() {
//!         let rec = rec.unwrap();
//!         match rec {
//!             Record::Snv { pos, alt, .. } => { let _ = (pos, alt); }
//!             Record::Complex { pos, ref_, alts, .. } => {
//!                 let _ = (pos, ref_, alts);
//!             }
//!         }
//!     }
//! }
//! ```

mod reader;

pub use reader::{ChromEntry, Record, RgvcfReader, RgvcfRecordIter};

pub const MAGIC_HEADER: &[u8; 4] = b"RGVC";
pub const MAGIC_FOOTER: &[u8; 4] = b"CFGR";
pub const VERSION: u8 = 1;
pub const DIGEST_LEN: usize = 32;

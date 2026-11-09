//! `rgblastn`: blastn's flags (single-dash spelling) plus the store and selection flags.
//! The binary normalises `-flag` to `--flag` before clap sees it, so `RgblastnArgs` uses
//! clap's long names; `TryFrom<&RgblastnArgs>` re-checks every rule for callers that skip
//! clap. blastn flags rgblastn does not implement are refused with `Error::unsupported`.

use super::SelectionArgs;
use crate::api::rgblastn::{QuerySource, RgblastnRequest, rgblastn, write_tabular};
use crate::engine::options::{DustSpec, Strand, Task};
use crate::engine::search::SearchParams;
use crate::engine::subject::DbStatsOverride;
use crate::error::Error;
use crate::progress::{ProgressEvent, describe};
use crate::select::Selection;
use clap::Parser;
use std::path::PathBuf;

/// blastn flags rgblastn knows about but does not implement; spelled as the user types them.
/// rgblastn covers the common searches; anything else needs a database and stock blastn.
pub const UNSUPPORTED_FLAGS: &[&str] = &[
    "reward",
    "penalty",
    "gapopen",
    "gapextend",
    "xdrop_ungap",
    "xdrop_gap",
    "xdrop_gap_final",
    "template_type",
    "template_length",
    "qcov_hsp_perc",
    "culling_limit",
    "best_hit_overhang",
    "best_hit_score_edge",
    "subject_besthit",
    "lcase_masking",
    "parse_deflines",
    "dbsize",
    "taxids",
    "negative_taxids",
    "taxidlist",
    "negative_taxidlist",
    "no_taxid_expansion",
    "db",
    "subject",
    "subject_loc",
    "query_loc",
    "sum_stats",
    "use_index",
    "index_name",
    "db_soft_mask",
    "db_hard_mask",
    "gilist",
    "negative_gilist",
    "seqidlist",
    "negative_seqidlist",
    "mt_mode",
    "window_masker_taxid",
    "window_masker_db",
    "ungapped",
    "no_greedy",
    "min_raw_gapped_score",
    "off_diagonal_range",
    "num_descriptions",
    "num_alignments",
    "line_length",
    "html",
    "show_gis",
    "import_search_strategy",
    "export_search_strategy",
    "remote",
    "entrez_query",
    "searchsp",
    "window_size",
    "sorthits",
    "sorthsps",
    "max_hsps_per_subject",
    "filtering_db",
    "ipglist",
    "negative_ipglist",
];

#[derive(Parser, Debug, Clone, PartialEq)]
#[command(
    name = "rgblastn",
    version,
    about = "Nucleotide BLAST that searches a RefgetStore directly (NCBI BLASTN 2.17.0+ algorithms)",
    after_help = "rgblastn covers common searches (tasks megablast, blastn and blastn-short; -outfmt 6). For any other blastn option, build a database with rgmakeblastdb and use blastn.",
    disable_help_flag = true
)]
pub struct RgblastnArgs {
    /// Local store directory, or an http(s) store URL
    #[arg(long, value_name = "PATH|URL")]
    pub store: String,
    #[command(flatten)]
    pub selection: SelectionArgs,
    /// Query FASTA file
    #[arg(long, value_name = "FILE")]
    pub query: Option<PathBuf>,
    /// megablast, blastn or blastn-short
    #[arg(long, value_name = "TASK", default_value = "megablast")]
    pub task: String,
    /// Expectation value threshold; default 10 (1000 for blastn-short), as blastn
    #[arg(long)]
    pub evalue: Option<f64>,
    /// 6, optionally followed by field names (blastn's grammar)
    #[arg(long, value_name = "FMT", default_value = "6")]
    pub outfmt: String,
    #[arg(long, value_name = "N", default_value_t = 500, value_parser = clap::value_parser!(u32).range(1..))]
    pub max_target_seqs: u32,
    #[arg(long, value_name = "N", value_parser = clap::value_parser!(u32).range(1..))]
    pub max_hsps: Option<u32>,
    #[arg(long, value_name = "N", value_parser = clap::value_parser!(u32).range(4..))]
    pub word_size: Option<u32>,
    /// `yes`, `no`, or `level window linker`
    #[arg(long, value_name = "SPEC")]
    pub dust: Option<String>,
    #[arg(long, value_name = "BOOL")]
    pub soft_masking: Option<String>,
    #[arg(long, value_name = "STRAND", default_value = "both")]
    pub strand: String,
    #[arg(long, default_value_t = 0.0)]
    pub perc_identity: f64,
    /// Search threads (blastn's -num_threads; the database is split, as -mt_mode 0)
    #[arg(long, value_name = "N", default_value_t = 1, value_parser = clap::value_parser!(usize))]
    pub num_threads: usize,
    /// STATS_NSEQ of a wrapper alias (with --stats-totlen): replaces both statistics terms
    #[arg(long, value_name = "N", value_parser = clap::value_parser!(u32).range(1..))]
    pub stats_nseq: Option<u32>,
    #[arg(long, value_name = "N", value_parser = clap::value_parser!(u64).range(1..))]
    pub stats_totlen: Option<u64>,
    /// Output file (default stdout)
    #[arg(long, value_name = "FILE")]
    pub out: Option<PathBuf>,
    #[command(flatten)]
    pub fetch: super::FetchFlags,
    /// Progress to stderr
    #[arg(long)]
    pub verbose: bool,
    #[arg(long, short = 'h')]
    pub help: bool,
}

impl Default for RgblastnArgs {
    fn default() -> Self {
        RgblastnArgs::try_parse_from(["rgblastn", "--store", "", "--all"]).expect("defaults parse")
    }
}

fn invalid(msg: impl Into<String>) -> Error {
    Error::InvalidOptions(msg.into())
}

impl TryFrom<&RgblastnArgs> for RgblastnRequest {
    type Error = Error;
    fn try_from(a: &RgblastnArgs) -> Result<Self, Error> {
        let query = a
            .query
            .clone()
            .ok_or_else(|| invalid("-query is required"))?;
        let task = match Task::parse(&a.task) {
            Some(t) => t,
            // blastn's other nucleotide tasks
            None if ["dc-megablast", "blastn-fast", "rmblastn"].contains(&a.task.as_str()) => {
                return Err(Error::unsupported(&format!("-task {}", a.task)));
            }
            None => return Err(invalid(format!("unknown -task {}", a.task))),
        };
        let strand = match a.strand.as_str() {
            "both" => Strand::Both,
            "plus" => Strand::Plus,
            "minus" => Strand::Minus,
            s => {
                return Err(invalid(format!(
                    "-strand must be both, plus or minus, not {s}"
                )));
            }
        };
        let dust = match a.dust.as_deref() {
            None => None,
            Some("no") => Some(None),
            Some("yes") => Some(Some(DustSpec::default())),
            Some(spec) => {
                let v: Vec<i32> = spec
                    .split_whitespace()
                    .filter_map(|t| t.parse().ok())
                    .collect();
                if v.len() != 3 {
                    return Err(invalid("-dust takes yes, no, or 'level window linker'"));
                }
                Some(Some(DustSpec {
                    level: v[0],
                    window: v[1],
                    linker: v[2],
                }))
            }
        };
        let soft_masking = match a.soft_masking.as_deref() {
            None => true,
            Some("true") | Some("T") | Some("yes") => true,
            Some("false") | Some("F") | Some("no") => false,
            Some(s) => {
                return Err(invalid(format!(
                    "-soft_masking must be true or false, not {s}"
                )));
            }
        };
        if a.evalue.is_some_and(|e| e <= 0.0) {
            return Err(invalid("-evalue must be positive"));
        }
        if a.max_target_seqs < 1 {
            return Err(invalid("-max_target_seqs must be at least 1"));
        }
        if let Some(w) = a.word_size
            && w < 4
        {
            return Err(invalid(
                "Word-size must be 4 or greater for nucleotide comparison",
            ));
        }
        let stats = match (a.stats_nseq, a.stats_totlen) {
            (None, None) => None,
            (Some(n), Some(l)) => Some(DbStatsOverride {
                num_seqs: n,
                total_length: l,
            }),
            _ => {
                return Err(invalid(
                    "--stats-nseq and --stats-totlen must be given together",
                ));
            }
        };
        if a.num_threads == 0 {
            return Err(invalid("-num_threads must be at least 1"));
        }
        a.fetch.check(&a.store)?;
        let selection = Selection::try_from(&a.selection)?;
        let params = SearchParams {
            task,
            evalue: a.evalue,
            word_size: a.word_size,
            dust,
            soft_masking,
            strand,
            perc_identity: a.perc_identity,
            max_target_seqs: a.max_target_seqs,
            max_hsps: a.max_hsps,
            stats,
            threads: a.num_threads,
            ..SearchParams::default()
        };
        // the engine validates the values (SearchParams::to_options) when the search runs;
        // validate here too so callers that skip clap get the same error before opening a store
        params.to_options()?;
        let mut req = RgblastnRequest::new(a.store.clone(), selection, QuerySource::Path(query));
        req.cache_dir = a.fetch.cache_dir.clone();
        req.params = params;
        req.outfmt = a.outfmt.clone();
        req.remote = a.fetch.remote();
        Ok(req)
    }
}

/// blastn spells long flags with one dash and underscores (`-max_target_seqs`); clap wants
/// two dashes and hyphens. Rewrite `-flag_name` to `--flag-name` unless it is a negative
/// number or a single letter.
pub fn normalize_argv<I: IntoIterator<Item = String>>(args: I) -> Vec<String> {
    args.into_iter()
        .map(|a| {
            if a.len() > 2
                && a.starts_with('-')
                && !a.starts_with("--")
                && !a[1..].chars().next().unwrap().is_ascii_digit()
                && !a[1..].starts_with('.')
            {
                format!("--{}", a[1..].replace('_', "-"))
            } else {
                a
            }
        })
        .collect()
}

/// Which raw argument names an unsupported blastn flag, if any.
pub fn unsupported_flag(args: &[String]) -> Option<String> {
    for a in args {
        let name = a.trim_start_matches('-');
        let name = name.split('=').next().unwrap_or(name);
        if a.starts_with('-') && UNSUPPORTED_FLAGS.contains(&name) {
            return Some(format!("-{name}"));
        }
    }
    None
}

/// The exit code and stderr line for an error, in blastn's words where blastn has them.
pub fn report_error(e: &Error) -> i32 {
    match e {
        Error::InvalidOptions(_) => {
            eprintln!("rgblastn: {e}");
            2
        }
        Error::Cancelled => {
            eprintln!("rgblastn: {e}");
            130
        }
        _ => {
            eprintln!("rgblastn: {e}");
            1
        }
    }
}

pub fn rgblastn_main() -> i32 {
    let raw: Vec<String> = std::env::args().collect();
    if let Some(f) = unsupported_flag(&raw[1..]) {
        eprintln!("{}", Error::unsupported(&f));
        return 2;
    }
    let argv = normalize_argv(raw);
    if argv.iter().any(|a| a == "--help" || a == "-h") {
        use clap::CommandFactory;
        let mut cmd = RgblastnArgs::command();
        let _ = cmd.print_long_help();
        println!();
        return 0;
    }
    let args = match RgblastnArgs::try_parse_from(&argv) {
        Ok(a) => a,
        Err(e) => {
            let _ = e.print();
            return 2;
        }
    };
    let req = match RgblastnRequest::try_from(&args) {
        Ok(r) => r,
        Err(e) => return report_error(&e),
    };
    let spec = match crate::format::tabular::TabularSpec::parse_outfmt(&req.outfmt) {
        Ok(s) => s,
        Err(e) => return report_error(&e),
    };
    for tok in &spec.ignored {
        eprintln!("rgblastn: -outfmt: '{tok}' is not a format specifier; ignored (as blastn does)");
    }
    let verbose = args.verbose;
    let mut progress = |e: &ProgressEvent| -> bool {
        if verbose {
            match e {
                ProgressEvent::Opened { .. }
                | ProgressEvent::SearchStarted { .. }
                | ProgressEvent::SearchProgress { .. }
                | ProgressEvent::Traceback { .. }
                | ProgressEvent::SearchFinished(_)
                | ProgressEvent::SearchWarning(_) => {
                    eprintln!("{}", describe(e))
                }
                _ => {}
            }
        } else if let ProgressEvent::SearchWarning(_) = e {
            eprintln!("{}", describe(e));
        }
        true
    };
    let outcome = match rgblastn(&req, &mut progress) {
        Ok(o) => o,
        Err(e) => return report_error(&e),
    };
    let mut out: Box<dyn std::io::Write> = match &args.out {
        Some(p) => match std::fs::File::create(p) {
            Ok(f) => Box::new(std::io::BufWriter::new(f)),
            Err(e) => {
                eprintln!("rgblastn: {}: {e}", p.display());
                return 1;
            }
        },
        None => Box::new(std::io::BufWriter::new(std::io::stdout().lock())),
    };
    if let Err(e) = write_tabular(&outcome, &mut out) {
        return report_error(&e);
    }
    if out.flush().is_err() {
        return 1;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(extra: &[&str]) -> Result<RgblastnArgs, clap::Error> {
        let mut v: Vec<String> = ["rgblastn", "--store", "s", "--all", "-query", "q"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        v.extend(extra.iter().map(|s| s.to_string()));
        RgblastnArgs::try_parse_from(normalize_argv(v))
    }

    #[test]
    fn normalize_and_defaults() {
        let v = normalize_argv(
            [
                "rgblastn",
                "-query",
                "q.fa",
                "-evalue",
                "1e-5",
                "--store",
                "s",
                "-h",
                "-5",
                "-max_target_seqs",
                "7",
            ]
            .map(String::from),
        );
        assert_eq!(
            v,
            [
                "rgblastn",
                "--query",
                "q.fa",
                "--evalue",
                "1e-5",
                "--store",
                "s",
                "-h",
                "-5",
                "--max-target-seqs",
                "7",
            ]
        );
        // every blastn flag the struct declares parses in blastn's spelling
        let a = parse(&[
            "-max_target_seqs",
            "7",
            "-word_size",
            "11",
            "-soft_masking",
            "false",
            "-perc_identity",
            "90",
            "-max_hsps",
            "3",
            "-strand",
            "plus",
            "-dust",
            "no",
            "-evalue",
            "1e-5",
            "-num_threads",
            "2",
        ])
        .unwrap();
        let req = RgblastnRequest::try_from(&a).unwrap();
        let p = &req.params;
        assert_eq!(
            (p.max_target_seqs, p.word_size, p.max_hsps, p.strand, p.dust),
            (7, Some(11), Some(3), Strand::Plus, Some(None))
        );
        assert_eq!(
            (p.soft_masking, p.perc_identity, p.threads),
            (false, 90.0, 2)
        );
        let d = RgblastnArgs::default();
        assert_eq!(
            (
                d.evalue,
                d.max_target_seqs,
                d.task.as_str(),
                d.outfmt.as_str()
            ),
            (None, 500, "megablast", "6")
        );
    }

    #[test]
    fn removed_blastn_flags_are_refused_by_name() {
        assert_eq!(unsupported_flag(&["-db".to_string()]), Some("-db".into()));
        for f in [
            "-dbsize",
            "-taxids",
            "-template_type",
            "-culling_limit",
            "-best_hit_overhang",
            "-subject_besthit",
            "-lcase_masking",
            "-parse_deflines",
            "-reward",
            "-gapopen",
            "-xdrop_gap",
            "-qcov_hsp_perc",
        ] {
            assert_eq!(unsupported_flag(&[f.to_string()]), Some(f.into()), "{f}");
        }
        let e = Error::unsupported("-dbsize");
        assert!(e.to_string().contains("rgmakeblastdb"), "{e}");
        for t in ["dc-megablast", "blastn-fast"] {
            let a = parse(&["-task", t]).unwrap();
            assert!(matches!(
                RgblastnRequest::try_from(&a),
                Err(Error::InvalidOptions(m)) if m.contains("not supported by rgblastn")
            ));
        }
        let a = parse(&["-task", "bogus"]).unwrap();
        assert!(matches!(
            RgblastnRequest::try_from(&a),
            Err(Error::InvalidOptions(_))
        ));
        // rgblastn's own removed flags are plain parse errors
        assert!(parse(&["--dump-stats"]).is_err());
        assert!(parse(&["--batch-mbp", "4"]).is_err());
    }

    #[test]
    fn stats_override_needs_both_values() {
        let a = parse(&["--stats-nseq", "4", "--stats-totlen", "3088286401"]).unwrap();
        assert_eq!(
            RgblastnRequest::try_from(&a).unwrap().params.stats,
            Some(DbStatsOverride {
                num_seqs: 4,
                total_length: 3_088_286_401
            })
        );
        let a = parse(&["--stats-nseq", "4"]).unwrap();
        assert!(RgblastnRequest::try_from(&a).is_err());
    }
}

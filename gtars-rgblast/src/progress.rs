//! Progress events and the one wording for them. Library code never prints: everything a user
//! sees goes through these events (or the returned summary), and every surface (the binary,
//! gtars-cli, Python, R, the refget CLI) prints [`describe`]'s lines.

use crate::blastdb::BuildSummary;
use crate::blastdb::preflight::PreflightSummary;
use serde::{Deserialize, Serialize};
use std::time::{Duration, Instant};

/// Events are rare (a few per second at most), so their size does not matter.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ProgressEvent {
    Opened {
        store: String,
        mode: String,
        remote: bool,
    },
    Planned(PreflightSummary),
    Sequences {
        done: u64,
        total: u64,
        bases_done: u64,
        bases_total: u64,
        elapsed_s: f64,
    },
    /// A remote store: the selection's sequences are being pulled into the cache.
    Pulling {
        sequences: u64,
        bases: u64,
    },
    /// Counting N bases of sequences the cache does not know (announced with `done` = 0).
    NScan {
        done: u64,
        total: u64,
        bytes_read: u64,
        bytes_total: u64,
    },
    /// "lmdb", "oid files", "sidecars", "njs", "views"
    Finishing {
        stage: String,
    },
    Warning(String),
    Finished(Box<BuildSummary>),
    /// `rgblast`: the engine is starting, with its final arguments.
    EngineStart {
        stage: String,
        engine: String,
        args: Vec<String>,
        /// `-num_threads` was added (not given).
        threads_added: bool,
    },
    /// `rgblast`: the heartbeat while an engine runs (at least once a second).
    EngineRunning {
        stage: String,
        elapsed_s: f64,
    },
    /// `rgblast`: one line the engine wrote to stderr.
    EngineStderr {
        line: String,
    },
    /// `rgblast origin`: a stage starts.
    Stage {
        stage: String,
        queries: usize,
        subjects_gbp: f64,
    },
    /// rgblastn: the search starts.
    SearchStarted {
        queries: u32,
        subjects: u32,
        bases: u64,
    },
    /// rgblastn: preliminary-stage progress (at least once per `progress_interval_ms`).
    SearchProgress {
        subjects_done: u32,
        subjects_total: u32,
        bases_done: u64,
        hits: u64,
        elapsed_s: f64,
    },
    /// rgblastn: traceback progress, in subjects with hits.
    Traceback {
        done: u32,
        total: u32,
    },
    /// rgblastn: a notice the binary prints (for example the `-taxids` warning).
    SearchWarning(String),
    /// rgblastn: done.
    SearchFinished(crate::engine::search::SearchSummary),
}

/// At most one event per interval.
pub struct Throttle {
    interval: Duration,
    last: Option<Instant>,
}

impl Throttle {
    pub fn new(interval_ms: u64) -> Self {
        Throttle {
            interval: Duration::from_millis(interval_ms),
            last: None,
        }
    }

    /// True (and restart the clock) if an interval has passed since the last `true`.
    /// The first call starts the clock and returns false.
    pub fn ready(&mut self) -> bool {
        let now = Instant::now();
        match self.last {
            None => {
                self.last = Some(now);
                false
            }
            Some(t) if now.duration_since(t) >= self.interval => {
                self.last = Some(now);
                true
            }
            _ => false,
        }
    }
}

/// 27721 -> "27,721".
pub fn commas(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::new();
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// Bytes, decimal units: "10.48 GB", "52.3 MB", "812 B".
pub fn gb(bytes: u64) -> String {
    let b = bytes as f64;
    if b >= 1e9 {
        format!("{:.2} GB", b / 1e9)
    } else if b >= 1e6 {
        format!("{:.1} MB", b / 1e6)
    } else if b >= 1e3 {
        format!("{:.1} kB", b / 1e3)
    } else {
        format!("{bytes} B")
    }
}

/// Bases: "41.81 Gbp", "3.10 Mbp", "12.0 kbp", "800 bp".
pub fn bp(n: u64) -> String {
    let b = n as f64;
    if b >= 1e9 {
        format!("{:.2} Gbp", b / 1e9)
    } else if b >= 1e6 {
        format!("{:.2} Mbp", b / 1e6)
    } else if b >= 1e3 {
        format!("{:.1} kbp", b / 1e3)
    } else {
        format!("{n} bp")
    }
}

fn file_name(p: &std::path::Path) -> String {
    p.file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| p.display().to_string())
}

/// The line(s) `rgmakeblastdb` prints for an event.
pub fn describe(e: &ProgressEvent) -> String {
    match e {
        ProgressEvent::Opened {
            store,
            mode,
            remote,
        } => {
            format!(
                "rgmakeblastdb: opened {store} ({mode}, {})",
                if *remote { "remote" } else { "local" }
            )
        }
        ProgressEvent::Planned(s) => {
            let skipped = format!(
                "  skipped: {} oversize, {} non-nucleotide, {} empty; sanitized {} ids{}",
                commas(s.skipped_oversize),
                commas(s.skipped_non_nucleotide),
                commas(s.skipped_empty),
                commas(s.sanitized),
                if s.sanitized > 0 {
                    format!(" ({}.sanitized.tsv)", file_name(&s.out))
                } else {
                    String::new()
                }
            );
            let free = match s.free_bytes {
                Some(f) => format!("{} free", gb(f)),
                None => "free space not checked".to_string(),
            };
            let selection = if s.selection.is_empty() {
                String::new()
            } else {
                format!("\n  selection: {}", s.selection)
            };
            format!(
                "rgmakeblastdb: {} ({}, {}): {} collections, {} names{selection}\n  database: {} sequences, {}, {} deflines, {} ids\n{skipped}\n  output: {}: {} exact + {} allowance, about {} volume{}; {free}",
                s.store,
                s.mode,
                if s.remote { "remote" } else { "local" },
                commas(s.collections),
                commas(s.names),
                commas(s.subjects),
                bp(s.total_bases),
                commas(s.deflines),
                commas(s.ids),
                s.out.display(),
                gb(s.exact_bytes),
                gb(s.allowance_bytes),
                s.predicted_volumes,
                if s.predicted_volumes == 1 { "" } else { "s" },
            )
        }
        ProgressEvent::Sequences {
            done,
            total,
            bases_done,
            bases_total,
            elapsed_s,
        } => {
            let pct = if *bases_total > 0 {
                100.0 * *bases_done as f64 / *bases_total as f64
            } else {
                100.0
            };
            let rate = *bases_done as f64 / elapsed_s.max(1e-9);
            let eta = if rate > 0.0 {
                bases_total.saturating_sub(*bases_done) as f64 / rate
            } else {
                0.0
            };
            format!(
                "sequences {}/{}  {}/{}  {:.0}%  {}/s  eta {:.0} s",
                commas(*done),
                commas(*total),
                bp(*bases_done),
                bp(*bases_total),
                pct,
                bp(rate as u64),
                eta
            )
        }
        ProgressEvent::Pulling { sequences, bases } => format!(
            "pulling {} sequences ({}) into the cache",
            commas(*sequences),
            bp(*bases)
        ),
        ProgressEvent::NScan {
            done,
            total,
            bytes_read,
            bytes_total,
        } => {
            if *done == 0 && *bytes_read == 0 {
                format!(
                    "computing N counts for {} sequences, {} to read",
                    commas(*total),
                    gb(*bytes_total)
                )
            } else {
                format!(
                    "N counts {}/{}  {}/{}",
                    commas(*done),
                    commas(*total),
                    gb(*bytes_read),
                    gb(*bytes_total)
                )
            }
        }
        ProgressEvent::Finishing { stage } => format!("finishing: {stage}"),
        ProgressEvent::Warning(w) => format!("warning: {w}"),
        ProgressEvent::EngineStart {
            stage,
            engine,
            args,
            threads_added,
        } => format!(
            "rgblast: {}{engine} {}{}",
            if stage == "search" {
                String::new()
            } else {
                format!("{stage}: ")
            },
            shell_join(args),
            if *threads_added {
                " (-num_threads added)"
            } else {
                ""
            }
        ),
        ProgressEvent::EngineRunning { stage, elapsed_s } => {
            format!("rgblast: {stage} running ({elapsed_s:.0} s)")
        }
        ProgressEvent::EngineStderr { line } => format!("  engine: {line}"),
        ProgressEvent::Stage {
            stage,
            queries,
            subjects_gbp,
        } => format!(
            "origin: {stage}: {} queries against {subjects_gbp:.2} Gbp",
            commas(*queries as u64)
        ),
        ProgressEvent::SearchStarted {
            queries,
            subjects,
            bases,
        } => format!(
            "searching {} queries against {} subjects ({})",
            commas(*queries as u64),
            commas(*subjects as u64),
            bp(*bases)
        ),
        ProgressEvent::SearchProgress {
            subjects_done,
            subjects_total,
            bases_done,
            hits,
            elapsed_s,
        } => format!(
            "  {}/{} subjects, {} scanned, {} subjects with hits, {elapsed_s:.1} s",
            commas(*subjects_done as u64),
            commas(*subjects_total as u64),
            bp(*bases_done),
            commas(*hits)
        ),
        ProgressEvent::Traceback { done, total } => {
            format!(
                "  traceback {}/{} subjects",
                commas(*done as u64),
                commas(*total as u64)
            )
        }
        ProgressEvent::SearchWarning(w) => w.clone(),
        ProgressEvent::SearchFinished(s) => format!(
            "done: {} hits, {} HSPs for {} queries in {:.1} s (search {:.1} s, traceback {:.1} s)",
            commas(s.hits),
            commas(s.hsps),
            commas(s.queries as u64),
            s.elapsed_s,
            s.prelim_s,
            s.traceback_s
        ),
        ProgressEvent::Finished(s) => format!(
            "wrote {}: {} sequences, {}, {} volume{}, {} in {:.1} s\nsearch it with: blastn -db {} -query QUERY.fa",
            s.preflight.out.display(),
            commas(s.preflight.subjects),
            bp(s.preflight.total_bases),
            s.volumes,
            if s.volumes == 1 { "" } else { "s" },
            gb(s.bytes_written),
            s.elapsed_s,
            s.preflight.out.display()
        ),
    }
}

/// Arguments for display, quoting the ones with spaces or shell characters.
pub fn shell_join(args: &[String]) -> String {
    args.iter()
        .map(|a| {
            if a.is_empty() || a.chars().any(|c| " \t\"'$|;&<>*?()".contains(c)) {
                format!("\"{}\"", a.replace('"', "\\\""))
            } else {
                a.clone()
            }
        })
        .collect::<Vec<_>>()
        .join(" ")
}

/// The hint printed after a build that wrote a real taxid: how to filter exactly, and the
/// `-taxids` trap (without `-no_taxid_expansion` and NCBI's `taxonomy4blast.sqlite3`, blastn
/// warns and searches the whole database).
pub fn taxid_hint(taxa: &[crate::blastdb::TaxonCount]) -> Option<String> {
    let first = taxa.first()?;
    let n = taxa.len();
    Some(format!(
        "taxids written for {n} {}. Filter exactly with: blastn -taxids {} -no_taxid_expansion\n\
         Without -no_taxid_expansion, blastn needs NCBI's taxonomy4blast.sqlite3; if it is missing,\n\
         blastn only warns and searches the whole database. For descendant taxa and sscinames,\n\
         download https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz (66 MB: taxdb.btd, taxdb.bti,\n\
         taxonomy4blast.sqlite3) into a directory on BLASTDB, or run `update_blastdb.pl --decompress taxdb` there.",
        if n == 1 { "taxon" } else { "taxa" },
        first.taxid
    ))
}

/// One event as a flat JSON object for the bindings: `event` (the variant in snake case),
/// the variant's fields (a plain value goes under `value`) and `message`, the line the
/// command line prints.
pub fn event_json(ev: &ProgressEvent) -> serde_json::Value {
    let mut m = serde_json::Map::new();
    match serde_json::to_value(ev).unwrap_or_default() {
        serde_json::Value::String(name) => {
            m.insert("event".into(), name.into());
        }
        serde_json::Value::Object(o) => {
            for (name, body) in o {
                m.insert("event".into(), name.into());
                match body {
                    serde_json::Value::Object(fields) => m.extend(fields),
                    other => {
                        m.insert("value".into(), other);
                    }
                }
            }
        }
        _ => {}
    }
    m.insert("message".into(), describe(ev).into());
    serde_json::Value::Object(m)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blastdb::{BuildSummary, Timings};

    fn summary() -> PreflightSummary {
        PreflightSummary {
            store: "/path/jungle".into(),
            store_label: "jungle".into(),
            mode: "Encoded".into(),
            remote: false,
            collections: 79,
            names: 27_721,
            subjects: 6_359,
            deflines: 34_080,
            ids: 209_845,
            total_bases: 41_810_552_907,
            max_len: 248_956_422,
            sanitized: 225,
            exact_bytes: 10_480_000_000,
            allowance_bytes: 520_000_000,
            predicted_volumes: 4,
            free_bytes: Some(812_000_000_000),
            out: "/out/jungle".into(),
            ..Default::default()
        }
    }

    #[test]
    fn describe_every_event() {
        let e = ProgressEvent::Opened {
            store: "/path/jungle".into(),
            mode: "Encoded".into(),
            remote: false,
        };
        assert_eq!(
            describe(&e),
            "rgmakeblastdb: opened /path/jungle (Encoded, local)"
        );
        assert_eq!(
            describe(&ProgressEvent::Planned(summary())),
            "rgmakeblastdb: /path/jungle (Encoded, local): 79 collections, 27,721 names\n  database: 6,359 sequences, 41.81 Gbp, 34,080 deflines, 209,845 ids\n  skipped: 0 oversize, 0 non-nucleotide, 0 empty; sanitized 225 ids (jungle.sanitized.tsv)\n  output: /out/jungle: 10.48 GB exact + 520.0 MB allowance, about 4 volumes; 812.00 GB free"
        );
        let e = ProgressEvent::Sequences {
            done: 1234,
            total: 6359,
            bases_done: 12_300_000_000,
            bases_total: 41_800_000_000,
            elapsed_s: 10.0,
        };
        assert_eq!(
            describe(&e),
            "sequences 1,234/6,359  12.30 Gbp/41.80 Gbp  29%  1.23 Gbp/s  eta 24 s"
        );
        assert_eq!(
            describe(&ProgressEvent::Pulling {
                sequences: 2,
                bases: 2000
            }),
            "pulling 2 sequences (2.0 kbp) into the cache"
        );
        assert_eq!(
            describe(&ProgressEvent::Finishing {
                stage: "lmdb".into()
            }),
            "finishing: lmdb"
        );
        assert_eq!(
            describe(&ProgressEvent::Warning("skipped 1 empty sequences".into())),
            "warning: skipped 1 empty sequences"
        );
        let s = BuildSummary {
            preflight: summary(),
            taxa: vec![],
            files: vec![],
            sidecars: vec![],
            bytes_written: 10_740_000_000,
            ambiguous_bases: 0,
            volumes: 4,
            warnings: vec![],
            timings: Timings::default(),
            network: None,
            elapsed_s: 16.3,
        };
        assert_eq!(
            describe(&ProgressEvent::Finished(Box::new(s))),
            "wrote /out/jungle: 6,359 sequences, 41.81 Gbp, 4 volumes, 10.74 GB in 16.3 s\nsearch it with: blastn -db /out/jungle -query QUERY.fa"
        );
    }

    #[test]
    fn number_formats() {
        assert_eq!(commas(0), "0");
        assert_eq!(commas(999), "999");
        assert_eq!(commas(1000), "1,000");
        assert_eq!(commas(1_234_567), "1,234,567");
    }
}

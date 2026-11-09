//! `--engine rgblastn` (the default): the in-process engine, so bindings get it without
//! another install.
//!
//! `run` turns the pass-through blastn arguments into an `RgblastnRequest` with the same
//! clap struct and `TryFrom` the `rgblastn` binary uses, searches the `ResolvedSelection` it
//! was given (`subjects` is applied as a `--sequences-file` subject filter), takes
//! `RunSpec.stats` as the `STATS_NSEQ`/`STATS_TOTLEN` override, and writes the tabular output
//! where `EngineOutput` says. Flags rgblastn does not implement come back as
//! `Error::InvalidOptions` with rgblastn's wording.

use super::{BlastRunner, EngineInfo, EngineOutput, RunOutcome, RunSpec, SearchTarget};
use crate::api::rgblastn::{RgblastnRequest, run_resolved, write_tabular};
use crate::cli::RgblastnArgs;
use crate::cli::rgblastn::normalize_argv;
use crate::engine::subject::DbStatsOverride;
use crate::error::Error;
use crate::progress::ProgressEvent;
use crate::store::OpenedStore;
use clap::Parser;
use std::path::PathBuf;

pub struct RgblastnRunner;

impl RgblastnRunner {
    /// The request for a run's pass-through arguments.
    fn request(
        spec: &RunSpec,
        location: String,
        subjects: Option<&Vec<String>>,
    ) -> Result<(RgblastnRequest, Option<PathBuf>), Error> {
        let args = spec.args.clone();
        if let Some(f) = crate::cli::rgblastn::unsupported_flag(&args) {
            return Err(Error::unsupported(&f));
        }
        let mut argv = vec![
            "rgblastn".to_string(),
            "--store".to_string(),
            location,
            "--all".to_string(),
        ];
        if let Some(q) = &spec.query {
            argv.push("--query".into());
            argv.push(q.display().to_string());
        }
        argv.extend(normalize_argv(args));
        let a = RgblastnArgs::try_parse_from(&argv)
            .map_err(|e| Error::InvalidOptions(format!("rgblastn: {}", e.to_string().trim())))?;
        let mut req = RgblastnRequest::try_from(&a)?;
        if let Some(d) = subjects {
            req.selection.sequences = Some(d.clone());
        }
        if let Some(s) = &spec.stats {
            req.params.stats = Some(DbStatsOverride {
                num_seqs: s.num_seqs,
                total_length: s.total_length,
            });
        }
        Ok((req, a.out.clone()))
    }
}

impl BlastRunner for RgblastnRunner {
    fn describe(&self) -> String {
        format!("rgblastn {} (in process)", env!("CARGO_PKG_VERSION"))
    }

    fn info(&self) -> EngineInfo {
        EngineInfo {
            name: "rgblastn".into(),
            version: env!("CARGO_PKG_VERSION").into(),
            via: "in process".into(),
        }
    }

    fn run(
        &self,
        spec: &RunSpec,
        store: Option<&mut OpenedStore>,
        progress: &mut dyn FnMut(&ProgressEvent) -> bool,
    ) -> Result<RunOutcome, Error> {
        let (
            SearchTarget::Store {
                selection,
                subjects,
            },
            Some(store),
        ) = (&spec.target, store)
        else {
            return Err(Error::InvalidOptions(
                "rgblastn searches a store: give --store PATH|URL".into(),
            ));
        };
        let subjects = subjects.as_deref();
        let (req, out_flag) = Self::request(spec, store.location_string(), subjects)?;
        if !progress(&ProgressEvent::EngineStart {
            stage: spec.stage.clone(),
            engine: self.describe(),
            args: spec.args.clone(),
            threads_added: spec.threads_added,
        }) {
            return Err(Error::Cancelled);
        }
        let mut stderr_tail: Vec<String> = Vec::new();
        let mut sel = selection.clone();
        if let Some(d) = subjects {
            sel.sequences = Some(d.clone());
        }
        let outcome = {
            let mut relay = |e: &ProgressEvent| -> bool {
                if let ProgressEvent::SearchWarning(w) = e {
                    stderr_tail.push(w.clone());
                    if stderr_tail.len() > 50 {
                        stderr_tail.remove(0);
                    }
                    return progress(&ProgressEvent::EngineStderr { line: w.clone() });
                }
                progress(e)
            };
            run_resolved(store, &sel, &req, &mut relay)?
        };
        let mut captured = None;
        let target = match (&spec.output, out_flag) {
            (EngineOutput::File(p), _) => Some(p.clone()),
            (EngineOutput::Inherit, Some(p)) => Some(p),
            _ => None,
        };
        match (&spec.output, target) {
            (EngineOutput::Capture, _) => {
                let mut buf = Vec::new();
                write_tabular(&outcome, &mut buf)?;
                captured = Some(String::from_utf8_lossy(&buf).into_owned());
            }
            (_, Some(p)) => {
                let f = std::fs::File::create(&p).map_err(|e| Error::io(&p, e))?;
                let mut w = std::io::BufWriter::new(f);
                write_tabular(&outcome, &mut w)?;
                std::io::Write::flush(&mut w).map_err(|e| Error::io(&p, e))?;
            }
            (_, None) => {
                let stdout = std::io::stdout();
                let mut w = std::io::BufWriter::new(stdout.lock());
                write_tabular(&outcome, &mut w)?;
                std::io::Write::flush(&mut w).map_err(|e| Error::io("<stdout>", e))?;
            }
        }
        Ok(RunOutcome {
            exit_code: 0,
            db: None,
            output: captured,
            engine_stderr: stderr_tail,
            args: spec.args.clone(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn added_threads_reach_rgblastn() {
        let mut args = vec!["-task".to_string(), "megablast".to_string()];
        let added = crate::run::passthrough::ensure_threads(&mut args).unwrap();
        let spec = RunSpec {
            stage: "search".into(),
            target: SearchTarget::Db(PathBuf::new()),
            args,
            stats: None,
            query: Some(PathBuf::from("q.fa")),
            output: EngineOutput::Capture,
            threads_added: true,
        };
        let (req, _) = RgblastnRunner::request(&spec, "store".into(), None).unwrap();
        assert_eq!(req.params.threads, added);
    }
}

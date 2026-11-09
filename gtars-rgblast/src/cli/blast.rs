//! `BlastCommand`: every verb behind one enum, so `rgmakeblastdb`, `rgblast`, `gtars blast` and
//! the bindings (`run_json`) mount the same parsers and the same runners.

use super::rgblast::{
    INTERRUPTED, fail, install_sigint, run_origin, run_report, run_reps, run_search, usage,
};
use super::{
    BuildArgs, NstatsArgs, ORIGIN_AFTER_HELP, OriginArgs, ReportArgs, RepsArgs, SearchArgs,
    VerifyArgs, ViewArgs,
};
use crate::api::{self, ViewRequest};
use crate::progress::{self, describe};
use crate::{BuildRequest, NstatsRequest, ProgressEvent, VerifyRequest};
use std::sync::atomic::Ordering;

/// All the verbs. `rgmakeblastdb` takes the first four, `rgblast` the rest.
#[derive(clap::Subcommand, Debug, Clone, PartialEq)]
pub enum BlastCommand {
    /// Write a BLAST version 5 nucleotide database from a store
    Build(BuildArgs),
    /// Check a database against its store: every id resolves, every md5 matches
    Verify(VerifyArgs),
    /// Write subset views of a database (.nal + .bsl) with exact statistics, no rebuild
    View(ViewArgs),
    /// Count N bases per collection (fills the N-count cache) to choose --max-n-fraction
    Nstats(NstatsArgs),
    /// Search the selection with rgblastn, or a database you built with blastn (--engine blastn)
    Search(SearchArgs),
    /// Best hit per genome, grouped by species, one line per query
    Report(ReportArgs),
    /// Where did these sequences come from: representatives first, then fan-out
    #[command(after_help = ORIGIN_AFTER_HELP)]
    Origin(Box<OriginArgs>),
    /// The representative of each species and the staged-search cost prediction
    Reps(RepsArgs),
}

impl BlastCommand {
    /// Run one command and return its exit code: 0 success, 1 runtime error, 2 usage error or
    /// missing engine, 130 cancelled, the engine's own code for `search`. A closed output pipe
    /// exits 0. `prog` prefixes messages ("gtars blast", "rgmakeblastdb", "rgblast").
    pub fn run_code(self, prog: &str) -> i32 {
        install_sigint();
        match &self {
            BlastCommand::Build(a) => run_build(prog, a),
            BlastCommand::Verify(a) => run_verify(prog, a),
            BlastCommand::View(a) => run_view(prog, a),
            BlastCommand::Nstats(a) => run_nstats(prog, a),
            BlastCommand::Search(a) => run_search(prog, a),
            BlastCommand::Report(a) => run_report(prog, a),
            BlastCommand::Origin(a) => run_origin(prog, a),
            BlastCommand::Reps(a) => run_reps(prog, a),
        }
    }
}

/// The `main` of a binary that mounts some of the verbs (`rgmakeblastdb`, `rgblast`): parse the
/// process arguments, run the verb and return its exit code.
pub fn main_with_verbs(prog: &'static str, about: &'static str, verbs: &[&str]) -> i32 {
    use clap::{FromArgMatches, Subcommand};
    let all = BlastCommand::augment_subcommands(clap::Command::new(prog));
    let cmd = clap::Command::new(prog)
        .version(env!("CARGO_PKG_VERSION"))
        .about(about)
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommands(
            all.get_subcommands()
                .filter(|s| verbs.contains(&s.get_name()))
                .cloned(),
        );
    match BlastCommand::from_arg_matches(&cmd.get_matches()) {
        Ok(c) => c.run_code(prog),
        Err(e) => e.exit(),
    }
}

/// The progress printer of build.
pub fn db_printer(quiet: bool, verbose: bool) -> impl FnMut(&ProgressEvent) -> bool {
    move |e: &ProgressEvent| -> bool {
        let show = match e {
            ProgressEvent::Warning(_) => true,
            ProgressEvent::Opened { .. } | ProgressEvent::Finishing { .. } => verbose,
            _ => !quiet,
        };
        if show {
            eprintln!("{}", describe(e));
        }
        !INTERRUPTED.load(Ordering::SeqCst)
    }
}

fn run_build(prog: &str, a: &BuildArgs) -> i32 {
    let req = match BuildRequest::try_from(a) {
        Ok(r) => r,
        Err(e) => return usage(prog, &e),
    };
    if a.dry_run {
        return match api::preflight(&req) {
            Ok(s) => {
                eprintln!("{}", describe(&ProgressEvent::Planned(s)));
                0
            }
            Err(e) => fail(prog, &e),
        };
    }
    let (quiet, verbose) = (a.quiet, a.verbose);
    let mut cb = db_printer(quiet, verbose);
    match api::build(&req, &mut cb) {
        Ok(s) => {
            if !quiet && let Some(h) = progress::taxid_hint(&s.taxa) {
                eprintln!("{h}");
            }
            if verbose {
                let t = &s.timings;
                eprintln!(
                    "time: plan {:.2}s  sequences {:.2}s  finish {:.2}s  total {:.2}s",
                    t.plan_s, t.sequences_s, t.finish_s, s.elapsed_s
                );
                eprintln!(
                    "cpu-sum: read {:.2}s  convert {:.2}s  write (writer thread) {:.2}s; {} ambiguous bases",
                    t.cpu_read_s,
                    t.cpu_convert_s,
                    t.write_s,
                    progress::commas(s.ambiguous_bases)
                );
                if let Some(n) = &s.network {
                    eprintln!("{}", n.line());
                }
                for f in s.files.iter().chain(s.sidecars.iter()) {
                    let len = std::fs::metadata(f).map(|m| m.len()).unwrap_or(0);
                    eprintln!("  {}  {}", f.display(), progress::commas(len));
                }
            }
            0
        }
        Err(e) => fail(prog, &e),
    }
}

fn run_verify(prog: &str, a: &VerifyArgs) -> i32 {
    let req = match VerifyRequest::try_from(a) {
        Ok(r) => r,
        Err(e) => return usage(prog, &e),
    };
    let mut cb = |e: &ProgressEvent| -> bool {
        if matches!(e, ProgressEvent::Sequences { .. }) {
            eprintln!("{}", describe(e));
        }
        !INTERRUPTED.load(Ordering::SeqCst)
    };
    match api::verify(&req, &mut cb) {
        Ok(r) => {
            eprintln!(
                "verify {}: {} ids checked, {} failed; {} sequences checked, {} md5 mismatches; totals {}",
                req.db.display(),
                progress::commas(r.ids_checked),
                progress::commas(r.ids_failed),
                progress::commas(r.sequences_checked),
                progress::commas(r.md5_mismatches),
                if r.totals_ok {
                    "match .njs"
                } else {
                    "DO NOT match .njs"
                }
            );
            for ex in &r.examples {
                eprintln!("  {ex}");
            }
            if r.ok() {
                eprintln!("OK");
                0
            } else {
                eprintln!("FAILED");
                1
            }
        }
        Err(e) => fail(prog, &e),
    }
}

fn run_view(prog: &str, a: &ViewArgs) -> i32 {
    let req = match ViewRequest::try_from(a) {
        Ok(r) => r,
        Err(e) => return usage(prog, &e),
    };
    match api::view(&req, &mut |_| !INTERRUPTED.load(Ordering::SeqCst)) {
        Ok(s) => {
            for v in &s.views {
                eprintln!(
                    "view {}: {} sequences, {} bases ({})",
                    v.name,
                    progress::commas(v.nseq),
                    progress::commas(v.length),
                    v.nal.display()
                );
            }
            eprintln!(
                "wrote {} view{} in {} (index {})\nsearch one with: blastn -db {}/VIEW -query QUERY.fa",
                s.views.len(),
                if s.views.len() == 1 { "" } else { "s" },
                s.out_dir.display(),
                s.index.display(),
                s.out_dir.display()
            );
            0
        }
        Err(e) => fail(prog, &e),
    }
}

fn run_nstats(prog: &str, a: &NstatsArgs) -> i32 {
    let req = match NstatsRequest::try_from(a) {
        Ok(r) => r,
        Err(e) => return usage(prog, &e),
    };
    let mut cb = |e: &ProgressEvent| -> bool {
        if matches!(e, ProgressEvent::NScan { .. } | ProgressEvent::Warning(_)) {
            eprintln!("{}", describe(e));
        }
        !INTERRUPTED.load(Ordering::SeqCst)
    };
    match api::nstats(&req, &mut cb) {
        Ok(s) => {
            println!("#collection\tdisplay\tspecies\trows\tbases\tn\tn_fraction");
            for r in &s.collections {
                println!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    r.digest,
                    r.display,
                    r.species.as_deref().unwrap_or(""),
                    r.rows,
                    r.bases,
                    r.n,
                    r.n_fraction.map(|f| format!("{f:.4}")).unwrap_or_default()
                );
            }
            eprintln!(
                "{} collections, {} sequences counted{}",
                s.collections.len(),
                progress::commas(s.sequences),
                if s.excluded.is_empty() {
                    String::new()
                } else {
                    format!("; {} excluded by the N filter", s.excluded.len())
                }
            );
            0
        }
        Err(e) => fail(prog, &e),
    }
}

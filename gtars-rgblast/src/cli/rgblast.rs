//! The `rgblast` command line: parses with the library's argument structs, calls
//! `gtars_rgblast::api`, and renders results and progress lines.

use super::origin::{OriginArgs, RepsArgs};
use super::report::{FormatArg, ReportArgs, ReportOutputArgs};
use super::search::SearchArgs;
use crate::api::{self, OriginRequest, ReportRequest, RepsRequest, SearchRequest};
use crate::error::Error;
use crate::progress::{ProgressEvent, describe};
use crate::report::output;
use crate::report::rank::Report;
use crate::report::taxonomy::TAXDUMP_HINT;
use std::io::Write;
use std::sync::atomic::{AtomicBool, Ordering};

pub(super) static INTERRUPTED: AtomicBool = AtomicBool::new(false);

extern "C" fn on_sigint(_: libc::c_int) {
    INTERRUPTED.store(true, Ordering::SeqCst);
}

/// Make Ctrl-C cancel through `printer` (removes partial output, exit 130).
pub fn install_sigint() {
    #[cfg(unix)]
    unsafe {
        libc::signal(libc::SIGINT, on_sigint as *const () as libc::sighandler_t);
    }
}

pub(super) fn usage(prog: &str, e: &Error) -> i32 {
    eprintln!("{prog}: error: {e}");
    eprintln!("\nFor more information, try '--help'.");
    2
}

/// A closed output pipe (`| head`) is not an error.
fn is_broken_pipe(e: &Error) -> bool {
    matches!(e, Error::Io { source, .. } if source.kind() == std::io::ErrorKind::BrokenPipe)
}

pub(super) fn fail(prog: &str, e: &Error) -> i32 {
    if is_broken_pipe(e) {
        return 0;
    }
    match e {
        Error::Cancelled => {
            eprintln!("{prog}: cancelled; partial output removed");
            130
        }
        Error::EngineNotFound { hint } => {
            eprintln!("{hint}");
            2
        }
        Error::InvalidOptions(_) => usage(prog, e),
        _ => {
            eprintln!("{prog}: error: {e}");
            if let Some(h) = e.hint() {
                eprintln!("  hint: {h}");
            }
            1
        }
    }
}

/// Progress lines: engine lines, warnings, search progress; the heartbeat is silent.
pub fn printer(quiet: bool) -> impl FnMut(&ProgressEvent) -> bool {
    move |e: &ProgressEvent| {
        match e {
            ProgressEvent::EngineStderr { line } => eprintln!("{line}"),
            ProgressEvent::EngineRunning { .. }
            | ProgressEvent::Opened { .. }
            | ProgressEvent::Finishing { .. } => {}
            ProgressEvent::Warning(_) => eprintln!("{}", describe(e)),
            _ if !quiet => eprintln!("{}", describe(e)),
            _ => {}
        }
        !INTERRUPTED.load(Ordering::SeqCst)
    }
}

fn write_out(path: Option<&std::path::Path>, text: &str) -> Result<(), Error> {
    match path {
        Some(p) => std::fs::write(p, text).map_err(|e| Error::io(p, e)),
        None => {
            let mut o = std::io::stdout().lock();
            o.write_all(text.as_bytes())
                .and_then(|_| o.flush())
                .map_err(|e| Error::io("<stdout>", e))
        }
    }
}

/// The summary in `format`.
pub(crate) fn report_text(r: &Report, format: FormatArg) -> String {
    match format {
        FormatArg::Tsv => output::tsv(r),
        FormatArg::Json => output::json(r),
    }
}

/// Write the summary and detail table as asked.
pub fn render_report(prog: &str, r: &Report, out: &ReportOutputArgs) -> Result<(), Error> {
    write_out(out.output.as_deref(), &report_text(r, out.format))?;
    if let Some(p) = &out.detail {
        write_out(Some(p), &output::detail(r, out.detail_top))?;
    }
    for w in &r.meta.warnings {
        eprintln!("{prog}: warning: {w}");
    }
    if r.meta.taxonomy.is_none() && r.queries.iter().any(|q| q.status == "multi") {
        eprintln!("{prog}: some queries tie across taxa ('multi'); {TAXDUMP_HINT}");
    }
    Ok(())
}

pub fn run_search(prog: &str, a: &SearchArgs) -> i32 {
    let req = match SearchRequest::try_from(a) {
        Ok(r) => r,
        Err(e) => return usage(prog, &e),
    };
    match api::search(&req, &mut printer(a.quiet)) {
        Ok(o) => o.exit_code,
        Err(e) => fail(prog, &e),
    }
}

pub fn run_report(prog: &str, a: &ReportArgs) -> i32 {
    let req = match ReportRequest::try_from(a) {
        Ok(r) => r,
        Err(e) => return usage(prog, &e),
    };
    match api::report(&req).and_then(|r| render_report(prog, &r, &a.out)) {
        Ok(()) => 0,
        Err(e) => fail(prog, &e),
    }
}

pub fn run_origin(prog: &str, a: &OriginArgs) -> i32 {
    let req = match OriginRequest::try_from(a) {
        Ok(r) => r,
        Err(e) => return usage(prog, &e),
    };
    match api::origin(&req, &mut printer(a.quiet)) {
        Ok(o) => {
            for n in &o.meta.notes {
                eprintln!("{n}");
            }
            match render_report(prog, &o.report, &a.out) {
                Ok(()) => 0,
                Err(e) => fail(prog, &e),
            }
        }
        Err(e) => fail(prog, &e),
    }
}

pub fn run_reps(prog: &str, a: &RepsArgs) -> i32 {
    let req = match RepsRequest::try_from(a) {
        Ok(r) => r,
        Err(e) => return usage(prog, &e),
    };
    match api::reps(&req).and_then(|t| write_out(None, &t.render())) {
        Ok(()) => 0,
        Err(e) => fail(prog, &e),
    }
}

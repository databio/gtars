//! `run_json`: run one verb from an argv and return its result as JSON. The Python and R
//! bindings call this after turning their keyword options into flags, so every front end
//! shares the clap parsers and their error messages.

use super::{BlastCommand, ReportOutputArgs};
use crate::api::{self, OriginRequest, ReportRequest, RepsRequest, SearchRequest, ViewRequest};
use crate::error::Error;
use crate::progress::ProgressEvent;
use crate::report::{Report, output};
use crate::run::EngineOutput;
use crate::run::passthrough::value;
use crate::{BuildRequest, NstatsRequest, VerifyRequest};
use clap::Parser;
use serde_json::{Map, Value};

#[derive(Parser)]
#[command(name = "blast", no_binary_name = true, subcommand_required = true)]
struct Argv {
    #[command(subcommand)]
    cmd: BlastCommand,
}

fn json<T: serde::Serialize>(v: &T) -> Result<Value, Error> {
    serde_json::to_value(v).map_err(|e| Error::Store(anyhow::anyhow!("serializing result: {e}")))
}

/// Run one verb. `argv` starts with the verb (`["build", "--store", "s/", ...]`), as after
/// `gtars blast`. A help request (`--help`) returns the help text as a JSON string; any other
/// parse error is `Error::InvalidOptions` with clap's message.
///
/// Results: the verb's summary (`build`, `verify` with `ok`, `view`, `nstats`, `reps`);
/// `build --dry-run` the preflight summary; `search` the run metadata plus `output` (the
/// engine's output when no `-out` was given) and `engine_stderr`; `report` `{summary, report}`
/// and `origin` `{summary, report, meta}`, after writing any `--output` and `--detail` files.
pub fn run_json(
    argv: &[String],
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<Value, Error> {
    let cmd = match Argv::try_parse_from(argv) {
        Ok(a) => a.cmd,
        Err(e) => {
            use clap::error::ErrorKind as K;
            let text = e.render().to_string();
            return match e.kind() {
                K::DisplayHelp | K::DisplayVersion => Ok(Value::String(text)),
                _ => Err(Error::InvalidOptions(
                    text.trim_start_matches("error: ").trim_end().to_string(),
                )),
            };
        }
    };
    match cmd {
        BlastCommand::Build(a) => {
            let req = BuildRequest::try_from(&a)?;
            if a.dry_run {
                json(&api::preflight(&req)?)
            } else {
                json(&api::build(&req, progress)?)
            }
        }
        BlastCommand::Verify(a) => {
            let r = api::verify(&VerifyRequest::try_from(&a)?, progress)?;
            let mut j = json(&r)?;
            if let Some(m) = j.as_object_mut() {
                m.insert("ok".into(), r.ok().into());
            }
            Ok(j)
        }
        BlastCommand::View(a) => json(&api::view(&ViewRequest::try_from(&a)?, progress)?),
        BlastCommand::Nstats(a) => json(&api::nstats(&NstatsRequest::try_from(&a)?, progress)?),
        BlastCommand::Search(a) => {
            let mut req = SearchRequest::try_from(&a)?;
            if value(&a.blastn_args, "-out").is_none() {
                req.output = EngineOutput::Capture;
            }
            let o = api::search(&req, progress)?;
            let mut j = json(&o)?;
            if let Some(m) = j.as_object_mut() {
                m.insert("output".into(), json(&o.output)?);
                m.insert("engine_stderr".into(), json(&o.engine_stderr)?);
            }
            Ok(j)
        }
        BlastCommand::Report(a) => {
            let r = api::report(&ReportRequest::try_from(&a)?)?;
            write_report(&r, &a.out)?;
            Ok(report_json(&r)?.into())
        }
        BlastCommand::Origin(a) => {
            let o = api::origin(&OriginRequest::try_from(a.as_ref())?, progress)?;
            write_report(&o.report, &a.out)?;
            let mut m = report_json(&o.report)?;
            m.insert("meta".into(), json(&o.meta)?);
            Ok(m.into())
        }
        BlastCommand::Reps(a) => json(&api::reps(&RepsRequest::try_from(&a)?)?),
    }
}

/// Run one verb from an argv as the command line does (messages on stderr, the summary or
/// report on stdout, clap's help and errors) and return the exit code; `prog` prefixes the
/// messages. The refget CLI passes `refget blast` through to this.
pub fn main_argv(prog: &str, argv: &[String]) -> i32 {
    match Argv::try_parse_from(argv) {
        Ok(a) => a.cmd.run_code(prog),
        Err(e) => {
            let _ = e.print();
            e.exit_code()
        }
    }
}

fn report_json(r: &Report) -> Result<Map<String, Value>, Error> {
    let mut m = Map::new();
    m.insert("summary".into(), json(&r.summary_rows())?);
    m.insert("report".into(), json(r)?);
    Ok(m)
}

/// The files a report or origin run asked for; the summary goes nowhere without `--output`
/// (it is in the returned JSON).
fn write_report(r: &Report, out: &ReportOutputArgs) -> Result<(), Error> {
    let write =
        |p: &std::path::Path, text: String| std::fs::write(p, text).map_err(|e| Error::io(p, e));
    if let Some(p) = &out.output {
        write(p, super::rgblast::report_text(r, out.format))?;
    }
    if let Some(p) = &out.detail {
        write(p, output::detail(r, out.detail_top))?;
    }
    Ok(())
}

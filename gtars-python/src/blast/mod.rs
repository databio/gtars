//! `gtars.blast`: BLAST over a RefgetStore. `_run(argv, progress)` runs one `gtars blast`
//! verb and returns its result as Python objects; the wrappers in `gtars/blast/__init__.py`
//! turn keyword options into flags. `_main(argv)` runs a verb as the command line does (the
//! refget CLI passes `refget blast` through to it).

mod convert;
mod errors;

use convert::json_to_py;
use errors::to_pyerr;
use gtars_rgblast::{Error, ProgressEvent};
use pyo3::prelude::*;

/// What `progress=` asked for.
enum Progress {
    Off,
    Print,
    Call(Py<PyAny>),
}

impl Progress {
    fn from_py(obj: Option<Bound<'_, PyAny>>) -> PyResult<Progress> {
        match obj {
            None => Ok(Progress::Off),
            Some(o) if o.is_none() => Ok(Progress::Off),
            Some(o) => {
                if let Ok(b) = o.extract::<bool>() {
                    if o.is_instance_of::<pyo3::types::PyBool>() {
                        return Ok(if b { Progress::Print } else { Progress::Off });
                    }
                }
                if o.is_callable() {
                    Ok(Progress::Call(o.unbind()))
                } else {
                    Err(errors::invalid(
                        o.py(),
                        "progress must be None, True or a callable",
                    ))
                }
            }
        }
    }

    /// Runs on the calling thread with the GIL taken again: check for Ctrl-C, then print the
    /// event's message or call the callback with the event dict.
    fn tick(&self, py: Python<'_>, ev: &ProgressEvent) -> PyResult<()> {
        py.check_signals()?;
        match self {
            Progress::Off => Ok(()),
            Progress::Print => {
                let line = format!("{}\n", gtars_rgblast::describe(ev));
                py.import("sys")?
                    .getattr("stderr")?
                    .call_method1("write", (line,))?;
                Ok(())
            }
            Progress::Call(cb) => {
                cb.call1(py, (json_to_py(py, &gtars_rgblast::event_json(ev))?,))?;
                Ok(())
            }
        }
    }
}

/// Run a library call with the GIL released. The library calls the callback on this thread;
/// the callback takes the GIL back, checks signals and runs the Python progress. A Python
/// error (or Ctrl-C) is stored and cancels the call; it is raised afterwards, even when the
/// library finished before it saw the cancel.
fn run<T: Send>(
    py: Python<'_>,
    progress: Option<Bound<'_, PyAny>>,
    f: impl FnOnce(&mut dyn FnMut(&ProgressEvent) -> bool) -> Result<T, Error> + Send,
) -> PyResult<T> {
    let progress = Progress::from_py(progress)?;
    let (result, stored) = py.detach(|| {
        let mut stored: Option<PyErr> = None;
        let mut cb = |ev: &ProgressEvent| -> bool {
            Python::attach(|py| match progress.tick(py, ev) {
                Ok(()) => true,
                Err(e) => {
                    if stored.is_none() {
                        stored = Some(e);
                    }
                    false
                }
            })
        };
        let r = f(&mut cb);
        (r, stored)
    });
    if let Some(e) = stored {
        return Err(e);
    }
    result.map_err(|e| to_pyerr(py, e, None))
}

/// Run one verb (`argv` starts with the verb) and return its result. A non-zero engine exit
/// status raises `EngineError`.
#[pyfunction]
#[pyo3(signature = (argv, progress=None))]
fn _run(
    py: Python<'_>,
    argv: Vec<String>,
    progress: Option<Bound<'_, PyAny>>,
) -> PyResult<Py<PyAny>> {
    let v = run(py, progress, |cb| gtars_rgblast::cli::run_json(&argv, cb))?;
    if let Some(code) = v
        .get("exit_code")
        .and_then(|c| c.as_i64())
        .filter(|c| *c != 0)
    {
        let tail: Vec<&str> = v
            .get("engine_stderr")
            .and_then(|s| s.as_array())
            .map(|a| a.iter().filter_map(|l| l.as_str()).collect())
            .unwrap_or_default();
        let mut msg = format!("the search engine exited with status {code}");
        if !tail.is_empty() {
            msg = format!("{msg}:\n{}", tail.join("\n"));
        }
        return Err(errors::engine_error(py, code as i32, msg));
    }
    json_to_py(py, &v)
}

/// Run one verb as `gtars blast` does and return its exit code.
#[pyfunction]
fn _main(py: Python<'_>, prog: String, argv: Vec<String>) -> i32 {
    py.detach(|| gtars_rgblast::cli::main_argv(&prog, &argv))
}

#[pymodule]
pub fn blast(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_run, m)?)?;
    m.add_function(wrap_pyfunction!(_main, m)?)?;
    m.add(
        "REPORT_FIELDS",
        gtars_rgblast::run::outfmt::REPORT_FIELDS.to_vec(),
    )?;
    m.add(
        "SUMMARY_COLUMNS",
        gtars_rgblast::report::rank::SUMMARY_COLUMNS.to_vec(),
    )?;
    m.add(
        "EngineNotFoundError",
        py.get_type::<errors::EngineNotFoundError>(),
    )?;
    m.add("EngineError", py.get_type::<errors::EngineError>())?;
    Ok(())
}

//! `gtars_rgblast::Error` to Python exceptions, by `ErrorKind`. Every exception gets a `kind`
//! attribute ("not_found", "invalid_input", ...); messages are the library's `Display` text.

use gtars_rgblast::{Error, ErrorKind};
use pyo3::create_exception;
use pyo3::exceptions::*;
use pyo3::prelude::*;

create_exception!(
    gtars.blast,
    EngineNotFoundError,
    PyFileNotFoundError,
    "No BLAST engine was found; the message says how to install one."
);
create_exception!(
    gtars.blast,
    EngineError,
    PyRuntimeError,
    "The engine exited with a non-zero status (`returncode`)."
);

pub fn kind_name(k: ErrorKind) -> &'static str {
    match k {
        ErrorKind::NotFound => "not_found",
        ErrorKind::InvalidInput => "invalid_input",
        ErrorKind::Conflict => "conflict",
        ErrorKind::Io => "io",
        ErrorKind::Network => "network",
        ErrorKind::Data => "data",
        ErrorKind::Cancelled => "cancelled",
    }
}

fn tag(py: Python<'_>, e: PyErr, kind: &str) -> PyErr {
    let _ = e.value(py).setattr("kind", kind);
    e
}

/// Map a library error. `stored` is the error a progress callback raised (or Ctrl-C); it wins
/// whenever the build was cancelled.
pub fn to_pyerr(py: Python<'_>, e: Error, stored: Option<PyErr>) -> PyErr {
    if let Some(s) = stored {
        return s;
    }
    let kind = e.kind();
    let name = kind_name(kind);
    let msg = e.to_string();
    let err = match (&e, kind) {
        (Error::EngineNotFound { hint }, _) => EngineNotFoundError::new_err(hint.clone()),
        (_, ErrorKind::NotFound) => PyFileNotFoundError::new_err(msg),
        (_, ErrorKind::InvalidInput) => PyValueError::new_err(msg),
        (Error::OutputExists { .. }, _) => PyFileExistsError::new_err(msg),
        (_, ErrorKind::Conflict) => PyValueError::new_err(msg),
        (Error::InsufficientSpace { .. }, _) => PyOSError::new_err((libc_enospc(), msg)),
        (Error::Io { source, .. }, _) if source.kind() == std::io::ErrorKind::BrokenPipe => {
            PyBrokenPipeError::new_err(msg)
        }
        (_, ErrorKind::Io) => PyOSError::new_err(msg),
        (_, ErrorKind::Network) => PyConnectionError::new_err(msg),
        (_, ErrorKind::Data) => PyRuntimeError::new_err(msg),
        (_, ErrorKind::Cancelled) => PyKeyboardInterrupt::new_err("cancelled"),
    };
    tag(py, err, name)
}

fn libc_enospc() -> i32 {
    28
}

/// A non-zero engine exit status (used by `search` and `origin`).
pub fn engine_error(py: Python<'_>, returncode: i32, msg: String) -> PyErr {
    let e = EngineError::new_err(msg);
    let _ = e.value(py).setattr("returncode", returncode);
    tag(py, e, "engine_failed")
}

/// A binding-side argument problem (the same class a library `InvalidOptions` gives).
pub fn invalid(py: Python<'_>, msg: impl Into<String>) -> PyErr {
    tag(py, PyValueError::new_err(msg.into()), "invalid_input")
}

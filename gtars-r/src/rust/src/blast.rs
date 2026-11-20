//! `blast_*`: BLAST over a RefgetStore. The R wrappers turn their options into a `gtars blast`
//! argv; `blast_start_rs` runs it on a worker thread (`gtars_rgblast::cli::run_json`) and R
//! drains its progress events on the main thread (`blast_poll_rs`). Results come back as
//! `list(ok = TRUE, value = ...)` or `list(ok = FALSE, class = c(...), message = "...")`; the
//! R side raises the condition.

use extendr_api::prelude::*;
use gtars_rgblast::{Error, ErrorKind, ProgressEvent};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{sync_channel, Receiver, RecvTimeoutError};
use std::sync::{Arc, Mutex};
use std::thread::JoinHandle;
use std::time::Duration;

// ---------------------------------------------------------------- results and conditions

fn ok(value: Robj) -> Robj {
    list!(ok = true, value = value).into()
}

fn fail(classes: &[&str], message: &str) -> Robj {
    list!(ok = false, class = classes.to_vec(), message = message).into()
}

fn kind_class(k: ErrorKind) -> &'static str {
    match k {
        ErrorKind::NotFound => "gtars_rgblast_not_found",
        ErrorKind::InvalidInput => "gtars_rgblast_invalid_input",
        ErrorKind::Conflict => "gtars_rgblast_conflict",
        ErrorKind::Io => "gtars_rgblast_io",
        ErrorKind::Network => "gtars_rgblast_network",
        ErrorKind::Data => "gtars_rgblast_data",
        ErrorKind::Cancelled => "gtars_rgblast_cancelled",
    }
}

/// A library error as a condition description. The classes are `kind_class`, plus a specific
/// one for a missing engine.
fn error_robj(e: &Error) -> Robj {
    let mut classes = vec![kind_class(e.kind())];
    if let Error::EngineNotFound { .. } = e {
        classes.push("gtars_rgblast_engine_not_found");
    }
    let msg = match e {
        Error::EngineNotFound { hint } => hint.clone(),
        _ => e.to_string(),
    };
    fail(&classes, &msg)
}

// ---------------------------------------------------------------- JSON to R

/// null becomes NULL, every number a double (counts such as 41.8 Gbp overflow R integers), a
/// homogeneous scalar array an atomic vector and an object a named list.
pub(crate) fn json_to_robj(v: &serde_json::Value) -> Robj {
    use serde_json::Value as V;
    match v {
        V::Null => ().into(),
        V::Bool(b) => (*b).into(),
        V::Number(n) => n.as_f64().unwrap_or(f64::NAN).into(),
        V::String(s) => s.as_str().into(),
        V::Array(a) => {
            if !a.is_empty() && a.iter().all(|x| x.is_number()) {
                let v: Vec<f64> = a.iter().map(|x| x.as_f64().unwrap_or(f64::NAN)).collect();
                v.into()
            } else if !a.is_empty() && a.iter().all(|x| x.is_string()) {
                let v: Vec<String> = a
                    .iter()
                    .map(|x| x.as_str().unwrap_or("").to_string())
                    .collect();
                v.into()
            } else if !a.is_empty() && a.iter().all(|x| x.is_boolean()) {
                let v: Vec<Rbool> = a
                    .iter()
                    .map(|x| Rbool::from(x.as_bool().unwrap_or(false)))
                    .collect();
                v.into()
            } else {
                List::from_values(a.iter().map(json_to_robj)).into()
            }
        }
        V::Object(o) => {
            let names: Vec<&str> = o.keys().map(|k| k.as_str()).collect();
            let vals: Vec<Robj> = o.values().map(json_to_robj).collect();
            List::from_names_and_values(names, vals).unwrap().into()
        }
    }
}

// ---------------------------------------------------------------- jobs

/// A library call, ready to run on the worker thread.
type Task = Box<
    dyn FnOnce(&mut dyn FnMut(&ProgressEvent) -> bool) -> Result<serde_json::Value, Error> + Send,
>;

/// What the worker sends back: progress events, then the outcome.
enum Msg {
    Event(Box<ProgressEvent>),
    Done(Result<serde_json::Value, Error>),
}

/// A running library call. The call runs on its own thread and R drains its events from the
/// main thread (`blast_poll_rs`), so R functions run in an ordinary R frame: condition
/// handlers such as `suppressMessages()` and `tryCatch()` see what they raise, and an R error
/// or Ctrl-C in the progress function cancels the job (`blast_cancel_rs`).
struct Job {
    /// Dropped on cancel, so a worker blocked on a full channel wakes up and stops.
    rx: Mutex<Option<Receiver<Msg>>>,
    cancel: Arc<AtomicBool>,
    thread: Mutex<Option<JoinHandle<()>>>,
}

fn start(task: Task) -> Robj {
    // A one-slot channel keeps the worker at most an event ahead of R, so a cancel arrives
    // while the work is still running.
    let (tx, rx) = sync_channel::<Msg>(1);
    let cancel = Arc::new(AtomicBool::new(false));
    let flag = Arc::clone(&cancel);
    let thread = std::thread::spawn(move || {
        let tx_events = tx.clone();
        let mut cb = move |ev: &ProgressEvent| -> bool {
            tx_events.send(Msg::Event(Box::new(ev.clone()))).is_ok() && !flag.load(Ordering::SeqCst)
        };
        let r = task(&mut cb);
        let _ = tx.send(Msg::Done(r));
    });
    let job = Job {
        rx: Mutex::new(Some(rx)),
        cancel,
        thread: Mutex::new(Some(thread)),
    };
    ok(ExternalPtr::new(job).into())
}

fn done_robj(r: Result<serde_json::Value, Error>) -> Robj {
    match r {
        Ok(v) => ok(json_to_robj(&v)),
        Err(e) => error_robj(&e),
    }
}

/// Wait up to `wait_ms` for the job's next message: list(kind = "event", event = <list>),
/// list(kind = "done", result = <ok/fail list>) or list(kind = "none").
/// @noRd
#[extendr]
pub fn blast_poll_rs(job: Robj, wait_ms: f64) -> Robj {
    let Ok(job) = ExternalPtr::<Job>::try_from(&job) else {
        return list!(
            kind = "done",
            result = fail(&["gtars_rgblast_data"], "not a job")
        )
        .into();
    };
    let wait = Duration::from_millis(wait_ms.max(0.0) as u64);
    let msg = match job.rx.lock().unwrap().as_ref() {
        Some(rx) => rx.recv_timeout(wait),
        None => Err(RecvTimeoutError::Disconnected),
    };
    match msg {
        Ok(Msg::Event(ev)) => list!(
            kind = "event",
            event = json_to_robj(&gtars_rgblast::event_json(&ev))
        )
        .into(),
        Ok(Msg::Done(r)) => {
            if let Some(t) = job.thread.lock().unwrap().take() {
                let _ = t.join();
            }
            list!(kind = "done", result = done_robj(r)).into()
        }
        Err(RecvTimeoutError::Timeout) => list!(kind = "none").into(),
        Err(RecvTimeoutError::Disconnected) => list!(
            kind = "done",
            result = fail(
                &["gtars_rgblast_data"],
                "the BLAST worker stopped without a result"
            )
        )
        .into(),
    }
}

/// Cancel a job and wait for it to clean up (the library removes partial output).
/// @noRd
#[extendr]
pub fn blast_cancel_rs(job: Robj) {
    if let Ok(job) = ExternalPtr::<Job>::try_from(&job) {
        job.cancel.store(true, Ordering::SeqCst);
        drop(job.rx.lock().unwrap().take());
        if let Some(t) = job.thread.lock().unwrap().take() {
            let _ = t.join();
        }
    }
}

/// The columns of the TSV summary, in order.
/// @noRd
#[extendr]
pub fn blast_summary_columns_rs() -> Vec<String> {
    gtars_rgblast::report::rank::SUMMARY_COLUMNS
        .iter()
        .map(|s| s.to_string())
        .collect()
}

/// The blastn `-outfmt 6/7` fields that `blast_search()` appends for `blast_report()`.
/// @noRd
#[extendr]
pub fn blast_report_fields_rs() -> Vec<String> {
    gtars_rgblast::run::outfmt::REPORT_FIELDS
        .iter()
        .map(|s| s.to_string())
        .collect()
}

/// Start one verb on a worker thread. `argv` starts with the verb, as after `gtars blast`.
/// Returns list(ok = TRUE, value = <job>); option errors come back from `blast_poll_rs`.
/// @param argv The verb and its flags
/// @noRd
#[extendr]
pub fn blast_start_rs(argv: Vec<String>) -> Robj {
    start(Box::new(move |cb| gtars_rgblast::cli::run_json(&argv, cb)))
}

extendr_module! {
    mod blast;
    fn blast_start_rs;
    fn blast_poll_rs;
    fn blast_cancel_rs;
    fn blast_report_fields_rs;
    fn blast_summary_columns_rs;
}

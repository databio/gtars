//! The blastn runner (`--engine blastn --db`) with a fake blastn (a committed shell script that
//! sleeps): the heartbeat, cancel, a non-zero exit returned as a value, and a missing engine.

mod common;

use gtars_rgblast::api::{SearchRequest, search};
use gtars_rgblast::run::{EngineChoice, EngineOutput};
use gtars_rgblast::select::Selection;
use gtars_rgblast::{Error, ProgressEvent};
use std::path::{Path, PathBuf};
use std::time::Instant;

fn fake() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/fixtures/origin/fake_blastn.sh")
}

fn req(dir: &Path, blastn: PathBuf) -> SearchRequest {
    SearchRequest {
        selection: Selection::all(),
        engine: EngineChoice::Blastn {
            db: dir.join("db"),
            blastn: Some(blastn),
        },
        extra_fields: true,
        args: vec!["-query".into(), "q.fa".into(), "-outfmt".into(), "6".into()],
        output: EngineOutput::Capture,
        ..Default::default()
    }
}

#[test]
fn heartbeat_and_nonzero_exit() {
    let t = common::tempdir();
    let mut beats = 0usize;
    let mut stderr = Vec::new();
    let mut started = None;
    let o = search(&req(t.path(), fake()), &mut |e| {
        match e {
            ProgressEvent::EngineStart { args, .. } => {
                started = Some(Instant::now());
                assert!(
                    args.contains(&"6 std qlen score sallseqid salltitles staxids".to_string())
                );
            }
            ProgressEvent::EngineRunning { .. } => beats += 1,
            ProgressEvent::EngineStderr { line } => stderr.push(line.clone()),
            _ => {}
        }
        true
    })
    .unwrap();
    let secs = started.unwrap().elapsed().as_secs_f64();
    assert!(secs >= 4.5, "the fake sleeps 5 s ({secs})");
    assert!(
        beats as f64 >= secs.floor(),
        "{beats} heartbeats in {secs:.1} s"
    );
    assert_eq!(o.exit_code, 3, "a non-zero exit is a value, not an error");
    assert_eq!(stderr, ["fake blastn: sleeping"]);
    assert_eq!(o.engine_stderr, stderr);
    assert_eq!(o.db.unwrap(), t.path().join("db"));
    assert!(o.store.is_none());
}

#[test]
fn cancel_kills_within_a_second() {
    let t = common::tempdir();
    let mut cancelled_at = None;
    let r = search(&req(t.path(), fake()), &mut |e| {
        if matches!(e, ProgressEvent::EngineRunning { .. }) {
            cancelled_at = Some(Instant::now());
            return false;
        }
        true
    });
    assert!(matches!(r, Err(Error::Cancelled)), "{r:?}");
    assert!(cancelled_at.unwrap().elapsed().as_secs_f64() < 1.0);
}

#[test]
fn missing_engine_has_the_install_hint() {
    let t = common::tempdir();
    let e = search(&req(t.path(), "/nonexistent/blastn".into()), &mut |_| true).unwrap_err();
    assert!(matches!(e, Error::EngineNotFound { .. }));
    assert!(
        e.to_string().contains("Install BLAST+ 2.10 or newer"),
        "{e}"
    );
    assert_eq!(e.kind(), gtars_rgblast::ErrorKind::NotFound);
}

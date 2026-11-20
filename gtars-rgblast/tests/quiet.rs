//! The `rgblast` library calls print nothing: stdout and stderr (the file descriptors) are
//! redirected to files around `search`, `report`, `origin` and `reps`. One
//! test in its own binary, so the harness prints nothing meanwhile.

mod common;

use common::origin::{ORIGIN_FIXTURES, origin_store};
use gtars_rgblast::api::{self, OriginRequest, ReportRequest, RepsRequest, SearchRequest};
use gtars_rgblast::origin::{ModeChoice, OriginOptions};
use gtars_rgblast::report::output;
use gtars_rgblast::report::parse::DEFAULT_FIELDS;
use gtars_rgblast::run::{EngineChoice, EngineOutput};
use gtars_rgblast::select::Selection;
use std::path::Path;

struct Redirect {
    saved: [i32; 2],
}

impl Redirect {
    fn to(dir: &Path) -> Redirect {
        let mut saved = [0; 2];
        for (i, fd) in [1, 2].into_iter().enumerate() {
            let f = std::fs::File::create(dir.join(format!("fd{fd}"))).unwrap();
            unsafe {
                saved[i] = libc::dup(fd);
                libc::dup2(std::os::fd::AsRawFd::as_raw_fd(&f), fd);
            }
        }
        Redirect { saved }
    }
}

impl Drop for Redirect {
    fn drop(&mut self) {
        for (i, fd) in [1, 2].into_iter().enumerate() {
            unsafe {
                libc::dup2(self.saved[i], fd);
                libc::close(self.saved[i]);
            }
        }
    }
}

#[test]
fn library_calls_print_nothing() {
    let f = origin_store();
    let t = common::tempdir();
    let cache = t.path().join("cache");
    let fake =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/fixtures/origin/fake_blastn.sh");
    unsafe {
        std::env::set_var("FAKE_BLASTN_SLEEP", "0");
        std::env::set_var("FAKE_BLASTN_EXIT", "0");
    }
    let engine = EngineChoice::Blastn {
        db: t.path().join("db"),
        blastn: Some(fake),
    };
    let out_dir = t.path().join("out");
    std::fs::create_dir_all(&out_dir).unwrap();
    {
        let _r = Redirect::to(&out_dir);
        let s = api::search(
            &SearchRequest {
                engine,
                extra_fields: true,
                args: vec![
                    "-query".into(),
                    f.queries.display().to_string(),
                    "-outfmt".into(),
                    "6".into(),
                ],
                output: EngineOutput::Capture,
                ..Default::default()
            },
            &mut |_| true,
        )
        .unwrap();
        assert_eq!(s.exit_code, 0);
        let r = api::report(&ReportRequest {
            hits: vec![Path::new(ORIGIN_FIXTURES).join("small_jungle_hits.tsv")],
            db: Some(Path::new(ORIGIN_FIXTURES).join("small_jungle")),
            fields: DEFAULT_FIELDS.into(),
            ..Default::default()
        })
        .unwrap();
        let o = api::origin(
            &OriginRequest {
                store: f.store.display().to_string(),
                selection: Selection::all(),
                engine: EngineChoice::Rgblastn,
                cache_dir: Some(cache.clone()),
                query: f.queries.clone(),
                args: Vec::new(),
                reps: Vec::new(),
                options: OriginOptions {
                    mode: ModeChoice::Staged,
                    ..Default::default()
                },
                taxdump: None,
                keep_hits: None,
                fetch: Default::default(),
            },
            &mut |_| true,
        )
        .unwrap();
        let reps = api::reps(&RepsRequest {
            store: f.store.display().to_string(),
            selection: Selection::all(),
            rep_namespace: "reference".into(),
            cache_dir: Some(cache),
            ..Default::default()
        })
        .unwrap();
        // checked after the redirect ends
        drop(_r);
        assert_eq!(r.queries.len(), 4);
        assert_eq!(o.report.queries.len(), 8);
        assert_eq!(reps.groups.len(), 3);
        assert!(output::tsv(&r).contains("\t37\t1\t79/79"));
    }
    for fd in [1, 2] {
        let got = std::fs::read_to_string(out_dir.join(format!("fd{fd}"))).unwrap();
        assert_eq!(got, "", "fd {fd} got output");
    }
}

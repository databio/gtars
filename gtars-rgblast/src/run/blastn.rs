//! Stock blastn (`--engine blastn`) on a database the user built, searched as given.

use super::{BlastRunner, EngineInfo, EngineOutput, RunOutcome, RunSpec, SearchTarget};
use crate::error::Error;
use crate::progress::ProgressEvent;
use crate::store::OpenedStore;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::mpsc;
use std::time::{Duration, Instant};

pub const TESTED_VERSION: &str = "2.17.0";
pub const NOT_FOUND_HINT: &str = "rgblast: blastn not found. Install BLAST+ 2.10 or newer (conda install -c bioconda blast, or apt install ncbi-blast+) or give its path with --blastn PATH; or search the store with the default engine, rgblastn.";

pub struct BlastnRunner {
    pub path: PathBuf,
    /// "2.17.0+"
    pub version: String,
    note: std::sync::Mutex<Option<String>>,
}

fn on_path(prog: &str) -> Option<PathBuf> {
    let path = std::env::var_os("PATH")?;
    std::env::split_paths(&path)
        .map(|d| d.join(prog))
        .find(|p| p.is_file())
}

/// "blastn: 2.17.0+" -> (2, 17, "2.17.0+").
pub fn parse_version(out: &str) -> Option<(u32, u32, String)> {
    let line = out
        .lines()
        .find(|l| l.trim_start().starts_with("blastn:"))?;
    let v = line.split_whitespace().nth(1)?.to_string();
    let mut it = v.trim_end_matches('+').split('.');
    let major = it.next()?.parse().ok()?;
    let minor = it.next()?.parse().ok()?;
    Some((major, minor, v))
}

fn not_found() -> Error {
    Error::EngineNotFound {
        hint: NOT_FOUND_HINT.into(),
    }
}

impl BlastnRunner {
    /// `--blastn PATH`, else `$GTARS_RGBLAST_BLASTN`, else `blastn` on PATH; then the version
    /// check (2.10 or newer).
    pub fn resolve(blastn: Option<&Path>) -> Result<Self, Error> {
        let path = match blastn {
            Some(p) => p.to_path_buf(),
            None => match std::env::var("GTARS_RGBLAST_BLASTN") {
                Ok(p) if !p.is_empty() => PathBuf::from(p),
                _ => on_path("blastn").ok_or_else(not_found)?,
            },
        };
        let out = match Command::new(&path)
            .arg("-version")
            .stdin(Stdio::null())
            .output()
        {
            Ok(o) => o,
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Err(not_found()),
            Err(e) => return Err(Error::Engine(format!("could not run blastn -version: {e}"))),
        };
        let text = String::from_utf8_lossy(&out.stdout);
        let Some((major, minor, version)) = parse_version(&text) else {
            return Err(Error::Engine(format!(
                "{}: 'blastn -version' did not print a blastn version (exit {}): {}",
                path.display(),
                out.status,
                String::from_utf8_lossy(&out.stderr).trim()
            )));
        };
        if (major, minor) < (2, 10) {
            return Err(Error::Engine(format!(
                "blastn {version} is too old: BLAST version 5 databases need 2.10 or newer"
            )));
        }
        let note = (version.trim_end_matches('+') != TESTED_VERSION)
            .then(|| format!("rgblast: tested with blastn {TESTED_VERSION}; this is {version}"));
        Ok(BlastnRunner {
            path,
            version,
            note: std::sync::Mutex::new(note),
        })
    }
}

/// The result of one engine process.
pub struct Executed {
    pub exit_code: i32,
    pub output: Option<String>,
    pub stderr_tail: Vec<String>,
}

/// Run `cmd`, forwarding stderr lines and a heartbeat (every 500 ms) to `progress`. Returning
/// `false` from it kills the process and gives `Error::Cancelled`.
pub fn execute(
    mut cmd: Command,
    output: &EngineOutput,
    stdin: bool,
    stage: &str,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<Executed, Error> {
    cmd.stdin(if stdin {
        Stdio::inherit()
    } else {
        Stdio::null()
    });
    cmd.stderr(Stdio::piped());
    match output {
        EngineOutput::Inherit => {
            cmd.stdout(Stdio::inherit());
        }
        EngineOutput::File(p) => {
            let f = std::fs::File::create(p).map_err(|e| Error::io(p, e))?;
            cmd.stdout(Stdio::from(f));
        }
        EngineOutput::Capture => {
            cmd.stdout(Stdio::piped());
        }
    }
    let mut child = cmd.spawn().map_err(|e| {
        if e.kind() == std::io::ErrorKind::NotFound {
            not_found()
        } else {
            Error::Engine(format!("could not start the engine: {e}"))
        }
    })?;
    let (tx, rx) = mpsc::channel::<String>();
    let err = child.stderr.take().unwrap();
    let err_thread = std::thread::spawn(move || {
        for l in BufReader::new(err).lines().map_while(Result::ok) {
            if tx.send(l).is_err() {
                break;
            }
        }
    });
    let out_thread = child.stdout.take().map(|mut o| {
        std::thread::spawn(move || {
            let mut s = String::new();
            let _ = o.read_to_string(&mut s);
            s
        })
    });
    let t0 = Instant::now();
    let mut last_beat = Instant::now();
    let mut tail: Vec<String> = Vec::new();
    let mut push = |l: String, progress: &mut dyn FnMut(&ProgressEvent) -> bool| -> bool {
        let ok = progress(&ProgressEvent::EngineStderr { line: l.clone() });
        tail.push(l);
        if tail.len() > 50 {
            tail.remove(0);
        }
        ok
    };
    let status = loop {
        let mut cancel = false;
        while let Ok(l) = rx.try_recv() {
            cancel |= !push(l, progress);
        }
        if let Some(st) = child.try_wait().map_err(|e| Error::Engine(e.to_string()))? {
            break st;
        }
        if last_beat.elapsed() >= Duration::from_millis(500) {
            last_beat = Instant::now();
            cancel |= !progress(&ProgressEvent::EngineRunning {
                stage: stage.to_string(),
                elapsed_s: t0.elapsed().as_secs_f64(),
            });
        }
        if cancel {
            let _ = child.kill();
            let _ = child.wait();
            return Err(Error::Cancelled);
        }
        std::thread::sleep(Duration::from_millis(250));
    };
    let _ = err_thread.join();
    while let Ok(l) = rx.try_recv() {
        push(l, progress);
    }
    let out = out_thread.and_then(|t| t.join().ok());
    #[cfg(unix)]
    let code = {
        use std::os::unix::process::ExitStatusExt;
        status
            .code()
            .unwrap_or_else(|| 128 + status.signal().unwrap_or(0))
    };
    #[cfg(not(unix))]
    let code = status.code().unwrap_or(-1);
    Ok(Executed {
        exit_code: code,
        output: out,
        stderr_tail: tail,
    })
}

pub(crate) fn abs(p: &Path) -> PathBuf {
    if p.is_absolute() {
        p.to_path_buf()
    } else {
        std::env::current_dir()
            .map(|c| c.join(p))
            .unwrap_or(p.to_path_buf())
    }
}

impl BlastRunner for BlastnRunner {
    fn describe(&self) -> String {
        format!("blastn {} ({})", self.version, self.path.display())
    }

    fn info(&self) -> EngineInfo {
        EngineInfo {
            name: "blastn".into(),
            version: self.version.clone(),
            via: self.path.display().to_string(),
        }
    }

    fn run(
        &self,
        spec: &RunSpec,
        _store: Option<&mut OpenedStore>,
        progress: &mut dyn FnMut(&ProgressEvent) -> bool,
    ) -> Result<RunOutcome, Error> {
        let SearchTarget::Db(db) = &spec.target else {
            return Err(Error::InvalidOptions(
                "blastn searches a database: give --db BASE".into(),
            ));
        };
        if spec.stats.is_some() {
            return Err(Error::InvalidOptions(
                "blastn searches the database as built, with its own statistics".into(),
            ));
        }
        if let Some(n) = self.note.lock().unwrap().take()
            && !progress(&ProgressEvent::Warning(n))
        {
            return Err(Error::Cancelled);
        }
        let db = abs(db);
        let mut args: Vec<String> = vec!["-db".into(), db.display().to_string()];
        if let Some(q) = &spec.query {
            args.push("-query".into());
            args.push(abs(q).display().to_string());
        }
        args.extend(spec.args.iter().cloned());
        let stdin = matches!(super::passthrough::value(&args, "-query"), None | Some("-"));
        let mut cmd = Command::new(&self.path);
        cmd.args(&args);
        if !progress(&ProgressEvent::EngineStart {
            stage: spec.stage.clone(),
            engine: self.describe(),
            args: spec.args.clone(),
            threads_added: spec.threads_added,
        }) {
            return Err(Error::Cancelled);
        }
        let ex = match execute(cmd, &spec.output, stdin, &spec.stage, progress) {
            Ok(ex) => ex,
            Err(e) => {
                if let EngineOutput::File(p) = &spec.output {
                    let _ = std::fs::remove_file(p);
                }
                return Err(e);
            }
        };
        Ok(RunOutcome {
            exit_code: ex.exit_code,
            db: Some(db),
            output: ex.output,
            engine_stderr: ex.stderr_tail,
            args,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_parse() {
        let v = parse_version("blastn: 2.17.0+\n Package: blast 2.17.0, build Jun  4 2025\n");
        assert_eq!(v, Some((2, 17, "2.17.0+".into())));
        assert_eq!(parse_version("blastn: 2.9.0+").unwrap().1, 9);
        assert!(parse_version("nope").is_none());
    }
}

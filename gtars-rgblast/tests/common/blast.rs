//! BLAST+ for the feature tests: binaries from `GTARS_RGBLAST_BIN` (a directory) or `PATH`,
//! else docker image `GTARS_RGBLAST_DOCKER` (default `ncbi/blast:2.17.0`). With a feature on
//! and no BLAST+ found, tests fail with a message saying how to provide it.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

pub const DEFAULT_IMAGE: &str = "ncbi/blast:2.17.0";

pub enum Blast {
    Bin(PathBuf),
    Docker(String),
}

fn on_path(prog: &str) -> Option<PathBuf> {
    let path = std::env::var_os("PATH")?;
    std::env::split_paths(&path)
        .map(|d| d.join(prog))
        .find(|p| p.is_file())
}

impl Blast {
    pub fn find() -> Blast {
        if let Some(d) = std::env::var_os("GTARS_RGBLAST_BIN") {
            let d = PathBuf::from(d);
            assert!(
                d.join("makeblastdb").is_file(),
                "GTARS_RGBLAST_BIN={} has no makeblastdb",
                d.display()
            );
            return Blast::Bin(d);
        }
        if let Some(p) = on_path("makeblastdb") {
            return Blast::Bin(p.parent().unwrap().to_path_buf());
        }
        let image = std::env::var("GTARS_RGBLAST_DOCKER").unwrap_or_else(|_| DEFAULT_IMAGE.into());
        let ok = Command::new("docker")
            .args(["image", "inspect", &image])
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false);
        assert!(
            ok,
            "BLAST+ not found: put makeblastdb/blastn/blastdbcmd on PATH, set GTARS_RGBLAST_BIN to their \
             directory, or make docker image {image} available (docker pull {image}; or set GTARS_RGBLAST_DOCKER)"
        );
        Blast::Docker(image)
    }

    /// Run a BLAST+ program in `workdir`; file arguments must be relative to it.
    pub fn run(&self, workdir: &Path, prog: &str, args: &[&str]) -> Output {
        self.run_env(workdir, prog, args, &[])
    }

    /// `run` with environment variables set for the program (`docker run -e`), e.g. blastn's
    /// `BATCH_SIZE`.
    pub fn run_env(
        &self,
        workdir: &Path,
        prog: &str,
        args: &[&str],
        env: &[(&str, &str)],
    ) -> Output {
        let out = match self {
            Blast::Bin(d) => Command::new(d.join(prog))
                .args(args)
                .envs(env.iter().copied())
                .current_dir(workdir)
                .output(),
            Blast::Docker(image) => {
                let (uid, gid) = unsafe { (libc::getuid(), libc::getgid()) };
                let mut c = Command::new("docker");
                c.args(["run", "--rm", "-u", &format!("{uid}:{gid}"), "-v"])
                    .arg(format!("{}:/w", workdir.display()))
                    .args(["-w", "/w"]);
                for (k, v) in env {
                    c.arg("-e").arg(format!("{k}={v}"));
                }
                c.args([image.as_str(), prog]).args(args).output()
            }
        };
        out.unwrap_or_else(|e| panic!("could not run {prog}: {e}"))
    }

    /// `run_env` and require success; returns stdout.
    pub fn ok_env(
        &self,
        workdir: &Path,
        prog: &str,
        args: &[&str],
        env: &[(&str, &str)],
    ) -> String {
        let o = self.run_env(workdir, prog, args, env);
        assert!(
            o.status.success(),
            "{prog} {args:?} {env:?} failed: {}\n{}",
            o.status,
            String::from_utf8_lossy(&o.stderr)
        );
        String::from_utf8(o.stdout).unwrap()
    }

    /// Run and require success; returns stdout.
    pub fn ok(&self, workdir: &Path, prog: &str, args: &[&str]) -> String {
        let o = self.run(workdir, prog, args);
        assert!(
            o.status.success(),
            "{prog} {args:?} failed: {}\n{}",
            o.status,
            String::from_utf8_lossy(&o.stderr)
        );
        String::from_utf8(o.stdout).unwrap()
    }

    /// A blastn for `--blastn PATH` (the runner): the binary, or for docker a test-only shim
    /// script written to `dir` that runs the image with the test directories (the test temp
    /// root, the system temp dir, this crate and the working directory) mounted at their own
    /// paths.
    pub fn blastn_path(&self, dir: &Path) -> PathBuf {
        match self {
            Blast::Bin(d) => d.join("blastn"),
            Blast::Docker(image) => {
                let mut dirs: Vec<PathBuf> = Vec::new();
                for d in [
                    PathBuf::from(env!("CARGO_TARGET_TMPDIR")),
                    std::env::temp_dir(),
                    PathBuf::from(env!("CARGO_MANIFEST_DIR")),
                    std::env::current_dir().unwrap(),
                ] {
                    let d = d.canonicalize().unwrap_or(d);
                    if !dirs.contains(&d) {
                        dirs.push(d);
                    }
                }
                let mounts: String = dirs
                    .iter()
                    .map(|d| format!(" -v '{0}:{0}'", d.display()))
                    .collect();
                let (uid, gid) = unsafe { (libc::getuid(), libc::getgid()) };
                let shim = dir.join("blastn");
                std::fs::write(
                    &shim,
                    format!(
                        "#!/bin/sh\n# test-only: blastn in docker\nexec docker run --rm -i -u {uid}:{gid}{mounts} -w \"$PWD\" {image} blastn \"$@\"\n"
                    ),
                )
                .unwrap();
                use std::os::unix::fs::PermissionsExt;
                std::fs::set_permissions(&shim, std::fs::Permissions::from_mode(0o755)).unwrap();
                shim
            }
        }
    }

    /// `makeblastdb -version` first line's version number, e.g. "2.17.0".
    pub fn version(&self, workdir: &Path) -> String {
        let s = self.ok(workdir, "makeblastdb", &["-version"]);
        s.split_whitespace()
            .nth(1)
            .unwrap_or("")
            .trim_end_matches('+')
            .to_string()
    }
}

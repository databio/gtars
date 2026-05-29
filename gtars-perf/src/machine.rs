//! Machine + gtars-commit metadata for run records. Uses only std + `/proc`
//! reads and a couple of subprocesses (git, rustc) to avoid extra dependencies.

use std::path::Path;
use std::process::Command;

use anyhow::Result;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Clone)]
pub struct RunMeta {
    pub timestamp_utc: String,
    pub gtars_commit: String,
    pub gtars_commit_date: String,
    pub gtars_dirty: bool,
    pub host: String,
    pub cpu_model: String,
    pub logical_cpus: usize,
    pub total_ram_mb: u64,
    pub rustc_version: String,
    pub profile: String,
}

pub fn collect_run_meta(repo: &Path) -> Result<RunMeta> {
    Ok(RunMeta {
        timestamp_utc: now_utc_rfc3339(),
        gtars_commit: git(repo, &["rev-parse", "--short", "HEAD"]).unwrap_or_else(|| "unknown".into()),
        gtars_commit_date: git(repo, &["show", "-s", "--format=%cI", "HEAD"]).unwrap_or_default(),
        gtars_dirty: git(repo, &["status", "--porcelain"]).map(|s| !s.is_empty()).unwrap_or(false),
        host: hostname(),
        cpu_model: cpu_model(),
        logical_cpus: std::thread::available_parallelism().map(|n| n.get()).unwrap_or(0),
        total_ram_mb: total_ram_mb(),
        rustc_version: rustc_version(),
        profile: if cfg!(debug_assertions) { "debug".into() } else { "release".into() },
    })
}

/// Seconds-resolution UTC timestamp without any date dependency. Computes Y-M-D
/// H:M:S from the Unix epoch (proleptic Gregorian).
fn now_utc_rfc3339() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let days = secs / 86_400;
    let rem = secs % 86_400;
    let (h, mi, s) = (rem / 3600, (rem % 3600) / 60, rem % 60);
    let (y, mo, d) = civil_from_days(days as i64);
    format!("{y:04}-{mo:02}-{d:02}T{h:02}:{mi:02}:{s:02}Z")
}

/// Howard Hinnant's days->civil algorithm.
fn civil_from_days(z: i64) -> (i64, u32, u32) {
    let z = z + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    (if m <= 2 { y + 1 } else { y }, m, d)
}

fn git(repo: &Path, args: &[&str]) -> Option<String> {
    let out = Command::new("git").arg("-C").arg(repo).args(args).output().ok()?;
    if !out.status.success() {
        return None;
    }
    Some(String::from_utf8_lossy(&out.stdout).trim().to_string())
}

fn hostname() -> String {
    std::fs::read_to_string("/proc/sys/kernel/hostname")
        .map(|s| s.trim().to_string())
        .or_else(|_| std::env::var("HOSTNAME"))
        .unwrap_or_else(|_| "unknown".into())
}

fn cpu_model() -> String {
    let info = std::fs::read_to_string("/proc/cpuinfo").unwrap_or_default();
    for line in info.lines() {
        if let Some(rest) = line.strip_prefix("model name") {
            if let Some(v) = rest.split(':').nth(1) {
                return v.trim().to_string();
            }
        }
    }
    "unknown".into()
}

fn total_ram_mb() -> u64 {
    let info = std::fs::read_to_string("/proc/meminfo").unwrap_or_default();
    for line in info.lines() {
        if let Some(rest) = line.strip_prefix("MemTotal:") {
            if let Some(kb) = rest.split_whitespace().next().and_then(|s| s.parse::<u64>().ok()) {
                return kb / 1024;
            }
        }
    }
    0
}

fn rustc_version() -> String {
    Command::new("rustc")
        .arg("--version")
        .output()
        .ok()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_else(|| "unknown".into())
}

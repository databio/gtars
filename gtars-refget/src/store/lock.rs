//! Exclusive writer lock for a store directory.
//!
//! # Why a lockfile and not `flock`
//!
//! These stores live on Lustre and NFS on a shared HPC filesystem. `flock(2)` on
//! Lustre only works cluster-wide when the filesystem is mounted with the `flock`
//! option; mounted with `localflock` it silently degrades to node-local locking,
//! and mounted with neither it can fail with `ENOLCK`. Node-local locking that
//! *looks* like it works is a worse failure mode than no lock at all, because it
//! is invisible — two SLURM jobs on two nodes would both believe they hold it.
//!
//! `O_CREAT|O_EXCL` file creation is atomic on NFSv3 (guarded CREATE), NFSv4,
//! Lustre, and Windows, needs no mount options, and matches the precedent already
//! in this crate (`transcripts::builder::BuildLock`). So: a lockfile.
//!
//! # Scope
//!
//! WRITERS ONLY. Readers never acquire this — reads must never block on a shared
//! filesystem, and after atomic publish (see [`super::atomic`]) a reader always
//! observes a complete file. The lock guards a *commit*, not a handle: it is held
//! for the seconds it takes to merge and publish the indexes, not for the hours a
//! FASTA import takes.
//!
//! # Staleness
//!
//! A holder can die without releasing (SIGKILL, node failure, OOM). PID checks are
//! meaningless across nodes, so liveness is decided in two tiers:
//!
//! - **Definitive:** the payload's `hostname` equals ours and `/proc/<pid>` does
//!   not exist. That process is gone; break immediately.
//! - **Presumed:** the payload's `heartbeat_at` AND the file's mtime are both
//!   older than `stale_after` (default 120s = 8 missed heartbeats). Break, but
//!   warn loudly naming host, pid, and age.
//!
//! Breaking is itself serialized. Contenders do NOT `remove_file` then re-create
//! (two of them would both observe staleness, both unlink, both create, both
//! proceed). Nor is a bare `rename(2)` of the lock enough: rename acts on the
//! current pathname, not the inode a contender inspected, so if A steals and
//! re-acquires before B's rename runs, B renames A's fresh lock away and both
//! proceed. Instead a contender first claims `.rgstore.lock.break` with
//! `O_EXCL`; only the claimant re-reads the lock, confirms it is still the very
//! instance it judged stale (pid, host, `started_at`), and renames it to
//! `.rgstore.lock.stale.<pid>.<nanos>`. Everyone else retries acquisition. A
//! break marker whose claimant died is itself stolen once its mtime passes
//! `BREAK_STALE_AFTER`.

use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::thread::JoinHandle;
use std::time::{Duration, SystemTime, UNIX_EPOCH};

use anyhow::{Context, Result, anyhow};
use serde::{Deserialize, Serialize};

use super::atomic::atomic_write_bytes;

/// Name of the lockfile inside a store directory.
pub const LOCK_FILENAME: &str = ".rgstore.lock";

/// Prefix of a stolen (stale) lockfile awaiting cleanup.
pub const STALE_LOCK_PREFIX: &str = ".rgstore.lock.stale.";

/// Marker claimed (`O_EXCL`) by the one contender allowed to break a stale lock.
pub const BREAK_FILENAME: &str = ".rgstore.lock.break";

/// A break marker older than this belongs to a claimant that died mid-break and
/// may be taken over. Breaking is a handful of syscalls, so this is generous.
const BREAK_STALE_AFTER: Duration = Duration::from_secs(60);

/// How often the heartbeat thread refreshes the payload.
const HEARTBEAT_INTERVAL: Duration = Duration::from_secs(15);

/// Granularity at which the heartbeat thread checks whether it should stop, so
/// releasing the lock does not wait a full heartbeat interval.
const HEARTBEAT_TICK: Duration = Duration::from_millis(100);

const DEFAULT_TIMEOUT_SECS: u64 = 30 * 60;
const DEFAULT_STALE_AFTER_SECS: u64 = 120;

const BACKOFF_START: Duration = Duration::from_millis(50);
const BACKOFF_MAX: Duration = Duration::from_secs(5);

/// The JSON payload written into the lockfile.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LockInfo {
    pub pid: u32,
    pub hostname: String,
    /// RFC 3339 timestamp of acquisition.
    pub started_at: String,
    /// RFC 3339 timestamp of the last heartbeat refresh.
    pub heartbeat_at: String,
    /// Unix epoch seconds of the last heartbeat, for age arithmetic that does
    /// not depend on parsing RFC 3339.
    pub heartbeat_epoch: u64,
    /// What the holder is doing, e.g. `"write_index_files"`. Shown in errors.
    pub operation: String,
    pub gtars_version: String,
}

impl LockInfo {
    fn new(operation: &str) -> Self {
        let now = chrono::Utc::now().to_rfc3339();
        LockInfo {
            pid: std::process::id(),
            hostname: hostname(),
            started_at: now.clone(),
            heartbeat_at: now,
            heartbeat_epoch: epoch_secs(),
            operation: operation.to_string(),
            gtars_version: env!("CARGO_PKG_VERSION").to_string(),
        }
    }

    /// Human-readable holder description for error messages.
    pub fn describe(&self) -> String {
        format!(
            "pid {} on {} (operation={}, gtars {}, since {})",
            self.pid, self.hostname, self.operation, self.gtars_version, self.started_at
        )
    }
}

/// Tunables for lock acquisition.
///
/// Defaults come from the environment when set:
/// `GTARS_STORE_LOCK_TIMEOUT` and `GTARS_STORE_LOCK_STALE_AFTER`, both in
/// seconds. A timeout of 0 means "wait forever".
#[derive(Debug, Clone)]
pub struct LockOptions {
    /// How long to block waiting for the lock before erroring. `None` = forever.
    pub timeout: Option<Duration>,
    /// How old a heartbeat must be before the lock is presumed stale.
    pub stale_after: Duration,
    /// Break any existing lock immediately, without waiting. Operator escape
    /// hatch (`--force-unlock`); do not set this by default.
    pub force: bool,
}

impl Default for LockOptions {
    fn default() -> Self {
        LockOptions {
            timeout: env_secs("GTARS_STORE_LOCK_TIMEOUT")
                .unwrap_or(Some(Duration::from_secs(DEFAULT_TIMEOUT_SECS))),
            stale_after: env_secs("GTARS_STORE_LOCK_STALE_AFTER")
                .flatten()
                .unwrap_or(Duration::from_secs(DEFAULT_STALE_AFTER_SECS)),
            force: false,
        }
    }
}

impl LockOptions {
    /// Never block: fail immediately if the lock is held by a live writer.
    pub fn non_blocking() -> Self {
        LockOptions {
            timeout: Some(Duration::ZERO),
            ..Default::default()
        }
    }

    #[must_use]
    pub fn timeout_secs(mut self, secs: u64) -> Self {
        self.timeout = if secs == 0 {
            None
        } else {
            Some(Duration::from_secs(secs))
        };
        self
    }

    #[must_use]
    pub fn force(mut self, force: bool) -> Self {
        self.force = force;
        self
    }
}

/// An acquired exclusive writer lock. Released on drop.
pub struct StoreLock {
    path: PathBuf,
    info: LockInfo,
    stop: Arc<AtomicBool>,
    heartbeat: Option<JoinHandle<()>>,
}

impl std::fmt::Debug for StoreLock {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("StoreLock")
            .field("path", &self.path)
            .field("pid", &self.info.pid)
            .field("operation", &self.info.operation)
            .finish()
    }
}

impl StoreLock {
    /// Path of the lockfile for a store directory.
    pub fn lock_path(store_dir: &Path) -> PathBuf {
        store_dir.join(LOCK_FILENAME)
    }

    /// Acquire the exclusive writer lock, blocking with bounded exponential
    /// backoff up to `opts.timeout`.
    pub fn acquire(store_dir: &Path, operation: &str, opts: LockOptions) -> Result<Self> {
        let path = Self::lock_path(store_dir);
        fs::create_dir_all(store_dir)
            .with_context(|| format!("creating store directory {:?} for locking", store_dir))?;

        if opts.force {
            steal(&path)?;
        }

        let deadline = opts.timeout.map(|t| SystemTime::now() + t);
        let mut backoff = BACKOFF_START;

        loop {
            match try_create(&path, operation) {
                Ok(Some(info)) => {
                    let stop = Arc::new(AtomicBool::new(false));
                    let heartbeat = spawn_heartbeat(path.clone(), info.clone(), stop.clone());
                    return Ok(StoreLock {
                        path,
                        info,
                        stop,
                        heartbeat: Some(heartbeat),
                    });
                }
                Ok(None) => {} // someone holds it; fall through to evaluate
                Err(e) => return Err(e),
            }

            // Held. Decide whether the holder is alive.
            match read_info(&path) {
                Ok(Some(holder)) => {
                    if let Some(reason) = staleness(&path, &holder, opts.stale_after) {
                        eprintln!(
                            "warning: breaking stale RefgetStore lock {:?}: {} ({})",
                            path,
                            holder.describe(),
                            reason
                        );
                        // Serialized through the break marker; whoever does
                        // not get it just retries and finds the lock gone.
                        let _ = break_stale(&path, &holder);
                        continue;
                    }
                    if deadline.is_some_and(|d| SystemTime::now() >= d) {
                        return Err(anyhow!(
                            "timed out waiting for the RefgetStore write lock at {:?}, \
                             held by {}. If that writer is gone, clear it with \
                             `gtars refget lock-status <store> --force-unlock`, or raise \
                             GTARS_STORE_LOCK_TIMEOUT.",
                            path,
                            holder.describe()
                        ));
                    }
                }
                // Unreadable/corrupt payload, or the holder released between our
                // create and our read. Treat an unparseable lock as stale only
                // after `stale_after`, using mtime alone.
                Ok(None) | Err(_) => {
                    if mtime_age(&path).is_none_or(|age| age > opts.stale_after) {
                        eprintln!(
                            "warning: breaking unreadable RefgetStore lock {:?} (no valid payload)",
                            path
                        );
                        let _ = break_unreadable(&path, opts.stale_after);
                        continue;
                    }
                    if deadline.is_some_and(|d| SystemTime::now() >= d) {
                        return Err(anyhow!(
                            "timed out waiting for the RefgetStore write lock at {:?} \
                             (lockfile payload unreadable)",
                            path
                        ));
                    }
                }
            }

            std::thread::sleep(backoff);
            backoff = std::cmp::min(backoff * 2, BACKOFF_MAX);
        }
    }

    /// Try to acquire without blocking. Returns `Ok(None)` if held by a live writer.
    pub fn try_acquire(store_dir: &Path, operation: &str) -> Result<Option<Self>> {
        match Self::acquire(store_dir, operation, LockOptions::non_blocking()) {
            Ok(lock) => Ok(Some(lock)),
            Err(_) => Ok(None),
        }
    }

    /// The payload this process wrote.
    pub fn info(&self) -> &LockInfo {
        &self.info
    }
}

impl Drop for StoreLock {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::Relaxed);
        if let Some(handle) = self.heartbeat.take() {
            let _ = handle.join();
        }
        // Only remove the lockfile if it is still OURS. If another process broke
        // it as stale and re-acquired, the file now belongs to them and unlinking
        // it would hand the store to a third writer.
        if let Ok(Some(current)) = read_info(&self.path) {
            if current.pid != self.info.pid || current.hostname != self.info.hostname {
                return;
            }
        }
        let _ = fs::remove_file(&self.path);
    }
}

// =========================================================================
// Free functions (operator entry points)
// =========================================================================

/// Report the current holder of a store's write lock, if any.
pub fn lock_status(store_dir: &Path) -> Result<Option<LockInfo>> {
    read_info(&StoreLock::lock_path(store_dir))
}

/// Forcibly clear a store's write lock. Operator escape hatch: use only when the
/// holder is known to be gone. Also sweeps `.rgstore.lock.stale.*` leftovers.
pub fn force_unlock(store_dir: &Path) -> Result<bool> {
    let path = StoreLock::lock_path(store_dir);
    let existed = path.exists();
    if existed {
        steal(&path)?;
    }
    // Sweep any stolen-lock leftovers (and a dead break marker) from earlier breaks.
    if let Ok(entries) = fs::read_dir(store_dir) {
        for entry in entries.filter_map(|e| e.ok()) {
            if entry
                .file_name()
                .to_str()
                .is_some_and(|n| n.starts_with(STALE_LOCK_PREFIX) || n == BREAK_FILENAME)
            {
                let _ = fs::remove_file(entry.path());
            }
        }
    }
    Ok(existed)
}

// =========================================================================
// Internals
// =========================================================================

/// Attempt an `O_CREAT|O_EXCL` create. `Ok(None)` means someone else holds it.
fn try_create(path: &Path, operation: &str) -> Result<Option<LockInfo>> {
    match OpenOptions::new().write(true).create_new(true).open(path) {
        Ok(mut file) => {
            let info = LockInfo::new(operation);
            let json = serde_json::to_vec_pretty(&info)?;
            // Write the payload directly into the exclusively-created file. It
            // must NOT be published by rename here: the rename would clobber a
            // lock another contender legitimately holds.
            file.write_all(&json)?;
            file.sync_all()?;
            Ok(Some(info))
        }
        Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => Ok(None),
        Err(e) => Err(e).with_context(|| format!("creating lockfile {:?}", path)),
    }
}

fn read_info(path: &Path) -> Result<Option<LockInfo>> {
    match fs::read_to_string(path) {
        Ok(s) => Ok(serde_json::from_str(&s).ok()),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(None),
        Err(e) => Err(e).with_context(|| format!("reading lockfile {:?}", path)),
    }
}

/// Take a held lock away by rename.
///
/// Callers breaking a STALE lock must hold the break claim (see
/// [`break_stale`]); the rename alone is not race-safe, because it acts on
/// whatever currently sits at `path`. Do NOT replace this with `remove_file` +
/// `create_new` either: that lets two contenders both unlink, both create, and
/// both proceed. The unguarded callers are the operator escape hatches
/// (`--force-unlock`), which by definition assert no live holder exists.
fn steal(path: &Path) -> Result<()> {
    let dir = path.parent().unwrap_or_else(|| Path::new("."));
    let stolen = dir.join(format!(
        "{}{}.{}",
        STALE_LOCK_PREFIX,
        std::process::id(),
        nanos()
    ));
    match fs::rename(path, &stolen) {
        Ok(()) => {
            let _ = fs::remove_file(&stolen);
            Ok(())
        }
        // Lost the steal race, or the holder released first. Either way the lock
        // is no longer the file we saw; the caller retries.
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(e) => Err(e).with_context(|| format!("stealing stale lockfile {:?}", path)),
    }
}

/// Exclusive right to break the lock at `path`. Removes the marker on drop.
struct BreakClaim(PathBuf);

impl Drop for BreakClaim {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.0);
    }
}

/// Claim the break marker next to `path`. `Ok(None)` means another contender
/// holds it (and is, or was, breaking the lock): back off and retry acquisition.
///
/// A marker left by a claimant that died is taken over by rename, which again
/// succeeds for exactly one contender; that one then claims afresh.
fn claim_break(path: &Path) -> Result<Option<BreakClaim>> {
    let dir = path.parent().unwrap_or_else(|| Path::new("."));
    let marker = dir.join(BREAK_FILENAME);
    match OpenOptions::new().write(true).create_new(true).open(&marker) {
        Ok(mut f) => {
            let _ = write!(f, "{} {}", std::process::id(), hostname());
            Ok(Some(BreakClaim(marker)))
        }
        Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {
            if mtime_age(&marker).is_some_and(|age| age > BREAK_STALE_AFTER) {
                let dead = dir.join(format!("{}break.{}.{}", STALE_LOCK_PREFIX, std::process::id(), nanos()));
                if fs::rename(&marker, &dead).is_ok() {
                    let _ = fs::remove_file(&dead);
                }
            }
            Ok(None)
        }
        Err(e) => Err(e).with_context(|| format!("creating break marker {:?}", marker)),
    }
}

/// Break the lock at `path`, but only if it is still the exact instance judged
/// stale. Under the claim nobody else steals, so the only way the lock can be a
/// different instance is that the stale holder released and a live writer
/// acquired in between; in that case leave it alone.
fn break_stale(path: &Path, judged: &LockInfo) -> Result<()> {
    let Some(_claim) = claim_break(path)? else {
        return Ok(());
    };
    match read_info(path)? {
        Some(current)
            if current.pid == judged.pid
                && current.hostname == judged.hostname
                && current.started_at == judged.started_at =>
        {
            steal(path)
        }
        // Gone, or already a different holder. Nothing to break.
        _ => Ok(()),
    }
}

/// Break a lock with no readable payload, re-checking its age under the claim so
/// a marker holder never removes a lock that was rewritten in the meantime.
fn break_unreadable(path: &Path, stale_after: Duration) -> Result<()> {
    let Some(_claim) = claim_break(path)? else {
        return Ok(());
    };
    match read_info(path)? {
        Some(_) => Ok(()), // readable now: a live writer owns it
        None => {
            if mtime_age(path).is_none_or(|age| age > stale_after) {
                steal(path)
            } else {
                Ok(())
            }
        }
    }
}

/// Returns `Some(reason)` if the lock should be considered stale.
fn staleness(path: &Path, holder: &LockInfo, stale_after: Duration) -> Option<String> {
    // Tier 1 (definitive): same host, dead pid.
    if holder.hostname == hostname() && !pid_alive(holder.pid) {
        return Some(format!("pid {} is not running on this host", holder.pid));
    }

    // Tier 2 (presumed): heartbeat AND mtime both older than the threshold. Both
    // are required — a clock skew between nodes could make one alone lie.
    let now = epoch_secs();
    let heartbeat_age = Duration::from_secs(now.saturating_sub(holder.heartbeat_epoch));
    let file_age = mtime_age(path)?;
    if heartbeat_age > stale_after && file_age > stale_after {
        return Some(format!(
            "no heartbeat for {}s (threshold {}s)",
            heartbeat_age.as_secs(),
            stale_after.as_secs()
        ));
    }
    None
}

/// `None` if the file is gone or its mtime is unreadable.
fn mtime_age(path: &Path) -> Option<Duration> {
    let modified = fs::metadata(path).ok()?.modified().ok()?;
    SystemTime::now().duration_since(modified).ok()
}

#[cfg(target_os = "linux")]
fn pid_alive(pid: u32) -> bool {
    Path::new(&format!("/proc/{}", pid)).exists()
}

/// Without `/proc` we cannot cheaply prove a pid is dead, so never claim tier-1
/// staleness; the heartbeat path still handles genuinely dead holders.
#[cfg(not(target_os = "linux"))]
fn pid_alive(_pid: u32) -> bool {
    true
}

fn hostname() -> String {
    for path in ["/proc/sys/kernel/hostname", "/etc/hostname"] {
        if let Ok(s) = fs::read_to_string(path) {
            let s = s.trim();
            if !s.is_empty() {
                return s.to_string();
            }
        }
    }
    for var in ["HOSTNAME", "COMPUTERNAME"] {
        if let Ok(s) = std::env::var(var) {
            if !s.is_empty() {
                return s;
            }
        }
    }
    "unknown-host".to_string()
}

fn epoch_secs() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0)
}

fn nanos() -> u128 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0)
}

/// Parse a duration (seconds) from an env var. Outer `None` = unset, inner
/// `None` = explicitly "wait forever" (`0`).
fn env_secs(var: &str) -> Option<Option<Duration>> {
    let raw = std::env::var(var).ok()?;
    let secs: u64 = raw.trim().parse().ok()?;
    Some(if secs == 0 {
        None
    } else {
        Some(Duration::from_secs(secs))
    })
}

/// Refresh `heartbeat_at` every [`HEARTBEAT_INTERVAL`] so other writers can tell
/// a live holder from a dead one. The refresh is an atomic publish, so a reader
/// never sees a half-written payload — but it is skipped if the file has stopped
/// being ours (another writer broke the lock), to avoid clobbering them.
fn spawn_heartbeat(path: PathBuf, info: LockInfo, stop: Arc<AtomicBool>) -> JoinHandle<()> {
    std::thread::spawn(move || {
        let mut elapsed = Duration::ZERO;
        while !stop.load(Ordering::Relaxed) {
            std::thread::sleep(HEARTBEAT_TICK);
            elapsed += HEARTBEAT_TICK;
            if elapsed < HEARTBEAT_INTERVAL {
                continue;
            }
            elapsed = Duration::ZERO;

            match read_info(&path) {
                Ok(Some(current))
                    if current.pid == info.pid && current.hostname == info.hostname => {}
                // Gone or taken over: stop refreshing rather than fight for it.
                _ => return,
            }

            let mut refreshed = info.clone();
            refreshed.heartbeat_at = chrono::Utc::now().to_rfc3339();
            refreshed.heartbeat_epoch = epoch_secs();
            if let Ok(json) = serde_json::to_vec_pretty(&refreshed) {
                let _ = atomic_write_bytes(&path, &json);
            }
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn acquire_and_release() {
        let dir = tempfile::tempdir().unwrap();
        let path = StoreLock::lock_path(dir.path());
        {
            let lock = StoreLock::acquire(dir.path(), "test", LockOptions::default()).unwrap();
            assert!(path.exists());
            assert_eq!(lock.info().operation, "test");
            assert!(lock_status(dir.path()).unwrap().is_some());
        }
        assert!(!path.exists());
        assert!(lock_status(dir.path()).unwrap().is_none());
    }

    #[test]
    fn second_acquire_times_out_while_first_is_held() {
        let dir = tempfile::tempdir().unwrap();
        let _held = StoreLock::acquire(dir.path(), "first", LockOptions::default()).unwrap();

        let err = StoreLock::acquire(
            dir.path(),
            "second",
            LockOptions::default().timeout_secs(1),
        );
        assert!(err.is_err());
        let msg = err.unwrap_err().to_string();
        assert!(msg.contains("timed out"), "unexpected error: {}", msg);
        assert!(msg.contains("first"), "error should name the holder: {}", msg);
    }

    #[test]
    fn stale_lock_from_dead_pid_on_this_host_is_broken() {
        let dir = tempfile::tempdir().unwrap();
        let path = StoreLock::lock_path(dir.path());

        // A pid that cannot be running: max_pid+1 is not allocatable.
        let dead_pid = 4_194_305u32;
        let mut info = LockInfo::new("abandoned");
        info.pid = dead_pid;
        info.hostname = hostname();
        fs::write(&path, serde_json::to_vec_pretty(&info).unwrap()).unwrap();

        // Even with a fresh heartbeat, tier-1 liveness breaks it immediately.
        let lock = StoreLock::acquire(dir.path(), "recovering", LockOptions::default()).unwrap();
        assert_eq!(lock.info().pid, std::process::id());
    }

    #[test]
    fn old_heartbeat_from_another_host_is_presumed_stale() {
        let dir = tempfile::tempdir().unwrap();
        let path = StoreLock::lock_path(dir.path());

        let mut info = LockInfo::new("remote");
        info.hostname = "some-other-node".to_string();
        info.heartbeat_epoch = epoch_secs().saturating_sub(10_000);
        fs::write(&path, serde_json::to_vec_pretty(&info).unwrap()).unwrap();
        // Backdate the file so the mtime half of the test also passes.
        let backdated = SystemTime::now() - Duration::from_secs(10_000);
        let f = fs::File::open(&path).unwrap();
        f.set_modified(backdated).unwrap();
        drop(f);

        let opts = LockOptions {
            stale_after: Duration::from_secs(120),
            ..LockOptions::default()
        };
        let lock = StoreLock::acquire(dir.path(), "taking-over", opts).unwrap();
        assert_eq!(lock.info().hostname, hostname());
    }

    #[test]
    fn live_remote_holder_is_not_broken() {
        let dir = tempfile::tempdir().unwrap();
        let path = StoreLock::lock_path(dir.path());

        let mut info = LockInfo::new("remote");
        info.hostname = "some-other-node".to_string();
        fs::write(&path, serde_json::to_vec_pretty(&info).unwrap()).unwrap();

        let result = StoreLock::acquire(
            dir.path(),
            "contender",
            LockOptions::default().timeout_secs(1),
        );
        assert!(result.is_err(), "a live remote holder must not be broken");
    }

    #[test]
    fn force_unlock_clears_a_held_lock() {
        let dir = tempfile::tempdir().unwrap();
        let mut info = LockInfo::new("remote");
        info.hostname = "some-other-node".to_string();
        fs::write(
            StoreLock::lock_path(dir.path()),
            serde_json::to_vec_pretty(&info).unwrap(),
        )
        .unwrap();

        assert!(force_unlock(dir.path()).unwrap());
        assert!(lock_status(dir.path()).unwrap().is_none());
        // Idempotent.
        assert!(!force_unlock(dir.path()).unwrap());
    }

    #[test]
    fn drop_does_not_remove_a_lock_another_writer_took_over() {
        let dir = tempfile::tempdir().unwrap();
        let path = StoreLock::lock_path(dir.path());
        let lock = StoreLock::acquire(dir.path(), "original", LockOptions::default()).unwrap();

        // Simulate a takeover: someone broke our lock and wrote their own.
        let mut other = LockInfo::new("takeover");
        other.pid = lock.info().pid + 1;
        fs::write(&path, serde_json::to_vec_pretty(&other).unwrap()).unwrap();

        drop(lock);
        let remaining = lock_status(dir.path()).unwrap().expect("lock was clobbered");
        assert_eq!(remaining.operation, "takeover");
    }
}

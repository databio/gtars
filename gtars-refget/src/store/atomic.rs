//! Atomic publish: write-to-temp, fsync, `rename(2)`.
//!
//! Every whole-file store artifact (`sequences.rgsi`, `collections.rgci`,
//! `rgstore.json`, alias TSVs, per-collection `.rgsi`, FHR sidecars, `.reftx`)
//! goes through [`atomic_write`]. A reader therefore only ever observes a
//! complete file: the destination either still holds the previous contents or
//! flips wholesale to the new ones. There is no torn intermediate, no truncated
//! index after a crash, and no need for readers to take any lock.
//!
//! The discipline is:
//!
//! 1. Create a temp file in the SAME directory as the destination (so the final
//!    rename is same-filesystem and therefore atomic).
//! 2. Let the caller write the full contents into it.
//! 3. `sync_all()` the temp file — bytes and metadata hit the platter BEFORE
//!    anything references them.
//! 4. `rename(2)` onto the destination.
//! 5. `sync_all()` the CONTAINING DIRECTORY, so the rename itself is durable.
//!    (Step 5 is what the older `.reftx` builder omitted: without it a crash can
//!    lose the directory entry update even though the file contents were synced.)

use std::fs::{self, File, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

/// Atomically publish `dest`, with contents produced by `f`.
///
/// `f` receives a buffered writer over a temp file in `dest`'s directory. If `f`
/// returns an error the temp file is discarded and `dest` is left untouched.
///
/// The parent directory is created if missing.
pub(crate) fn atomic_write<F>(dest: &Path, f: F) -> Result<()>
where
    F: FnOnce(&mut dyn Write) -> Result<()>,
{
    let dir = parent_dir(dest);
    fs::create_dir_all(&dir)
        .with_context(|| format!("creating directory {:?} for atomic write", dir))?;

    let (mut file, tmp_path) = create_temp_in(&dir)?;

    // Write the payload through a BufWriter, then unwrap it so we can fsync the
    // underlying File (BufWriter::flush alone does not reach the disk).
    let write_result = (|| -> Result<()> {
        {
            let mut buf = std::io::BufWriter::new(&mut file);
            f(&mut buf)?;
            buf.flush()?;
        }
        file.sync_all()
            .with_context(|| format!("syncing temp file {:?}", tmp_path))?;
        Ok(())
    })();

    if let Err(e) = write_result {
        drop(file);
        let _ = fs::remove_file(&tmp_path);
        return Err(e);
    }
    drop(file);

    fs::rename(&tmp_path, dest).with_context(|| {
        format!(
            "atomically publishing {:?} (temp file {:?})",
            dest, tmp_path
        )
    })?;

    sync_dir(&dir);

    Ok(())
}

/// Convenience wrapper: atomically publish a byte payload.
pub(crate) fn atomic_write_bytes(dest: &Path, bytes: &[u8]) -> Result<()> {
    atomic_write(dest, |w| {
        w.write_all(bytes)?;
        Ok(())
    })
}

/// Publish a content-addressed payload (a `.seq` file) via temp-file + rename,
/// WITHOUT an fsync.
///
/// Per-item payloads get the rename discipline but not the durability step:
/// a transcriptome import writes hundreds of thousands of tiny files, and an
/// fsync per file would dominate its runtime. The rename is what matters for
/// correctness under concurrency -- a reader never sees a half-written file,
/// and a slow writer re-publishing the same digest replaces a live file with
/// identical bytes instead of truncating it in place. Durability of the
/// payloads across a power loss is provided by the index commit that
/// references them, which does fsync, so the window is the same as before.
pub(crate) fn publish_bytes_nosync(dest: &Path, bytes: &[u8]) -> Result<()> {
    let dir = parent_dir(dest);
    let (mut file, tmp_path) = create_temp_in(&dir)?;
    if let Err(e) = file.write_all(bytes) {
        drop(file);
        let _ = fs::remove_file(&tmp_path);
        return Err(e).with_context(|| format!("writing temp file {:?}", tmp_path));
    }
    drop(file);
    fs::rename(&tmp_path, dest).with_context(|| {
        format!("publishing {:?} (temp file {:?})", dest, tmp_path)
    })?;
    Ok(())
}

/// The directory a temp file must be created in for the rename to be atomic.
fn parent_dir(dest: &Path) -> PathBuf {
    match dest.parent().filter(|p| !p.as_os_str().is_empty()) {
        Some(p) => p.to_path_buf(),
        None => PathBuf::from("."),
    }
}

/// Create a fresh temp file in `dir`, returning the open handle and its path.
///
/// Uses a `create_new` loop rather than `tempfile::NamedTempFile` so this works
/// identically with and without the `filesystem` feature (`tempfile` is
/// feature-gated, but `write_index_files` is not).
///
/// The name is deliberately `.rgstore.tmp.*`: anything that mirrors a store
/// directory (`aws s3 sync`, integrity checkers) excludes that prefix along with
/// the lock files.
fn create_temp_in(dir: &Path) -> Result<(File, PathBuf)> {
    let pid = std::process::id();
    for attempt in 0..1024u32 {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        let candidate = dir.join(format!(".rgstore.tmp.{}.{}.{}", pid, nanos, attempt));
        match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&candidate)
        {
            Ok(f) => return Ok((f, candidate)),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(e) => {
                return Err(e).with_context(|| format!("creating temp file in {:?}", dir));
            }
        }
    }
    anyhow::bail!("could not create a temp file in {:?} after 1024 attempts", dir)
}

/// `fsync` a directory so a rename into it is durable.
///
/// Best-effort: not every platform lets you open a directory for this (Windows
/// notably), and failing to fsync the directory costs durability across a power
/// loss, not correctness for concurrent readers. Never fail a publish over it.
fn sync_dir(dir: &Path) {
    if let Ok(handle) = File::open(dir) {
        let _ = handle.sync_all();
    }
}

/// True if `name` is a transient artifact of the atomic-publish or locking
/// machinery, and therefore must be ignored by anything enumerating or
/// mirroring a store directory.
pub fn is_transient_store_file(name: &str) -> bool {
    name.starts_with(".rgstore.tmp.")
        || name == ".rgstore.lock"
        || name == ".rgstore.lock.break"
        || name.starts_with(".rgstore.lock.stale.")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn publishes_full_contents() {
        let dir = tempfile::tempdir().unwrap();
        let dest = dir.path().join("out.txt");
        atomic_write(&dest, |w| {
            w.write_all(b"hello")?;
            Ok(())
        })
        .unwrap();
        assert_eq!(fs::read_to_string(&dest).unwrap(), "hello");
    }

    #[test]
    fn failed_write_leaves_destination_untouched_and_no_temp() {
        let dir = tempfile::tempdir().unwrap();
        let dest = dir.path().join("out.txt");
        atomic_write_bytes(&dest, b"original").unwrap();

        let err = atomic_write(&dest, |w| {
            w.write_all(b"partial")?;
            anyhow::bail!("boom")
        });
        assert!(err.is_err());
        assert_eq!(fs::read_to_string(&dest).unwrap(), "original");

        let leftovers: Vec<_> = fs::read_dir(dir.path())
            .unwrap()
            .filter_map(|e| e.ok())
            .filter(|e| {
                e.file_name()
                    .to_str()
                    .map_or(false, |n| n.starts_with(".rgstore.tmp."))
            })
            .collect();
        assert!(leftovers.is_empty(), "temp file was not cleaned up");
    }

    #[test]
    fn creates_missing_parent_directory() {
        let dir = tempfile::tempdir().unwrap();
        let dest = dir.path().join("nested/deeper/out.txt");
        atomic_write_bytes(&dest, b"x").unwrap();
        assert_eq!(fs::read_to_string(&dest).unwrap(), "x");
    }

    #[test]
    fn transient_file_detection() {
        assert!(is_transient_store_file(".rgstore.lock"));
        assert!(is_transient_store_file(".rgstore.lock.stale.123.456"));
        assert!(is_transient_store_file(".rgstore.tmp.1.2.0"));
        assert!(!is_transient_store_file("rgstore.json"));
        assert!(!is_transient_store_file("sequences.rgsi"));
    }
}

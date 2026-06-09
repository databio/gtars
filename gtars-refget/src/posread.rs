//! Crate-internal positioned-read helper.
//!
//! Provides a single cursorless byte-window read over a shared `&File`.
//!
//! ## Platform strategy
//!
//! On Unix the read is implemented via `pread` (`FileExt::read_exact_at`): the
//! call takes an explicit byte offset and does NOT advance the file cursor, so it
//! is safe to call concurrently on a shared `&File` (the basis for the cached
//! `Arc<File>` in the sequence store and the single `Arc<File>` in the pread
//! transcript backend).
//!
//! On non-Unix platforms there is no `pread` equivalent in stable std, so the
//! fallback clones the handle via `try_clone()` (giving it an **independent**
//! cursor), seeks to `byte_start`, and calls `read_exact`. Cloning is slightly
//! more expensive than `pread` but avoids cursor races on a shared handle.
//!
//! ## Why this module is ungated
//!
//! The sequence-store caller in `store/readonly.rs` is declared without a feature
//! gate and compiles in ALL configurations (including `--no-default-features`).
//! Gating this module behind `filesystem` would break that build, so the module
//! must be ungated. The transcript-backend caller in `transcripts/mmap.rs` IS
//! gated behind `filesystem`, so the function is always reachable from at least
//! one call site regardless of feature combination.

use std::fs::File;

use anyhow::Result;

/// Read exactly `len` bytes starting at byte offset `byte_start` from `file`.
///
/// The read is cursorless on Unix (`pread`) and therefore safe to call
/// concurrently on a shared `&File`. On non-Unix a cloned handle is used so
/// that concurrent readers do not race on a single seek position.
pub(crate) fn read_exact_window(file: &File, byte_start: usize, len: usize) -> Result<Vec<u8>> {
    let mut buf = vec![0u8; len];
    #[cfg(unix)]
    {
        use std::os::unix::fs::FileExt;
        file.read_exact_at(&mut buf, byte_start as u64)?;
    }
    #[cfg(not(unix))]
    {
        use std::io::{Read, Seek, SeekFrom};
        // try_clone gives an independent cursor so concurrent reads on the
        // shared handle do not race on a single seek position.
        let mut handle = file.try_clone()?;
        handle.seek(SeekFrom::Start(byte_start as u64))?;
        handle.read_exact(&mut buf)?;
    }
    Ok(buf)
}

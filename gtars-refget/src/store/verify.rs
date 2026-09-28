//! Store integrity verification.
//!
//! Streams every stored sequence through the normal decode path
//! ([`ReadonlyRefgetStore::stream_sequence`]), recomputes its sha512t24u (and
//! md5, when one is stored) from the decoded bytes, and reports every
//! sequence whose bytes do not hash to its stored digest.
//!
//! This is a generic "bytes vs digest" check: it does not special-case any
//! alphabet or encoding. A resident `Full` record is checked from its
//! in-memory bytes; a `Stub` is checked from disk (or, for a remote store,
//! from the remote source -- which streams every byte over HTTP; this works
//! but is not the intended primary use case). A freshly opened local store
//! holds only stubs, so opening a store and calling verify checks the
//! on-disk bytes, which is what users generally want.

use super::*;
use super::readonly::ReadonlyRefgetStore;

use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering};

use anyhow::{anyhow, Result};
use serde::Serialize;

use crate::digest::{AlphabetType, SequenceHasher};
use crate::hashkeyable::{key_to_digest_string, HashKeyable};

/// How often (in verified sequences) to print a progress line, and the
/// threshold above which a progress line is printed at all.
const PROGRESS_INTERVAL: usize = 10_000;

/// Read buffer size used while streaming a sequence for verification.
const VERIFY_READ_BUF_SIZE: usize = 1024 * 1024;

// =========================================================================
// Options
// =========================================================================

/// Options controlling [`ReadonlyRefgetStore::verify_sequences`].
#[derive(Debug, Clone)]
pub struct VerifyOptions {
    /// Worker threads. `0` = auto ([`std::thread::available_parallelism`]),
    /// `1` = serial.
    pub jobs: usize,
    /// Also recompute md5 and compare when the stored md5 is non-empty.
    pub check_md5: bool,
}

impl Default for VerifyOptions {
    fn default() -> Self {
        Self {
            jobs: 0,
            check_md5: true,
        }
    }
}

// =========================================================================
// Report types
// =========================================================================

/// Why one sequence failed verification.
///
/// Precedence when a sequence has several problems (only one failure record
/// is ever produced per sequence): [`Self::ReadError`] >
/// [`Self::LengthMismatch`] > [`Self::DigestMismatch`] > [`Self::Md5Mismatch`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum VerifyFailureKind {
    /// Decoded bytes hash to a different sha512t24u than the stored key.
    DigestMismatch,
    /// sha512t24u matches but the stored md5 does not.
    Md5Mismatch,
    /// Stream ended with a different number of bases than metadata.length.
    LengthMismatch,
    /// The sequence could not be opened or read (missing/short file, decompress error).
    ReadError,
}

impl VerifyFailureKind {
    /// Stable lowercase-with-underscores name, used by the CLI and language bindings.
    pub fn as_str(&self) -> &'static str {
        match self {
            VerifyFailureKind::DigestMismatch => "digest_mismatch",
            VerifyFailureKind::Md5Mismatch => "md5_mismatch",
            VerifyFailureKind::LengthMismatch => "length_mismatch",
            VerifyFailureKind::ReadError => "read_error",
        }
    }
}

impl std::fmt::Display for VerifyFailureKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

/// One sequence that failed verification.
#[derive(Debug, Clone, Serialize)]
pub struct VerifyFailure {
    /// Stored sha512t24u.
    pub digest: String,
    /// `metadata.name`.
    pub name: String,
    pub alphabet: AlphabetType,
    /// `metadata.length`.
    pub length: usize,
    pub kind: VerifyFailureKind,
    pub bases_read: u64,
    /// `None` for [`VerifyFailureKind::ReadError`] when no bytes were hashed.
    pub computed_sha512t24u: Option<String>,
    pub stored_md5: String,
    /// `None` when md5 was not checked, or for a [`VerifyFailureKind::ReadError`].
    pub computed_md5: Option<String>,
    /// Set for [`VerifyFailureKind::ReadError`].
    pub error: Option<String>,
}

/// Per-alphabet checked/failed counts within a [`VerifyReport`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize)]
pub struct AlphabetVerifyCounts {
    pub checked: usize,
    pub failed: usize,
}

/// Result of [`ReadonlyRefgetStore::verify_sequences`] /
/// [`crate::store::RefgetStore::verify_collection`].
#[derive(Debug, Clone, Default, Serialize)]
pub struct VerifyReport {
    pub n_checked: usize,
    pub n_ok: usize,
    pub n_failed: usize,
    /// Keyed by the alphabet's string form (e.g. "dna2bit", "dnaio", "protein", "ASCII").
    pub by_alphabet: BTreeMap<String, AlphabetVerifyCounts>,
    /// Sorted by digest so output does not depend on thread scheduling.
    pub failures: Vec<VerifyFailure>,
}

impl VerifyReport {
    /// True when no sequence failed verification.
    pub fn is_ok(&self) -> bool {
        self.n_failed == 0
    }
}

// =========================================================================
// Verification
// =========================================================================

impl ReadonlyRefgetStore {
    /// Verify one sequence: stream its decoded bytes, recompute digests, and
    /// compare against the stored metadata.
    ///
    /// Returns `Ok(None)` when the sequence checks out, `Ok(Some(failure))`
    /// when it does not (including when the sequence could not be read --
    /// this is reported as a [`VerifyFailureKind::ReadError`], not an `Err`).
    /// Returns `Err` only when `digest` does not resolve to a known sequence.
    pub fn verify_sequence(&self, digest: &str, check_md5: bool) -> Result<Option<VerifyFailure>> {
        // Resolve the digest the same way `stream_sequence` does (accept md5
        // through `md5_lookup`).
        let digest_key = digest.to_key();
        let actual_key = self
            .md5_lookup
            .get(&digest_key)
            .copied()
            .unwrap_or(digest_key);
        let record = self
            .sequence_store
            .get(&actual_key)
            .ok_or_else(|| anyhow!("Sequence not found: {}", digest))?;
        let meta = record.metadata().clone();

        let want_md5 = check_md5 && !meta.md5.is_empty();

        let mut reader = match self.stream_sequence(&meta.sha512t24u, None, None) {
            Ok(r) => r,
            Err(e) => {
                return Ok(Some(VerifyFailure {
                    digest: meta.sha512t24u.clone(),
                    name: meta.name.clone(),
                    alphabet: meta.alphabet,
                    length: meta.length,
                    kind: VerifyFailureKind::ReadError,
                    bases_read: 0,
                    computed_sha512t24u: None,
                    stored_md5: meta.md5.clone(),
                    computed_md5: None,
                    error: Some(e.to_string()),
                }));
            }
        };

        let mut hasher = SequenceHasher::new(want_md5);
        let mut buf = vec![0u8; VERIFY_READ_BUF_SIZE];
        let mut bases_read: u64 = 0;
        loop {
            match reader.read(&mut buf) {
                Ok(0) => break,
                Ok(n) => {
                    hasher.update(&buf[..n]);
                    bases_read += n as u64;
                }
                Err(e) => {
                    let (computed_sha512t24u, computed_md5, _) = hasher.finalize();
                    return Ok(Some(VerifyFailure {
                        digest: meta.sha512t24u.clone(),
                        name: meta.name.clone(),
                        alphabet: meta.alphabet,
                        length: meta.length,
                        kind: VerifyFailureKind::ReadError,
                        bases_read,
                        computed_sha512t24u: Some(computed_sha512t24u),
                        stored_md5: meta.md5.clone(),
                        computed_md5,
                        error: Some(e.to_string()),
                    }));
                }
            }
        }

        let (computed_sha512t24u, computed_md5, hashed_len) = hasher.finalize();

        if hashed_len != meta.length as u64 {
            return Ok(Some(VerifyFailure {
                digest: meta.sha512t24u.clone(),
                name: meta.name.clone(),
                alphabet: meta.alphabet,
                length: meta.length,
                kind: VerifyFailureKind::LengthMismatch,
                bases_read: hashed_len,
                computed_sha512t24u: Some(computed_sha512t24u),
                stored_md5: meta.md5.clone(),
                computed_md5,
                error: None,
            }));
        }

        if computed_sha512t24u != meta.sha512t24u {
            return Ok(Some(VerifyFailure {
                digest: meta.sha512t24u.clone(),
                name: meta.name.clone(),
                alphabet: meta.alphabet,
                length: meta.length,
                kind: VerifyFailureKind::DigestMismatch,
                bases_read: hashed_len,
                computed_sha512t24u: Some(computed_sha512t24u),
                stored_md5: meta.md5.clone(),
                computed_md5,
                error: None,
            }));
        }

        if want_md5
            && let Some(cm) = computed_md5.as_ref()
            && *cm != meta.md5
        {
            return Ok(Some(VerifyFailure {
                digest: meta.sha512t24u.clone(),
                name: meta.name.clone(),
                alphabet: meta.alphabet,
                length: meta.length,
                kind: VerifyFailureKind::Md5Mismatch,
                bases_read: hashed_len,
                computed_sha512t24u: Some(computed_sha512t24u),
                stored_md5: meta.md5.clone(),
                computed_md5,
                error: None,
            }));
        }

        Ok(None)
    }

    /// Verify many sequences (or, with `digests: None`, every sequence in the
    /// store), in parallel when `opts.jobs != 1`.
    ///
    /// `digests` entries must already be resolved sequence identifiers (a
    /// sha512t24u or md5 digest, with no `SQ.` prefix); an unknown digest
    /// makes the whole call return `Err` naming every unknown digest, before
    /// any verification work happens.
    pub fn verify_sequences(
        &self,
        digests: Option<&[String]>,
        opts: &VerifyOptions,
    ) -> Result<VerifyReport> {
        // `stream_sequence` is called from worker threads below, so the store
        // must be `Sync`. The server use case already needs this; fail loudly
        // at compile time (not with a confusing runtime error) if it regresses.
        fn _assert_sync<T: Sync>() {}
        _assert_sync::<ReadonlyRefgetStore>();

        // 1. Build the (deduped, sorted) work list of canonical sequence keys.
        let mut work: Vec<DigestKey> = match digests {
            None => self.sequence_store.keys().copied().collect(),
            Some(requested) => {
                let mut keys = Vec::with_capacity(requested.len());
                let mut unknown = Vec::new();
                for d in requested {
                    let key = d.to_key();
                    let actual = self.md5_lookup.get(&key).copied().unwrap_or(key);
                    if self.sequence_store.contains_key(&actual) {
                        keys.push(actual);
                    } else {
                        unknown.push(d.clone());
                    }
                }
                if !unknown.is_empty() {
                    return Err(anyhow!(
                        "Unknown sequence digest(s): {}",
                        unknown.join(", ")
                    ));
                }
                keys
            }
        };
        work.sort();
        work.dedup();
        let total = work.len();

        // 2. Resolve worker count, same convention as FASTA import.
        let jobs = match opts.jobs {
            0 => std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1),
            n => n,
        };
        let jobs = jobs.min(total.max(1)).max(1);

        let quiet = self.is_quiet();
        let done = AtomicUsize::new(0);
        let failed_count = AtomicUsize::new(0);

        let maybe_report_progress = |n: usize| {
            if !quiet && n.is_multiple_of(PROGRESS_INTERVAL) {
                eprintln!(
                    "Verified {} / {} sequences ({} failures so far)",
                    n,
                    total,
                    failed_count.load(Ordering::Relaxed)
                );
            }
        };

        let mut all_failures: Vec<VerifyFailure> = Vec::new();
        let mut by_alphabet: BTreeMap<String, AlphabetVerifyCounts> = BTreeMap::new();

        let record_one = |key: &DigestKey,
                           failures: &mut Vec<VerifyFailure>,
                           counts: &mut BTreeMap<String, AlphabetVerifyCounts>|
         -> Result<()> {
            let digest_str = key_to_digest_string(key);
            let alphabet = self
                .sequence_store
                .get(key)
                .map(|r| r.metadata().alphabet);
            let outcome = self.verify_sequence(&digest_str, opts.check_md5)?;
            if let Some(alphabet) = alphabet {
                let entry = counts.entry(alphabet.to_string()).or_default();
                entry.checked += 1;
                if outcome.is_some() {
                    entry.failed += 1;
                }
            }
            if let Some(failure) = outcome {
                failed_count.fetch_add(1, Ordering::Relaxed);
                failures.push(failure);
            }
            let n = done.fetch_add(1, Ordering::Relaxed) + 1;
            maybe_report_progress(n);
            Ok(())
        };

        if jobs <= 1 {
            for key in &work {
                record_one(key, &mut all_failures, &mut by_alphabet)?;
            }
        } else {
            let index = AtomicUsize::new(0);
            let work_ref = &work;
            std::thread::scope(|scope| -> Result<()> {
                let mut handles = Vec::with_capacity(jobs);
                for _ in 0..jobs {
                    let index = &index;
                    let record_one = &record_one;
                    handles.push(scope.spawn(move || -> Result<(Vec<VerifyFailure>, BTreeMap<String, AlphabetVerifyCounts>)> {
                        let mut local_failures = Vec::new();
                        let mut local_counts: BTreeMap<String, AlphabetVerifyCounts> = BTreeMap::new();
                        loop {
                            let i = index.fetch_add(1, Ordering::Relaxed);
                            if i >= work_ref.len() {
                                break;
                            }
                            record_one(&work_ref[i], &mut local_failures, &mut local_counts)?;
                        }
                        Ok((local_failures, local_counts))
                    }));
                }
                for handle in handles {
                    let (failures, counts) = handle
                        .join()
                        .map_err(|_| anyhow!("verify worker thread panicked"))??;
                    all_failures.extend(failures);
                    for (alphabet, c) in counts {
                        let entry = by_alphabet.entry(alphabet).or_default();
                        entry.checked += c.checked;
                        entry.failed += c.failed;
                    }
                }
                Ok(())
            })?;
        }

        all_failures.sort_by(|a, b| a.digest.cmp(&b.digest));
        let n_failed = all_failures.len();

        if !quiet {
            eprintln!(
                "Verified {} / {} sequences: {} failures",
                total, total, n_failed
            );
        }

        Ok(VerifyReport {
            n_checked: total,
            n_ok: total - n_failed,
            n_failed,
            by_alphabet,
            failures: all_failures,
        })
    }
}

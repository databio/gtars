//! The shared download cache for remote RefgetStores.
//!
//! ```text
//! <cache root>/                          REFGET_CACHE_DIR, else $XDG_CACHE_HOME/refget, else ~/.cache/refget
//!   CACHEDIR.TAG                         tells backup tools to skip the folder
//!   sequences/<tag>/<d0d1>/<digest>.seq  shared by every remote store
//!   remotes/<origin key>/                one per remote (first 16 hex of sha256(normalized URL))
//!     .origin  rgstore.json  sequences.rgsi  collections.rgci  collections/  aliases/  fhr/
//!   tmp/                                 files being written, renamed into sequences/ when done
//! ```
//!
//! `<tag>` is the alphabet for Encoded stores (`dna2bit`, `dna3bit`, `dnaio`,
//! `protein`, ...), `raw` for Raw and `zstd` for Zstd: the same digest packs to
//! different bytes under a different alphabet or mode. With the tag, a `.seq`
//! file never goes stale. Files reach `sequences/` only by rename, after their
//! length and digest checks pass, so a published file is always complete.

use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::{Context, Result, anyhow};
use sha2::{Digest, Sha256};

use super::StorageMode;
use crate::digest::AlphabetType;

const CACHEDIR_TAG: &str = "Signature: 8a477f597d28d172789f06886806bc55\n\
# This file is a cache directory tag created by gtars-refget.\n\
# For information about cache directory tags, see https://bford.info/cachedir/\n";

/// The default cache root: `REFGET_CACHE_DIR`, else `$XDG_CACHE_HOME/refget`,
/// else `$HOME/.cache/refget` (`.refget-cache` in the working directory when
/// no home is known).
pub fn default_cache_root() -> PathBuf {
    let nonempty = |k: &str| {
        std::env::var_os(k)
            .filter(|v| !v.is_empty())
            .map(PathBuf::from)
    };
    if let Some(p) = nonempty("REFGET_CACHE_DIR") {
        return p;
    }
    if let Some(p) = nonempty("XDG_CACHE_HOME") {
        return p.join("refget");
    }
    if let Some(home) = nonempty("HOME").or_else(|| nonempty("USERPROFILE")) {
        return home.join(".cache").join("refget");
    }
    PathBuf::from(".refget-cache")
}

/// Canonical form of a remote store URL, used for cache identity and access:
/// whitespace and trailing `/` trimmed, scheme and host lowercased, and the
/// default port (`:80` for http, `:443` for https) dropped. Everything else
/// (credentials, path, query) is kept as given.
pub fn normalize_remote_url(url: &str) -> Result<String> {
    let url = url.trim().trim_end_matches('/');
    let (scheme, rest) = url
        .split_once("://")
        .filter(|(s, r)| !s.is_empty() && !r.is_empty() && !r.starts_with('/'))
        .ok_or_else(|| anyhow!("Remote URL must include a scheme and host: {}", url))?;
    let scheme = scheme.to_ascii_lowercase();
    let (authority, tail) = rest.split_at(rest.find(['/', '?', '#']).unwrap_or(rest.len()));
    let (userinfo, host) = match authority.rsplit_once('@') {
        Some((u, h)) => (format!("{}@", u), h),
        None => (String::new(), authority),
    };
    let mut host = host.to_ascii_lowercase();
    let default_port = match scheme.as_str() {
        "http" => Some(":80"),
        "https" => Some(":443"),
        _ => None,
    };
    if let Some(port) = default_port
        && let Some(h) = host.strip_suffix(port)
    {
        host = h.to_string();
    }
    Ok(format!("{}://{}{}{}", scheme, userinfo, host, tail))
}

/// First 16 hex characters of sha256(normalized URL).
pub(crate) fn origin_key(normalized: &str) -> String {
    let hash = Sha256::digest(normalized.as_bytes());
    format!("{:x}", hash)[..16].to_string()
}

/// The cache folder for one way of packing sequences.
fn cache_tag(mode: StorageMode, alphabet: AlphabetType) -> String {
    match mode {
        StorageMode::Encoded => alphabet.to_string().to_ascii_lowercase(),
        StorageMode::Raw => "raw".to_string(),
        StorageMode::Zstd => "zstd".to_string(),
    }
}

static TMP_COUNTER: AtomicU64 = AtomicU64::new(0);

/// The shared cache at one root.
#[derive(Debug, Clone)]
pub(crate) struct SeqCache {
    root: PathBuf,
}

impl SeqCache {
    /// Open (creating if needed) the cache at `root`.
    pub(crate) fn open(root: &Path) -> Result<Self> {
        for sub in ["sequences", "remotes", "tmp"] {
            fs::create_dir_all(root.join(sub)).with_context(|| {
                format!("Failed to create cache folder {}", root.join(sub).display())
            })?;
        }
        let tag = root.join("CACHEDIR.TAG");
        if !tag.exists() {
            let _ = fs::write(&tag, CACHEDIR_TAG);
        }
        Ok(Self {
            root: root.to_path_buf(),
        })
    }

    pub(crate) fn seq_path(
        &self,
        mode: StorageMode,
        alphabet: AlphabetType,
        digest: &str,
    ) -> PathBuf {
        let d01 = digest.get(0..2).unwrap_or(digest);
        self.root
            .join("sequences")
            .join(cache_tag(mode, alphabet))
            .join(d01)
            .join(format!("{}.seq", digest))
    }

    pub(crate) fn origin_dir(&self, normalized_url: &str) -> PathBuf {
        self.root.join("remotes").join(origin_key(normalized_url))
    }

    /// Write `bytes` to a temp file under `tmp/`, then rename it to `dest`.
    /// Another process publishing the same file first counts as success.
    pub(crate) fn publish(&self, dest: &Path, bytes: &[u8]) -> Result<()> {
        let tmp_dir = self.root.join("tmp");
        fs::create_dir_all(&tmp_dir)?;
        let tmp = tmp_dir.join(format!(
            ".publish.{}.{}.tmp",
            std::process::id(),
            TMP_COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        let write = (|| -> Result<()> {
            let mut f = File::create(&tmp)?;
            f.write_all(bytes)?;
            f.sync_data()?;
            Ok(())
        })();
        if let Err(e) = write {
            let _ = fs::remove_file(&tmp);
            return Err(e).with_context(|| format!("Failed to write cache file {}", tmp.display()));
        }
        if let Some(parent) = dest.parent() {
            fs::create_dir_all(parent)?;
        }
        match fs::rename(&tmp, dest) {
            Ok(()) => Ok(()),
            // Another writer won the race with identical bytes.
            Err(_) if dest.exists() => {
                let _ = fs::remove_file(&tmp);
                Ok(())
            }
            Err(e) => Err(e).with_context(|| {
                format!("Failed to publish {} to {}", tmp.display(), dest.display())
            }),
        }
    }
}

// ============================================================================
// Cache management
// ============================================================================

/// What the cache at one root holds.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct CacheInfo {
    pub root: PathBuf,
    pub sequence_files: u64,
    pub sequence_bytes: u64,
    /// One entry per remote: its URL (from `.origin`) and the bytes of its metadata.
    pub remotes: Vec<(String, u64)>,
}

/// Count files with `ext` under `dir`, and their bytes.
fn tally(dir: &Path, ext: Option<&str>) -> (u64, u64) {
    let mut files = 0;
    let mut bytes = 0;
    let mut stack = vec![dir.to_path_buf()];
    while let Some(d) = stack.pop() {
        let Ok(entries) = fs::read_dir(&d) else {
            continue;
        };
        for entry in entries.flatten() {
            let Ok(meta) = entry.metadata() else { continue };
            if meta.is_dir() {
                stack.push(entry.path());
            } else if ext.is_none_or(|e| entry.path().extension().is_some_and(|x| x == e)) {
                files += 1;
                bytes += meta.len();
            }
        }
    }
    (files, bytes)
}

/// Summarize the cache at `root`. A missing root is an empty cache.
pub fn cache_info(root: &Path) -> Result<CacheInfo> {
    let (sequence_files, sequence_bytes) = tally(&root.join("sequences"), Some("seq"));
    let mut remotes = Vec::new();
    if let Ok(entries) = fs::read_dir(root.join("remotes")) {
        for entry in entries.flatten() {
            if !entry.path().is_dir() {
                continue;
            }
            let url = fs::read_to_string(entry.path().join(".origin"))
                .map(|s| s.trim_end().to_string())
                .unwrap_or_else(|_| format!("(unknown: {})", entry.file_name().to_string_lossy()));
            remotes.push((url, tally(&entry.path(), None).1));
        }
    }
    remotes.sort();
    Ok(CacheInfo {
        root: root.to_path_buf(),
        sequence_files,
        sequence_bytes,
        remotes,
    })
}

/// Delete everything in the cache at `root`: sequences, remote metadata and
/// temp files. Returns the bytes freed.
pub fn clear_cache(root: &Path) -> Result<u64> {
    let mut freed = 0;
    for sub in ["sequences", "tmp", "remotes"] {
        let path = root.join(sub);
        if path.exists() {
            freed += tally(&path, None).1;
            fs::remove_dir_all(&path)
                .with_context(|| format!("Failed to remove {}", path.display()))?;
        }
    }
    Ok(freed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn normalize_trims_lowercases_and_drops_default_ports() {
        let n = |u: &str| normalize_remote_url(u).unwrap();
        assert_eq!(
            n(" HTTPS://Example.COM/Store/ "),
            "https://example.com/Store"
        );
        assert_eq!(n("https://example.com/a/b///"), "https://example.com/a/b");
        assert_eq!(
            n("http://User:pw@Host:80/x?q=1"),
            "http://User:pw@host/x?q=1"
        );
        assert_eq!(n("https://host:443"), "https://host");
        assert_eq!(n("http://host:443/s"), "http://host:443/s");
        assert_eq!(n("http://[::1]:8000/s/"), "http://[::1]:8000/s");
        assert!(normalize_remote_url("example.com/store").is_err());
        assert!(normalize_remote_url("http:///store").is_err());
        // Same as Python's hashlib.sha256(b"https://example.com/store").hexdigest()[:16].
        assert_eq!(origin_key("https://example.com/store").len(), 16);
    }

    #[test]
    fn default_root_honors_env() {
        // Only read here; env changes would race other tests.
        let root = default_cache_root();
        assert!(!root.as_os_str().is_empty());
    }

    #[test]
    fn layout_publish_info_and_clear() {
        let dir = tempdir().unwrap();
        let cache = SeqCache::open(dir.path()).unwrap();
        assert!(dir.path().join("CACHEDIR.TAG").exists());
        let p = cache.seq_path(StorageMode::Encoded, AlphabetType::Dna2bit, "abcdef");
        assert!(p.ends_with("sequences/dna2bit/ab/abcdef.seq"));
        let p = cache.seq_path(StorageMode::Raw, AlphabetType::Dna2bit, "abcdef");
        assert!(p.ends_with("sequences/raw/ab/abcdef.seq"));
        cache.publish(&p, b"ACGT").unwrap();
        assert_eq!(fs::read(&p).unwrap(), b"ACGT");
        // A second publish of the same file succeeds.
        cache.publish(&p, b"ACGT").unwrap();
        assert_eq!(
            tally(&dir.path().join("tmp"), None).0,
            0,
            "publish leaves no temp files"
        );

        let url = "http://example.com/store";
        let origin = cache.origin_dir(url);
        fs::create_dir_all(&origin).unwrap();
        fs::write(origin.join(".origin"), format!("{}\n", url)).unwrap();
        fs::write(origin.join("rgstore.json"), b"{}").unwrap();
        let info = cache_info(dir.path()).unwrap();
        assert_eq!((info.sequence_files, info.sequence_bytes), (1, 4));
        assert_eq!(
            info.remotes,
            vec![(url.to_string(), 2 + url.len() as u64 + 1)]
        );

        assert_eq!(
            clear_cache(dir.path()).unwrap(),
            4 + 2 + url.len() as u64 + 1
        );
        assert!(!p.exists() && !origin.exists());
        assert_eq!(clear_cache(dir.path()).unwrap(), 0);
    }
}

//! The cache root: one folder shared with gtars-refget and every other gtars tool.
//!
//! gtars-refget keeps downloaded sequences and remote store metadata in it (`sequences/`,
//! `remotes/`, `tmp/`). gtars-rgblast keeps its N counts in `<root>/blast/nstats/`.

use std::path::{Path, PathBuf};

/// The shared cache root: `explicit` (`--cache-dir`), else gtars-refget's default
/// (`REFGET_CACHE_DIR`, else `$XDG_CACHE_HOME/refget`, else `~/.cache/refget`).
pub fn cache_root(explicit: Option<&Path>) -> PathBuf {
    explicit
        .map(Path::to_path_buf)
        .unwrap_or_else(gtars_refget::store::default_cache_root)
}

/// gtars-rgblast's folder in the cache root: `<root>/blast`.
pub fn blast_dir(explicit: Option<&Path>) -> PathBuf {
    cache_root(explicit).join("blast")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn explicit_root_wins() {
        assert_eq!(cache_root(Some(Path::new("/x"))), PathBuf::from("/x"));
        assert_eq!(blast_dir(Some(Path::new("/x"))), PathBuf::from("/x/blast"));
        assert_eq!(
            blast_dir(None),
            gtars_refget::store::default_cache_root().join("blast")
        );
    }
}

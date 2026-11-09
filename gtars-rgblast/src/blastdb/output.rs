//! Output safety: never silently replace a database, never leave a half-written one.

use crate::error::Error;
use std::cell::RefCell;
use std::path::{Path, PathBuf};

/// Database file extensions at the base path.
pub const BASE_EXTS: [&str; 11] = [
    "nin", "nhr", "nsq", "ndb", "nal", "njs", "nos", "not", "ntf", "nto", "nog",
];
/// Volume file extensions (`<base>.<NN>.<ext>`).
pub const VOLUME_EXTS: [&str; 3] = ["nin", "nhr", "nsq"];
/// Sidecars beyond plan 1's five.
pub const EXTRA_SIDECARS: [&str; 2] = ["oids.tsv", "manifest.json"];
/// The default view folder, `<base>.views/` (owned when it holds `views.tsv`).
pub const VIEWS_DIR: &str = "views";
/// Sidecar suffixes (`<base>.<suffix>`).
pub const SIDECARS: [&str; 5] = [
    "ids.tsv",
    "collections.tsv",
    "members.tsv",
    "sanitized.tsv",
    "skipped.tsv",
];

/// `<base>.<ext>`, appending (never replacing) an extension.
pub fn with_ext(base: &Path, ext: &str) -> PathBuf {
    let mut s = base.as_os_str().to_owned();
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

/// Every file a database at `base` could own.
pub struct OutputSet {
    pub base: PathBuf,
}

/// The directory a database base path lives in (`.` for a bare name).
pub(crate) fn db_dir(base: &Path) -> PathBuf {
    match base.parent() {
        Some(p) if !p.as_os_str().is_empty() => p.to_path_buf(),
        _ => PathBuf::from("."),
    }
}

impl OutputSet {
    pub fn new(base: &Path) -> Self {
        OutputSet {
            base: base.to_path_buf(),
        }
    }

    fn dir(&self) -> PathBuf {
        db_dir(&self.base)
    }

    /// True if `name` (a file name in the output directory) belongs to this database.
    pub fn owns(&self, name: &str) -> bool {
        let Some(stem) = self.base.file_name().and_then(|s| s.to_str()) else {
            return false;
        };
        let Some(rest) = name.strip_prefix(stem).and_then(|r| r.strip_prefix('.')) else {
            return false;
        };
        if BASE_EXTS.contains(&rest) || SIDECARS.contains(&rest) || EXTRA_SIDECARS.contains(&rest) {
            return true;
        }
        let digits = |x: &str| !x.is_empty() && x.bytes().all(|b| b.is_ascii_digit());
        // <NN>.<ext>, any number of digits
        matches!(rest.split_once('.'), Some((nn, ext)) if digits(nn) && VOLUME_EXTS.contains(&ext))
    }

    /// `<base>.views/`, when it exists and holds a `views.tsv` (rgmakeblastdb made it).
    pub fn views_dir(&self) -> Option<PathBuf> {
        let d = with_ext(&self.base, VIEWS_DIR);
        d.join("views.tsv").is_file().then_some(d)
    }

    /// Existing files of the set, sorted.
    pub fn existing(&self) -> Result<Vec<PathBuf>, Error> {
        let dir = self.dir();
        let rd = match std::fs::read_dir(&dir) {
            Ok(rd) => rd,
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(Vec::new()),
            Err(e) => return Err(Error::io(&dir, e)),
        };
        let mut out = Vec::new();
        for e in rd {
            let e = e.map_err(|e| Error::io(&dir, e))?;
            if let Some(n) = e.file_name().to_str()
                && self.owns(n)
                && e.file_type().map(|t| !t.is_dir()).unwrap_or(true)
            {
                out.push(dir.join(n));
            }
        }
        if let Some(v) = self.views_dir() {
            out.push(v);
        }
        out.sort();
        Ok(out)
    }

    /// Fail with `OutputExists` if any file of the set exists and `force` is false; with
    /// `force`, delete them (faster than overwriting on ext4, where truncating or renaming
    /// over a file forces writeback).
    pub fn check(&self, force: bool) -> Result<(), Error> {
        let existing = self.existing()?;
        if let Some(first) = existing.first()
            && !force
        {
            return Err(Error::OutputExists {
                path: first.clone(),
            });
        }
        Ok(())
    }

    /// Delete every existing file of the set.
    pub fn remove_existing(&self) -> Result<(), Error> {
        for p in self.existing()? {
            if p.is_dir() {
                std::fs::remove_dir_all(&p).map_err(|e| Error::io(&p, e))?;
            } else {
                std::fs::remove_file(&p).map_err(|e| Error::io(&p, e))?;
            }
        }
        Ok(())
    }

    /// Create the output directory if it is missing.
    pub fn create_dir(&self) -> Result<(), Error> {
        let d = self.dir();
        std::fs::create_dir_all(&d).map_err(|e| Error::io(&d, e))
    }
}

/// Records every file the writer creates; unless [`PartialOutput::commit`] is called, dropping
/// it deletes them, so an error, a cancel or a full disk leaves no half-written database.
#[derive(Default)]
pub struct PartialOutput {
    files: RefCell<Vec<PathBuf>>,
    committed: bool,
}

impl PartialOutput {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn record(&self, p: impl Into<PathBuf>) {
        self.files.borrow_mut().push(p.into());
    }

    pub fn commit(mut self) {
        self.committed = true;
    }
}

impl Drop for PartialOutput {
    fn drop(&mut self) {
        if !self.committed {
            for p in self.files.borrow().iter() {
                let _ = std::fs::remove_file(p);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_set_membership() {
        let s = OutputSet::new(Path::new("/x/db"));
        for n in [
            "db.nin",
            "db.ndb",
            "db.njs",
            "db.00.nsq",
            "db.123.nhr",
            "db.ids.tsv",
            "db.skipped.tsv",
            "db.oids.tsv",
            "db.manifest.json",
        ] {
            assert!(s.owns(n), "{n}");
        }
        for n in [
            "db",
            "db.fa",
            "db.00.ndb",
            "db..nin",
            "dbx.nin",
            "other.nin",
            "db.a1.nin",
            "db.views",
        ] {
            assert!(!s.owns(n), "{n}");
        }
    }
}

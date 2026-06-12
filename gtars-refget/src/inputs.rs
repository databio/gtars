//! Reusable expansion of multi-FASTA input specifications.
//!
//! Both the CLI (`gtars refget build`) and the Python binding share this single
//! code path so that the resulting list of FASTA files — and therefore the
//! resulting store — is identical regardless of how the inputs were specified.
//!
//! See [`expand_fasta_inputs`] for the ordering contract.

use std::path::{Path, PathBuf};

use anyhow::{anyhow, Context};

/// Recognized FASTA extensions (case-insensitive). Used to filter directory
/// contents. Plain explicit file paths are NOT required to match these.
pub const FASTA_EXTENSIONS: &[&str] = &["fa", "fasta", "fna", "fa.gz", "fasta.gz", "fna.gz"];

/// A bundle of user-supplied FASTA input specifications.
///
/// Entries in `paths` may be literal file paths, glob patterns, or directories.
/// `file_list` is an optional file-of-filenames (fofn).
#[derive(Debug, Default, Clone)]
pub struct FastaInputs {
    /// Literal paths to files, glob patterns, or directories, in the order
    /// the caller supplied them (e.g. CLI positional args).
    pub paths: Vec<PathBuf>,
    /// Optional fofn: a file containing one path per line.
    pub file_list: Option<PathBuf>,
}

/// Returns true if the entry's file name (case-insensitive) ends with one of
/// the recognized FASTA extensions.
fn has_fasta_extension(path: &Path) -> bool {
    let name = match path.file_name().and_then(|n| n.to_str()) {
        Some(n) => n.to_ascii_lowercase(),
        None => return false,
    };
    FASTA_EXTENSIONS
        .iter()
        .any(|ext| name.ends_with(&format!(".{ext}")))
}

/// Returns true if the literal string looks like a glob pattern (contains any
/// of the glob metacharacters `*`, `?`, `[`).
fn looks_like_glob(s: &str) -> bool {
    s.contains('*') || s.contains('?') || s.contains('[')
}

/// Expand a single working-list entry (path, glob, or directory) into zero or
/// more concrete file paths, appending to `out`.
fn expand_entry(entry: &Path, out: &mut Vec<PathBuf>) -> Result<(), anyhow::Error> {
    let entry_str = entry.to_string_lossy();

    if entry.is_dir() {
        // Directory: keep top-level files with a FASTA extension, sorted.
        let mut matches: Vec<PathBuf> = std::fs::read_dir(entry)
            .with_context(|| format!("Failed to read directory '{}'", entry.display()))?
            .filter_map(|e| e.ok().map(|e| e.path()))
            .filter(|p| p.is_file() && has_fasta_extension(p))
            .collect();
        matches.sort();
        out.extend(matches);
        Ok(())
    } else if looks_like_glob(&entry_str) {
        // Glob pattern: expand, keep files only, sorted; empty match is an error.
        let paths = glob::glob(&entry_str)
            .with_context(|| format!("Invalid glob pattern '{}'", entry_str))?;
        let mut matches: Vec<PathBuf> = Vec::new();
        for p in paths {
            let p = p.with_context(|| {
                format!("Error reading path matched by glob '{}'", entry_str)
            })?;
            if p.is_file() {
                matches.push(p);
            }
        }
        if matches.is_empty() {
            return Err(anyhow!("Glob pattern '{}' matched no files", entry_str));
        }
        matches.sort();
        out.extend(matches);
        Ok(())
    } else {
        // Plain file path: must exist. No extension requirement.
        if !entry.exists() {
            return Err(anyhow!("Input path does not exist: '{}'", entry.display()));
        }
        out.push(entry.to_path_buf());
        Ok(())
    }
}

/// Read a fofn (file-of-filenames) and return its meaningful lines as paths.
///
/// Each line is trimmed; blank lines and lines beginning with `#` are skipped.
/// Relative paths are resolved relative to the current working directory (not
/// the fofn's location).
fn read_file_list(file_list: &Path) -> Result<Vec<PathBuf>, anyhow::Error> {
    let contents = std::fs::read_to_string(file_list)
        .with_context(|| format!("Failed to read file list '{}'", file_list.display()))?;
    Ok(contents
        .lines()
        .map(|l| l.trim())
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(PathBuf::from)
        .collect())
}

/// Expand a [`FastaInputs`] bundle into a de-duplicated, deterministically
/// ordered list of FASTA file paths.
///
/// # Ordering contract
///
/// Order is: explicit paths first (in argv order), then `file_list` lines (in
/// file order); within each directory or glob the matches are sorted
/// lexicographically; duplicates are removed keeping the first occurrence. This
/// order is what feeds collection insertion and `name_lookup`, so it is stable
/// and identical across the CLI and Python entry points.
///
/// # Expansion rules
///
/// - **Directory**: top-level (non-recursive) files whose name ends with a
///   [`FASTA_EXTENSIONS`] entry (case-insensitive), sorted lexicographically.
/// - **Glob pattern** (contains `*`, `?`, or `[`): expanded with the `glob`
///   crate, keeping files only, sorted lexicographically. A glob matching zero
///   files is an error.
/// - **Plain file path**: used as-is (no extension requirement). A
///   non-existent path is an error.
///
/// A `file_list` line may itself be a glob or directory (so a fofn can mix
/// forms).
///
/// # Errors
///
/// Returns an error if: a glob matches nothing, an explicit path does not
/// exist, the fofn is unreadable, a glob pattern is malformed, or the final
/// expanded list is empty.
pub fn expand_fasta_inputs(inputs: &FastaInputs) -> Result<Vec<PathBuf>, anyhow::Error> {
    // 1. Build the ordered working list: explicit paths first, then fofn lines.
    let mut working: Vec<PathBuf> = inputs.paths.clone();
    if let Some(file_list) = &inputs.file_list {
        working.extend(read_file_list(file_list)?);
    }

    // 2. Expand each working-list entry.
    let mut expanded: Vec<PathBuf> = Vec::new();
    for entry in &working {
        expand_entry(entry, &mut expanded)?;
    }

    // 3. De-duplicate by canonical path, preserving first-seen order.
    let mut seen: std::collections::HashSet<PathBuf> = std::collections::HashSet::new();
    let mut result: Vec<PathBuf> = Vec::with_capacity(expanded.len());
    for path in expanded {
        let key = std::fs::canonicalize(&path).unwrap_or_else(|_| path.clone());
        if seen.insert(key) {
            result.push(path);
        }
    }

    // 4. The final list must not be empty.
    if result.is_empty() {
        return Err(anyhow!("no FASTA inputs found after expansion"));
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;
    use tempfile::tempdir;

    fn touch(dir: &Path, name: &str) -> PathBuf {
        let p = dir.join(name);
        File::create(&p).unwrap();
        p
    }

    #[test]
    fn fofn_parsing_skips_blanks_and_comments() {
        let dir = tempdir().unwrap();
        let a = touch(dir.path(), "a.fa");
        let b = touch(dir.path(), "b.fa");

        let fofn = dir.path().join("list.txt");
        let mut f = File::create(&fofn).unwrap();
        writeln!(f, "# a comment").unwrap();
        writeln!(f, "").unwrap();
        writeln!(f, "  {}  ", a.display()).unwrap();
        writeln!(f, "   ").unwrap();
        writeln!(f, "{}", b.display()).unwrap();
        writeln!(f, "# trailing comment").unwrap();
        drop(f);

        let inputs = FastaInputs {
            paths: vec![],
            file_list: Some(fofn),
        };
        let out = expand_fasta_inputs(&inputs).unwrap();
        assert_eq!(out, vec![a, b]);
    }

    #[test]
    fn glob_expansion_is_lexicographically_ordered_and_deterministic() {
        let dir = tempdir().unwrap();
        // Create out of lexical order.
        touch(dir.path(), "c.fa");
        touch(dir.path(), "a.fa");
        touch(dir.path(), "b.fa");

        let pattern = dir.path().join("*.fa");
        let inputs = FastaInputs {
            paths: vec![pattern.clone()],
            file_list: None,
        };
        let out1 = expand_fasta_inputs(&inputs).unwrap();
        let names: Vec<String> = out1
            .iter()
            .map(|p| p.file_name().unwrap().to_str().unwrap().to_string())
            .collect();
        assert_eq!(names, vec!["a.fa", "b.fa", "c.fa"]);

        // Repeated calls are identical.
        let out2 = expand_fasta_inputs(&inputs).unwrap();
        assert_eq!(out1, out2);
    }

    #[test]
    fn directory_input_filters_by_extension_and_sorts() {
        let dir = tempdir().unwrap();
        touch(dir.path(), "z.fa");
        touch(dir.path(), "a.fasta.gz");
        touch(dir.path(), "m.fna");
        touch(dir.path(), "decoy.txt");
        touch(dir.path(), "readme.md");

        let inputs = FastaInputs {
            paths: vec![dir.path().to_path_buf()],
            file_list: None,
        };
        let out = expand_fasta_inputs(&inputs).unwrap();
        let names: Vec<String> = out
            .iter()
            .map(|p| p.file_name().unwrap().to_str().unwrap().to_string())
            .collect();
        assert_eq!(names, vec!["a.fasta.gz", "m.fna", "z.fa"]);
    }

    #[test]
    fn dedup_preserves_first_seen_order() {
        let dir = tempdir().unwrap();
        let x = touch(dir.path(), "x.fa");
        let y = touch(dir.path(), "y.fa");

        // x explicit, then the directory (which contains x and y), then x again.
        let inputs = FastaInputs {
            paths: vec![x.clone(), dir.path().to_path_buf(), x.clone()],
            file_list: None,
        };
        let out = expand_fasta_inputs(&inputs).unwrap();
        // x should appear once, at its first position; y appears once (from dir).
        let names: Vec<String> = out
            .iter()
            .map(|p| p.file_name().unwrap().to_str().unwrap().to_string())
            .collect();
        assert_eq!(names, vec!["x.fa", "y.fa"]);
        assert_eq!(out[0], x);
        assert!(out.contains(&y));
    }

    #[test]
    fn glob_matching_nothing_is_error() {
        let dir = tempdir().unwrap();
        let pattern = dir.path().join("*.fa");
        let inputs = FastaInputs {
            paths: vec![pattern],
            file_list: None,
        };
        assert!(expand_fasta_inputs(&inputs).is_err());
    }

    #[test]
    fn nonexistent_explicit_path_is_error() {
        let dir = tempdir().unwrap();
        let missing = dir.path().join("nope.fa");
        let inputs = FastaInputs {
            paths: vec![missing],
            file_list: None,
        };
        assert!(expand_fasta_inputs(&inputs).is_err());
    }

    #[test]
    fn empty_result_is_error() {
        let inputs = FastaInputs::default();
        assert!(expand_fasta_inputs(&inputs).is_err());
    }

    #[test]
    fn composition_explicit_then_fofn_with_glob() {
        let dir = tempdir().unwrap();
        let explicit = touch(dir.path(), "0explicit.fa");
        // glob targets: created out of order
        touch(dir.path(), "g_c.fa");
        touch(dir.path(), "g_a.fa");
        touch(dir.path(), "g_b.fa");

        // fofn contains a glob line.
        let fofn = dir.path().join("list.txt");
        let mut f = File::create(&fofn).unwrap();
        writeln!(f, "{}", dir.path().join("g_*.fa").display()).unwrap();
        drop(f);

        let inputs = FastaInputs {
            paths: vec![explicit.clone()],
            file_list: Some(fofn),
        };
        let out = expand_fasta_inputs(&inputs).unwrap();
        let names: Vec<String> = out
            .iter()
            .map(|p| p.file_name().unwrap().to_str().unwrap().to_string())
            .collect();
        // 0explicit.fa would also match g_*.fa? No (prefix g_). But note list.txt
        // itself is excluded by the glob (only g_*.fa). Explicit first, then
        // glob matches sorted.
        assert_eq!(names, vec!["0explicit.fa", "g_a.fa", "g_b.fa", "g_c.fa"]);
    }
}

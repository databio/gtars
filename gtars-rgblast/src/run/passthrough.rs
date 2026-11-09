//! The engine arguments after our own flags: which are refused, the `-taxids` guard and the
//! default `-num_threads`.

use crate::error::Error;
use std::path::{Path, PathBuf};

/// Flags `search` refuses, with the reason.
pub const SEARCH_REJECTED: &[(&str, &str)] = &[
    (
        "-db",
        "give a database you built as --db BASE, with --engine blastn",
    ),
    (
        "-subject",
        "search runs against the store or --db; for a FASTA subject run blastn directly",
    ),
    (
        "-subject_loc",
        "search runs against the store or --db; for a FASTA subject run blastn directly",
    ),
    ("-remote", "search runs locally"),
];

/// Flags `origin` owns (on top of `search`'s), with the reason.
pub const ORIGIN_OWNED: &[(&str, &str)] = &[
    (
        "-outfmt",
        "origin always writes the tabular fields its report needs",
    ),
    ("-out", "origin writes the report to stdout or --output"),
    (
        "-dbsize",
        "origin searches every stage under the full selection's statistics",
    ),
    ("-seqidlist", "origin chooses the subjects of each stage"),
    (
        "-negative_seqidlist",
        "origin chooses the subjects of each stage",
    ),
    ("-taxids", "select taxa with --taxon"),
    ("-taxidlist", "select taxa with --taxon"),
    ("-negative_taxids", "select taxa with --taxon"),
    ("-negative_taxidlist", "select taxa with --taxon"),
];

pub const TAXIDS_GUARD: &str = "-taxids without taxonomy4blast.sqlite3 makes blastn 2.17 search the whole database (it only warns). Add -no_taxid_expansion to filter on these exact taxids, or select with --taxon.";

pub fn has(args: &[String], flag: &str) -> bool {
    args.iter().any(|a| a == flag)
}

/// The value after `flag`.
pub fn value<'a>(args: &'a [String], flag: &str) -> Option<&'a str> {
    let i = args.iter().position(|a| a == flag)?;
    args.get(i + 1).map(|s| s.as_str())
}

fn refuse(rules: &[(&str, &str)], args: &[String]) -> Result<(), Error> {
    for (flag, why) in rules {
        if has(args, flag) {
            return Err(Error::InvalidOptions(format!(
                "{flag} is not allowed: {why}"
            )));
        }
    }
    Ok(())
}

/// True when `taxonomy4blast.sqlite3` is in a `$BLASTDB` directory, the working directory or
/// `db_dir`.
pub fn taxonomy4blast_present(db_dir: Option<&Path>) -> bool {
    let mut dirs: Vec<PathBuf> = std::env::var_os("BLASTDB")
        .map(|v| std::env::split_paths(&v).collect())
        .unwrap_or_default();
    dirs.push(PathBuf::from("."));
    if let Some(d) = db_dir {
        dirs.push(d.to_path_buf());
    }
    dirs.iter()
        .any(|d| d.join("taxonomy4blast.sqlite3").is_file())
}

/// `search`'s rules: the refused flags and the `-taxids` guard.
pub fn check_search(args: &[String], taxonomy4blast: bool) -> Result<(), Error> {
    refuse(SEARCH_REJECTED, args)?;
    if (has(args, "-taxids") || has(args, "-taxidlist"))
        && !has(args, "-no_taxid_expansion")
        && !taxonomy4blast
    {
        return Err(Error::InvalidOptions(TAXIDS_GUARD.into()));
    }
    Ok(())
}

/// `origin`'s rules: `search`'s refused flags plus the ones origin owns.
pub fn check_origin(args: &[String]) -> Result<(), Error> {
    refuse(SEARCH_REJECTED, args)?;
    refuse(ORIGIN_OWNED, args)
}

/// Add `-num_threads <cores>` when absent; returns the number added.
pub fn ensure_threads(args: &mut Vec<String>) -> Option<usize> {
    if has(args, "-num_threads") {
        return None;
    }
    let n = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    args.push("-num_threads".into());
    args.push(n.to_string());
    Some(n)
}

/// `-max_target_seqs`, else blastn's default of 500.
pub fn max_target_seqs(args: &[String]) -> usize {
    value(args, "-max_target_seqs")
        .and_then(|v| v.parse().ok())
        .unwrap_or(500)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn a(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn refused_flags() {
        for (f, _) in SEARCH_REJECTED {
            let e = check_search(&a(&["-query", "q", f, "x"]), false).unwrap_err();
            assert!(e.to_string().starts_with(&format!("{f} is not allowed")));
        }
        let e = check_search(&a(&["-db", "x"]), false).unwrap_err();
        assert!(e.to_string().contains("--db BASE, with --engine blastn"));
        assert!(
            check_search(
                &a(&["-query", "q", "-seqidlist", "s", "-dbsize", "9"]),
                false
            )
            .is_ok()
        );
        for (f, _) in ORIGIN_OWNED {
            assert!(check_origin(&a(&[f, "x"])).is_err(), "{f}");
        }
        assert!(check_origin(&a(&["-task", "blastn", "-evalue", "1e-5"])).is_ok());
    }

    #[test]
    fn taxids_guard() {
        let e = check_search(&a(&["-taxids", "9606"]), false).unwrap_err();
        assert_eq!(e.to_string(), TAXIDS_GUARD);
        assert!(check_search(&a(&["-taxidlist", "t.txt"]), false).is_err());
        assert!(check_search(&a(&["-taxids", "9606", "-no_taxid_expansion"]), false).is_ok());
        assert!(check_search(&a(&["-taxids", "9606"]), true).is_ok());
    }

    #[test]
    fn threads_and_mts() {
        let mut v = a(&["-query", "q"]);
        assert!(ensure_threads(&mut v).is_some());
        assert_eq!(v[2], "-num_threads");
        let mut v = a(&["-num_threads", "2"]);
        assert_eq!(ensure_threads(&mut v), None);
        assert_eq!(max_target_seqs(&a(&["-max_target_seqs", "20"])), 20);
        assert_eq!(max_target_seqs(&[]), 500);
    }
}

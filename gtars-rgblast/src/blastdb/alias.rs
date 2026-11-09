//! Native BLAST alias files (`.nal`) and binary seqid lists (`.bsl`), for views. No BLAST+
//! install is needed, and unlike `blastdb_aliastool -seqidlist` the `.nal` carries the exact
//! subset `LENGTH` (the tool writes SeqDB's approximate one).
//!
//! `.bsl`, as `blastdb_aliastool -seqid_file_in` writes it: u8 `0`; u64 LE file size; u64 LE id
//! count; u32 LE title length + title; u8 date length + date; u64 LE database length (0 when no
//! database info is pinned, the only form written here); then the ids, normalized and sorted
//! bytewise, each as a u8 length + bytes, or `0xFF` + u32 LE length when 255 bytes or longer.

use super::date::BlastDate;
use crate::error::Error;
use std::path::Path;

/// What a view `.nal` says.
#[derive(Debug, Clone, Default)]
pub struct NalSpec<'a> {
    pub title: &'a str,
    /// Database or volume names, relative to the `.nal`.
    pub dblist: &'a [String],
    /// `SEQIDLIST <file>`.
    pub seqidlist: Option<&'a str>,
    /// `NSEQ`.
    pub nseq: Option<u64>,
    /// `LENGTH`, exact.
    pub length: Option<u64>,
    /// `STATS_NSEQ` and `STATS_TOTLEN` (views with `--full-db-stats`).
    pub full_stats: Option<(u64, u64)>,
    /// makeblastdb date text.
    pub date: &'a str,
}

fn quote(e: &str) -> String {
    if e.contains(' ') {
        format!("\"{e}\"")
    } else {
        e.to_string()
    }
}

/// The text of a `.nal`.
pub fn nal_text(s: &NalSpec) -> String {
    let mut t = format!(
        "#\n# Alias file created by rgmakeblastdb {} on {}\n#\nTITLE {}\nDBLIST {}\n",
        env!("CARGO_PKG_VERSION"),
        s.date,
        s.title,
        s.dblist
            .iter()
            .map(|e| quote(e))
            .collect::<Vec<_>>()
            .join(" ")
    );
    if let Some(f) = s.seqidlist {
        t += &format!("SEQIDLIST {}\n", quote(f));
    }
    if let Some(n) = s.nseq {
        t += &format!("NSEQ {n}\n");
    }
    if let Some(l) = s.length {
        t += &format!("LENGTH {l}\n");
    }
    if let Some((n, l)) = s.full_stats {
        t += &format!("STATS_NSEQ {n}\nSTATS_TOTLEN {l}\n");
    }
    t
}

pub(crate) fn write_atomic(path: &Path, bytes: &[u8]) -> Result<(), Error> {
    let mut tmp = path.as_os_str().to_owned();
    tmp.push(format!(".tmp{}", std::process::id()));
    let tmp = std::path::PathBuf::from(tmp);
    std::fs::write(&tmp, bytes).map_err(|e| Error::io(&tmp, e))?;
    std::fs::rename(&tmp, path).map_err(|e| {
        let _ = std::fs::remove_file(&tmp);
        Error::io(path, e)
    })
}

/// The id spelling SeqDB matches a seqid list against: `gnl|db|tag` -> `db:tag`, `lcl|x` ->
/// `x`; `SQ.<digest>` and `db:tag` stay as they are.
pub fn normalize_id(id: &str) -> String {
    if let Some(r) = id.strip_prefix("gnl|")
        && let Some((db, tag)) = r.split_once('|')
    {
        return format!("{db}:{tag}");
    }
    if let Some(r) = id.strip_prefix("lcl|") {
        return r.to_string();
    }
    id.to_string()
}

/// `.bsl` bytes for these ids, title and date (no database info).
pub fn bsl_bytes(ids: &[String], title: &str, date: &str) -> Vec<u8> {
    let mut v: Vec<String> = ids.iter().map(|i| normalize_id(i.trim())).collect();
    v.retain(|i| !i.is_empty());
    v.sort();
    v.dedup();
    let mut b = vec![0u8];
    b.extend_from_slice(&0u64.to_le_bytes()); // file size, filled in below
    b.extend_from_slice(&(v.len() as u64).to_le_bytes());
    b.extend_from_slice(&(title.len() as u32).to_le_bytes());
    b.extend_from_slice(title.as_bytes());
    b.push(date.len().min(255) as u8);
    b.extend_from_slice(&date.as_bytes()[..date.len().min(255)]);
    b.extend_from_slice(&0u64.to_le_bytes()); // database length: no database info
    for id in &v {
        if id.len() >= 255 {
            b.push(0xFF);
            b.extend_from_slice(&(id.len() as u32).to_le_bytes());
        } else {
            b.push(id.len() as u8);
        }
        b.extend_from_slice(id.as_bytes());
    }
    let n = b.len() as u64;
    b[1..9].copy_from_slice(&n.to_le_bytes());
    b
}

/// Write a `.bsl` with the file stem as title and today's date (`SOURCE_DATE_EPOCH` when set).
/// Accepts `SQ.<digest>`, `lcl|x`, `gnl|db|tag` and `db:tag` spellings.
pub fn write_bsl(ids: &[String], path: &Path) -> Result<(), Error> {
    let date = BlastDate::resolve(None)?.blast();
    write_bsl_dated(ids, path, &date)
}

/// [`write_bsl`] with an explicit date text.
pub fn write_bsl_dated(ids: &[String], path: &Path, date: &str) -> Result<(), Error> {
    let title = path
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    write_atomic(path, &bsl_bytes(ids, &title, date))
}

/// A parsed `.bsl`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bsl {
    pub title: String,
    pub date: String,
    pub db_length: u64,
    pub ids: Vec<String>,
}

pub fn read_bsl(path: &Path) -> Result<Bsl, Error> {
    let b = std::fs::read(path).map_err(|e| Error::io(path, e))?;
    parse_bsl(&b).map_err(|m| Error::Store(anyhow::anyhow!("{}: bad .bsl: {m}", path.display())))
}

pub fn parse_bsl(b: &[u8]) -> Result<Bsl, String> {
    let mut p = 0usize;
    let take = |p: &mut usize, n: usize| -> Result<&[u8], String> {
        let s = b.get(*p..*p + n).ok_or("truncated")?;
        *p += n;
        Ok(s)
    };
    if take(&mut p, 1)?[0] != 0 {
        return Err("not a binary seqid list".into());
    }
    let size = u64::from_le_bytes(take(&mut p, 8)?.try_into().unwrap());
    if size != b.len() as u64 {
        return Err(format!("size field {size}, file {}", b.len()));
    }
    let n = u64::from_le_bytes(take(&mut p, 8)?.try_into().unwrap());
    let tl = u32::from_le_bytes(take(&mut p, 4)?.try_into().unwrap()) as usize;
    let title = String::from_utf8_lossy(take(&mut p, tl)?).into_owned();
    let dl = take(&mut p, 1)?[0] as usize;
    let date = String::from_utf8_lossy(take(&mut p, dl)?).into_owned();
    let db_length = u64::from_le_bytes(take(&mut p, 8)?.try_into().unwrap());
    if db_length != 0 {
        // database date and volume names follow; skip them
        let l = take(&mut p, 1)?[0] as usize;
        take(&mut p, l)?;
        let l = u32::from_le_bytes(take(&mut p, 4)?.try_into().unwrap()) as usize;
        take(&mut p, l)?;
    }
    let mut ids = Vec::with_capacity(n as usize);
    for _ in 0..n {
        let mut l = take(&mut p, 1)?[0] as usize;
        if l == 0xFF {
            l = u32::from_le_bytes(take(&mut p, 4)?.try_into().unwrap()) as usize;
        }
        ids.push(String::from_utf8_lossy(take(&mut p, l)?).into_owned());
    }
    if p != b.len() {
        return Err(format!("{} trailing bytes", b.len() - p));
    }
    Ok(Bsl {
        title,
        date,
        db_length,
        ids,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalize_spellings() {
        assert_eq!(normalize_id("gnl|asm|chr1"), "asm:chr1");
        assert_eq!(normalize_id("gnl|-abc|7"), "-abc:7");
        assert_eq!(normalize_id("lcl|zz"), "zz");
        assert_eq!(normalize_id("SQ.abc"), "SQ.abc");
        assert_eq!(normalize_id("db:tag"), "db:tag");
    }

    #[test]
    fn bsl_roundtrip() {
        let long = format!("x:{}", "y".repeat(298));
        let ids: Vec<String> = vec![
            "gnl|b|2".into(),
            "SQ.Zz".into(),
            "a:1".into(),
            "lcl|q".into(),
            "a:1".into(),
            long.clone(),
        ];
        let b = bsl_bytes(&ids, "t", "Nov 8, 2026  06:13 PM");
        let p = parse_bsl(&b).unwrap();
        assert_eq!(p.title, "t");
        assert_eq!(p.date, "Nov 8, 2026  06:13 PM");
        assert_eq!(p.db_length, 0);
        assert_eq!(p.ids, vec!["SQ.Zz", "a:1", "b:2", "q", long.as_str()]);
        assert_eq!(b[b.len() - 300 - 5], 0xFF);
    }

    #[test]
    fn nal_lines() {
        let t = nal_text(&NalSpec {
            title: "db: Homo sapiens hg38",
            dblist: &["../db".into(), "with space".into()],
            seqidlist: Some("v.bsl"),
            nseq: Some(2),
            length: Some(12),
            full_stats: Some((3, 16)),
            date: "Nov 8, 2026  12:51 PM",
        });
        assert!(t.starts_with("#\n# Alias file created by rgmakeblastdb "));
        assert!(t.ends_with(
            " on Nov 8, 2026  12:51 PM\n#\nTITLE db: Homo sapiens hg38\nDBLIST ../db \"with space\"\nSEQIDLIST v.bsl\nNSEQ 2\nLENGTH 12\nSTATS_NSEQ 3\nSTATS_TOTLEN 16\n"
        ), "{t}");
        let t = nal_text(&NalSpec {
            title: "x",
            dblist: &["db.00".into(), "db.01".into()],
            date: "d",
            ..Default::default()
        });
        assert!(t.ends_with("TITLE x\nDBLIST db.00 db.01\n"), "{t}");
    }
}

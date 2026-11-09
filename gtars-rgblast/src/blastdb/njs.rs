//! `.njs`: the Blast-db-metadata JSON makeblastdb writes by reopening the finished database
//! (seqdb.cpp `GetDBMetaData`), reproduced from what we wrote.

use super::output::with_ext;
use crate::error::Error;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

/// Keys in Blast-db-metadata order.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Njs {
    pub version: String,
    pub dbname: String,
    pub dbtype: String,
    #[serde(rename = "db-version")]
    pub db_version: u32,
    pub description: String,
    #[serde(rename = "number-of-letters")]
    pub number_of_letters: u64,
    #[serde(rename = "number-of-sequences")]
    pub number_of_sequences: u64,
    #[serde(rename = "last-updated")]
    pub last_updated: String,
    #[serde(rename = "number-of-volumes")]
    pub number_of_volumes: u64,
    /// Distinct taxids in `.nto`, 0 included; omitted when 0 is the only one, as makeblastdb does.
    #[serde(
        rename = "number-of-taxids",
        skip_serializing_if = "Option::is_none",
        default
    )]
    pub number_of_taxids: Option<u64>,
    #[serde(rename = "bytes-total")]
    pub bytes_total: u64,
    #[serde(rename = "bytes-to-cache")]
    pub bytes_to_cache: u64,
    pub files: Vec<String>,
}

pub struct NjsInput<'a> {
    pub base: &'a Path,
    pub title: &'a str,
    pub letters: u64,
    pub sequences: u64,
    /// ISO date, minute precision.
    pub last_updated: String,
    pub volumes: u64,
    /// Every database file (volumes, `.nal`, LMDB and OID files; never `.njs` or sidecars).
    pub files: &'a [PathBuf],
    /// Distinct taxids in `.nto` (0 included).
    pub taxids: &'a [i32],
}

pub fn build(i: &NjsInput) -> Result<Njs, Error> {
    let mut total = 0u64;
    let mut cache = 0u64;
    let mut names = Vec::with_capacity(i.files.len());
    for f in i.files {
        let len = std::fs::metadata(f).map_err(|e| Error::io(f, e))?.len();
        total += len;
        let name = super::writer::file_name(f);
        if name.ends_with(".nin") || name.ends_with(".nsq") {
            cache += len;
        }
        names.push(name);
    }
    names.sort();
    Ok(Njs {
        version: "1.2".into(),
        dbname: super::writer::file_name(i.base),
        dbtype: "Nucleotide".into(),
        db_version: 5,
        description: i.title.into(),
        number_of_letters: i.letters,
        number_of_sequences: i.sequences,
        last_updated: i.last_updated.clone(),
        number_of_volumes: i.volumes,
        number_of_taxids: (i.taxids.iter().any(|&t| t != 0)).then_some(i.taxids.len() as u64),
        bytes_total: total,
        bytes_to_cache: cache,
        files: names,
    })
}

/// `serde_json::to_string_pretty` plus a trailing newline, as makeblastdb writes it.
pub fn to_text(n: &Njs) -> String {
    let mut s = serde_json::to_string_pretty(n).expect("plain struct serializes");
    s.push('\n');
    s
}

/// Write `<base>.njs` (call it last, after every other file's size is final).
pub fn write(i: &NjsInput) -> Result<PathBuf, Error> {
    let n = build(i)?;
    let p = with_ext(i.base, "njs");
    std::fs::write(&p, to_text(&n)).map_err(|e| Error::io(&p, e))?;
    Ok(p)
}

pub fn read(path: &Path) -> Result<Njs, Error> {
    let t = std::fs::read_to_string(path).map_err(|e| Error::io(path, e))?;
    serde_json::from_str(&t)
        .map_err(|e| Error::Store(anyhow::anyhow!("{}: bad .njs: {e}", path.display())))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn njs_field_order_and_format() {
        // The reference is makeblastdb 2.17's hg38sub.njs.
        let n = Njs {
            version: "1.2".into(),
            dbname: "hg38sub".into(),
            dbtype: "Nucleotide".into(),
            db_version: 5,
            description: "hg38sub.fa".into(),
            number_of_letters: 352290346,
            number_of_sequences: 6,
            last_updated: "2026-11-08T12:59:00".into(),
            number_of_volumes: 1,
            number_of_taxids: None,
            bytes_total: 88206258,
            bytes_to_cache: 88155919,
            files: [
                "ndb", "nhr", "nin", "nog", "nos", "not", "nsq", "ntf", "nto",
            ]
            .iter()
            .map(|e| format!("hg38sub.{e}"))
            .collect(),
        };
        let want = r#"{
  "version": "1.2",
  "dbname": "hg38sub",
  "dbtype": "Nucleotide",
  "db-version": 5,
  "description": "hg38sub.fa",
  "number-of-letters": 352290346,
  "number-of-sequences": 6,
  "last-updated": "2026-11-08T12:59:00",
  "number-of-volumes": 1,
  "bytes-total": 88206258,
  "bytes-to-cache": 88155919,
  "files": [
    "hg38sub.ndb",
    "hg38sub.nhr",
    "hg38sub.nin",
    "hg38sub.nog",
    "hg38sub.nos",
    "hg38sub.not",
    "hg38sub.nsq",
    "hg38sub.ntf",
    "hg38sub.nto"
  ]
}
"#;
        assert_eq!(to_text(&n), want);
        let mut t = n.clone();
        t.number_of_taxids = Some(3);
        let s = to_text(&t);
        assert!(
            s.contains("\"number-of-volumes\": 1,\n  \"number-of-taxids\": 3,\n  \"bytes-total\"")
        );
    }
}

// The LMDB writer and reader in this module are a port of parts of `mdb.c` from LMDB
// (liblmdb 0.9.70, the mdb.master branch vendored by the lmdb-master-sys 0.2.6 crate), which
// carries these notices:
//
//   Copyright 2011-2021 Howard Chu, Symas Corp.
//   All rights reserved.
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted only as authorized by the OpenLDAP
//   Public License.
//
//   A copy of this license is available in the file LICENSE in the
//   top-level directory of the distribution or, alternatively, at
//   <http://www.OpenLDAP.org/license.html>.
//
//   This code is derived from btree.c written by Martin Hedenfalk.
//
//   Copyright (c) 2009, 2010 Martin Hedenfalk <martin@bzero.se>
//
//   Permission to use, copy, modify, and distribute this software for any
//   purpose with or without fee is hereby granted, provided that the above
//   copyright notice and this permission notice appear in all copies.
//
//   THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
//   WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
//   MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
//   ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
//   WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
//   ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
//   OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
//
// The OpenLDAP Public License, Version 2.8, is in LICENSE-LMDB at the root of this crate.

//! The version 5 LMDB files: `.ndb` (acc2oid, volinfo, volname) and `.ntf` (taxid2offset);
//! writing, lookup and a canonical dump for comparisons.
//!
//! The files are written by [`write`], a pure-Rust port of the LMDB write path that gives the
//! bytes liblmdb gives for the same calls (one new environment, one write transaction), and
//! read by [`read`]. [`page`] has the file format both share.

mod page;
mod read;
mod write;

use crate::error::Error;
use crate::plan::{Defline, SeqId};
use page::{MDB_DUPFIXED, MDB_DUPSORT, MDB_INTEGERKEY};
use read::Reader;
use std::path::{Path, PathBuf};

/// An error of the LMDB writer or reader.
#[derive(Debug)]
pub enum LmdbError {
    /// `MDB_MAP_FULL`: the map size is too small (handled by growing the map).
    MapFull,
    /// `MDB_BAD_VALSIZE`: a key (or duplicate value) of 0 or more than 511 bytes.
    BadValSize(usize, usize),
    /// Something LMDB supports that these files never need (and this port leaves out).
    Unsupported(String),
    /// The file is not a valid LMDB file, or an internal invariant failed.
    Corrupted(String),
    Io(std::io::Error),
}

impl std::fmt::Display for LmdbError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LmdbError::MapFull => write!(f, "MDB_MAP_FULL: environment map size reached"),
            LmdbError::BadValSize(k, d) => write!(
                f,
                "MDB_BAD_VALSIZE: unsupported size of key ({k} bytes) or value ({d} bytes)"
            ),
            LmdbError::Unsupported(s) => write!(f, "not supported by this LMDB writer: {s}"),
            LmdbError::Corrupted(s) => write!(f, "MDB_CORRUPTED: {s}"),
            LmdbError::Io(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for LmdbError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            LmdbError::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for LmdbError {
    fn from(e: std::io::Error) -> Self {
        LmdbError::Io(e)
    }
}

/// An [`LmdbError`] about `path` as a crate [`Error`] (I/O errors as [`Error::Io`]).
fn lmdb_err(path: &Path) -> impl Fn(LmdbError) -> Error + '_ {
    move |e| match e {
        LmdbError::Io(io) => Error::io(path, io),
        e => Error::lmdb(path, e),
    }
}

/// The prototype's map-size estimate for `.ndb`: Σ(key length + 24) + 64 × volumes.
pub fn ndb_estimate(key_bytes_plus_24: u64, n_volumes: usize) -> u64 {
    key_bytes_plus_24 + 64 * n_volumes as u64
}

/// One acc2oid key of a [`KeyTable`]: a span of its byte arena, the OID, and whether the key
/// also goes to `.nos` (the original spelling, not the upper-case copy).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KeyEntry {
    pub start: u64,
    pub len: u16,
    pub nos: bool,
    pub oid: u32,
}

/// Every acc2oid key of a database in one byte arena, in the order they were added (OID order)
/// until [`sort_dedup`](Self::sort_dedup) puts them in LMDB order.
#[derive(Debug, Default)]
pub struct KeyTable {
    bytes: Vec<u8>,
    entries: Vec<KeyEntry>,
}

impl KeyTable {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn entries(&self) -> &[KeyEntry] {
        &self.entries
    }

    #[inline]
    pub fn key(&self, e: &KeyEntry) -> &[u8] {
        &self.bytes[e.start as usize..e.start as usize + e.len as usize]
    }

    fn push(&mut self, start: usize, end: usize, nos: bool, oid: u32) -> Result<(), Error> {
        let len: u16 = (end - start).try_into().map_err(|_| {
            Error::InvalidOptions(format!(
                "a sequence id of {} bytes is longer than the 65,535 bytes of a lookup key",
                end - start
            ))
        })?;
        self.entries.push(KeyEntry {
            start: start as u64,
            len,
            nos,
            oid,
        });
        Ok(())
    }

    /// Add the keys of one id ([`SeqId::lmdb_keys`]: an upper-case copy when it differs, then
    /// the original spelling, which also goes to `.nos`).
    pub fn push_id(&mut self, id: &SeqId, oid: u32) -> Result<(), Error> {
        let start = self.bytes.len();
        id.write_lmdb_key(&mut self.bytes);
        let end = self.bytes.len();
        let k = &self.bytes[start..end];
        if k.is_ascii() {
            if k.iter().any(u8::is_ascii_lowercase) {
                self.bytes.extend_from_within(start..end);
                self.bytes[end..].make_ascii_uppercase();
            }
        } else {
            let orig = std::str::from_utf8(k).expect("ids are strings");
            let u = orig.to_uppercase();
            if u != orig {
                self.bytes.extend_from_slice(u.as_bytes());
            }
        }
        let copy_end = self.bytes.len();
        if copy_end > end {
            self.push(end, copy_end, false, oid)?;
        }
        self.push(start, end, true, oid)
    }

    /// Add the keys of every id of every defline of one OID.
    pub fn push_deflines(&mut self, deflines: &[Defline], oid: u32) -> Result<(), Error> {
        for d in deflines {
            for id in &d.ids {
                self.push_id(id, oid)?;
            }
        }
        Ok(())
    }

    /// The `.nos` ids (original spellings) as (OID, id), in the order added.
    pub fn nos_ids(&self) -> impl Iterator<Item = (u32, &[u8])> {
        self.entries
            .iter()
            .filter(|e| e.nos)
            .map(|e| (e.oid, self.key(e)))
    }

    /// Sort by (key, OID) and drop repeated (key, OID) pairs: the order acc2oid is filled in.
    pub fn sort_dedup(&mut self) {
        let bytes = &self.bytes;
        let key = |e: &KeyEntry| &bytes[e.start as usize..e.start as usize + e.len as usize];
        self.entries
            .sort_unstable_by(|a, b| key(a).cmp(key(b)).then(a.oid.cmp(&b.oid)));
        self.entries
            .dedup_by(|a, b| a.oid == b.oid && key(a) == key(b));
    }

    /// Σ(key length + 24) over the entries.
    pub fn bytes_plus_24(&self) -> u64 {
        self.entries.iter().map(|e| e.len as u64 + 24).sum()
    }
}

/// Create `path` (replacing any file) as a new environment of `map` bytes and run `fill` in
/// one write transaction. On MDB_MAP_FULL the transaction is aborted, the map doubled with
/// `mdb_env_set_mapsize` and `fill` run again in the same environment (the file is never
/// deleted and rewritten). Returns the final map size.
fn write_env(
    path: &Path,
    mut map: u64,
    fill: impl Fn(&mut write::Env) -> Result<(), LmdbError>,
) -> Result<u64, Error> {
    let e = lmdb_err(path);
    let _ = std::fs::remove_file(path);
    let mut env = write::Env::create(path, map).map_err(&e)?;
    loop {
        env.begin();
        match fill(&mut env).and_then(|()| env.commit()) {
            Ok(()) => return Ok(map),
            Err(LmdbError::MapFull) => {
                env.abort();
                map *= 2;
                env.set_mapsize(map);
            }
            Err(x) => {
                env.abort();
                return Err(e(x));
            }
        }
    }
}

/// The starting map size for `.ndb` from the exact key bytes.
fn ndb_map(keys: &KeyTable, n_volumes: usize) -> u64 {
    let est = ndb_estimate(keys.bytes_plus_24(), n_volumes);
    (est * 4 + (1 << 20)).max(10 << 20).div_ceil(4096) * 4096
}

/// Write `.ndb`: acc2oid (DUPSORT|DUPFIXED, id -> u32 OID), volinfo (u32 -> u32 #OIDs) and
/// volname (u32 -> volume file name), in one transaction. `keys` is sorted and deduplicated
/// first. On MDB_MAP_FULL the map grows and the transaction runs again in the same file.
pub fn write_ndb(path: &Path, vols: &[(String, u32)], keys: &mut KeyTable) -> Result<(), Error> {
    keys.sort_dedup();
    let map = ndb_map(keys, vols.len());
    write_ndb_with_map(path, vols, keys, map).map(|_| ())
}

/// [`write_ndb`] from a given starting map size (sorted `keys`); returns the final map size.
fn write_ndb_with_map(
    path: &Path,
    vols: &[(String, u32)],
    keys: &KeyTable,
    map: u64,
) -> Result<u64, Error> {
    write_env(path, map, |env| {
        if !keys.is_empty() {
            let acc = env.create_db("acc2oid", MDB_DUPSORT | MDB_DUPFIXED)?;
            for k in keys.entries() {
                env.put(acc, keys.key(k), &k.oid.to_le_bytes())?;
            }
        }
        let vi = env.create_db("volinfo", MDB_INTEGERKEY)?;
        let vn = env.create_db("volname", MDB_INTEGERKEY)?;
        for (i, (name, n)) in vols.iter().enumerate() {
            let k = (i as u32).to_ne_bytes();
            env.put(vn, &k, name.as_bytes())?;
            env.put(vi, &k, &n.to_ne_bytes())?;
        }
        Ok(())
    })
}

/// Write `.ntf`: LMDB taxid2offset (DUPSORT|DUPFIXED): i32 LE taxid -> u64 LE byte offset of
/// its count in `.nto`. Grows the map on MDB_MAP_FULL, as `.ndb` does.
pub fn write_ntf(path: &Path, offsets: &[(i32, u64)]) -> Result<(), Error> {
    let map = ((offsets.len() as u64 * 64) * 4 + (1 << 20))
        .max(10 << 20)
        .div_ceil(4096)
        * 4096;
    write_env(path, map, |env| {
        let db = env.create_db("taxid2offset", MDB_DUPSORT | MDB_DUPFIXED)?;
        for (t, off) in offsets {
            env.put(db, &t.to_le_bytes(), &off.to_le_bytes())?;
        }
        Ok(())
    })?;
    Ok(())
}

/// The `.ntf` entries (taxid, offset), sorted by taxid.
pub fn read_ntf(path: &Path) -> Result<Vec<(i32, u64)>, Error> {
    let e = lmdb_err(path);
    let r = Reader::open(path).map_err(&e)?;
    let Some(db) = r.table("taxid2offset").map_err(&e)? else {
        return Ok(Vec::new());
    };
    let mut out = Vec::new();
    r.each(&db, &mut |k, v| {
        out.push((
            i32::from_le_bytes(k.try_into().unwrap_or([0xFF; 4])),
            u64::from_le_bytes(v.try_into().unwrap_or([0xFF; 8])),
        ));
        Ok(())
    })
    .map_err(&e)?;
    out.sort();
    Ok(out)
}

/// Every table, key and value of an LMDB file as sorted text lines
/// (`table\tkey\tvalue`; acc2oid keys and volname values as text, the rest hex).
/// Two files with equal dumps have identical content even if their page layout differs.
pub fn dump(path: &Path) -> Result<Vec<String>, Error> {
    let e = lmdb_err(path);
    let r = Reader::open(path).map_err(&e)?;
    let hex = |b: &[u8]| b.iter().map(|x| format!("{x:02x}")).collect::<String>();
    let mut out = Vec::new();
    for (name, db) in r.tables().map_err(&e)? {
        let n = String::from_utf8_lossy(&name).into_owned();
        out.push(format!("#table\t{n}"));
        r.each(&db, &mut |k, v| {
            let ks = if n == "acc2oid" {
                String::from_utf8_lossy(k).into_owned()
            } else {
                hex(k)
            };
            let vs = if n == "volname" {
                String::from_utf8_lossy(v).into_owned()
            } else {
                hex(v)
            };
            out.push(format!("{n}\t{ks}\t{vs}"));
            Ok(())
        })
        .map_err(&e)?;
    }
    out.sort();
    Ok(out)
}

/// Read access to a `.ndb` file.
pub struct Ndb {
    path: PathBuf,
    reader: Reader,
}

impl Ndb {
    pub fn open(path: &Path) -> Result<Self, Error> {
        Ok(Ndb {
            path: path.to_path_buf(),
            reader: Reader::open(path).map_err(lmdb_err(path))?,
        })
    }

    /// Every (key, value) of table `name` (none if the table is absent).
    fn each(&self, name: &str, f: &mut dyn FnMut(&[u8], &[u8])) -> Result<(), Error> {
        let e = lmdb_err(&self.path);
        if let Some(db) = self.reader.table(name).map_err(&e)? {
            self.reader
                .each(&db, &mut |k, v| {
                    f(k, v);
                    Ok(())
                })
                .map_err(&e)?;
        }
        Ok(())
    }

    /// OIDs stored under an acc2oid key (exact bytes).
    pub fn lookup(&self, key: &str) -> Result<Vec<u32>, Error> {
        let e = lmdb_err(&self.path);
        let Some(db) = self.reader.table("acc2oid").map_err(&e)? else {
            return Ok(Vec::new());
        };
        Ok(self
            .reader
            .get(&db, key.as_bytes())
            .map_err(&e)?
            .into_iter()
            .map(|v| u32::from_le_bytes(v.try_into().unwrap_or([0xFF; 4])))
            .collect())
    }

    /// Every acc2oid key and OID.
    pub fn all_keys(&self) -> Result<Vec<(String, u32)>, Error> {
        let mut out = Vec::new();
        self.each("acc2oid", &mut |k, v| {
            out.push((
                String::from_utf8_lossy(k).into_owned(),
                u32::from_le_bytes(v.try_into().unwrap_or([0xFF; 4])),
            ))
        })?;
        Ok(out)
    }

    /// (volume name, OIDs) from volname and volinfo, in volume order.
    pub fn volumes(&self) -> Result<Vec<(String, u32)>, Error> {
        let mut counts = std::collections::HashMap::new();
        self.each("volinfo", &mut |k, v| {
            counts.insert(
                k.to_vec(),
                u32::from_ne_bytes(v.try_into().unwrap_or([0; 4])),
            );
        })?;
        let mut out = Vec::new();
        self.each("volname", &mut |k, name| {
            out.push((
                String::from_utf8_lossy(name).into_owned(),
                counts.get(k).copied().unwrap_or(0),
            ))
        })?;
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table(n: u32) -> KeyTable {
        let mut t = KeyTable::new();
        for oid in 0..n {
            t.push_id(&SeqId::Local(format!("SQ.digest{oid:07}abc")), oid)
                .unwrap();
            t.push_id(
                &SeqId::General {
                    db: "gencode".into(),
                    tag: format!("ENST{:011}.1", n - oid),
                },
                oid,
            )
            .unwrap();
        }
        t
    }

    #[test]
    fn key_table_matches_lmdb_keys() {
        let ids = [
            SeqId::Local("SQ.abC".into()),
            SeqId::Local("SQ.ABC".into()),
            SeqId::General {
                db: "x".into(),
                tag: "chr1".into(),
            },
            SeqId::General {
                db: "É".into(),
                tag: "straße".into(),
            },
        ];
        let mut t = KeyTable::new();
        let mut want: Vec<(Vec<u8>, u32, bool)> = Vec::new();
        for (oid, id) in ids.iter().enumerate() {
            t.push_id(id, oid as u32).unwrap();
            let mut k = Vec::new();
            id.lmdb_keys(&mut k);
            want.extend(k.into_iter().map(|(k, s)| (k.into_bytes(), oid as u32, s)));
        }
        let got: Vec<(Vec<u8>, u32, bool)> = t
            .entries()
            .iter()
            .map(|e| (t.key(e).to_vec(), e.oid, e.nos))
            .collect();
        assert_eq!(got, want);
        let nos: Vec<(u32, Vec<u8>)> = t.nos_ids().map(|(o, k)| (o, k.to_vec())).collect();
        assert_eq!(nos.len(), ids.len());
        // SQ.ABC: the copy of OID 0 and the original of OID 1 are distinct (key, OID) pairs
        t.sort_dedup();
        assert_eq!(t.len(), want.len());
        let keys: Vec<&[u8]> = t.entries().iter().map(|e| t.key(e)).collect();
        assert!(keys.windows(2).all(|w| w[0] <= w[1]));
        let long = SeqId::Local("x".repeat(70_000));
        assert!(matches!(
            KeyTable::new().push_id(&long, 0),
            Err(Error::InvalidOptions(_))
        ));
    }

    /// A map far too small grows in place (same file, same environment) and the content equals
    /// a write that never ran out of room.
    #[test]
    fn ndb_grows_the_map_instead_of_rewriting() {
        let d = tempfile::tempdir().unwrap();
        let vols = vec![("db".to_string(), 20_000u32)];
        let mut t = table(20_000);
        t.sort_dedup();
        let small = d.path().join("small.ndb");
        let end = write_ndb_with_map(&small, &vols, &t, 64 << 10).unwrap();
        assert!(end > 64 << 10, "the map grew: {end}");
        let big = d.path().join("big.ndb");
        let mut t2 = table(20_000);
        write_ndb(&big, &vols, &mut t2).unwrap();
        assert_eq!(dump(&small).unwrap(), dump(&big).unwrap());
        assert_eq!(
            Ndb::open(&small)
                .unwrap()
                .lookup("SQ.digest0000007abc")
                .unwrap(),
            vec![7]
        );
    }

    /// A key table from raw (key, OID) pairs.
    fn raw(pairs: &[(Vec<u8>, u32)]) -> KeyTable {
        let mut t = KeyTable::new();
        for (k, oid) in pairs {
            let start = t.bytes.len();
            t.bytes.extend_from_slice(k);
            t.push(start, t.bytes.len(), true, *oid).unwrap();
        }
        t
    }

    /// A deterministic pseudo-random sequence.
    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 33
        }
    }

    fn vols(n: usize) -> Vec<(String, u32)> {
        (0..n)
            .map(|i| (format!("db.{i:02}"), 1000 + i as u32))
            .collect()
    }

    fn hex(b: &[u8]) -> String {
        b.iter().map(|x| format!("{x:02x}")).collect()
    }

    fn sha(b: &[u8]) -> String {
        use sha2::Digest;
        hex(&sha2::Sha256::digest(b))
    }

    /// SHA-256 and final map size of the `.ndb` liblmdb wrote for each case of [`ndb_cases`]
    /// (through heed 0.22 / lmdb-master-sys 0.2.6, built with `MDB_IDL_LOGN` 8 like these
    /// tests), recorded before the port replaced it. They hold for 4096-byte pages on a
    /// little-endian machine.
    const LIBLMDB_NDB: &[(&str, &str, u64)] = &[
        (
            "dup-13",
            "8c17d96e7347773b3afb01cb58fedc98692631004b0f56df8d7659df41adfb6d",
            10_485_760,
        ),
        (
            "dup-2",
            "aea2605bca181e5ee08494cfb21fd95c1b3b1c9ecf0735d491cc8ae65fdb8c1b",
            10_485_760,
        ),
        (
            "dup-3",
            "b487829f4063e2e2acec86507d0be8830940b3a2cdeb61c3d5239904eb7710f6",
            10_485_760,
        ),
        (
            "dup-300",
            "e5bbc78f848bb16a1a3098c9f9458a092cd5f3b2e70e0a5fd593a70e7c395662",
            10_485_760,
        ),
        (
            "dup-3000",
            "41f79c4e8c2c77334d6a449c8a55075cbc0f1e981c06b5467c96b723e6c76ea7",
            10_485_760,
        ),
        (
            "dup-4",
            "ea1468be00262bb309821118f84b95ce0fb6f113aa981c42bc091d5c7de12176",
            10_485_760,
        ),
        (
            "dup-40",
            "35ec8d167eee464acc88fc4f5aeab8cd6dad9042c40f5a87c26893b81e716fc1",
            10_485_760,
        ),
        (
            "dup-5",
            "a6060da240b9aba9aa5ebb1ca32c81ea28eff15dc94901de8ce329f3d3aced8e",
            10_485_760,
        ),
        (
            "dup-6",
            "7908d5a4f244f8ec59b722799435a3eaf7571aaaf485792075217a098c540cab",
            10_485_760,
        ),
        (
            "dup-600",
            "a328ba7a0f214e9eee3f45af64f2442d9f4cff7ef902f5a547c8d217188fbfd2",
            10_485_760,
        ),
        (
            "dup-9",
            "b93433b37a10bab10bfc86b6d99c7116d0b119c7faf5c50f504031f0f9db02f5",
            10_485_760,
        ),
        (
            "grow",
            "126db9157f7e11c441bde6e4bb74cc5ff850576fa9625188866dffb1a6c9ffd9",
            4_194_304,
        ),
        (
            "mixed-1",
            "63a66a39eab99a6690deafd8679f7587470c0b19c39783ce7156fb70bf08e693",
            10_485_760,
        ),
        (
            "mixed-2",
            "5e0a31540ae0a77332b651cc8c846dda1b2466bb77e20eb7a481ce10e0d67e39",
            24_674_304,
        ),
        (
            "mixed-3",
            "e853380ce48546412dabec9b9259ee0286152fc5f0c053c6a71660516c3d09d5",
            10_485_760,
        ),
        (
            "mixed-4",
            "d0a8edfe52e6bcdeeb16db3f495f0bea606cd833ed1ac8cb7de9bec1ec30c65a",
            10_522_624,
        ),
        (
            "no-keys",
            "b24186bc51deb9e6c0bf1eef26e5cd21f4194b24fd3c12bc5a6e08bb445edac3",
            10_485_760,
        ),
        (
            "seq-1",
            "648c51498ac9422ade497df56d45b5682ec4ebc4ddd3f98fc148620e6247e75d",
            10_485_760,
        ),
        (
            "seq-2",
            "ddc1b0b85d2e2dc7d585c224bec29c6f19ef262dd4bc26115603c95e9ee16085",
            10_485_760,
        ),
        (
            "seq-20000",
            "6984ec0b8c1b55e4de5823bbad5c097e4e891551488f102e8c41f170672e8ad3",
            15_769_600,
        ),
        (
            "seq-3000",
            "e6c41015cb01ee9959d0377ebe896ec594cbf1bd81b90b68464a5db53098f54e",
            10_485_760,
        ),
        (
            "seq-50",
            "d5dfa9bb7cd35a20aef1ddf0f126af9778650c04741ff4a23a4f39d0149f55a8",
            10_485_760,
        ),
        (
            "spill",
            "43cf84a9e646de45e2d5d8a4bb8044f24415472dad2072198de45c64d83a239d",
            22_388_736,
        ),
        (
            "vols-1",
            "29c5977e3a5db752a113a6bef78aa8cb39ffbac6dae3cda67f4dcbf506349563",
            10_485_760,
        ),
        (
            "vols-400",
            "eed6d20a110967422dfa02299267d98c4f88f03f54adff91db66dd9e6ead866d",
            10_485_760,
        ),
        (
            "vols-5",
            "c1601516b83f0da5350cdf3ef7b1694b3b7a5c103b2bda1ea5dd1871f3b2168f",
            10_485_760,
        ),
    ];

    /// SHA-256 of the `.ntf` liblmdb wrote for each case of [`ntf_cases`] (as above).
    const LIBLMDB_NTF: &[(&str, &str)] = &[
        (
            "ntf-0",
            "24375a9730e3a4b3a81ce4c3edceaedc07582238a89df264caad70d8637ebe83",
        ),
        (
            "ntf-1",
            "84d3340f637716b8e9e4a6c5b76b659c1a0b9573b751134d1c02e061d9374eb9",
        ),
        (
            "ntf-1000",
            "402ef200ad3275909006e9d243218f51b5a0fbba96f8c5330f1b3a2c9242bc3b",
        ),
        (
            "ntf-100000",
            "3e9e000c8c773fa0b0b2c51693818bb201f506dfce29df7b12fa223a8af4b1ff",
        ),
        (
            "ntf-dups",
            "165b9caa0e35e08f7a27c4533e55248dabee338cb90665da0d89477c83f737d0",
        ),
        (
            "ntf-spill",
            "bba7ee1035a15a536846c85292b5cad43a5262568aeecd237bb175501c4a1621",
        ),
    ];

    /// Whether `file` has the page size and byte order the recorded hashes were made with.
    fn comparable(file: &[u8]) -> bool {
        cfg!(target_endian = "little") && page::r32(file, 16 + 24) == 4096
    }

    /// One `.ndb` input: volumes, keys (unsorted), and a starting map size.
    struct NdbCase {
        name: String,
        vols: Vec<(String, u32)>,
        keys: KeyTable,
        map: Option<u64>,
    }

    fn ndb_cases() -> Vec<NdbCase> {
        let case = |name: String, vols, keys, map| NdbCase {
            name,
            vols,
            keys,
            map,
        };
        let mut out = Vec::new();
        // sequential keys (an upper-case copy each), and no keys at all
        for n in [1, 2, 50, 3000, 20_000] {
            out.push(case(format!("seq-{n}"), vols(1), table(n), None));
        }
        out.push(case("no-keys".into(), vols(3), KeyTable::new(), None));
        // one key with k OIDs: sub-pages that grow, then sub-databases that split
        for k in [2u32, 3, 4, 5, 6, 9, 13, 40, 300, 600, 3000] {
            let mut pairs = vec![(b"a".to_vec(), 7), (b"z".to_vec(), 1)];
            pairs.extend((0..k).map(|i| (b"dup".to_vec(), i * 37)));
            out.push(case(format!("dup-{k}"), vols(1), raw(&pairs), None));
        }
        // random keys of 1 to 511 bytes, some with several OIDs
        for (seed, n, maxlen, maxdup) in [
            (1u64, 500usize, 40usize, 3u64),
            (2, 20_000, 120, 6),
            (3, 4000, 511, 2),
            (4, 60_000, 30, 1),
        ] {
            let mut r = Lcg(seed);
            let mut pairs = Vec::new();
            for _ in 0..n {
                let len = 1 + r.next() as usize % maxlen;
                let key: Vec<u8> = (0..len)
                    .map(|_| b"ACGT.:_|aZ\xff"[r.next() as usize % 11])
                    .collect();
                for _ in 0..1 + r.next() % maxdup {
                    pairs.push((key.clone(), (r.next() % 100_000) as u32));
                }
            }
            out.push(case(
                format!("mixed-{seed}"),
                vols(1 + n % 7),
                raw(&pairs),
                None,
            ));
        }
        // a map far too small: it grows in place
        out.push(case("grow".into(), vols(2), table(20_000), Some(64 << 10)));
        for n in [1, 5, 400] {
            let v = (0..n)
                .map(|i| (format!("{}{i}", "v".repeat(i % 200)), i as u32))
                .collect();
            out.push(case(format!("vols-{n}"), v, table(10), None));
        }
        // more dirty pages than the transaction holds (the spill), with a few busy keys
        let mut r = Lcg(5);
        let mut pairs = Vec::new();
        for i in 0..150_000u32 {
            pairs.push((format!("SQ.{:08x}", r.next()).into_bytes(), i));
            if i % 50 == 0 {
                pairs.push((format!("hot{}", i % 7).into_bytes(), i));
            }
        }
        out.push(case("spill".into(), vols(2), raw(&pairs), None));
        out
    }

    fn ntf_cases() -> Vec<(String, Vec<(i32, u64)>)> {
        let mut out = Vec::new();
        // random taxids, inserted in random order
        for (seed, n) in [(1u64, 0usize), (2, 1), (3, 1000), (4, 100_000)] {
            let mut r = Lcg(seed);
            let offs = (0..n)
                .map(|_| ((r.next() as i64 - (1 << 30)) as i32, r.next()))
                .collect();
            out.push((format!("ntf-{n}"), offs));
        }
        // repeated taxids: 8-byte duplicates
        out.push((
            "ntf-dups".into(),
            (0..5000i32).map(|i| (i % 37, i as u64 * 8)).collect(),
        ));
        // the spill, also while in sub-databases (taxids with hundreds of offsets)
        let mut r = Lcg(6);
        let offs = (0..300_000)
            .map(|_| {
                let t = if r.next().is_multiple_of(23) {
                    (r.next() % 20) as i32
                } else {
                    r.next() as i32
                };
                (t, r.next())
            })
            .collect();
        out.push(("ntf-spill".into(), offs));
        out
    }

    /// The dump an `.ndb` of these inputs must give.
    fn ndb_expected(vols: &[(String, u32)], keys: &KeyTable) -> Vec<String> {
        let mut out = vec!["#table\tvolinfo".to_string(), "#table\tvolname".to_string()];
        if !keys.is_empty() {
            out.push("#table\tacc2oid".into());
        }
        for e in keys.entries() {
            out.push(format!(
                "acc2oid\t{}\t{}",
                String::from_utf8_lossy(keys.key(e)),
                hex(&e.oid.to_le_bytes())
            ));
        }
        for (i, (name, n)) in vols.iter().enumerate() {
            let k = hex(&(i as u32).to_ne_bytes());
            out.push(format!("volinfo\t{k}\t{}", hex(&n.to_ne_bytes())));
            out.push(format!("volname\t{k}\t{name}"));
        }
        out.sort();
        out
    }

    #[test]
    fn ndb_bytes_match_liblmdb() {
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("x.ndb");
        for mut c in ndb_cases() {
            c.keys.sort_dedup();
            let map = c.map.unwrap_or_else(|| ndb_map(&c.keys, c.vols.len()));
            let end = write_ndb_with_map(&p, &c.vols, &c.keys, map).unwrap();
            let bytes = std::fs::read(&p).unwrap();
            let &(_, want, want_map) = LIBLMDB_NDB.iter().find(|r| r.0 == c.name).unwrap();
            if comparable(&bytes) {
                assert_eq!(sha(&bytes), want, "{}", c.name);
                assert_eq!(end, want_map, "{}", c.name);
            }
            assert_eq!(
                dump(&p).unwrap(),
                ndb_expected(&c.vols, &c.keys),
                "{}",
                c.name
            );
            let ndb = Ndb::open(&p).unwrap();
            assert_eq!(ndb.volumes().unwrap(), c.vols, "{}", c.name);
            assert_eq!(ndb.all_keys().unwrap().len(), c.keys.len(), "{}", c.name);
            for e in c.keys.entries().iter().step_by(97) {
                let k = String::from_utf8_lossy(c.keys.key(e)).into_owned();
                if k.as_bytes() == c.keys.key(e) {
                    assert!(ndb.lookup(&k).unwrap().contains(&e.oid), "{}: {k}", c.name);
                }
            }
            assert!(ndb.lookup("not a key").unwrap().is_empty());
        }
    }

    #[test]
    fn ntf_bytes_match_liblmdb() {
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("x.ntf");
        for (name, offs) in ntf_cases() {
            write_ntf(&p, &offs).unwrap();
            let bytes = std::fs::read(&p).unwrap();
            let &(_, want) = LIBLMDB_NTF.iter().find(|r| r.0 == name).unwrap();
            if comparable(&bytes) {
                assert_eq!(sha(&bytes), want, "{name}");
            }
            // LMDB order: taxids, then offsets, as bytes; a repeated pair is stored once
            let mut want: Vec<(i32, u64)> = offs.clone();
            want.sort_by_key(|x| (x.0.to_le_bytes(), x.1.to_le_bytes()));
            want.dedup();
            let mut lines: Vec<String> = want
                .iter()
                .map(|(t, o)| {
                    format!(
                        "taxid2offset\t{}\t{}",
                        hex(&t.to_le_bytes()),
                        hex(&o.to_le_bytes())
                    )
                })
                .collect();
            lines.push("#table\ttaxid2offset".into());
            lines.sort();
            assert_eq!(dump(&p).unwrap(), lines, "{name}");
            want.sort();
            assert_eq!(read_ntf(&p).unwrap(), want, "{name}");
        }
    }

    #[test]
    fn bad_key_sizes_fail() {
        let d = tempfile::tempdir().unwrap();
        for k in [Vec::new(), vec![b'x'; 512]] {
            let mut t = raw(&[(k, 0)]);
            assert!(write_ndb(&d.path().join("x.ndb"), &vols(1), &mut t).is_err());
        }
        let mut t = raw(&[(vec![b'x'; 511], 0)]);
        write_ndb(&d.path().join("y.ndb"), &vols(1), &mut t).unwrap();
    }
}

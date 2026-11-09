//! BLAST version 5 nucleotide volumes (`.nin`/`.nhr`/`.nsq`) and the alias file (`.nal`),
//! byte-identical to makeblastdb 2.17 (CWriteDB_Volume and friends).

use super::lmdb::KeyTable;
use super::output::{PartialOutput, with_ext};
use crate::error::Error;
use crate::plan::Defline;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

pub struct WriterOptions {
    pub base: PathBuf,
    pub title: String,
    pub max_file_size: u64,
    /// makeblastdb date text.
    pub date: String,
}

struct Volume {
    stem: PathBuf, // e.g. out/db.00
    nhr: BufWriter<File>,
    nsq: BufWriter<File>,
    hdr_off: Vec<u32>,
    seq_off: Vec<u32>,
    amb_off: Vec<u32>,
    hdr_size: u64,
    seq_size: u64,
    letters: u64,
    max_len: u32,
    idx_overhead: u64,
}

/// What `DbWriter::finish` leaves for the remaining files.
pub struct Written {
    /// Volume and alias files written, in writing order.
    pub files: Vec<PathBuf>,
    /// (volume file name, OIDs in it).
    pub volumes: Vec<(String, u32)>,
    /// Letters of each volume.
    pub volume_letters: Vec<u64>,
    /// acc2oid keys of every OID written, in OID order (`.nos` takes the original spellings).
    pub keys: KeyTable,
    pub n_oids: u32,
    pub total_letters: u64,
    pub max_len: u32,
}

pub struct DbWriter<'a> {
    o: WriterOptions,
    out: &'a PartialOutput,
    cur: Option<Volume>,
    done: Vec<(PathBuf, u32)>,
    done_letters: Vec<u64>,
    oid: u32,
    keys: KeyTable,
    /// The current subject's encoded header (reused).
    hdr: Vec<u8>,
    pub total_letters: u64,
    max_len: u32,
}

pub(crate) fn file_name(p: &Path) -> String {
    p.file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default()
}

fn be32(w: &mut Vec<u8>, v: u32) {
    w.extend_from_slice(&v.to_be_bytes());
}
fn bstr(w: &mut Vec<u8>, s: &[u8]) {
    be32(w, s.len() as u32);
    w.extend_from_slice(s);
}

/// Size of the fixed `.nin` header fields with an unpadded date of `date_len` bytes
/// (CWriteDB_IndexFile::x_Overhead, version 5).
pub fn nin_overhead(title: &str, lmdb_name: &str, date_len: usize) -> u64 {
    40 + title.len() as u64 + lmdb_name.len() as u64 + date_len as u64
}

/// The `.nin` size of a volume with `n` OIDs.
pub fn nin_size(title: &str, lmdb_name: &str, date: &str, n: u64) -> u64 {
    nin_overhead(title, lmdb_name, date.len()).div_ceil(8) * 8 + if n > 0 { 12 * n + 12 } else { 8 }
}

/// The volume-split rule (CWriteDB_Volume::WriteSequence's CanFit checks): would one more
/// sequence with `hdr` header bytes and `data` sequence+ambiguity bytes fit a volume that
/// already holds `n` OIDs, `hdr_size` header bytes and `seq_size` sequence bytes?
pub fn fits(
    max: u64,
    idx_overhead: u64,
    n: u64,
    hdr_size: u64,
    seq_size: u64,
    hdr: u64,
    data: u64,
) -> bool {
    let idx_ok = n == 0 || idx_overhead + 12 * n < max - 12;
    let hdr_ok = hdr_size == 0 || hdr_size + hdr <= max - 1;
    let seq_ok = seq_size <= 1 || seq_size + data <= max - 1;
    idx_ok && hdr_ok && seq_ok
}

/// The fields of a version 5 nucleotide `.nin`.
pub struct NinFields<'a> {
    pub volume: u32,
    pub title: &'a str,
    pub lmdb_name: &'a str,
    /// Date text without padding (NULs are added to make the header a multiple of 8 bytes).
    pub date: &'a str,
    pub letters: u64,
    pub max_len: u32,
    /// n + 1 header offsets.
    pub hdr_off: &'a [u32],
    /// n + 1 sequence offsets.
    pub seq_off: &'a [u32],
    /// n ambiguity offsets.
    pub amb_off: &'a [u32],
}

/// `.nin` bytes, as CWriteDB_IndexFile writes them.
pub fn nin_bytes(f: &NinFields) -> Vec<u8> {
    let n = f.amb_off.len() as u32;
    let mut date = f.date.as_bytes().to_vec();
    while !nin_overhead(f.title, f.lmdb_name, date.len()).is_multiple_of(8) {
        date.push(0);
    }
    let mut w = Vec::with_capacity(64 + 12 * (n as usize + 1));
    be32(&mut w, 5); // version
    be32(&mut w, 0); // nucleotide
    be32(&mut w, f.volume);
    bstr(&mut w, f.title.as_bytes());
    bstr(&mut w, f.lmdb_name.as_bytes());
    bstr(&mut w, &date);
    be32(&mut w, n);
    w.extend_from_slice(&f.letters.to_le_bytes());
    be32(&mut w, f.max_len);
    for x in f.hdr_off {
        be32(&mut w, *x);
    }
    for x in f.seq_off {
        be32(&mut w, *x);
    }
    for x in f.amb_off {
        be32(&mut w, *x);
    }
    if n > 0 {
        be32(&mut w, *f.seq_off.last().unwrap());
    }
    w
}

impl<'a> DbWriter<'a> {
    pub fn new(o: WriterOptions, out: &'a PartialOutput) -> Self {
        DbWriter {
            o,
            out,
            cur: None,
            done: Vec::new(),
            done_letters: Vec::new(),
            oid: 0,
            keys: KeyTable::new(),
            hdr: Vec::new(),
            total_letters: 0,
            max_len: 0,
        }
    }

    pub fn lmdb_name(&self) -> String {
        format!("{}.ndb", file_name(&self.o.base))
    }

    fn create(&self, p: PathBuf) -> Result<File, Error> {
        self.out.record(&p);
        File::create(&p).map_err(|e| Error::io(&p, e))
    }

    fn open_volume(&mut self) -> Result<(), Error> {
        let idx = self.done.len();
        let stem = PathBuf::from(format!("{}.{:02}", self.o.base.display(), idx));
        let p = with_ext(&stem, "nsq");
        let mut nsq = BufWriter::with_capacity(8 << 20, self.create(p.clone())?);
        nsq.write_all(&[0]).map_err(|e| Error::io(&p, e))?;
        let nhr = BufWriter::with_capacity(1 << 20, self.create(with_ext(&stem, "nhr"))?);
        let ov = nin_overhead(&self.o.title, &self.lmdb_name(), self.o.date.len()).div_ceil(8) * 8;
        self.cur = Some(Volume {
            stem,
            nhr,
            nsq,
            hdr_off: vec![0],
            seq_off: vec![1],
            amb_off: Vec::new(),
            hdr_size: 0,
            seq_size: 1,
            letters: 0,
            max_len: 0,
            idx_overhead: ov,
        });
        Ok(())
    }

    fn fits(&self, hdr: usize, data: usize) -> bool {
        let v = self.cur.as_ref().unwrap();
        fits(
            self.o.max_file_size,
            v.idx_overhead,
            v.amb_off.len() as u64,
            v.hdr_size,
            v.seq_size,
            hdr as u64,
            data as u64,
        )
    }

    /// Append one sequence (in OID order): its sequence and ambiguity bytes, and its deflines,
    /// which are encoded here (`.nhr`) and give its acc2oid keys.
    pub fn add(
        &mut self,
        len: u32,
        seq: &[u8],
        amb: &[u8],
        deflines: &[Defline],
    ) -> Result<(), Error> {
        let mut hdr = std::mem::take(&mut self.hdr);
        super::asn::encode(deflines, &mut hdr);
        let r = self.add_encoded(len, seq, amb, &hdr, deflines);
        self.hdr = hdr;
        r
    }

    fn add_encoded(
        &mut self,
        len: u32,
        seq: &[u8],
        amb: &[u8],
        hdr: &[u8],
        deflines: &[Defline],
    ) -> Result<(), Error> {
        if self.cur.is_none() {
            self.open_volume()?;
        }
        if !self.fits(hdr.len(), seq.len() + amb.len()) {
            self.close_volume()?;
            self.open_volume()?;
        }
        let v = self.cur.as_mut().unwrap();
        let werr = |p: PathBuf| move |e| Error::io(p, e);
        v.nhr
            .write_all(hdr)
            .map_err(werr(with_ext(&v.stem, "nhr")))?;
        v.hdr_size += hdr.len() as u64;
        v.nsq
            .write_all(seq)
            .map_err(werr(with_ext(&v.stem, "nsq")))?;
        let amb_start = v.seq_size + seq.len() as u64;
        v.nsq
            .write_all(amb)
            .map_err(werr(with_ext(&v.stem, "nsq")))?;
        v.seq_size = amb_start + amb.len() as u64;
        if v.seq_size > u32::MAX as u64 || v.hdr_size > u32::MAX as u64 {
            return Err(Error::InvalidOptions(
                "volume offset exceeds 4 bytes; lower --max-file-size".into(),
            ));
        }
        v.hdr_off.push(v.hdr_size as u32);
        v.amb_off.push(amb_start as u32);
        v.seq_off.push(v.seq_size as u32);
        v.letters += len as u64;
        v.max_len = v.max_len.max(len);
        self.total_letters += len as u64;
        self.max_len = self.max_len.max(len);
        self.keys.push_deflines(deflines, self.oid)?;
        self.oid += 1;
        Ok(())
    }

    fn close_volume(&mut self) -> Result<(), Error> {
        let Some(mut v) = self.cur.take() else {
            return Ok(());
        };
        v.nhr
            .flush()
            .map_err(|e| Error::io(with_ext(&v.stem, "nhr"), e))?;
        v.nsq
            .flush()
            .map_err(|e| Error::io(with_ext(&v.stem, "nsq"), e))?;
        let idx = self.done.len() as u32;
        let n = v.amb_off.len() as u32;
        let w = nin_bytes(&NinFields {
            volume: idx,
            title: &self.o.title,
            lmdb_name: &self.lmdb_name(),
            date: &self.o.date,
            letters: v.letters,
            max_len: v.max_len,
            hdr_off: &v.hdr_off,
            seq_off: &v.seq_off,
            amb_off: &v.amb_off,
        });
        debug_assert_eq!(n as usize + 1, v.hdr_off.len());
        let p = with_ext(&v.stem, "nin");
        self.out.record(&p);
        std::fs::write(&p, &w).map_err(|e| Error::io(&p, e))?;
        self.done.push((v.stem, n));
        self.done_letters.push(v.letters);
        Ok(())
    }

    /// makeblastdb's rename for more than 100 volumes: a volume whose number has fewer digits
    /// than the count gets zeros before its two-digit name (`db.07` becomes `db.007`, and
    /// `db.100` keeps its name).
    fn widen_volume_names(&mut self) -> Result<(), Error> {
        let digits = self.done.len().ilog10() + 1;
        let pad = "0".repeat(digits as usize - 2);
        for (i, (stem, _)) in self.done.iter_mut().enumerate() {
            if i > 0 && i.ilog10() + 1 == digits {
                continue;
            }
            let wide = PathBuf::from(format!("{}.{pad}{i:02}", self.o.base.display()));
            for e in ["nin", "nhr", "nsq"] {
                let (from, to) = (with_ext(stem, e), with_ext(&wide, e));
                self.out.record(&to);
                std::fs::rename(&from, &to).map_err(|e| Error::io(&to, e))?;
            }
            *stem = wide;
        }
        Ok(())
    }

    /// Close the last volume, rename a lone volume to the base name, widen the names of more
    /// than 100 volumes, write `.nal` for several volumes. The LMDB and OID files are written by the caller from the result.
    pub fn finish(mut self) -> Result<Written, Error> {
        if self.cur.is_none() && self.done.is_empty() {
            // an empty database still has one (empty) volume
            self.open_volume()?;
        }
        self.close_volume()?;
        let mut files = Vec::new();
        let exts = ["nin", "nhr", "nsq"];
        if self.done.len() == 1 {
            let stem = self.done[0].0.clone();
            for e in exts {
                let (from, to) = (with_ext(&stem, e), with_ext(&self.o.base, e));
                self.out.record(&to);
                std::fs::rename(&from, &to).map_err(|e| Error::io(&to, e))?;
                files.push(to);
            }
            self.done[0].0 = self.o.base.clone();
        } else {
            if self.done.len() > 100 {
                self.widen_volume_names()?;
            }
            for (s, _) in &self.done {
                for e in exts {
                    files.push(with_ext(s, e));
                }
            }
            let list: Vec<String> = self.done.iter().map(|(s, _)| file_name(s)).collect();
            let nal = format!(
                "#\n# Alias file created: {}\n#\nTITLE {}\nDBLIST {}\n",
                self.o.date,
                self.o.title,
                list.join(" ")
            );
            let p = with_ext(&self.o.base, "nal");
            self.out.record(&p);
            std::fs::write(&p, nal).map_err(|e| Error::io(&p, e))?;
            files.push(p);
        }
        Ok(Written {
            files,
            volumes: self.done.iter().map(|(s, n)| (file_name(s), *n)).collect(),
            volume_letters: self.done_letters.clone(),
            keys: std::mem::take(&mut self.keys),
            n_oids: self.oid,
            total_letters: self.total_letters,
            max_len: self.max_len,
        })
    }
}

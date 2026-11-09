//! Reading a BLAST version 5 nucleotide database back: `.nal`, `.nin`, `.nsq` (with IUPAC
//! restored from the ambiguity tables), `.nhr` (BER deflines) and the LMDB files. Used by
//! `verify`, the tests and the golden canary; it needs no BLAST+.

use super::output::with_ext;
use crate::convert::NA4_TO_IUPAC;
use crate::error::Error;
use crate::plan::{Defline, SeqId};
use std::fs::File;
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};

pub use super::lmdb::{Ndb, dump as dump_lmdb};

/// The parsed fields of a `.nin` header.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NinHeader {
    pub version: u32,
    pub dbtype: u32,
    pub volume: u32,
    pub title: String,
    pub lmdb_name: String,
    /// Date text with the NUL padding removed.
    pub date: String,
    pub n_oids: u32,
    pub letters: u64,
    pub max_len: u32,
}

pub struct Volume {
    pub stem: PathBuf,
    pub header: NinHeader,
    pub first_oid: u32,
    hdr_off: Vec<u32>,
    seq_off: Vec<u32>,
    amb_off: Vec<u32>,
    nhr: File,
    nsq: File,
}

struct Cur<'a> {
    b: &'a [u8],
    p: usize,
    path: &'a Path,
}

impl Cur<'_> {
    fn bad(&self, what: &str) -> Error {
        Error::Store(anyhow::anyhow!(
            "{}: malformed .nin ({what} at byte {})",
            self.path.display(),
            self.p
        ))
    }
    fn take(&mut self, n: usize) -> Result<&[u8], Error> {
        if self.p + n > self.b.len() {
            return Err(self.bad("truncated"));
        }
        let s = &self.b[self.p..self.p + n];
        self.p += n;
        Ok(s)
    }
    fn be32(&mut self) -> Result<u32, Error> {
        Ok(u32::from_be_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn str(&mut self) -> Result<String, Error> {
        let n = self.be32()? as usize;
        Ok(String::from_utf8_lossy(self.take(n)?)
            .trim_end_matches('\0')
            .to_string())
    }
}

/// Parse a `.nin` file: header and the three offset arrays.
/// A parsed `.nin`: header, header offsets, sequence offsets, ambiguity offsets.
pub type Nin = (NinHeader, Vec<u32>, Vec<u32>, Vec<u32>);

pub fn parse_nin(path: &Path) -> Result<Nin, Error> {
    let b = std::fs::read(path).map_err(|e| Error::io(path, e))?;
    let mut c = Cur { b: &b, p: 0, path };
    let version = c.be32()?;
    let dbtype = c.be32()?;
    if version != 5 {
        return Err(c.bad("not a version 5 database"));
    }
    let volume = c.be32()?;
    let title = c.str()?;
    let lmdb_name = c.str()?;
    let date = c.str()?;
    let n = c.be32()?;
    let letters = u64::from_le_bytes(c.take(8)?.try_into().unwrap());
    let max_len = c.be32()?;
    let mut arr = |k: usize| -> Result<Vec<u32>, Error> { (0..k).map(|_| c.be32()).collect() };
    let hdr = arr(n as usize + 1)?;
    let seq = arr(n as usize + 1)?;
    let amb = arr(n as usize)?;
    let h = NinHeader {
        version,
        dbtype,
        volume,
        title,
        lmdb_name,
        date,
        n_oids: n,
        letters,
        max_len,
    };
    Ok((h, hdr, seq, amb))
}

pub struct DbReader {
    pub base: PathBuf,
    pub volumes: Vec<Volume>,
}

fn open_file(p: &Path) -> Result<File, Error> {
    File::open(p).map_err(|e| Error::io(p, e))
}

fn read_at(f: &File, path: &Path, start: u64, end: u64) -> Result<Vec<u8>, Error> {
    let mut v = vec![0u8; end.saturating_sub(start) as usize];
    f.read_exact_at(&mut v, start)
        .map_err(|e| Error::io(path, e))?;
    Ok(v)
}

impl DbReader {
    /// Open the database at `base`: the volumes listed in `<base>.nal`, else `<base>.nin`.
    pub fn open(base: &Path) -> Result<Self, Error> {
        let dir = base.parent().map(Path::to_path_buf).unwrap_or_default();
        let nal = with_ext(base, "nal");
        let stems: Vec<PathBuf> = if nal.exists() {
            let t = std::fs::read_to_string(&nal).map_err(|e| Error::io(&nal, e))?;
            let list = t
                .lines()
                .find_map(|l| l.strip_prefix("DBLIST "))
                .ok_or_else(|| {
                    Error::Store(anyhow::anyhow!("{}: no DBLIST line", nal.display()))
                })?;
            list.split_whitespace().map(|v| dir.join(v)).collect()
        } else {
            vec![base.to_path_buf()]
        };
        let mut volumes = Vec::new();
        let mut first = 0u32;
        for stem in stems {
            let (header, hdr_off, seq_off, amb_off) = parse_nin(&with_ext(&stem, "nin"))?;
            let n = header.n_oids;
            volumes.push(Volume {
                nhr: open_file(&with_ext(&stem, "nhr"))?,
                nsq: open_file(&with_ext(&stem, "nsq"))?,
                stem,
                header,
                first_oid: first,
                hdr_off,
                seq_off,
                amb_off,
            });
            first += n;
        }
        Ok(DbReader {
            base: base.to_path_buf(),
            volumes,
        })
    }

    pub fn num_oids(&self) -> u32 {
        self.volumes.iter().map(|v| v.header.n_oids).sum()
    }

    pub fn total_letters(&self) -> u64 {
        self.volumes.iter().map(|v| v.header.letters).sum()
    }

    pub fn max_len(&self) -> u32 {
        self.volumes
            .iter()
            .map(|v| v.header.max_len)
            .max()
            .unwrap_or(0)
    }

    fn locate(&self, oid: u32) -> Result<(&Volume, usize), Error> {
        for v in &self.volumes {
            if oid >= v.first_oid && oid < v.first_oid + v.header.n_oids {
                return Ok((v, (oid - v.first_oid) as usize));
            }
        }
        Err(Error::Store(anyhow::anyhow!(
            "OID {oid} is not in {}",
            self.base.display()
        )))
    }

    /// The `.nsq` sequence bytes and ambiguity bytes of an OID.
    pub fn raw(&self, oid: u32) -> Result<(Vec<u8>, Vec<u8>), Error> {
        let (v, i) = self.locate(oid)?;
        let p = with_ext(&v.stem, "nsq");
        let (s, a, e) = (
            v.seq_off[i] as u64,
            v.amb_off[i] as u64,
            v.seq_off[i + 1] as u64,
        );
        Ok((read_at(&v.nsq, &p, s, a)?, read_at(&v.nsq, &p, a, e)?))
    }

    /// Size of an OID's ambiguity table in bytes.
    pub fn ambiguity_bytes(&self, oid: u32) -> Result<u64, Error> {
        let (v, i) = self.locate(oid)?;
        Ok((v.seq_off[i + 1] - v.amb_off[i]) as u64)
    }

    /// Sequence length, from the `.nsq` bytes (the low bits of the last byte hold len % 4).
    pub fn length(&self, oid: u32) -> Result<u32, Error> {
        let (v, i) = self.locate(oid)?;
        let (s, a) = (v.seq_off[i] as u64, v.amb_off[i] as u64);
        let last = read_at(&v.nsq, &with_ext(&v.stem, "nsq"), a - 1, a)?[0];
        Ok(((a - s - 1) * 4) as u32 + (last & 3) as u32)
    }

    /// The sequence as uppercase IUPAC: ncbi2na decoded, then the ambiguity runs overwritten
    /// with their ncbi4na letter.
    pub fn sequence(&self, oid: u32) -> Result<Vec<u8>, Error> {
        let (seq, amb) = self.raw(oid)?;
        let last = *seq
            .last()
            .ok_or_else(|| Error::Store(anyhow::anyhow!("OID {oid}: empty sequence record")))?;
        let len = (seq.len() - 1) * 4 + (last & 3) as usize;
        let mut out = Vec::with_capacity(len);
        const NA2: [u8; 4] = *b"ACGT";
        for p in 0..len {
            out.push(NA2[((seq[p / 4] >> (6 - 2 * (p % 4))) & 3) as usize]);
        }
        for (v, start, run) in ambiguity_runs(&amb)? {
            let (s, e) = (start as usize, start as usize + run as usize);
            if e > len {
                return Err(Error::Store(anyhow::anyhow!(
                    "OID {oid}: ambiguity run past the end"
                )));
            }
            out[s..e].fill(NA4_TO_IUPAC[v as usize & 15]);
        }
        Ok(out)
    }

    pub fn header_bytes(&self, oid: u32) -> Result<Vec<u8>, Error> {
        let (v, i) = self.locate(oid)?;
        read_at(
            &v.nhr,
            &with_ext(&v.stem, "nhr"),
            v.hdr_off[i] as u64,
            v.hdr_off[i + 1] as u64,
        )
    }

    pub fn deflines(&self, oid: u32) -> Result<Vec<Defline>, Error> {
        decode_deflines(&self.header_bytes(oid)?)
    }

    /// The `.ndb` of this database: the LMDB its volumes name (else `<base>.ndb`).
    pub fn ndb(&self) -> Result<Ndb, Error> {
        match self.volumes.first() {
            Some(v) if !v.header.lmdb_name.is_empty() => {
                let dir = v.stem.parent().map(Path::to_path_buf).unwrap_or_default();
                Ndb::open(&dir.join(&v.header.lmdb_name))
            }
            _ => Ndb::open(&with_ext(&self.base, "ndb")),
        }
    }
}

/// (ncbi4na value, start, length) runs of an ambiguity table, old or new format.
pub fn ambiguity_runs(amb: &[u8]) -> Result<Vec<(u8, u32, u32)>, Error> {
    if amb.is_empty() {
        return Ok(Vec::new());
    }
    let bad = || Error::Store(anyhow::anyhow!("malformed ambiguity table"));
    let word = u32::from_be_bytes(amb.get(..4).ok_or_else(bad)?.try_into().unwrap());
    let newf = word & 0x8000_0000 != 0;
    let words = (word & 0x7FFF_FFFF) as usize;
    if amb.len() != 4 + 4 * words {
        return Err(bad());
    }
    let mut out = Vec::new();
    if newf {
        for e in amb[4..].chunks_exact(8) {
            let v = e[0] >> 4;
            let run = (((e[0] & 15) as u32) << 8 | e[1] as u32) + 1;
            out.push((v, u32::from_be_bytes(e[4..8].try_into().unwrap()), run));
        }
    } else {
        for e in amb[4..].chunks_exact(4) {
            let v = e[0] >> 4;
            let run = (e[0] & 15) as u32 + 1;
            out.push((
                v,
                (e[1] as u32) << 16 | (e[2] as u32) << 8 | e[3] as u32,
                run,
            ));
        }
    }
    Ok(out)
}

// ---------------------------------------------------------------------------
// BER decoding of Blast-def-line-set (definite or indefinite lengths)
// ---------------------------------------------------------------------------

enum Node<'a> {
    Prim(u8, &'a [u8]),
    Cons(u8, Vec<Node<'a>>),
}

impl<'a> Node<'a> {
    fn tag(&self) -> u8 {
        match self {
            Node::Prim(t, _) | Node::Cons(t, _) => *t,
        }
    }
    fn children(&self) -> Result<&[Node<'a>], Error> {
        match self {
            Node::Cons(_, c) => Ok(c),
            Node::Prim(t, _) => Err(ber_err(&format!(
                "expected a constructed value, got tag 0x{t:02x}"
            ))),
        }
    }
    fn only(&self) -> Result<&Node<'a>, Error> {
        match self.children()? {
            [one] => Ok(one),
            _ => Err(ber_err("expected exactly one inner value")),
        }
    }
    fn string(&self) -> Result<String, Error> {
        match self {
            Node::Prim(0x1A, b) => Ok(String::from_utf8_lossy(b).into_owned()),
            _ => Err(ber_err("expected a VisibleString")),
        }
    }
    fn int(&self) -> Result<i64, Error> {
        match self {
            Node::Prim(0x02, b) if !b.is_empty() && b.len() <= 8 => {
                let mut v: i64 = if b[0] & 0x80 != 0 { -1 } else { 0 };
                for &x in *b {
                    v = (v << 8) | x as i64;
                }
                Ok(v)
            }
            _ => Err(ber_err("expected an INTEGER")),
        }
    }
}

fn ber_err(msg: &str) -> Error {
    Error::Store(anyhow::anyhow!("bad .nhr defline data: {msg}"))
}

fn parse_node<'a>(b: &'a [u8], p: &mut usize) -> Result<Node<'a>, Error> {
    let tag = *b.get(*p).ok_or_else(|| ber_err("truncated"))?;
    *p += 1;
    let l0 = *b.get(*p).ok_or_else(|| ber_err("truncated"))?;
    *p += 1;
    let constructed = tag & 0x20 != 0;
    if l0 == 0x80 {
        if !constructed {
            return Err(ber_err("indefinite length on a primitive value"));
        }
        let mut kids = Vec::new();
        loop {
            if b.get(*p..*p + 2) == Some(&[0, 0]) {
                *p += 2;
                return Ok(Node::Cons(tag, kids));
            }
            kids.push(parse_node(b, p)?);
        }
    }
    let len = if l0 & 0x80 == 0 {
        l0 as usize
    } else {
        let n = (l0 & 0x7F) as usize;
        let mut v = 0usize;
        for _ in 0..n {
            v = (v << 8) | *b.get(*p).ok_or_else(|| ber_err("truncated"))? as usize;
            *p += 1;
        }
        v
    };
    let end = *p + len;
    if end > b.len() {
        return Err(ber_err("truncated"));
    }
    if constructed {
        let mut kids = Vec::new();
        while *p < end {
            kids.push(parse_node(b, p)?);
        }
        Ok(Node::Cons(tag, kids))
    } else {
        let s = &b[*p..end];
        *p = end;
        Ok(Node::Prim(tag, s))
    }
}

fn object_id(n: &Node) -> Result<String, Error> {
    let inner = n.only()?;
    match inner.tag() {
        0xA0 => Ok(inner.only()?.int()?.to_string()),
        0xA1 => inner.only()?.string(),
        t => Err(ber_err(&format!("unexpected Object-id choice 0x{t:02x}"))),
    }
}

fn seq_id(n: &Node) -> Result<SeqId, Error> {
    match n.tag() {
        0xA0 => Ok(SeqId::Local(object_id(n)?)),
        0xAA => {
            let dbtag = n.only()?;
            let mut db = None;
            let mut tag = None;
            for c in dbtag.children()? {
                match c.tag() {
                    0xA0 => db = Some(c.only()?.string()?),
                    0xA1 => tag = Some(object_id(c)?),
                    _ => {}
                }
            }
            Ok(SeqId::General {
                db: db.ok_or_else(|| ber_err("Dbtag without db"))?,
                tag: tag.ok_or_else(|| ber_err("Dbtag without tag"))?,
            })
        }
        t => Err(ber_err(&format!(
            "unsupported Seq-id choice [{}]",
            t & 0x1F
        ))),
    }
}

/// Decode a Blast-def-line-set. Integer Object-ids come back as their decimal string.
pub fn decode_deflines(ber: &[u8]) -> Result<Vec<Defline>, Error> {
    let mut p = 0;
    let top = parse_node(ber, &mut p)?;
    let mut out = Vec::new();
    for d in top.children()? {
        let mut title = String::new();
        let mut ids = Vec::new();
        let mut taxid = None;
        for f in d.children()? {
            match f.tag() {
                0xA0 => title = f.only()?.string()?,
                0xA1 => {
                    for id in f.only()?.children()? {
                        ids.push(seq_id(id)?);
                    }
                }
                0xA2 => taxid = Some(f.only()?.int()? as i32),
                _ => {}
            }
        }
        out.push(Defline { ids, title, taxid });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blastdb::asn;

    #[test]
    fn encode_decode_roundtrip() {
        let long_title = "t".repeat(300);
        let d = vec![
            Defline {
                ids: vec![
                    SeqId::Local("SQ.abc".into()),
                    SeqId::Local("NC_000001.11".into()),
                    SeqId::Local("123456".into()),
                    SeqId::Local("0099".into()),
                ],
                title: "md5:x length:5 collections:1".into(),
                taxid: Some(0),
            },
            Defline {
                ids: vec![
                    SeqId::General {
                        db: "hg38".into(),
                        tag: "chr1".into(),
                    },
                    SeqId::General {
                        db: "hg38".into(),
                        tag: "22".into(),
                    },
                    SeqId::General {
                        db: "hg38".into(),
                        tag: "007".into(),
                    },
                    SeqId::General {
                        db: "-e70".into(),
                        tag: "2147483647".into(),
                    },
                ],
                title: long_title,
                taxid: Some(9606),
            },
            // the field omitted (an SQ. defline whose collections disagree)
            Defline {
                ids: vec![SeqId::Local("SQ.mixed".into())],
                title: "Homo sapiens, Mus musculus chr1; 2 assemblies".into(),
                taxid: None,
            },
        ];
        let mut b = Vec::new();
        asn::encode(&d, &mut b);
        assert_eq!(decode_deflines(&b).unwrap(), d);
    }

    #[test]
    fn ambiguity_tables_parse() {
        let old = [0, 0, 0, 1, (15 << 4) | 2, 0, 0, 4];
        assert_eq!(ambiguity_runs(&old).unwrap(), vec![(15, 4, 3)]);
        let new = [0x80, 0, 0, 2, (15 << 4) | 0xF, 0xFE, 0, 0, 0, 0, 0, 9];
        assert_eq!(ambiguity_runs(&new).unwrap(), vec![(15, 9, 4095)]);
    }
}

//! A read-only view of an LMDB file: the newest meta page, its named databases, and their
//! B+trees (branch and leaf pages, duplicate sub-pages and sub-databases, overflow values).
//! Every offset is checked, so a damaged file gives [`LmdbError::Corrupted`], not a panic.

use super::LmdbError;
use super::page::*;
use memmap2::Mmap;
use std::cmp::Ordering;
use std::fs::File;
use std::path::Path;

/// Deeper than any real tree; stops a corrupt file from looping.
const MAX_DEPTH: usize = 64;

fn bad(what: &str) -> LmdbError {
    LmdbError::Corrupted(what.to_string())
}

/// A leaf node: its key, flags, and value (data, a sub-page, or a sub-database record).
struct Node<'a> {
    key: &'a [u8],
    flags: u16,
    data: &'a [u8],
}

/// A visitor of (key, value) pairs.
pub type EachFn<'a> = dyn FnMut(&[u8], &[u8]) -> Result<(), LmdbError> + 'a;

pub struct Reader {
    map: Mmap,
    psize: usize,
    main: Db,
}

impl Reader {
    pub fn open(path: &Path) -> Result<Self, LmdbError> {
        let file = File::open(path)?;
        let len = file.metadata()?.len();
        if len < (PAGEHDRSZ + META_SIZE) as u64 {
            return Err(bad("file too short for an LMDB meta page"));
        }
        // SAFETY: the file is opened read-only and only read; the files this crate reads are
        // not modified while open.
        let map = unsafe { Mmap::map(&file)? };
        let meta = |off: usize| -> Option<(u64, Db, u32)> {
            let m = map.get(off + PAGEHDRSZ..off + PAGEHDRSZ + META_SIZE)?;
            if r32(m, 0) != MDB_MAGIC || r32(m, 4) != MDB_DATA_VERSION {
                return None;
            }
            let free = Db::from_bytes(&m[24..])?;
            Some((r64(m, 128), Db::from_bytes(&m[72..])?, free.pad))
        };
        let (t0, main0, psize) = meta(0).ok_or_else(|| bad("not an LMDB file (meta page 0)"))?;
        let psize = psize as usize;
        if !psize.is_power_of_two() || !(256..=MAX_PAGESIZE).contains(&psize) {
            return Err(bad("page size"));
        }
        let main = match meta(psize) {
            Some((t1, main1, _)) if t1 > t0 => main1,
            _ => main0,
        };
        Ok(Reader { map, psize, main })
    }

    fn bytes(&self) -> &[u8] {
        &self.map
    }

    fn page(&self, pgno: u64) -> Result<&[u8], LmdbError> {
        let off = (pgno as usize)
            .checked_mul(self.psize)
            .ok_or_else(|| bad("page number"))?;
        let p = self
            .bytes()
            .get(off..off + self.psize)
            .ok_or_else(|| bad("page beyond the end of the file"))?;
        if r64(p, 0) != pgno {
            return Err(bad("page number mismatch"));
        }
        Ok(p)
    }

    /// The key count of a (sub-)page at `b` in `p`, checked against its size `len`.
    fn count(p: &[u8], b: usize, len: usize) -> Result<usize, LmdbError> {
        if len < PAGEHDRSZ || b + len > p.len() {
            return Err(bad("page header"));
        }
        let (l, u) = (pg_lower(p, b), pg_upper(p, b));
        if l < PAGEHDRSZ || l > u || u > len {
            return Err(bad("page bounds"));
        }
        Ok(numkeys(p, b))
    }

    /// Node `i` of the whole page `p`.
    fn node<'a>(&'a self, p: &'a [u8], i: usize) -> Result<Node<'a>, LmdbError> {
        let o = ptr(p, 0, i);
        if o < PAGEHDRSZ || o + NODESIZE > p.len() {
            return Err(bad("node offset"));
        }
        let (flags, ks) = (n_flags(p, o), n_ksize(p, o));
        let k = o + NODESIZE;
        let key = p.get(k..k + ks).ok_or_else(|| bad("node key"))?;
        let data = if pg_flags(p, 0) & P_LEAF == 0 {
            &[][..]
        } else if flags & F_BIGDATA != 0 {
            let pg = p.get(k + ks..k + ks + 8).ok_or_else(|| bad("node data"))?;
            let ov = self.page(r64(pg, 0))?;
            let len = n_dsize(p, o);
            let start = r64(pg, 0) as usize * self.psize + PAGEHDRSZ;
            if pg_flags(ov, 0) & P_OVERFLOW == 0 {
                return Err(bad("overflow page"));
            }
            self.bytes()
                .get(start..start + len)
                .ok_or_else(|| bad("overflow data"))?
        } else {
            p.get(k + ks..k + ks + n_dsize(p, o))
                .ok_or_else(|| bad("node data"))?
        };
        Ok(Node { key, flags, data })
    }

    /// The main database's records: (name, record) for each named database.
    pub fn tables(&self) -> Result<Vec<(Vec<u8>, Db)>, LmdbError> {
        let mut out = Vec::new();
        self.walk(&self.main, 0, &mut |n: &Node| {
            if n.flags & (F_SUBDATA | F_DUPDATA) != F_SUBDATA {
                return Err(bad("main database entry is not a database"));
            }
            let db = Db::from_bytes(n.data).ok_or_else(|| bad("database record"))?;
            out.push((n.key.to_vec(), db));
            Ok(())
        })?;
        Ok(out)
    }

    /// The record of the named database `name`.
    pub fn table(&self, name: &str) -> Result<Option<Db>, LmdbError> {
        Ok(self
            .tables()?
            .into_iter()
            .find(|(n, _)| n == name.as_bytes())
            .map(|(_, d)| d))
    }

    /// Visit every leaf node of a tree in key order.
    fn walk(
        &self,
        db: &Db,
        depth: usize,
        f: &mut dyn FnMut(&Node) -> Result<(), LmdbError>,
    ) -> Result<(), LmdbError> {
        if db.root == P_INVALID {
            return Ok(());
        }
        self.walk_page(db.root, db, depth, f)
    }

    fn walk_page(
        &self,
        pgno: u64,
        db: &Db,
        depth: usize,
        f: &mut dyn FnMut(&Node) -> Result<(), LmdbError>,
    ) -> Result<(), LmdbError> {
        if depth > MAX_DEPTH {
            return Err(bad("tree too deep"));
        }
        let p = self.page(pgno)?;
        let n = Self::count(p, 0, self.psize)?;
        let fl = pg_flags(p, 0);
        if fl & P_BRANCH != 0 {
            for i in 0..n {
                let o = ptr(p, 0, i);
                if o + NODESIZE > p.len() {
                    return Err(bad("node offset"));
                }
                self.walk_page(n_pgno(p, o), db, depth + 1, f)?;
            }
        } else if fl & P_LEAF2 != 0 {
            let ks = db.pad as usize;
            for i in 0..n {
                let o = PAGEHDRSZ + i * ks;
                let key = p.get(o..o + ks).ok_or_else(|| bad("fixed-size key"))?;
                f(&Node {
                    key,
                    flags: 0,
                    data: &[],
                })?;
            }
        } else if fl & P_LEAF != 0 {
            for i in 0..n {
                f(&self.node(p, i)?)?;
            }
        } else {
            return Err(bad("page type"));
        }
        Ok(())
    }

    /// The values of a leaf node: its data, or each duplicate.
    fn values(
        &self,
        n: &Node,
        depth: usize,
        f: &mut dyn FnMut(&[u8]) -> Result<(), LmdbError>,
    ) -> Result<(), LmdbError> {
        if n.flags & F_DUPDATA == 0 {
            return f(n.data);
        }
        if n.flags & F_SUBDATA != 0 {
            let sub = Db::from_bytes(n.data).ok_or_else(|| bad("sub-database record"))?;
            return self.walk(&sub, depth + 1, &mut |d: &Node| f(d.key));
        }
        let sp = n.data;
        let k = Self::count(sp, 0, sp.len())?;
        if pg_flags(sp, 0) & P_LEAF2 != 0 {
            let ks = pg_pad(sp, 0) as usize;
            for i in 0..k {
                let o = PAGEHDRSZ + i * ks;
                f(sp.get(o..o + ks).ok_or_else(|| bad("sub-page key"))?)?;
            }
        } else {
            for i in 0..k {
                let o = ptr(sp, 0, i);
                let ks = sp
                    .get(o + 6..o + 8)
                    .map(|_| n_ksize(sp, o))
                    .ok_or_else(|| bad("sub-page node"))?;
                f(sp.get(o + NODESIZE..o + NODESIZE + ks)
                    .ok_or_else(|| bad("sub-page key"))?)?;
            }
        }
        Ok(())
    }

    /// Every (key, value) of a database in order, duplicates in their sorted order.
    pub fn each(&self, db: &Db, f: &mut EachFn) -> Result<(), LmdbError> {
        self.walk(db, 0, &mut |n: &Node| {
            self.values(n, 0, &mut |v| f(n.key, v))
        })
    }

    /// The values stored under `key` (none if absent).
    pub fn get(&self, db: &Db, key: &[u8]) -> Result<Vec<Vec<u8>>, LmdbError> {
        let mut out = Vec::new();
        if db.root == P_INVALID {
            return Ok(out);
        }
        let cmp = Cmp::key(db.flags);
        let mut pgno = db.root;
        for _ in 0..=MAX_DEPTH {
            let p = self.page(pgno)?;
            let n = Self::count(p, 0, self.psize)?;
            let fl = pg_flags(p, 0);
            let leaf = fl & P_LEAF != 0;
            let (mut low, mut high) = (if leaf { 0isize } else { 1 }, n as isize - 1);
            let (mut i, mut rc) = (0isize, Ordering::Equal);
            while low <= high {
                i = (low + high) >> 1;
                let nk = if fl & P_LEAF2 != 0 {
                    let o = PAGEHDRSZ + i as usize * db.pad as usize;
                    p.get(o..o + db.pad as usize)
                        .ok_or_else(|| bad("fixed-size key"))?
                } else {
                    let o = ptr(p, 0, i as usize);
                    let ks = p
                        .get(o..o + NODESIZE)
                        .map(|_| n_ksize(p, o))
                        .ok_or_else(|| bad("node"))?;
                    p.get(o + NODESIZE..o + NODESIZE + ks)
                        .ok_or_else(|| bad("node key"))?
                };
                rc = cmp.cmp(key, nk);
                match rc {
                    Ordering::Equal => break,
                    Ordering::Greater => low = i + 1,
                    Ordering::Less => high = i - 1,
                }
            }
            if leaf {
                if rc == Ordering::Equal && n > 0 && fl & P_LEAF2 == 0 {
                    let node = self.node(p, i as usize)?;
                    self.values(&node, 0, &mut |v| {
                        out.push(v.to_vec());
                        Ok(())
                    })?;
                }
                return Ok(out);
            }
            if n == 0 {
                return Err(bad("empty branch page"));
            }
            let c = match rc {
                Ordering::Greater if i as usize + 1 >= n => n - 1,
                Ordering::Greater => i as usize,
                Ordering::Less => (i as usize).saturating_sub(1),
                Ordering::Equal => i as usize,
            };
            let o = ptr(p, 0, c);
            if o + NODESIZE > p.len() {
                return Err(bad("node offset"));
            }
            pgno = n_pgno(p, o);
        }
        Err(bad("tree too deep"))
    }
}

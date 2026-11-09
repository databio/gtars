//! The LMDB write path: a new environment, one write transaction of puts, and its commit.
//! A port of the parts of `mdb.c` this crate exercises (see [`super`] for the copyright and
//! license): `mdb_env_open` of a new file, `mdb_txn_begin`/`commit`/`abort`, `mdb_dbi_open`
//! with `MDB_CREATE`, and `mdb_put` with flags 0, including `MDB_DUPSORT|MDB_DUPFIXED`
//! duplicates (sub-pages and sub-databases), `mdb_page_split`, and the dirty-page spill of
//! large transactions. The file it writes is byte-identical to what liblmdb writes for the
//! same calls, including the bytes LMDB leaves unset: LMDB reuses freed page buffers without
//! clearing them (`mdb_page_malloc`), so the unused middle of a page holds whatever that
//! buffer held before; this module keeps the same buffers in the same order to reproduce it.
//!
//! Not ported (refused with [`LmdbError::Unsupported`]): overflow pages (values larger than
//! about half a page), `MDB_DUPSORT` without `MDB_DUPFIXED`, reverse-order comparators, and
//! cursor operations beyond `put`. None of `.ndb`/`.ntf` needs them.

use super::LmdbError;
use super::page::*;
use std::cmp::Ordering;
use std::fs::File;
use std::io::{Seek, SeekFrom, Write};
use std::path::Path;

/// `MDB_MAXKEYSIZE` (the default build).
const MDB_MAXKEYSIZE: usize = 511;
/// `MDB_IDL_LOGN`, a liblmdb build option: 16 by default. Unit tests use 8 so that small
/// files exercise the dirty-page spill.
const MDB_IDL_LOGN: u32 = if cfg!(test) { 8 } else { 16 };
/// `MDB_IDL_UM_MAX`: the dirty pages a transaction may hold.
const MDB_IDL_UM_MAX: u64 = (1 << (MDB_IDL_LOGN + 1)) - 1;
/// `MDB_NOSUBDIR`, the low environment flag recorded in the meta pages.
const MDB_NOSUBDIR: u16 = 0x4000;
/// `NODE_ADD_FLAGS` that are node flags (`MDB_RESERVE`/`MDB_APPEND` are never used here).
const NODE_ADD_FLAGS: u16 = F_DUPDATA | F_SUBDATA;

/// The system page size (`sysconf(_SC_PAGESIZE)`), which LMDB uses for a new file.
fn os_page_size() -> usize {
    #[cfg(unix)]
    {
        // SAFETY: sysconf has no preconditions.
        let n = unsafe { libc::sysconf(libc::_SC_PAGESIZE) };
        if n > 0 {
            return n as usize;
        }
    }
    4096
}

/// A node's value: data of a leaf node, or the child page of a branch node.
#[derive(Clone, Copy)]
enum Val<'a> {
    Data(&'a [u8]),
    Pgno(u64),
}

/// A page of the transaction: dirty (in memory, to be written at commit) or spilled
/// (written early by [`Env::spill_check`]; this is its content on disk).
enum Slot {
    Dirty(Box<[u8]>),
    Spilled(Box<[u8]>),
}

/// A cursor stack: (page, index) from the root to a leaf.
type Stack = Vec<(u64, usize)>;

/// An LMDB environment open for writing, with at most one transaction at a time.
pub struct Env {
    file: File,
    psize: usize,
    /// `me_nodemax`: the largest node kept on its page.
    nodemax: usize,
    /// `me_mapsize`.
    mapsize: u64,
    /// The map size recorded in meta page 0 when the file was created.
    meta0_mapsize: u64,
    /// The meta `MDB_db` records the transaction starts from.
    meta_dbs: [Db; 2],
    /// `me_pbuf`: scratch page for building duplicate sub-pages; never cleared.
    pbuf: Vec<u8>,
    /// `me_dpages`: freed page buffers, reused last-in first-out without clearing.
    dpages: Vec<Box<[u8]>>,
    // The transaction.
    maxpg: u64,
    next_pgno: u64,
    slots: Vec<Slot>,
    dirty_room: u64,
    spills: bool,
    dbs: Vec<Db>,
    names: Vec<Vec<u8>>,
}

/// `mdb_leaf_size`.
fn leaf_size(nodemax: usize, k: usize, d: usize) -> usize {
    let mut sz = NODESIZE + k + d;
    if sz > nodemax {
        sz -= d - 8;
    }
    even(sz + 2)
}

/// `mdb_branch_size`.
fn branch_size(k: usize) -> usize {
    NODESIZE + k + 2
}

/// `mdb_node_search` on the (sub-)page at `b`: the index of the first entry not below `key`,
/// whether it equals `key`, and whether such an entry exists (the node pointer was non-NULL).
fn node_search(p: &[u8], b: usize, key: &[u8], cmp: Cmp, ks2: usize) -> (usize, bool, bool) {
    let n = numkeys(p, b);
    let fl = pg_flags(p, b);
    let mut low: isize = if fl & P_LEAF != 0 { 0 } else { 1 };
    let mut high = n as isize - 1;
    let mut i: isize = 0;
    let mut rc = Ordering::Equal;
    while low <= high {
        i = (low + high) >> 1;
        let nk = if fl & P_LEAF2 != 0 {
            let o = b + PAGEHDRSZ + i as usize * ks2;
            &p[o..o + ks2]
        } else {
            n_key(p, b + ptr(p, b, i as usize))
        };
        rc = cmp.cmp(key, nk);
        match rc {
            Ordering::Equal => break,
            Ordering::Greater => low = i + 1,
            Ordering::Less => high = i - 1,
        }
    }
    if rc == Ordering::Greater {
        i += 1;
    }
    let i = i as usize;
    (i, rc == Ordering::Equal && n > 0, i < n)
}

/// `mdb_node_add` for a `P_LEAF2` (sub-)page at `b`.
fn leaf2_add(p: &mut [u8], b: usize, indx: usize, key: &[u8], ks: usize) {
    let n = numkeys(p, b);
    let off = b + PAGEHDRSZ + indx * ks;
    if n > indx {
        p.copy_within(off..off + (n - indx) * ks, off + ks);
    }
    p[off..off + ks].copy_from_slice(&key[..ks]);
    let (l, u) = (pg_lower(p, b), pg_upper(p, b));
    set_lower(p, b, l + 2);
    set_upper(p, b, u.wrapping_sub(ks).wrapping_add(2));
}

/// `mdb_node_add` on a whole page.
fn node_add(
    p: &mut [u8],
    nodemax: usize,
    ks2: usize,
    indx: usize,
    key: &[u8],
    val: Val,
    flags: u16,
) -> Result<(), LmdbError> {
    let fl = pg_flags(p, 0);
    if fl & P_LEAF2 != 0 {
        leaf2_add(p, 0, indx, key, ks2);
        return Ok(());
    }
    let room = sizeleft(p, 0) as isize - 2;
    let mut node_size = NODESIZE + key.len();
    if fl & P_LEAF != 0 {
        let Val::Data(d) = val else {
            return Err(LmdbError::Corrupted("leaf node without data".into()));
        };
        if flags & F_BIGDATA != 0 || node_size + d.len() > nodemax {
            return Err(LmdbError::Unsupported(format!(
                "a value of {} bytes needs an overflow page",
                d.len()
            )));
        }
        node_size += d.len();
    }
    node_size = even(node_size);
    if node_size as isize > room {
        return Err(LmdbError::Corrupted("page full (MDB_PAGE_FULL)".into()));
    }
    let n = numkeys(p, 0);
    for i in (indx + 1..=n).rev() {
        let v = ptr(p, 0, i - 1);
        set_ptr(p, 0, i, v);
    }
    let ofs = pg_upper(p, 0) - node_size;
    set_ptr(p, 0, indx, ofs);
    set_upper(p, 0, ofs);
    let l = pg_lower(p, 0);
    set_lower(p, 0, l + 2);
    w16(p, ofs + 6, key.len() as u16);
    w16(p, ofs + 4, flags);
    match val {
        Val::Data(d) => w32(p, ofs, d.len() as u32),
        Val::Pgno(pg) => {
            w32(p, ofs, pg as u32);
            w16(p, ofs + 4, (pg >> 32) as u16);
        }
    }
    let k = ofs + NODESIZE;
    p[k..k + key.len()].copy_from_slice(key);
    if let (true, Val::Data(d)) = (fl & P_LEAF != 0, val) {
        let o = k + key.len();
        p[o..o + d.len()].copy_from_slice(d);
    }
    Ok(())
}

/// `mdb_node_del` on a whole page (`ks2`: key size of a `P_LEAF2` page).
fn node_del(p: &mut [u8], indx: usize, ks2: usize) {
    let n = numkeys(p, 0);
    let (l, u) = (pg_lower(p, 0), pg_upper(p, 0));
    if pg_flags(p, 0) & P_LEAF2 != 0 {
        let x = n - 1 - indx;
        let base = PAGEHDRSZ + indx * ks2;
        if x > 0 {
            p.copy_within(base + ks2..base + ks2 + x * ks2, base);
        }
        set_lower(p, 0, l - 2);
        set_upper(p, 0, u + ks2 - 2);
        return;
    }
    let pi = ptr(p, 0, indx);
    let mut sz = NODESIZE + n_ksize(p, pi);
    if pg_flags(p, 0) & P_LEAF != 0 {
        sz += if n_flags(p, pi) & F_BIGDATA != 0 {
            8
        } else {
            n_dsize(p, pi)
        };
    }
    let sz = even(sz);
    let mut j = 0;
    for i in 0..n {
        if i != indx {
            let v = ptr(p, 0, i);
            set_ptr(p, 0, j, if v < pi { v + sz } else { v });
            j += 1;
        }
    }
    p.copy_within(u..pi, u + sz);
    set_lower(p, 0, l - 2);
    set_upper(p, 0, u + sz);
}

/// `mdb_page_copy`: the used parts of a page (header and pointers, then nodes), 8-byte aligned.
fn page_copy(dst: &mut [u8], src: &[u8], psize: usize) {
    let (lower, upper) = (pg_lower(src, 0), pg_upper(src, 0));
    let unused = ((upper as u16).wrapping_sub(lower as u16) as usize) & !7;
    if unused != 0 && pg_flags(src, 0) & P_LEAF2 == 0 {
        let up = upper & !7;
        let lo = (lower + 7) & !7;
        dst[..lo].copy_from_slice(&src[..lo]);
        dst[up..psize].copy_from_slice(&src[up..psize]);
    } else {
        dst[..psize - unused].copy_from_slice(&src[..psize - unused]);
    }
}

/// Insert `key` into the `P_LEAF2` sub-page at `b` of `p` (an `_mdb_cursor_put` on a sub-page
/// cursor). Returns false when the item was already there (a no-op).
fn subpage_put(p: &mut [u8], b: usize, key: &[u8], cmp: Cmp, ks: usize) -> Result<bool, LmdbError> {
    let (ki, exact, _) = node_search(p, b, key, cmp, ks);
    if exact {
        let o = b + PAGEHDRSZ + ki * ks;
        p[o..o + ks].copy_from_slice(key);
        return Ok(false);
    }
    if sizeleft(p, b) < ks {
        return Err(LmdbError::Corrupted("duplicate sub-page full".into()));
    }
    leaf2_add(p, b, ki, key, ks);
    Ok(true)
}

impl Env {
    /// `mdb_env_create` + `mdb_env_set_mapsize` + `mdb_env_open` of a new, empty file with
    /// `MDB_NOSUBDIR | MDB_NOLOCK`: writes the two meta pages.
    pub fn create(path: &Path, mapsize: u64) -> Result<Self, LmdbError> {
        let psize = os_page_size().min(MAX_PAGESIZE);
        let mapsize = mapsize.max(2 * psize as u64);
        let mut free = Db::empty(MDB_NOSUBDIR | MDB_INTEGERKEY);
        free.pad = psize as u32;
        let main = Db::empty(0);
        let mut meta = Vec::with_capacity(META_SIZE);
        meta.extend_from_slice(&MDB_MAGIC.to_ne_bytes());
        meta.extend_from_slice(&MDB_DATA_VERSION.to_ne_bytes());
        meta.extend_from_slice(&0u64.to_ne_bytes()); // mm_address
        meta.extend_from_slice(&mapsize.to_ne_bytes());
        meta.extend_from_slice(&free.to_bytes());
        meta.extend_from_slice(&main.to_bytes());
        meta.extend_from_slice(&1u64.to_ne_bytes()); // mm_last_pg
        meta.extend_from_slice(&0u64.to_ne_bytes()); // mm_txnid
        let mut pages = vec![0u8; 2 * psize];
        for (i, pg) in pages.chunks_mut(psize).enumerate() {
            w64(pg, 0, i as u64);
            set_flags(pg, 0, P_META);
            pg[PAGEHDRSZ..PAGEHDRSZ + META_SIZE].copy_from_slice(&meta);
        }
        let mut file = File::create(path)?;
        file.write_all(&pages)?;
        Ok(Env {
            file,
            psize,
            nodemax: (((psize - PAGEHDRSZ) / 2) & !1) - 2,
            mapsize,
            meta0_mapsize: mapsize,
            meta_dbs: [free, main],
            pbuf: vec![0u8; psize],
            dpages: Vec::new(),
            maxpg: 0,
            next_pgno: 2,
            slots: Vec::new(),
            dirty_room: 0,
            spills: false,
            dbs: Vec::new(),
            names: Vec::new(),
        })
    }

    /// `mdb_env_set_mapsize` on the open environment (no transaction active).
    pub fn set_mapsize(&mut self, size: u64) {
        self.mapsize = size.max(2 * self.psize as u64);
    }

    /// `mdb_txn_begin` of a write transaction on the fresh file (meta txnid 0).
    pub fn begin(&mut self) {
        self.maxpg = self.mapsize / self.psize as u64;
        self.next_pgno = 2;
        self.slots.clear();
        self.dirty_room = MDB_IDL_UM_MAX;
        self.spills = false;
        self.dbs = self.meta_dbs.to_vec();
        self.names = vec![Vec::new(), Vec::new()];
    }

    /// `mdb_txn_abort`: every dirty page buffer goes back to the free list, in page order.
    pub fn abort(&mut self) {
        for s in std::mem::take(&mut self.slots) {
            if let Slot::Dirty(p) = s {
                self.dpages.push(p);
            }
        }
    }

    /// `mdb_dbi_open(txn, name, flags | MDB_CREATE)` for a name not in the file yet.
    pub fn create_db(&mut self, name: &str, flags: u16) -> Result<usize, LmdbError> {
        if flags & !SUPPORTED_DB_FLAGS != 0
            || (flags & MDB_DUPFIXED != 0 && flags & MDB_DUPSORT == 0)
        {
            return Err(LmdbError::Unsupported(format!("database flags {flags:#x}")));
        }
        let key = name.as_bytes();
        let main = self.dbs[MAIN_DBI];
        if main.root != P_INVALID && self.search(&main, Cmp::Memn, key).1 {
            return Err(LmdbError::Unsupported(format!("database {name} exists")));
        }
        self.cursor_put(MAIN_DBI, key, &Db::empty(flags).to_bytes(), F_SUBDATA)?;
        self.dbs.push(Db::empty(flags));
        self.names.push(key.to_vec());
        Ok(self.dbs.len() - 1)
    }

    /// `mdb_put(txn, dbi, key, data, 0)`.
    pub fn put(&mut self, dbi: usize, key: &[u8], data: &[u8]) -> Result<(), LmdbError> {
        self.cursor_put(dbi, key, data, 0)
    }

    /// `mdb_txn_commit`: store the database records in the main database, write the dirty
    /// pages, then meta page 1 (the first transaction's txnid is 1).
    pub fn commit(&mut self) -> Result<(), LmdbError> {
        for dbi in 2..self.dbs.len() {
            let rec = self.dbs[dbi].to_bytes();
            let name = self.names[dbi].clone();
            self.cursor_put(MAIN_DBI, &name, &rec, F_SUBDATA)?;
        }
        let psize = self.psize;
        let mut w = std::io::BufWriter::with_capacity(1 << 20, &self.file);
        w.seek(SeekFrom::Start(2 * psize as u64))?;
        let mut skip = 0i64;
        for s in &mut self.slots {
            match s {
                Slot::Dirty(p) => {
                    if skip != 0 {
                        w.seek(SeekFrom::Current(skip))?;
                        skip = 0;
                    }
                    let f = pg_flags(p, 0);
                    set_flags(p, 0, f & !P_DIRTY);
                    w.write_all(p)?;
                }
                Slot::Spilled(_) => skip += psize as i64,
            }
        }
        w.flush()?;
        drop(w);
        self.file.sync_data()?;
        let mut m = Vec::with_capacity(META_SIZE - 16);
        m.extend_from_slice(&self.meta0_mapsize.max(self.mapsize).to_ne_bytes());
        m.extend_from_slice(&self.dbs[FREE_DBI].to_bytes());
        m.extend_from_slice(&self.dbs[MAIN_DBI].to_bytes());
        m.extend_from_slice(&(self.next_pgno - 1).to_ne_bytes());
        m.extend_from_slice(&1u64.to_ne_bytes());
        self.file
            .seek(SeekFrom::Start((psize + PAGEHDRSZ + 16) as u64))?;
        self.file.write_all(&m)?;
        self.file.sync_all()?;
        for s in std::mem::take(&mut self.slots) {
            if let Slot::Dirty(p) = s {
                self.dpages.push(p);
            }
        }
        Ok(())
    }

    fn pg(&self, pgno: u64) -> &[u8] {
        match &self.slots[(pgno - 2) as usize] {
            Slot::Dirty(p) | Slot::Spilled(p) => p,
        }
    }

    /// A dirty page (every page this transaction writes to was touched first).
    fn pg_mut(&mut self, pgno: u64) -> &mut [u8] {
        match &mut self.slots[(pgno - 2) as usize] {
            Slot::Dirty(p) => p,
            Slot::Spilled(_) => unreachable!("page {pgno} written while spilled"),
        }
    }

    /// `mdb_page_malloc(txn, 1)`: a freed buffer as it is, or a new zeroed one.
    fn malloc(&mut self) -> Box<[u8]> {
        self.dpages
            .pop()
            .unwrap_or_else(|| vec![0u8; self.psize].into_boxed_slice())
    }

    /// `mdb_page_alloc(mc, 1)` on a file without free pages.
    fn alloc(&mut self) -> Result<u64, LmdbError> {
        if self.dirty_room == 0 {
            return Err(LmdbError::Corrupted(
                "transaction full (MDB_TXN_FULL)".into(),
            ));
        }
        let pgno = self.next_pgno;
        if pgno + 1 >= self.maxpg {
            return Err(LmdbError::MapFull);
        }
        let mut p = self.malloc();
        w64(&mut p, 0, pgno);
        self.slots.push(Slot::Dirty(p));
        self.next_pgno = pgno + 1;
        self.dirty_room -= 1;
        Ok(pgno)
    }

    /// `mdb_page_new`.
    fn page_new(&mut self, db: &mut Db, flags: u16) -> Result<u64, LmdbError> {
        let pgno = self.alloc()?;
        let psize = self.psize;
        let p = self.pg_mut(pgno);
        set_flags(p, 0, flags | P_DIRTY);
        set_lower(p, 0, PAGEHDRSZ);
        set_upper(p, 0, psize);
        if flags & P_BRANCH != 0 {
            db.branch_pages += 1;
        } else if flags & P_LEAF != 0 {
            db.leaf_pages += 1;
        }
        Ok(pgno)
    }

    /// `mdb_page_touch` of a page of this transaction: a spilled page becomes dirty again in
    /// a fresh buffer (`mdb_page_unspill`).
    fn touch(&mut self, pgno: u64) -> Result<(), LmdbError> {
        let i = (pgno - 2) as usize;
        if !self.spills || matches!(self.slots[i], Slot::Dirty(_)) {
            return Ok(());
        }
        if self.dirty_room == 0 {
            return Err(LmdbError::Corrupted(
                "transaction full (MDB_TXN_FULL)".into(),
            ));
        }
        let mut np = self.malloc();
        let Slot::Spilled(src) = &self.slots[i] else {
            unreachable!()
        };
        page_copy(&mut np, src, self.psize);
        let f = pg_flags(&np, 0);
        set_flags(&mut np, 0, f | P_DIRTY);
        self.slots[i] = Slot::Dirty(np);
        self.dirty_room -= 1;
        Ok(())
    }

    /// `mdb_page_spill` before a put on a main cursor: when the dirty list is nearly full,
    /// write out the newest dirty pages other than the cursor's and the database roots.
    fn spill_check(
        &mut self,
        dbi: usize,
        db: &Db,
        leafsize: usize,
        keep: &[u64],
    ) -> Result<(), LmdbError> {
        let mut i = db.depth as u64;
        if dbi >= 2 {
            i += self.dbs[MAIN_DBI].depth as u64;
        }
        i += ((leafsize + self.psize) / self.psize) as u64;
        i += i;
        if self.dirty_room > i {
            return Ok(());
        }
        let is_dirty = |s: &Self, pg: u64| matches!(s.slots[(pg - 2) as usize], Slot::Dirty(_));
        let mut kept: std::collections::HashSet<u64> = keep
            .iter()
            .copied()
            .filter(|&pg| is_dirty(self, pg))
            .collect();
        for d in self.dbs.iter() {
            if d.root != P_INVALID && is_dirty(self, d.root) {
                kept.insert(d.root);
            }
        }
        let dl: Vec<u64> = (0..self.slots.len())
            .filter(|&x| matches!(self.slots[x], Slot::Dirty(_)))
            .map(|x| x as u64 + 2)
            .collect();
        let mut need = i.max(MDB_IDL_UM_MAX / 8);
        let mut idx = dl.len();
        while idx > 0 && need > 0 {
            if !kept.contains(&dl[idx - 1]) {
                need -= 1;
            }
            idx -= 1;
        }
        for &pg in &dl[idx..] {
            if kept.contains(&pg) {
                continue;
            }
            let s = (pg - 2) as usize;
            let Slot::Dirty(mut p) =
                std::mem::replace(&mut self.slots[s], Slot::Spilled(Box::new([])))
            else {
                unreachable!()
            };
            let f = pg_flags(&p, 0);
            set_flags(&mut p, 0, f & !P_DIRTY);
            self.file.seek(SeekFrom::Start(pg * self.psize as u64))?;
            self.file.write_all(&p)?;
            self.slots[s] = Slot::Spilled(p.clone());
            self.dpages.push(p);
            self.dirty_room += 1;
        }
        self.spills = true;
        Ok(())
    }

    /// `mdb_page_search` + `mdb_node_search`: the cursor stack to `key` and whether it exists.
    fn search(&self, db: &Db, cmp: Cmp, key: &[u8]) -> (Stack, bool) {
        let mut path = Vec::with_capacity(db.depth as usize);
        let mut pg = db.root;
        loop {
            let p = self.pg(pg);
            if pg_flags(p, 0) & P_BRANCH != 0 {
                let n = numkeys(p, 0);
                let (ki, exact, found) = node_search(p, 0, key, cmp, 0);
                let i = if !found {
                    n - 1
                } else if !exact {
                    ki - 1
                } else {
                    ki
                };
                path.push((pg, i));
                pg = n_pgno(p, ptr(p, 0, i));
            } else {
                let (ki, exact, _) = node_search(p, 0, key, cmp, db.pad as usize);
                path.push((pg, ki));
                return (path, exact);
            }
        }
    }

    /// The pages from a tree's root to its first leaf (`mdb_page_search(MDB_PS_FIRST)`).
    fn first_path(&self, db: &Db) -> Vec<u64> {
        let mut out = Vec::new();
        let mut pg = db.root;
        while pg != P_INVALID {
            out.push(pg);
            let p = self.pg(pg);
            if pg_flags(p, 0) & P_BRANCH == 0 {
                break;
            }
            pg = n_pgno(p, ptr(p, 0, 0));
        }
        out
    }

    /// The page and node offset of an existing key.
    fn locate(&self, db: &Db, cmp: Cmp, key: &[u8]) -> Result<(u64, usize), LmdbError> {
        let (path, exact) = self.search(db, cmp, key);
        let &(pg, ki) = path.last().unwrap();
        if !exact {
            return Err(LmdbError::Corrupted("inserted key not found".into()));
        }
        Ok((pg, ptr(self.pg(pg), 0, ki)))
    }

    /// The `new_sub` step of `_mdb_cursor_put`: add the node at the cursor, splitting the
    /// page when it has no room.
    fn insert(
        &mut self,
        db: &mut Db,
        path: &mut Stack,
        key: &[u8],
        data: &[u8],
        nflags: u16,
    ) -> Result<(), LmdbError> {
        let (pg, ki) = *path.last().unwrap();
        let p = self.pg(pg);
        let nsize = if pg_flags(p, 0) & P_LEAF2 != 0 {
            key.len()
        } else {
            leaf_size(self.nodemax, key.len(), data.len())
        };
        if sizeleft(p, 0) < nsize {
            self.page_split(db, path, key, Val::Data(data), nflags)
        } else {
            let nodemax = self.nodemax;
            node_add(
                self.pg_mut(pg),
                nodemax,
                db.pad as usize,
                ki,
                key,
                Val::Data(data),
                nflags,
            )
        }
    }

    /// `_mdb_cursor_put` on a main cursor (a fresh cursor per call, as `mdb_put` uses).
    fn cursor_put(
        &mut self,
        dbi: usize,
        key: &[u8],
        data: &[u8],
        flags: u16,
    ) -> Result<(), LmdbError> {
        let mut db = self.dbs[dbi];
        if key.len().wrapping_sub(1) >= MDB_MAXKEYSIZE
            || (db.flags & MDB_DUPSORT != 0 && data.len() > MDB_MAXKEYSIZE)
        {
            return Err(LmdbError::BadValSize(key.len(), data.len()));
        }
        let kcmp = Cmp::key(db.flags);
        let found = (db.root != P_INVALID).then(|| self.search(&db, kcmp, key));
        // Check for room in the dirty list, keeping the cursor's pages (and those of its
        // sub-database cursor, which `mdb_cursor_set` put on the first duplicate).
        let mut keep: Vec<u64> = Vec::new();
        if let Some((path, exact)) = &found {
            keep.extend(path.iter().map(|x| x.0));
            if *exact && db.flags & MDB_DUPSORT != 0 {
                let &(lp, ki) = path.last().unwrap();
                let p = self.pg(lp);
                let o = ptr(p, 0, ki);
                if n_flags(p, o) & (F_DUPDATA | F_SUBDATA) == F_DUPDATA | F_SUBDATA {
                    let d = o + NODESIZE + n_ksize(p, o);
                    if let Some(sub) = Db::from_bytes(&p[d..d + DB_SIZE]) {
                        keep.extend(self.first_path(&sub));
                    }
                }
            }
        }
        self.spill_check(dbi, &db, NODESIZE + key.len() + data.len(), &keep)?;

        let Some((mut path, exact)) = found else {
            // A new database: its root leaf.
            let pg = self.page_new(&mut db, P_LEAF)?;
            db.root = pg;
            db.depth += 1;
            let mut path = vec![(pg, 0)];
            self.insert(&mut db, &mut path, key, data, flags & NODE_ADD_FLAGS)?;
            db.entries += 1;
            self.dbs[dbi] = db;
            return Ok(());
        };
        for &(pg, _) in &path {
            self.touch(pg)?;
        }
        if !exact {
            if db.flags & MDB_DUPSORT != 0 && NODESIZE + key.len() + data.len() > self.nodemax {
                return Err(LmdbError::Unsupported(
                    "a duplicate too large for a node".into(),
                ));
            }
            self.insert(&mut db, &mut path, key, data, flags & NODE_ADD_FLAGS)?;
            db.entries += 1;
            self.dbs[dbi] = db;
            return Ok(());
        }

        let (lp, ki) = *path.last().unwrap();
        let (o, lflags, lks, ldsize) = {
            let p = self.pg(lp);
            let o = ptr(p, 0, ki);
            (o, n_flags(p, o), n_ksize(p, o), n_dsize(p, o))
        };
        let dat = o + NODESIZE + lks;
        if db.flags & MDB_DUPSORT == 0 {
            // Replace the value ("current:").
            if (lflags ^ flags) & F_SUBDATA != 0 {
                return Err(LmdbError::Corrupted("MDB_INCOMPATIBLE".into()));
            }
            if lflags & F_BIGDATA != 0 {
                return Err(LmdbError::Unsupported("overflow pages".into()));
            }
            if ldsize == data.len() {
                self.pg_mut(lp)[dat..dat + ldsize].copy_from_slice(data);
            } else {
                node_del(self.pg_mut(lp), ki, 0);
                self.insert(&mut db, &mut path, key, data, flags & NODE_ADD_FLAGS)?;
            }
            self.dbs[dbi] = db;
            return Ok(());
        }

        // A duplicate of an existing key.
        if db.flags & MDB_DUPFIXED == 0 {
            return Err(LmdbError::Unsupported(
                "MDB_DUPSORT without MDB_DUPFIXED".into(),
            ));
        }
        let dcmp = Cmp::dup(db.flags);
        w64(&mut self.pbuf, 0, lp);
        if lflags & F_DUPDATA == 0 {
            // One value so far: equal means nothing to do, else make a sub-page of both.
            let old = self.pg(lp)[dat..dat + ldsize].to_vec();
            if dcmp.cmp(data, &old) == Ordering::Equal {
                self.pg_mut(lp)[dat..dat + ldsize].copy_from_slice(data);
                return Ok(());
            }
            if old.len() != data.len() {
                return Err(LmdbError::BadValSize(key.len(), data.len()));
            }
            let pb = &mut self.pbuf;
            pb[PAGEHDRSZ..PAGEHDRSZ + old.len()].copy_from_slice(&old);
            set_flags(pb, 0, P_LEAF | P_DIRTY | P_SUBP | P_LEAF2);
            set_lower(pb, 0, PAGEHDRSZ);
            set_pad(pb, 0, data.len() as u16);
            let xsize = PAGEHDRSZ + old.len() + data.len() + 2 * data.len();
            set_upper(pb, 0, xsize);
            if NODESIZE + lks + xsize > self.nodemax {
                return Err(LmdbError::Unsupported(
                    "a duplicate too large for a sub-page".into(),
                ));
            }
            let xdata = self.pbuf[..xsize].to_vec();
            node_del(self.pg_mut(lp), ki, 0);
            self.insert(&mut db, &mut path, key, &xdata, F_DUPDATA)?;
            let (pg, o) = self.locate(&db, kcmp, key)?;
            let b = o + NODESIZE + lks;
            let p = self.pg_mut(pg);
            subpage_put(p, b, &old, dcmp, data.len())?;
            subpage_put(p, b, data, dcmp, data.len())?;
            db.entries += 1;
            self.dbs[dbi] = db;
            return Ok(());
        }
        if lflags & F_SUBDATA != 0 {
            // A sub-database.
            let mut sub = Db::from_bytes(&self.pg(lp)[dat..dat + ldsize.min(DB_SIZE)])
                .ok_or_else(|| LmdbError::Corrupted("short sub-database record".into()))?;
            let inserted = self.subdb_put(&mut sub, dcmp, data)?;
            self.pg_mut(lp)[dat..dat + DB_SIZE].copy_from_slice(&sub.to_bytes());
            if inserted {
                db.entries += 1;
            }
            self.dbs[dbi] = db;
            return Ok(());
        }
        // A sub-page: add in place if it has room, else grow it by four items or, past
        // `me_nodemax`, move it to a sub-database.
        let (fpad, flower, fupper, fflags) = {
            let p = self.pg(lp);
            (
                pg_pad(p, dat) as usize,
                pg_lower(p, dat),
                pg_upper(p, dat),
                pg_flags(p, dat),
            )
        };
        if data.len() != fpad {
            return Err(LmdbError::BadValSize(key.len(), data.len()));
        }
        if sizeleft(self.pg(lp), dat) >= fpad {
            let p = self.pg_mut(lp);
            set_flags(p, dat, fflags | P_DIRTY);
            w64(p, dat, lp);
            if subpage_put(p, dat, data, dcmp, fpad)? {
                db.entries += 1;
            }
            self.dbs[dbi] = db;
            return Ok(());
        }
        let offset = fpad * 4;
        let xsize = ldsize + offset;
        let nkeys = numkeys(self.pg(lp), dat);
        let keys = self.pg(lp)[dat + PAGEHDRSZ..dat + PAGEHDRSZ + nkeys * fpad].to_vec();
        if NODESIZE + lks + xsize > self.nodemax {
            // prep_subDB
            let mut sub = Db {
                pad: fpad as u32,
                flags: MDB_DUPFIXED
                    | if db.flags & MDB_INTEGERDUP != 0 {
                        MDB_INTEGERKEY
                    } else {
                        0
                    },
                depth: 1,
                branch_pages: 0,
                leaf_pages: 1,
                overflow_pages: 0,
                entries: nkeys as u64,
                root: 0,
            };
            let mp = self.alloc()?;
            let off2 = self.psize - ldsize;
            sub.root = mp;
            let q = self.pg_mut(mp);
            set_flags(q, 0, (fflags & !P_SUBP) | P_DIRTY);
            set_pad(q, 0, fpad as u16);
            set_lower(q, 0, flower);
            set_upper(q, 0, fupper + off2);
            q[PAGEHDRSZ..PAGEHDRSZ + keys.len()].copy_from_slice(&keys);
            node_del(self.pg_mut(lp), ki, 0);
            self.insert(
                &mut db,
                &mut path,
                key,
                &sub.to_bytes(),
                F_DUPDATA | F_SUBDATA,
            )?;
            let (pg, o) = self.locate(&db, kcmp, key)?;
            let d = o + NODESIZE + lks;
            let inserted = self.subdb_put(&mut sub, dcmp, data)?;
            self.pg_mut(pg)[d..d + DB_SIZE].copy_from_slice(&sub.to_bytes());
            if inserted {
                db.entries += 1;
            }
            self.dbs[dbi] = db;
            return Ok(());
        }
        let pb = &mut self.pbuf;
        set_flags(pb, 0, fflags | P_DIRTY);
        set_pad(pb, 0, fpad as u16);
        set_lower(pb, 0, flower);
        set_upper(pb, 0, fupper + offset);
        pb[PAGEHDRSZ..PAGEHDRSZ + keys.len()].copy_from_slice(&keys);
        let xdata = self.pbuf[..xsize].to_vec();
        node_del(self.pg_mut(lp), ki, 0);
        self.insert(&mut db, &mut path, key, &xdata, F_DUPDATA)?;
        let (pg, o) = self.locate(&db, kcmp, key)?;
        if subpage_put(self.pg_mut(pg), o + NODESIZE + lks, data, dcmp, fpad)? {
            db.entries += 1;
        }
        self.dbs[dbi] = db;
        Ok(())
    }

    /// `_mdb_cursor_put` on a sub-database cursor: `key` is the duplicate, stored as a
    /// `P_LEAF2` key. Returns false when it was already there.
    fn subdb_put(&mut self, sub: &mut Db, cmp: Cmp, key: &[u8]) -> Result<bool, LmdbError> {
        if key.len().wrapping_sub(1) >= MDB_MAXKEYSIZE || key.len() != sub.pad as usize {
            return Err(LmdbError::BadValSize(key.len(), 0));
        }
        let (mut path, exact) = self.search(sub, cmp, key);
        for &(pg, _) in &path {
            self.touch(pg)?;
        }
        if exact {
            return Ok(false);
        }
        self.insert(sub, &mut path, key, &[], 0)?;
        sub.entries += 1;
        Ok(true)
    }

    /// `mdb_page_split`: split the page at the top of `path` to insert `newkey` at the
    /// cursor index, adding the separator to the parent (splitting it in turn if full).
    fn page_split(
        &mut self,
        db: &mut Db,
        path: &mut Stack,
        newkey: &[u8],
        newval: Val,
        nflags: u16,
    ) -> Result<(), LmdbError> {
        let top = path.len() - 1;
        let (mp, newindx) = path[top];
        let (mp_flags, mp_pad, nkeys) = {
            let p = self.pg(mp);
            (pg_flags(p, 0), pg_pad(p, 0), numkeys(p, 0))
        };
        let rp = self.page_new(db, mp_flags)?;
        set_pad(self.pg_mut(rp), 0, mp_pad);
        let ptop = if top < 1 {
            let pp = self.page_new(db, P_BRANCH)?;
            path.insert(0, (pp, 0));
            db.root = pp;
            db.depth += 1;
            let nodemax = self.nodemax;
            node_add(self.pg_mut(pp), nodemax, 0, 0, &[], Val::Pgno(mp), 0)?;
            0
        } else {
            top - 1
        };
        let mut copy = None;
        let r = self.split_body(
            db,
            path,
            ptop,
            (mp, rp),
            mp_flags,
            nkeys,
            newindx,
            newkey,
            newval,
            nflags,
            &mut copy,
        );
        if let Some(c) = copy {
            self.dpages.push(c);
        }
        r
    }

    #[allow(clippy::too_many_arguments)]
    fn split_body(
        &mut self,
        db: &mut Db,
        path: &Stack,
        ptop: usize,
        (mp, rp): (u64, u64),
        mp_flags: u16,
        nkeys: usize,
        newindx: usize,
        newkey: &[u8],
        newval: Val,
        nflags: u16,
        copy: &mut Option<Box<[u8]>>,
    ) -> Result<(), LmdbError> {
        let psize = self.psize;
        let nodemax = self.nodemax;
        let is_leaf = mp_flags & P_LEAF != 0;
        let mut split_indx = nkeys.div_ceil(2);
        let sepkey: Vec<u8>;
        if mp_flags & P_LEAF2 != 0 {
            let ks = db.pad as usize;
            let rsize = (nkeys - split_indx) * ks;
            let lsize = (nkeys - split_indx) * 2;
            let split = PAGEHDRSZ + split_indx * ks;
            let old = self.pg(mp).to_vec();
            {
                let p = self.pg_mut(mp);
                let (l, u) = (pg_lower(p, 0), pg_upper(p, 0));
                set_lower(p, 0, l.wrapping_sub(lsize));
                set_upper(p, 0, u.wrapping_add(rsize).wrapping_sub(lsize));
            }
            {
                let q = self.pg_mut(rp);
                let (l, u) = (pg_lower(q, 0), pg_upper(q, 0));
                set_lower(q, 0, l + lsize);
                set_upper(q, 0, u.wrapping_sub(rsize).wrapping_add(lsize));
            }
            let mut sk = if newindx == split_indx {
                newkey.to_vec()
            } else {
                old[split..split + ks].to_vec()
            };
            if newindx < split_indx {
                let ins = PAGEHDRSZ + newindx * ks;
                let q = self.pg_mut(rp);
                q[PAGEHDRSZ..PAGEHDRSZ + rsize].copy_from_slice(&old[split..split + rsize]);
                sk = q[PAGEHDRSZ..PAGEHDRSZ + ks].to_vec();
                let p = self.pg_mut(mp);
                p.copy_within(ins..ins + (split_indx - newindx) * ks, ins + ks);
                p[ins..ins + ks].copy_from_slice(&newkey[..ks]);
                let (l, u) = (pg_lower(p, 0), pg_upper(p, 0));
                set_lower(p, 0, l + 2);
                set_upper(p, 0, u.wrapping_sub(ks).wrapping_add(2));
            } else {
                let x = newindx - split_indx;
                let q = self.pg_mut(rp);
                if x > 0 {
                    q[PAGEHDRSZ..PAGEHDRSZ + x * ks].copy_from_slice(&old[split..split + x * ks]);
                }
                let ins = PAGEHDRSZ + x * ks;
                q[ins..ins + ks].copy_from_slice(&newkey[..ks]);
                q[ins + ks..ins + ks + rsize - x * ks]
                    .copy_from_slice(&old[split + x * ks..split + rsize]);
                let (l, u) = (pg_lower(q, 0), pg_upper(q, 0));
                set_lower(q, 0, l + 2);
                set_upper(q, 0, u.wrapping_sub(ks).wrapping_add(2));
            }
            sepkey = sk;
        } else {
            let pmax = psize - PAGEHDRSZ;
            let keythresh = psize >> 7;
            let nsize = even(match (is_leaf, newval) {
                (true, Val::Data(d)) => leaf_size(nodemax, newkey.len(), d.len()),
                _ => branch_size(newkey.len()),
            });
            let mut c = self.malloc();
            w64(&mut c, 0, mp);
            set_flags(&mut c, 0, mp_flags);
            set_lower(&mut c, 0, PAGEHDRSZ);
            set_upper(&mut c, 0, psize);
            let p = self.pg(mp);
            let mut j = 0;
            for i in 0..nkeys {
                if i == newindx {
                    set_ptr(&mut c, 0, j, 0);
                    j += 1;
                }
                set_ptr(&mut c, 0, j, ptr(p, 0, i));
                j += 1;
            }
            if nkeys < keythresh || nsize > pmax / 16 || newindx >= nkeys {
                let mut psz = 0usize;
                let (mut i, j, k): (isize, isize, isize) =
                    if newindx <= split_indx || newindx >= nkeys {
                        let k = if newindx >= nkeys {
                            nkeys
                        } else {
                            split_indx + 1 + is_leaf as usize
                        };
                        (0, 1, k as isize)
                    } else {
                        (nkeys as isize, -1, split_indx as isize - 1)
                    };
                while i != k {
                    if i as usize == newindx {
                        psz += nsize;
                    } else {
                        let o = ptr(&c, 0, i as usize);
                        psz += NODESIZE + n_ksize(p, o) + 2;
                        if is_leaf {
                            psz += if n_flags(p, o) & F_BIGDATA != 0 {
                                8
                            } else {
                                n_dsize(p, o)
                            };
                        }
                        psz = even(psz);
                    }
                    if psz > pmax || i == k - j {
                        split_indx = (i + (j < 0) as isize) as usize;
                        break;
                    }
                    i += j;
                }
            }
            sepkey = if split_indx == newindx {
                newkey.to_vec()
            } else {
                n_key(p, ptr(&c, 0, split_indx)).to_vec()
            };
            *copy = Some(c);
        }

        // Add the separator to the parent, splitting the parent if it is full.
        let (ppg, pki) = path[ptop];
        if sizeleft(self.pg(ppg), 0) < branch_size(sepkey.len()) {
            let mut mn: Stack = path[..=ptop].to_vec();
            mn[ptop].1 = pki + 1;
            self.page_split(db, &mut mn, &sepkey, Val::Pgno(rp), 0)?;
        } else {
            node_add(
                self.pg_mut(ppg),
                nodemax,
                0,
                pki + 1,
                &sepkey,
                Val::Pgno(rp),
                0,
            )?;
        }
        if mp_flags & P_LEAF2 != 0 {
            return Ok(());
        }

        // Move the nodes from split_indx on to the right page, the rest through the copy.
        let c = copy.as_mut().unwrap();
        let mut i = split_indx;
        let mut j = 0;
        let mut to_copy = false;
        loop {
            let (k, d, pg, fl) = if i == newindx {
                match newval {
                    Val::Data(d) => (newkey.to_vec(), d.to_vec(), 0, nflags),
                    Val::Pgno(pg) => (newkey.to_vec(), Vec::new(), pg, nflags),
                }
            } else {
                let p = self.pg(mp);
                let o = ptr(c, 0, i);
                let fl = n_flags(p, o);
                if is_leaf && fl & F_BIGDATA != 0 {
                    return Err(LmdbError::Unsupported("overflow pages".into()));
                }
                let ks = n_ksize(p, o);
                let d = if is_leaf {
                    p[o + NODESIZE + ks..o + NODESIZE + ks + n_dsize(p, o)].to_vec()
                } else {
                    Vec::new()
                };
                (n_key(p, o).to_vec(), d, n_pgno(p, o), fl)
            };
            let key: &[u8] = if !is_leaf && j == 0 { &[] } else { &k };
            let val = if is_leaf {
                Val::Data(&d)
            } else {
                Val::Pgno(pg)
            };
            let dst: &mut [u8] = if to_copy { &mut c[..] } else { self.pg_mut(rp) };
            node_add(dst, nodemax, 0, j, key, val, fl)?;
            if i == nkeys {
                i = 0;
                j = 0;
                to_copy = true;
            } else {
                i += 1;
                j += 1;
            }
            if i == split_indx {
                break;
            }
        }
        // Back to the original page: pointers, bounds and nodes; its gap keeps the old bytes.
        let c: &[u8] = c;
        let (nk, up) = (numkeys(c, 0), pg_upper(c, 0));
        let from = ptr(c, 0, nk - 1);
        let p = self.pg_mut(mp);
        for x in 0..nk {
            set_ptr(p, 0, x, ptr(c, 0, x));
        }
        set_lower(p, 0, pg_lower(c, 0));
        set_upper(p, 0, up);
        p[from..from + (psize - up)].copy_from_slice(&c[from..from + (psize - up)]);
        Ok(())
    }
}

//! The LMDB file format as `mdb.c` lays it out: page and node headers, the `MDB_db` record and
//! the key comparators. Shared by the writer ([`super::write`]) and the reader
//! ([`super::read`]). See [`super`] for the copyright and license of the code this follows.
//!
//! LMDB writes its structures in native byte order with the field layout of a 64-bit build
//! (`pgno_t`, `txnid_t` and `mdb_size_t` are 8 bytes); this module does the same.

use std::cmp::Ordering;

/// `PAGEHDRSZ`: pgno (8), pad (2), flags (2), lower (2), upper (2).
pub const PAGEHDRSZ: usize = 16;
/// `NODESIZE`: lo (2), hi (2), flags (2), ksize (2).
pub const NODESIZE: usize = 8;
/// `P_INVALID`: no page (an empty tree).
pub const P_INVALID: u64 = u64::MAX;

pub const P_BRANCH: u16 = 0x01;
pub const P_LEAF: u16 = 0x02;
pub const P_OVERFLOW: u16 = 0x04;
pub const P_META: u16 = 0x08;
pub const P_DIRTY: u16 = 0x10;
pub const P_LEAF2: u16 = 0x20;
pub const P_SUBP: u16 = 0x40;

/// Node flags: data on an overflow page, a sub-database record, duplicate data.
pub const F_BIGDATA: u16 = 0x01;
pub const F_SUBDATA: u16 = 0x02;
pub const F_DUPDATA: u16 = 0x04;

/// Database flags (`mdb_dbi_open`).
pub const MDB_DUPSORT: u16 = 0x04;
pub const MDB_INTEGERKEY: u16 = 0x08;
pub const MDB_DUPFIXED: u16 = 0x10;
pub const MDB_INTEGERDUP: u16 = 0x20;
/// The flags this port compares keys and duplicates for; others (`MDB_REVERSEKEY`,
/// `MDB_REVERSEDUP`) are refused.
pub const SUPPORTED_DB_FLAGS: u16 = MDB_DUPSORT | MDB_INTEGERKEY | MDB_DUPFIXED | MDB_INTEGERDUP;

pub const MDB_MAGIC: u32 = 0xBEEF_C0DE;
pub const MDB_DATA_VERSION: u32 = 1;
/// `MAX_PAGESIZE` (with `PAGEBASE` 0).
pub const MAX_PAGESIZE: usize = 0x8000;
/// Size of an `MDB_db` record.
pub const DB_SIZE: usize = 48;
/// Size of `MDB_meta`: magic, version, address, mapsize, two `MDB_db`, last_pg, txnid.
pub const META_SIZE: usize = 24 + 2 * DB_SIZE + 16;
/// `FREE_DBI` and `MAIN_DBI`.
pub const FREE_DBI: usize = 0;
pub const MAIN_DBI: usize = 1;

#[inline]
pub fn r16(p: &[u8], o: usize) -> u16 {
    u16::from_ne_bytes([p[o], p[o + 1]])
}
#[inline]
pub fn r32(p: &[u8], o: usize) -> u32 {
    u32::from_ne_bytes(p[o..o + 4].try_into().unwrap())
}
#[inline]
pub fn r64(p: &[u8], o: usize) -> u64 {
    u64::from_ne_bytes(p[o..o + 8].try_into().unwrap())
}
#[inline]
pub fn w16(p: &mut [u8], o: usize, v: u16) {
    p[o..o + 2].copy_from_slice(&v.to_ne_bytes());
}
#[inline]
pub fn w32(p: &mut [u8], o: usize, v: u32) {
    p[o..o + 4].copy_from_slice(&v.to_ne_bytes());
}
#[inline]
pub fn w64(p: &mut [u8], o: usize, v: u64) {
    p[o..o + 8].copy_from_slice(&v.to_ne_bytes());
}

/// `EVEN(n)`.
#[inline]
pub fn even(n: usize) -> usize {
    (n + 1) & !1
}

// Page header fields of the page (or sub-page) starting at `b`.
#[inline]
pub fn pg_pad(p: &[u8], b: usize) -> u16 {
    r16(p, b + 8)
}
#[inline]
pub fn pg_flags(p: &[u8], b: usize) -> u16 {
    r16(p, b + 10)
}
#[inline]
pub fn pg_lower(p: &[u8], b: usize) -> usize {
    r16(p, b + 12) as usize
}
#[inline]
pub fn pg_upper(p: &[u8], b: usize) -> usize {
    r16(p, b + 14) as usize
}
#[inline]
pub fn set_pad(p: &mut [u8], b: usize, v: u16) {
    w16(p, b + 8, v)
}
#[inline]
pub fn set_flags(p: &mut [u8], b: usize, v: u16) {
    w16(p, b + 10, v)
}
/// `mp_lower`/`mp_upper` are `indx_t` (u16): stores wrap like the C fields.
#[inline]
pub fn set_lower(p: &mut [u8], b: usize, v: usize) {
    w16(p, b + 12, v as u16)
}
#[inline]
pub fn set_upper(p: &mut [u8], b: usize, v: usize) {
    w16(p, b + 14, v as u16)
}
/// `NUMKEYS`.
#[inline]
pub fn numkeys(p: &[u8], b: usize) -> usize {
    pg_lower(p, b).saturating_sub(PAGEHDRSZ) >> 1
}
/// `SIZELEFT`.
#[inline]
pub fn sizeleft(p: &[u8], b: usize) -> usize {
    (pg_upper(p, b) as u16).wrapping_sub(pg_lower(p, b) as u16) as usize
}
/// `mp_ptrs[i]`: the offset of node `i` from the start of its page.
#[inline]
pub fn ptr(p: &[u8], b: usize, i: usize) -> usize {
    r16(p, b + PAGEHDRSZ + 2 * i) as usize
}
#[inline]
pub fn set_ptr(p: &mut [u8], b: usize, i: usize, v: usize) {
    w16(p, b + PAGEHDRSZ + 2 * i, v as u16)
}

// Node fields of the node at absolute offset `o`.
/// `NODEDSZ` (mn_lo | mn_hi << 16; the two halves are in native order, so this is a native u32).
#[inline]
pub fn n_dsize(p: &[u8], o: usize) -> usize {
    r32(p, o) as usize
}
#[inline]
pub fn n_flags(p: &[u8], o: usize) -> u16 {
    r16(p, o + 4)
}
#[inline]
pub fn n_ksize(p: &[u8], o: usize) -> usize {
    r16(p, o + 6) as usize
}
/// `NODEPGNO` of a branch node (mn_lo | mn_hi << 16 | mn_flags << 32).
#[inline]
pub fn n_pgno(p: &[u8], o: usize) -> u64 {
    r32(p, o) as u64 | (n_flags(p, o) as u64) << 32
}
#[inline]
pub fn n_key(p: &[u8], o: usize) -> &[u8] {
    &p[o + NODESIZE..o + NODESIZE + n_ksize(p, o)]
}

/// An `MDB_db` record: a B+tree's root and counts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Db {
    pub pad: u32,
    pub flags: u16,
    pub depth: u16,
    pub branch_pages: u64,
    pub leaf_pages: u64,
    pub overflow_pages: u64,
    pub entries: u64,
    pub root: u64,
}

impl Db {
    /// An empty tree (`md_root` = `P_INVALID`), as `mdb_dbi_open` creates it.
    pub fn empty(flags: u16) -> Self {
        Db {
            pad: 0,
            flags,
            depth: 0,
            branch_pages: 0,
            leaf_pages: 0,
            overflow_pages: 0,
            entries: 0,
            root: P_INVALID,
        }
    }

    pub fn to_bytes(self) -> [u8; DB_SIZE] {
        let mut b = [0u8; DB_SIZE];
        w32(&mut b, 0, self.pad);
        w16(&mut b, 4, self.flags);
        w16(&mut b, 6, self.depth);
        w64(&mut b, 8, self.branch_pages);
        w64(&mut b, 16, self.leaf_pages);
        w64(&mut b, 24, self.overflow_pages);
        w64(&mut b, 32, self.entries);
        w64(&mut b, 40, self.root);
        b
    }

    pub fn from_bytes(b: &[u8]) -> Option<Self> {
        (b.len() >= DB_SIZE).then(|| Db {
            pad: r32(b, 0),
            flags: r16(b, 4),
            depth: r16(b, 6),
            branch_pages: r64(b, 8),
            leaf_pages: r64(b, 16),
            overflow_pages: r64(b, 24),
            entries: r64(b, 32),
            root: r64(b, 40),
        })
    }
}

/// A key comparator: `mdb_cmp_memn` (bytes, then length) or the unsigned-integer comparators
/// (`mdb_cmp_int`, `mdb_cmp_cint`, `mdb_cmp_long`) for keys of one native-endian integer size.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Cmp {
    Memn,
    Int,
}

impl Cmp {
    /// `md_cmp` of a database (`mdb_default_cmp`).
    pub fn key(flags: u16) -> Self {
        if flags & MDB_INTEGERKEY != 0 {
            Cmp::Int
        } else {
            Cmp::Memn
        }
    }

    /// `md_dcmp` of a `MDB_DUPSORT` database.
    pub fn dup(flags: u16) -> Self {
        if flags & MDB_INTEGERDUP != 0 {
            Cmp::Int
        } else {
            Cmp::Memn
        }
    }

    #[inline]
    pub fn cmp(self, a: &[u8], b: &[u8]) -> Ordering {
        match self {
            Cmp::Memn => a.cmp(b),
            // Same-size native-endian unsigned integers: most significant byte first.
            Cmp::Int if cfg!(target_endian = "little") => a.iter().rev().cmp(b.iter().rev()),
            Cmp::Int => a.cmp(b),
        }
    }
}

//! Columnar store catalog: every sequence, membership and collection fact the
//! store keeps in RAM, each stored once.
//!
//! The old catalog was an object graph: an owned `SequenceRecord` per sequence
//! in a global map, a second owned copy per collection member, a third copy of
//! each name in a per-collection name map, digests stored as zero-padded text,
//! and every copy a separate heap allocation inside a hash table. This module
//! replaces it with plain columns and `u32` row numbers.
//!
//! Three row spaces:
//! - **sequence rows** ([`SeqRow`]): one per unique sequence (by sha512t24u).
//!   Raw digests, length, alphabet, and the store-wide label (first-seen name,
//!   which is what `sequences.rgsi` publishes).
//! - **member rows** ([`MemberRow`]): one per (collection, member). The member's
//!   sequence row and its per-collection name. A collection's members are one
//!   contiguous range of member rows, in FASTA order.
//! - **collection rows** ([`CollRow`]): one per collection, with its level-1
//!   metadata and whether its member list is loaded.
//!
//! Lookups are `hashbrown::HashTable`s of row numbers, hashed from the column
//! they index, so no key is stored twice. Descriptions and FAI records are rare
//! on the read path (FAI is never persisted), so they live in sparse maps.
//!
//! Deletes tombstone rows; [`Catalog::compact`] renumbers once dead rows pile up.
//! Row numbers are internal: callers outside the store go through digests.

use std::hash::BuildHasher;

use anyhow::{Result, anyhow};
use hashbrown::{HashMap, HashTable};

use crate::digest::raw::{Digest32, RawMd5, RawSha, md5_to_raw, raw_to_md5, raw_to_sha, sha_to_raw};
use crate::digest::{AlphabetType, FaiMetadata, SequenceCollectionMetadata, SequenceMetadata};

pub(crate) type SeqRow = u32;
pub(crate) type MemberRow = u32;
pub(crate) type CollRow = u32;

/// Fraction of dead rows above which [`Catalog::needs_compaction`] says yes.
const COMPACT_DEAD_FRACTION: f64 = 0.25;

/// A packed bit vector (one bit per row).
#[derive(Debug, Default, Clone)]
struct Bits {
    words: Vec<u64>,
    len: usize,
}

impl Bits {
    fn push(&mut self, v: bool) {
        if self.len % 64 == 0 {
            self.words.push(0);
        }
        self.len += 1;
        self.set(self.len - 1, v);
    }
    fn get(&self, i: usize) -> bool {
        (self.words[i / 64] >> (i % 64)) & 1 == 1
    }
    fn set(&mut self, i: usize, v: bool) {
        let w = &mut self.words[i / 64];
        if v {
            *w |= 1 << (i % 64);
        } else {
            *w &= !(1 << (i % 64));
        }
    }
    fn shrink_to_fit(&mut self) {
        self.words.shrink_to_fit();
    }
}

/// Drop the zero padding of a `DigestKey` passed as a digest.
fn trim_nul(d: &[u8]) -> &[u8] {
    let end = d.iter().rposition(|&b| b != 0).map_or(0, |i| i + 1);
    &d[..end]
}

fn hash_sha(raw: &RawSha) -> u64 {
    // sha512t24u bytes are uniformly distributed already.
    u64::from_le_bytes(raw[..8].try_into().unwrap())
}

fn hash_md5(raw: &RawMd5) -> u64 {
    u64::from_le_bytes(raw[..8].try_into().unwrap())
}

/// The fixed-size facts about one sequence that the byte-read paths need.
/// `Copy`, no heap: cheap to pull out of the catalog and hand to I/O code.
#[derive(Clone, Copy, Debug)]
pub(crate) struct SeqCore {
    pub(crate) sha512t24u: Digest32,
    pub(crate) length: usize,
    pub(crate) alphabet: AlphabetType,
}

/// Borrowed columns of one sequence, as handed to [`Catalog::add_sequence_parts`].
pub(crate) struct SeqParts<'a> {
    pub(crate) name: &'a str,
    pub(crate) description: Option<&'a str>,
    pub(crate) length: usize,
    pub(crate) alphabet: AlphabetType,
    pub(crate) sha512t24u: &'a str,
    pub(crate) md5: &'a str,
}

/// One member of a collection, as handed to [`Catalog::set_members`].
pub(crate) struct NewMember<'a> {
    pub(crate) seq: SeqRow,
    pub(crate) name: &'a str,
    pub(crate) description: Option<&'a str>,
    pub(crate) fai: Option<&'a FaiMetadata>,
}

#[derive(Debug, Default)]
pub(crate) struct Catalog {
    // --- one row per unique sequence ---
    sha: Vec<RawSha>,
    md5: Vec<RawMd5>,
    length: Vec<u64>,
    alphabet: Vec<AlphabetType>,
    seq_name_off: Vec<u32>,
    seq_name_len: Vec<u32>,
    seq_desc: hashbrown::HashMap<SeqRow, Box<str>>,
    seq_fai: hashbrown::HashMap<SeqRow, FaiMetadata>,
    seq_refcount: Vec<u32>,
    seq_dead: Bits,
    n_dead_seqs: usize,

    // --- one row per (collection, member) ---
    member_seq: Vec<SeqRow>,
    /// Owning collection of each member (needed to rehash `by_name`).
    member_coll: Vec<CollRow>,
    member_name_off: Vec<u32>,
    member_name_len: Vec<u32>,
    member_desc: hashbrown::HashMap<MemberRow, Box<str>>,
    member_fai: hashbrown::HashMap<MemberRow, FaiMetadata>,
    member_dead: Bits,
    n_dead_members: usize,

    /// One shared arena of all sequence and member names. A member whose name
    /// equals its sequence's label shares the label's span.
    names: String,

    // --- one row per collection ---
    coll_meta: Vec<SequenceCollectionMetadata>,
    coll_start: Vec<MemberRow>,
    coll_len: Vec<u32>,
    coll_loaded: Bits,
    coll_dead: Bits,
    n_dead_colls: usize,

    // --- lookups ---
    by_sha: HashTable<SeqRow>,
    by_md5: HashTable<SeqRow>,
    by_coll: HashTable<CollRow>,
    /// Key = (collection row, member name). Points at the LAST member with that
    /// name, matching the old `IndexMap::insert` overwrite semantics.
    by_name: HashTable<MemberRow>,
    hasher: hashbrown::DefaultHashBuilder,
}

impl Catalog {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    // =====================================================================
    // Sizes
    // =====================================================================

    /// Number of live sequences.
    pub(crate) fn n_sequences(&self) -> usize {
        self.sha.len() - self.n_dead_seqs
    }

    /// Number of live collections.
    pub(crate) fn n_collections(&self) -> usize {
        self.coll_meta.len() - self.n_dead_colls
    }

    // =====================================================================
    // Names arena
    // =====================================================================

    fn push_name(&mut self, name: &str) -> Result<(u32, u32)> {
        let off = u32::try_from(self.names.len())
            .map_err(|_| anyhow!("catalog name arena exceeds 4 GiB"))?;
        let len = u32::try_from(name.len()).map_err(|_| anyhow!("sequence name too long"))?;
        self.names.push_str(name);
        Ok((off, len))
    }

    fn name_at(&self, off: u32, len: u32) -> &str {
        &self.names[off as usize..(off + len) as usize]
    }

    // =====================================================================
    // Sequences
    // =====================================================================

    /// Find a live sequence by sha512t24u text (bare or `SQ.`-prefixed).
    pub(crate) fn find_sha(&self, digest: impl AsRef<[u8]>) -> Option<SeqRow> {
        let raw = sha_to_raw(trim_nul(digest.as_ref()))?;
        self.find_sha_raw(&raw)
    }

    pub(crate) fn find_sha_raw(&self, raw: &RawSha) -> Option<SeqRow> {
        self.by_sha
            .find(hash_sha(raw), |&r| &self.sha[r as usize] == raw)
            .copied()
    }

    /// Find a live sequence by MD5 hex text.
    pub(crate) fn find_md5(&self, digest: impl AsRef<[u8]>) -> Option<SeqRow> {
        let raw = md5_to_raw(trim_nul(digest.as_ref()))?;
        self.by_md5
            .find(hash_md5(&raw), |&r| self.md5[r as usize] == raw)
            .copied()
    }

    /// Resolve either digest kind: MD5 first (as the old `md5_lookup` fallback
    /// did), then sha512t24u.
    pub(crate) fn resolve(&self, digest: impl AsRef<[u8]>) -> Option<SeqRow> {
        let d = digest.as_ref();
        self.find_md5(d).or_else(|| self.find_sha(d))
    }

    /// Insert a sequence, deduplicating by sha512t24u. Returns its row and
    /// whether it was newly inserted. A dedup hit leaves the existing row
    /// untouched.
    pub(crate) fn add_sequence(&mut self, meta: &SequenceMetadata) -> Result<(SeqRow, bool)> {
        let (row, new) = self.add_sequence_parts(SeqParts {
            name: &meta.name,
            description: meta.description.as_deref(),
            length: meta.length,
            alphabet: meta.alphabet,
            sha512t24u: &meta.sha512t24u,
            md5: &meta.md5,
        })?;
        if new {
            if let Some(f) = &meta.fai {
                self.seq_fai.insert(row, f.clone());
            }
        }
        Ok((row, new))
    }

    /// [`add_sequence`](Self::add_sequence) from borrowed columns, so index
    /// parsers append without building a `SequenceMetadata` per line.
    pub(crate) fn add_sequence_parts(&mut self, p: SeqParts<'_>) -> Result<(SeqRow, bool)> {
        let sha = sha_to_raw(p.sha512t24u)
            .ok_or_else(|| anyhow!("invalid sha512t24u digest: {:?}", p.sha512t24u))?;
        if let Some(row) = self.find_sha_raw(&sha) {
            return Ok((row, false));
        }
        let md5 = md5_to_raw(p.md5).ok_or_else(|| anyhow!("invalid md5 digest: {:?}", p.md5))?;
        let row = u32::try_from(self.sha.len()).map_err(|_| anyhow!("too many sequences"))?;
        let (off, len) = self.push_name(p.name)?;
        self.sha.push(sha);
        self.md5.push(md5);
        self.length.push(p.length as u64);
        self.alphabet.push(p.alphabet);
        self.seq_name_off.push(off);
        self.seq_name_len.push(len);
        self.seq_refcount.push(0);
        self.seq_dead.push(false);
        if let Some(d) = p.description {
            self.seq_desc.insert(row, d.into());
        }
        self.index_seq(row);
        Ok((row, true))
    }

    fn index_seq(&mut self, row: SeqRow) {
        let Self { sha, md5, by_sha, by_md5, .. } = self;
        by_sha.insert_unique(hash_sha(&sha[row as usize]), row, |&r| hash_sha(&sha[r as usize]));
        // An MD5 shared by two sha rows cannot happen for real data; keep the
        // first, like a HashMap insert keyed by md5 would keep the last. Either
        // way lookups stay well defined.
        let m = md5[row as usize];
        if by_md5.find(hash_md5(&m), |&r| md5[r as usize] == m).is_none() {
            by_md5.insert_unique(hash_md5(&m), row, |&r| hash_md5(&md5[r as usize]));
        }
    }

    /// Replace a sequence's store-wide label (name, description, FAI). Used by
    /// forced re-adds, which overwrote the old global record.
    pub(crate) fn set_sequence_label(&mut self, row: SeqRow, meta: &SequenceMetadata) -> Result<()> {
        let (off, len) = self.push_name(&meta.name)?;
        self.seq_name_off[row as usize] = off;
        self.seq_name_len[row as usize] = len;
        match &meta.description {
            Some(d) => {
                self.seq_desc.insert(row, d.as_str().into());
            }
            None => {
                self.seq_desc.remove(&row);
            }
        }
        match &meta.fai {
            Some(f) => {
                self.seq_fai.insert(row, f.clone());
            }
            None => {
                self.seq_fai.remove(&row);
            }
        }
        Ok(())
    }

    /// Tombstone a sequence. The caller is responsible for not removing a
    /// sequence a live member still points at.
    pub(crate) fn remove_sequence(&mut self, row: SeqRow) {
        let i = row as usize;
        if self.seq_dead.get(i) {
            return;
        }
        let Self { sha, md5, by_sha, by_md5, .. } = self;
        if let Ok(e) = by_sha.find_entry(hash_sha(&sha[i]), |&r| r == row) {
            e.remove();
        }
        if let Ok(e) = by_md5.find_entry(hash_md5(&md5[i]), |&r| r == row) {
            e.remove();
        }
        self.seq_desc.remove(&row);
        self.seq_fai.remove(&row);
        self.seq_dead.set(i, true);
        self.n_dead_seqs += 1;
    }

    pub(crate) fn seq(&self, row: SeqRow) -> SeqView<'_> {
        SeqView { cat: self, row }
    }

    /// Iterate live sequence rows, in insertion order.
    pub(crate) fn seq_rows(&self) -> impl Iterator<Item = SeqRow> + '_ {
        (0..self.sha.len() as u32).filter(move |&r| !self.seq_dead.get(r as usize))
    }

    #[cfg(test)]
    pub(crate) fn seq_refcount(&self, row: SeqRow) -> u32 {
        self.seq_refcount[row as usize]
    }

    // =====================================================================
    // Collections
    // =====================================================================

    fn hash_coll(&self, digest: &[u8]) -> u64 {
        self.hasher.hash_one(digest)
    }

    pub(crate) fn find_coll(&self, digest: impl AsRef<[u8]>) -> Option<CollRow> {
        let d = trim_nul(digest.as_ref());
        self.by_coll
            .find(self.hash_coll(d), |&r| self.coll_meta[r as usize].digest.as_bytes() == d)
            .copied()
    }

    /// Insert a collection's metadata (member list not loaded), or replace the
    /// metadata of an existing one. Returns its row.
    pub(crate) fn upsert_collection(&mut self, meta: SequenceCollectionMetadata) -> Result<CollRow> {
        if let Some(row) = self.find_coll(meta.digest.as_bytes()) {
            self.coll_meta[row as usize] = meta;
            return Ok(row);
        }
        let row = u32::try_from(self.coll_meta.len()).map_err(|_| anyhow!("too many collections"))?;
        let h = self.hash_coll(meta.digest.as_bytes());
        self.coll_meta.push(meta);
        self.coll_start.push(0);
        self.coll_len.push(0);
        self.coll_loaded.push(false);
        self.coll_dead.push(false);
        let Self { coll_meta, by_coll, hasher, .. } = self;
        by_coll.insert_unique(h, row, |&r| hasher.hash_one(coll_meta[r as usize].digest.as_bytes()));
        Ok(row)
    }

    /// Load a collection's members, in order, as one new contiguous range. Any
    /// previously loaded members of this collection are tombstoned first.
    pub(crate) fn set_members<'a>(
        &mut self,
        coll: CollRow,
        members: impl IntoIterator<Item = NewMember<'a>>,
    ) -> Result<()> {
        self.start_members(coll);
        for m in members {
            self.push_member(coll, m)?;
        }
        Ok(())
    }

    /// Begin (re)loading a collection's members: drop any loaded ones and
    /// start an empty range at the end of the member columns. Follow with
    /// [`push_member`](Self::push_member) calls for this collection only.
    pub(crate) fn start_members(&mut self, coll: CollRow) {
        if self.coll_loaded.get(coll as usize) {
            self.drop_members(coll);
        }
        self.coll_start[coll as usize] = self.member_seq.len() as u32;
        self.coll_len[coll as usize] = 0;
        self.coll_loaded.set(coll as usize, true);
    }

    /// Append one member to the collection opened by `start_members`.
    pub(crate) fn push_member(&mut self, coll: CollRow, m: NewMember<'_>) -> Result<()> {
        let row = u32::try_from(self.member_seq.len()).map_err(|_| anyhow!("too many members"))?;
        debug_assert_eq!(
            row,
            self.coll_start[coll as usize] + self.coll_len[coll as usize],
            "members of one collection must be pushed contiguously"
        );
        let s = m.seq as usize;
        let (off, len) = if self.name_at(self.seq_name_off[s], self.seq_name_len[s]) == m.name {
            (self.seq_name_off[s], self.seq_name_len[s])
        } else {
            self.push_name(m.name)?
        };
        self.member_seq.push(m.seq);
        self.member_coll.push(coll);
        self.member_name_off.push(off);
        self.member_name_len.push(len);
        self.member_dead.push(false);
        if let Some(d) = m.description {
            self.member_desc.insert(row, d.into());
        }
        if let Some(f) = m.fai {
            self.member_fai.insert(row, f.clone());
        }
        self.seq_refcount[s] += 1;
        self.coll_len[coll as usize] += 1;
        self.index_member(coll, row);
        Ok(())
    }

    fn hash_member(&self, coll: CollRow, name: &str) -> u64 {
        self.hasher.hash_one((coll, name))
    }

    fn member_name(&self, row: MemberRow) -> &str {
        self.name_at(self.member_name_off[row as usize], self.member_name_len[row as usize])
    }

    fn index_member(&mut self, coll: CollRow, row: MemberRow) {
        let name = self.member_name(row);
        let h = self.hash_member(coll, name);
        let (start, end) = self.coll_range(coll);
        let in_coll = |r: MemberRow| r >= start && r < end;
        let existing = self
            .by_name
            .find(h, |&r| in_coll(r) && self.member_name(r) == name)
            .copied();
        match existing {
            Some(old) => {
                // Duplicate name: point at the later member.
                let slot = self
                    .by_name
                    .find_mut(h, |&r| r == old)
                    .expect("entry just found");
                *slot = row;
            }
            None => {
                let Self { by_name, names, member_name_off, member_name_len, member_coll, hasher, .. } = self;
                by_name.insert_unique(h, row, |&r| {
                    let i = r as usize;
                    let name = &names[member_name_off[i] as usize
                        ..(member_name_off[i] + member_name_len[i]) as usize];
                    let c = member_coll[i];
                    hasher.hash_one((c, name))
                });
            }
        }
    }

    fn drop_members(&mut self, coll: CollRow) {
        let (start, end) = self.coll_range(coll);
        for row in start..end {
            let name = self.member_name(row);
            let h = self.hash_member(coll, name);
            if let Ok(e) = self.by_name.find_entry(h, |&r| r == row) {
                e.remove();
            }
            self.member_desc.remove(&row);
            self.member_fai.remove(&row);
            let s = self.member_seq[row as usize] as usize;
            self.seq_refcount[s] = self.seq_refcount[s].saturating_sub(1);
            self.member_dead.set(row as usize, true);
            self.n_dead_members += 1;
        }
        self.coll_len[coll as usize] = 0;
        self.coll_loaded.set(coll as usize, false);
    }

    /// Tombstone a collection and its members. Sequences are left alone: which
    /// ones are orphans is the store's (disk-aware) decision.
    pub(crate) fn remove_collection(&mut self, coll: CollRow) {
        if self.coll_dead.get(coll as usize) {
            return;
        }
        if self.coll_loaded.get(coll as usize) {
            self.drop_members(coll);
        }
        let h = self.hash_coll(self.coll_meta[coll as usize].digest.as_bytes());
        if let Ok(e) = self.by_coll.find_entry(h, |&r| r == coll) {
            e.remove();
        }
        self.coll_dead.set(coll as usize, true);
        self.n_dead_colls += 1;
    }

    /// Mark a collection's members unloaded (tombstoning them), e.g. after a
    /// failed load.
    pub(crate) fn unload_members(&mut self, coll: CollRow) {
        if self.coll_loaded.get(coll as usize) {
            self.drop_members(coll);
        }
    }

    pub(crate) fn set_coll_n_sequences(&mut self, coll: CollRow, n: usize) {
        self.coll_meta[coll as usize].n_sequences = n;
    }

    fn coll_range(&self, coll: CollRow) -> (MemberRow, MemberRow) {
        let s = self.coll_start[coll as usize];
        (s, s + self.coll_len[coll as usize])
    }

    pub(crate) fn coll_meta(&self, coll: CollRow) -> &SequenceCollectionMetadata {
        &self.coll_meta[coll as usize]
    }

    pub(crate) fn is_loaded(&self, coll: CollRow) -> bool {
        self.coll_loaded.get(coll as usize)
    }

    /// Iterate live collection rows, in insertion order.
    pub(crate) fn coll_rows(&self) -> impl Iterator<Item = CollRow> + '_ {
        (0..self.coll_meta.len() as u32).filter(move |&r| !self.coll_dead.get(r as usize))
    }

    /// Members of a loaded collection, in FASTA order. Empty if not loaded.
    pub(crate) fn members(&self, coll: CollRow) -> impl ExactSizeIterator<Item = MemberView<'_>> + '_ {
        let (s, e) = self.coll_range(coll);
        (s..e).map(move |row| MemberView { cat: self, row })
    }

    pub(crate) fn n_members(&self, coll: CollRow) -> usize {
        self.coll_len[coll as usize] as usize
    }

    /// Find the member of `coll` named `name`.
    pub(crate) fn member_by_name(&self, coll: CollRow, name: &str) -> Option<MemberView<'_>> {
        let (start, end) = self.coll_range(coll);
        self.by_name
            .find(self.hash_member(coll, name), |&r| {
                r >= start && r < end && self.member_name(r) == name
            })
            .map(|&row| MemberView { cat: self, row })
    }

    // =====================================================================
    // Maintenance
    // =====================================================================

    /// Release spare capacity after a bulk load.
    pub(crate) fn shrink_to_fit(&mut self) {
        self.sha.shrink_to_fit();
        self.md5.shrink_to_fit();
        self.length.shrink_to_fit();
        self.alphabet.shrink_to_fit();
        self.seq_name_off.shrink_to_fit();
        self.seq_name_len.shrink_to_fit();
        self.seq_desc.shrink_to_fit();
        self.seq_fai.shrink_to_fit();
        self.seq_refcount.shrink_to_fit();
        self.seq_dead.shrink_to_fit();
        self.member_seq.shrink_to_fit();
        self.member_coll.shrink_to_fit();
        self.member_name_off.shrink_to_fit();
        self.member_name_len.shrink_to_fit();
        self.member_desc.shrink_to_fit();
        self.member_fai.shrink_to_fit();
        self.member_dead.shrink_to_fit();
        self.names.shrink_to_fit();
        self.coll_meta.shrink_to_fit();
        self.coll_start.shrink_to_fit();
        self.coll_len.shrink_to_fit();
        self.coll_loaded.shrink_to_fit();
        self.coll_dead.shrink_to_fit();
        // The tables keep their load factor; shrinking them would rehash for
        // little gain right before the first lookups.
    }

    /// Whether enough rows are dead that [`compact`](Self::compact) is worth it.
    pub(crate) fn needs_compaction(&self) -> bool {
        let dead = (self.n_dead_seqs + self.n_dead_members + self.n_dead_colls) as f64;
        let total = (self.sha.len() + self.member_seq.len() + self.coll_meta.len()) as f64;
        total > 0.0 && dead / total > COMPACT_DEAD_FRACTION
    }

    /// Rebuild without tombstoned rows. Returns the sequence-row remap
    /// (`old row -> new row`, `None` for dead rows) so the caller can rekey any
    /// side tables (resident bytes) it holds by sequence row.
    pub(crate) fn compact(&mut self) -> Vec<Option<SeqRow>> {
        let old = std::mem::take(self);
        let mut remap: Vec<Option<SeqRow>> = vec![None; old.sha.len()];
        for r in old.seq_rows() {
            let meta = old.seq(r).to_metadata();
            let (new, _) = self.add_sequence(&meta).expect("re-adding a valid row");
            remap[r as usize] = Some(new);
        }
        for c in old.coll_rows() {
            let new_c = self
                .upsert_collection(old.coll_meta(c).clone())
                .expect("re-adding a valid collection");
            if old.is_loaded(c) {
                let members: Vec<MemberView<'_>> = old.members(c).collect();
                self.set_members(
                    new_c,
                    members.iter().map(|m| NewMember {
                        seq: remap[m.seq_row() as usize].expect("member of a live sequence"),
                        name: m.name(),
                        description: m.description(),
                        fai: m.own_fai(),
                    }),
                )
                .expect("re-adding valid members");
            }
        }
        self.shrink_to_fit();
        remap
    }
}

// =========================================================================
// Pack location columns (packed layout)
// =========================================================================

const NO_PACK: u32 = u32::MAX;

/// Where each sequence row's bytes live in a packed store: pack (interned),
/// offset and length, indexed by [`SeqRow`]. Built from the
/// `sequences.pack.idx` sidecar; replaces keeping the sidecar's 100-byte text
/// rows in RAM and binary-searching them.
#[derive(Debug, Default)]
pub(crate) struct PackColumns {
    packs: Vec<std::sync::Arc<str>>,
    pack_id: Vec<u32>,
    offset: Vec<u64>,
    len: Vec<u64>,
}

impl PackColumns {
    /// Build from sidecar bytes. Rows for digests not in the catalog (another
    /// writer's sequences) are skipped: nothing here can ask for them.
    pub(crate) fn build(cat: &Catalog, sidecar: &[u8]) -> Result<Self> {
        let n = cat.sha.len();
        let mut out = PackColumns {
            packs: Vec::new(),
            pack_id: vec![NO_PACK; n],
            offset: vec![0; n],
            len: vec![0; n],
        };
        let mut intern: HashMap<String, u32> = HashMap::new();
        super::packidx::for_each_row(sidecar, |digest, pack, offset, len| {
            let Some(row) = cat.find_sha(digest) else { return Ok(()) };
            let id = match intern.get(pack) {
                Some(&id) => id,
                None => {
                    let id = out.packs.len() as u32;
                    out.packs.push(std::sync::Arc::from(pack));
                    intern.insert(pack.to_string(), id);
                    id
                }
            };
            let r = row as usize;
            out.pack_id[r] = id;
            out.offset[r] = offset;
            out.len[r] = len;
            Ok(())
        })?;
        Ok(out)
    }

    /// `(pack, offset, len)` of a sequence row, if it is in a pack.
    pub(crate) fn get(&self, row: SeqRow) -> Option<(std::sync::Arc<str>, u64, u64)> {
        let r = row as usize;
        let id = *self.pack_id.get(r)?;
        if id == NO_PACK {
            return None;
        }
        Some((self.packs[id as usize].clone(), self.offset[r], self.len[r]))
    }
}

// =========================================================================
// Views
// =========================================================================

/// A borrowed view of one sequence row.
#[derive(Clone, Copy)]
pub struct SeqView<'a> {
    cat: &'a Catalog,
    row: SeqRow,
}

impl<'a> SeqView<'a> {
    pub(crate) fn row(&self) -> SeqRow {
        self.row
    }
    pub fn sha512t24u(&self) -> Digest32 {
        raw_to_sha(&self.cat.sha[self.row as usize])
    }
    pub fn md5(&self) -> Digest32 {
        raw_to_md5(&self.cat.md5[self.row as usize])
    }
    pub fn length(&self) -> usize {
        self.cat.length[self.row as usize] as usize
    }
    pub fn alphabet(&self) -> AlphabetType {
        self.cat.alphabet[self.row as usize]
    }
    /// The store-wide label: the first name this sequence was added under.
    pub fn name(&self) -> &'a str {
        let i = self.row as usize;
        self.cat.name_at(self.cat.seq_name_off[i], self.cat.seq_name_len[i])
    }
    pub fn description(&self) -> Option<&'a str> {
        if self.cat.seq_desc.is_empty() {
            return None;
        }
        self.cat.seq_desc.get(&self.row).map(|d| &**d)
    }
    pub fn fai(&self) -> Option<&'a FaiMetadata> {
        if self.cat.seq_fai.is_empty() {
            return None;
        }
        self.cat.seq_fai.get(&self.row)
    }
    /// Stored size of this sequence's bytes under `mode` (see
    /// `SequenceMetadataExt::disk_size`).
    pub(crate) fn disk_size(&self, mode: &super::StorageMode) -> usize {
        let len = self.length();
        match mode {
            super::StorageMode::Raw | super::StorageMode::Zstd => len,
            super::StorageMode::Encoded => (len * self.alphabet().bits_per_symbol()).div_ceil(8),
        }
    }
    pub(crate) fn core(&self) -> SeqCore {
        SeqCore {
            sha512t24u: self.sha512t24u(),
            length: self.length(),
            alphabet: self.alphabet(),
        }
    }
    /// Owned metadata, for the public API and bindings.
    pub fn to_metadata(&self) -> SequenceMetadata {
        SequenceMetadata {
            name: self.name().to_string(),
            description: self.description().map(str::to_string),
            length: self.length(),
            sha512t24u: self.sha512t24u().into(),
            md5: self.md5().into(),
            alphabet: self.alphabet(),
            fai: self.fai().cloned(),
        }
    }
}

impl std::fmt::Debug for SeqView<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SeqView")
            .field("sha512t24u", &self.sha512t24u())
            .field("name", &self.name())
            .field("length", &self.length())
            .finish()
    }
}

/// A borrowed view of one collection member.
#[derive(Clone, Copy)]
pub struct MemberView<'a> {
    cat: &'a Catalog,
    row: MemberRow,
}

impl<'a> MemberView<'a> {
    pub(crate) fn seq_row(&self) -> SeqRow {
        self.cat.member_seq[self.row as usize]
    }
    /// The member's name in this collection.
    pub fn name(&self) -> &'a str {
        let i = self.row as usize;
        self.cat.name_at(self.cat.member_name_off[i], self.cat.member_name_len[i])
    }
    pub fn description(&self) -> Option<&'a str> {
        if self.cat.member_desc.is_empty() {
            return None;
        }
        self.cat.member_desc.get(&self.row).map(|d| &**d)
    }
    /// This member's own FAI record (from its FASTA), if any.
    pub(crate) fn own_fai(&self) -> Option<&'a FaiMetadata> {
        if self.cat.member_fai.is_empty() {
            return None;
        }
        self.cat.member_fai.get(&self.row)
    }
    /// FAI for this member, falling back to the sequence's store-wide record.
    pub fn fai(&self) -> Option<&'a FaiMetadata> {
        self.own_fai().or_else(|| self.seq().fai())
    }
    pub fn seq(&self) -> SeqView<'a> {
        SeqView { cat: self.cat, row: self.seq_row() }
    }
    /// Owned metadata carrying THIS collection's name and description.
    pub fn to_metadata(&self) -> SequenceMetadata {
        let s = self.seq();
        SequenceMetadata {
            name: self.name().to_string(),
            description: self.description().map(str::to_string),
            length: s.length(),
            sha512t24u: s.sha512t24u().into(),
            md5: s.md5().into(),
            alphabet: s.alphabet(),
            fai: self.fai().cloned(),
        }
    }
}

impl std::fmt::Debug for MemberView<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MemberView")
            .field("name", &self.name())
            .field("seq", &self.seq())
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::digest::{md5, sha512t24u};

    fn meta(name: &str, bases: &str) -> SequenceMetadata {
        SequenceMetadata {
            name: name.to_string(),
            description: None,
            length: bases.len(),
            sha512t24u: sha512t24u(bases),
            md5: md5(bases),
            alphabet: AlphabetType::Dna2bit,
            fai: None,
        }
    }

    fn coll(digest: &str, n: usize) -> SequenceCollectionMetadata {
        SequenceCollectionMetadata {
            digest: digest.to_string(),
            n_sequences: n,
            names_digest: String::new(),
            sequences_digest: String::new(),
            lengths_digest: String::new(),
            name_length_pairs_digest: None,
            sorted_name_length_pairs_digest: None,
            sorted_sequences_digest: None,
            file_path: None,
        }
    }

    fn add_coll(cat: &mut Catalog, digest: &str, members: &[(&str, &str)]) -> CollRow {
        let rows: Vec<(SeqRow, &str)> = members
            .iter()
            .map(|(name, bases)| (cat.add_sequence(&meta(name, bases)).unwrap().0, *name))
            .collect();
        let c = cat.upsert_collection(coll(digest, members.len())).unwrap();
        cat.set_members(
            c,
            rows.iter().map(|&(seq, name)| NewMember { seq, name, description: None, fai: None }),
        )
        .unwrap();
        c
    }

    #[test]
    fn sequences_dedup_and_resolve_by_either_digest() {
        let mut cat = Catalog::new();
        let (a, new_a) = cat.add_sequence(&meta("chr1", "ACGT")).unwrap();
        let (b, new_b) = cat.add_sequence(&meta("other", "ACGT")).unwrap();
        assert!(new_a && !new_b);
        assert_eq!(a, b);
        assert_eq!(cat.n_sequences(), 1);
        assert_eq!(cat.seq(a).name(), "chr1", "dedup keeps the first label");
        assert_eq!(cat.find_sha(sha512t24u("ACGT")), Some(a));
        assert_eq!(cat.find_sha(format!("SQ.{}", sha512t24u("ACGT"))), Some(a));
        assert_eq!(cat.resolve(md5("ACGT")), Some(a));
        assert_eq!(cat.resolve(sha512t24u("ACGT")), Some(a));
        assert_eq!(cat.resolve("nope"), None);
        let m = cat.seq(a).to_metadata();
        assert_eq!(m.sha512t24u, sha512t24u("ACGT"));
        assert_eq!(m.md5, md5("ACGT"));
        assert_eq!(m.length, 4);
    }

    #[test]
    fn invalid_digests_are_refused() {
        let mut cat = Catalog::new();
        let mut m = meta("x", "ACGT");
        m.sha512t24u = "not-a-digest".into();
        assert!(cat.add_sequence(&m).is_err());
    }

    #[test]
    fn same_name_in_two_collections_and_same_sequence_under_two_names() {
        let mut cat = Catalog::new();
        let c1 = add_coll(&mut cat, "c1", &[("chr1", "AAAA"), ("chr2", "CCCC")]);
        let c2 = add_coll(&mut cat, "c2", &[("1", "AAAA"), ("chr2", "GGGG")]);
        assert_eq!(cat.n_sequences(), 3);
        // chr2 means different sequences in each collection.
        let m1 = cat.member_by_name(c1, "chr2").unwrap();
        let m2 = cat.member_by_name(c2, "chr2").unwrap();
        assert_eq!(m1.seq().sha512t24u(), sha512t24u("CCCC"));
        assert_eq!(m2.seq().sha512t24u(), sha512t24u("GGGG"));
        // AAAA is chr1 in c1 and 1 in c2; the global label is the first name.
        let a = cat.member_by_name(c2, "1").unwrap();
        assert_eq!(a.seq().name(), "chr1");
        assert_eq!(a.to_metadata().name, "1");
        assert!(cat.member_by_name(c1, "1").is_none());
        let names: Vec<&str> = cat.members(c2).map(|m| m.name()).collect();
        assert_eq!(names, vec!["1", "chr2"]);
        assert_eq!(cat.seq_refcount(a.seq_row()), 2);
    }

    #[test]
    fn duplicate_names_point_at_the_last_member() {
        let mut cat = Catalog::new();
        let c = add_coll(&mut cat, "c", &[("dup", "AAAA"), ("dup", "CCCC")]);
        assert_eq!(cat.n_members(c), 2);
        let m = cat.member_by_name(c, "dup").unwrap();
        assert_eq!(m.seq().sha512t24u(), sha512t24u("CCCC"));
    }

    #[test]
    fn delete_then_re_add() {
        let mut cat = Catalog::new();
        let c1 = add_coll(&mut cat, "c1", &[("chr1", "AAAA"), ("chr2", "CCCC")]);
        let _c2 = add_coll(&mut cat, "c2", &[("chr1", "AAAA")]);
        cat.remove_collection(c1);
        assert!(cat.find_coll("c1").is_none());
        assert_eq!(cat.n_collections(), 1);
        let cccc = cat.find_sha(sha512t24u("CCCC")).unwrap();
        assert_eq!(cat.seq_refcount(cccc), 0);
        cat.remove_sequence(cccc);
        assert!(cat.find_sha(sha512t24u("CCCC")).is_none());
        assert!(cat.find_md5(md5("CCCC")).is_none());
        assert_eq!(cat.n_sequences(), 1);

        let c1 = add_coll(&mut cat, "c1", &[("chr1", "AAAA"), ("chr2", "CCCC")]);
        assert_eq!(cat.member_by_name(c1, "chr2").unwrap().seq().length(), 4);
        assert_eq!(cat.n_sequences(), 2);
        assert_eq!(cat.n_collections(), 2);
    }

    #[test]
    fn reloading_members_replaces_them() {
        let mut cat = Catalog::new();
        let c = add_coll(&mut cat, "c", &[("a", "AAAA")]);
        let s = cat.add_sequence(&meta("b", "CCCC")).unwrap().0;
        cat.set_members(c, [NewMember { seq: s, name: "b", description: Some("desc"), fai: None }])
            .unwrap();
        assert!(cat.member_by_name(c, "a").is_none());
        let b = cat.member_by_name(c, "b").unwrap();
        assert_eq!(b.description(), Some("desc"));
        assert_eq!(cat.n_members(c), 1);
    }

    #[test]
    fn compaction_preserves_everything_live() {
        let mut cat = Catalog::new();
        let c1 = add_coll(&mut cat, "c1", &[("x", "AAAA"), ("y", "CCCC")]);
        let c2 = add_coll(&mut cat, "c2", &[("z", "GGGG")]);
        let stub = cat.upsert_collection(coll("c3", 5)).unwrap();
        cat.remove_collection(c1);
        for b in ["AAAA", "CCCC"] {
            let r = cat.find_sha(sha512t24u(b)).unwrap();
            cat.remove_sequence(r);
        }
        assert!(cat.needs_compaction());
        let _ = (c2, stub);
        let remap = cat.compact();
        assert_eq!(remap.iter().filter(|r| r.is_some()).count(), 1);
        assert!(!cat.needs_compaction());
        let c2 = cat.find_coll("c2").unwrap();
        assert_eq!(cat.member_by_name(c2, "z").unwrap().seq().sha512t24u(), sha512t24u("GGGG"));
        let c3 = cat.find_coll("c3").unwrap();
        assert!(!cat.is_loaded(c3));
        assert_eq!(cat.coll_meta(c3).n_sequences, 5);
        assert_eq!(cat.n_sequences(), 1);
    }

    #[test]
    fn many_members_survive_table_growth() {
        let mut cat = Catalog::new();
        let n = 5000;
        let names: Vec<String> = (0..n).map(|i| format!("s{i}")).collect();
        let bases: Vec<String> = (0..n).map(|i| format!("ACGT{i:08b}").replace('0', "A").replace('1', "C")).collect();
        let mut rows = Vec::new();
        for i in 0..n {
            rows.push(cat.add_sequence(&meta(&names[i], &bases[i])).unwrap().0);
        }
        let ca = cat.upsert_collection(coll("a", n)).unwrap();
        cat.set_members(ca, (0..n).map(|i| NewMember { seq: rows[i], name: &names[i], description: None, fai: None }))
            .unwrap();
        let cb = cat.upsert_collection(coll("b", n)).unwrap();
        cat.set_members(cb, (0..n).map(|i| NewMember { seq: rows[n - 1 - i], name: &names[i], description: None, fai: None }))
            .unwrap();
        for i in (0..n).step_by(97) {
            assert_eq!(cat.member_by_name(ca, &names[i]).unwrap().seq_row(), rows[i]);
            assert_eq!(cat.member_by_name(cb, &names[i]).unwrap().seq_row(), rows[n - 1 - i]);
        }
    }
}

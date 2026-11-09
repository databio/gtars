//! `MappedFiles`: the stored bytes of one subject, borrowed from a mapping where possible.
//!
//! A local standalone store maps its own `.seq` files on demand through a bounded cache. A
//! remote store is pulled into the shared cache first (`open`) and then mapped from there.

use super::OpenedStore;
use crate::error::Error;
use crate::plan::SubjectPlan;
use memmap2::Mmap;
use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

/// Stored bytes of one subject: a shared mapping (kept alive by the `Arc`) or an owned copy.
/// (The plan wrote `Cow<'_, [u8]>`; a mapping evicted from the cache would leave a borrowed
/// slice dangling, so the handle owns its mapping instead.)
pub enum SubjectBytes {
    Mapped(Arc<Mmap>),
    Owned(Vec<u8>),
}

impl SubjectBytes {
    #[inline]
    pub fn as_slice(&self) -> &[u8] {
        match self {
            SubjectBytes::Mapped(m) => &m[..],
            SubjectBytes::Owned(v) => v,
        }
    }
}

/// Local standalone store: one `.seq` file per subject, memory-mapped on demand.
pub struct MappedFiles {
    paths: Vec<PathBuf>,
    expected: Vec<Option<u64>>,
    digests: Vec<String>,
    cache: Mutex<MapCache>,
}

struct MapCache {
    maps: HashMap<u32, Arc<Mmap>>,
    order: std::collections::VecDeque<u32>,
    limit: usize,
}

impl MappedFiles {
    /// Map cache bound (files); only the traceback and the aligned-sequence output map
    /// chromosome-sized files (`read_at` serves the first search stage), and they re-read only
    /// the subjects with hits.
    pub const DEFAULT_CACHE: usize = 64;
    /// Files below this size are read with one `read` call into an owned buffer instead of
    /// being mapped: on a transcriptome store (384k files of a few hundred bytes) the
    /// mmap/page-fault/munmap round trip was about 60 µs per subject and 95% of a megablast
    /// search; a read is a few µs. Chromosome-sized files are mapped for windows (page-cache
    /// backed, borrowed, not copied) and read piece by piece for the first search stage.
    pub const MMAP_MIN: u64 = 1 << 20;

    pub fn new(store: &OpenedStore, plan: &SubjectPlan) -> Result<MappedFiles, Error> {
        let mut paths = Vec::with_capacity(plan.subjects.len());
        let mut expected = Vec::with_capacity(plan.subjects.len());
        let mut digests = Vec::with_capacity(plan.subjects.len());
        for s in &plan.subjects {
            let p = store
                .store
                .sequence_file_path(&s.digest)
                .ok_or_else(|| Error::sequence(&s.digest, "store has no local sequence path"))?;
            paths.push(p);
            expected.push(s.stored_len);
            digests.push(s.digest.clone());
        }
        Ok(MappedFiles {
            paths,
            expected,
            digests,
            cache: Mutex::new(MapCache {
                maps: HashMap::new(),
                order: std::collections::VecDeque::new(),
                limit: Self::DEFAULT_CACHE,
            }),
        })
    }

    fn map(&self, oid: u32) -> Result<Arc<Mmap>, Error> {
        let mut c = self.cache.lock().unwrap();
        if let Some(m) = c.maps.get(&oid) {
            return Ok(m.clone());
        }
        let p = &self.paths[oid as usize];
        let f = std::fs::File::open(p).map_err(|e| Error::io(p, e))?;
        // SAFETY: the store is read-only while a search runs; a concurrent writer would be a
        // data race the caller must avoid (same contract as every other reader of the store).
        let m = unsafe { Mmap::map(&f) }.map_err(|e| Error::io(p, e))?;
        if let Some(n) = self.expected[oid as usize]
            && m.len() as u64 != n
        {
            return Err(Error::sequence(
                &self.digests[oid as usize],
                format!("{} is {} bytes, expected {n}", p.display(), m.len()),
            ));
        }
        let m = Arc::new(m);
        if c.maps.len() >= c.limit
            && let Some(old) = c.order.pop_front()
        {
            c.maps.remove(&old);
        }
        c.maps.insert(oid, m.clone());
        c.order.push_back(oid);
        Ok(m)
    }

    /// Stored bytes of one subject (packed, raw or zstd as the store keeps them).
    pub fn bytes(&self, oid: u32) -> Result<SubjectBytes, Error> {
        match self.expected[oid as usize] {
            Some(n) if n < Self::MMAP_MIN => {
                use std::io::Read;
                let p = &self.paths[oid as usize];
                let mut f = std::fs::File::open(p).map_err(|e| Error::io(p, e))?;
                // exactly the indexed length, with no size query (a `statx` per file was a
                // tenth of the search); a shorter file is an error, extra bytes are ignored
                let mut v = vec![0u8; n as usize];
                f.read_exact(&mut v).map_err(|e| {
                    Error::sequence(
                        &self.digests[oid as usize],
                        format!(
                            "{} is shorter than the indexed {n} bytes ({e})",
                            p.display()
                        ),
                    )
                })?;
                Ok(SubjectBytes::Owned(v))
            }
            _ => Ok(SubjectBytes::Mapped(self.map(oid)?)),
        }
    }

    /// A subject's stored file, for `read_at`.
    pub fn open(&self, oid: u32) -> Result<File, Error> {
        let p = &self.paths[oid as usize];
        File::open(p).map_err(|e| Error::io(p, e))
    }

    /// Bytes `off..off + out.len()` of a subject's stored file (`file` from `open`). The search
    /// reads a long subject this way, piece by piece, rather than through its mapping: on GPFS
    /// every 4 KB page of a mapping costs a page fault (about 5 µs, more with many threads),
    /// while a read is one copy.
    pub fn read_at(
        &self,
        oid: u32,
        file: &mut File,
        off: u64,
        out: &mut [u8],
    ) -> Result<(), Error> {
        file.seek(SeekFrom::Start(off))
            .and_then(|_| file.read_exact(out))
            .map_err(|e| Error::io(&self.paths[oid as usize], e))
    }
}

/// The access path for a store: `MappedFiles` over the store's own files, or, for a remote
/// store, over the cache after the subjects are pulled into it ([`OpenedStore::pull`]).
pub fn open(store: &mut OpenedStore, plan: &SubjectPlan) -> Result<MappedFiles, Error> {
    let digests: Vec<String> = plan.subjects.iter().map(|s| s.digest.clone()).collect();
    store.pull(&digests)?;
    MappedFiles::new(store, plan)
}

//! `MappedFiles`: the stored bytes of one subject, borrowed from a mapping where possible.
//!
//! A subject's stored record is a `.seq` file (standalone store) or a span of a pack file
//! (packed store, `RefgetStore::on_disk_packed`); gtars-refget's `stored_span` says which. A
//! pack holds each record verbatim, so a span is read or mapped like a `.seq` file at an
//! offset, through a bounded map cache. A remote store is pulled into the shared cache first
//! (`open`), as `.seq` files whatever the remote's layout, and then mapped from there.

use super::OpenedStore;
use crate::error::Error;
use crate::plan::SubjectPlan;
use memmap2::{Mmap, MmapOptions};
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::os::unix::fs::FileExt;
use std::path::PathBuf;
use std::sync::{Arc, Mutex, OnceLock};

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

/// Where one subject's record sits: `files[file]` at `offset`; `len` is the exact length of
/// a pack span, `None` for a whole `.seq` file.
struct Span {
    file: u32,
    offset: u64,
    len: Option<u64>,
}

impl Span {
    fn in_pack(&self) -> bool {
        self.len.is_some()
    }
}

/// Every subject's stored record (a `.seq` file or a pack span), memory-mapped or read on
/// demand.
pub struct MappedFiles {
    /// One per `.seq` file, or one per pack.
    files: Vec<PathBuf>,
    /// An open handle per pack, opened on first use and shared by the search threads. A
    /// `.seq` file's slot stays empty: a standalone transcriptome store would otherwise hold
    /// hundreds of thousands of open files.
    shared: Vec<OnceLock<File>>,
    /// By oid.
    spans: Vec<Span>,
    /// The plan's stored length, else the pack span's.
    expected: Vec<Option<u64>>,
    digests: Vec<String>,
    cache: Mutex<MapCache>,
}

struct MapCache {
    maps: HashMap<u32, Arc<Mmap>>,
    order: std::collections::VecDeque<u32>,
    limit: usize,
}

/// A handle `MappedFiles::read_at` reads one subject through: its own `.seq` file, or the
/// shared handle of its pack, with the record's offset in it.
pub struct SubjectReader<'a> {
    file: ReaderFile<'a>,
    base: u64,
}

enum ReaderFile<'a> {
    Owned(File),
    Shared(&'a File),
}

impl SubjectReader<'_> {
    fn file(&self) -> &File {
        match &self.file {
            ReaderFile::Owned(f) => f,
            ReaderFile::Shared(f) => f,
        }
    }
}

impl MappedFiles {
    /// Map cache bound (subjects); only the traceback and the aligned-sequence output map
    /// chromosome-sized records (`read_at` serves the first search stage), and they re-read
    /// only the subjects with hits.
    pub const DEFAULT_CACHE: usize = 64;
    /// Records below this size are read with one `read` call into an owned buffer instead of
    /// being mapped: on a transcriptome store (384k files of a few hundred bytes) the
    /// mmap/page-fault/munmap round trip was about 60 µs per subject and 95% of a megablast
    /// search; a read is a few µs. Chromosome-sized records are mapped for windows (page-cache
    /// backed, borrowed, not copied) and read piece by piece for the first search stage.
    pub const MMAP_MIN: u64 = 1 << 20;

    pub fn new(store: &OpenedStore, plan: &SubjectPlan) -> Result<MappedFiles, Error> {
        let n = plan.subjects.len();
        let mut files = Vec::new();
        let mut packs: HashMap<PathBuf, u32> = HashMap::new();
        let mut spans = Vec::with_capacity(n);
        let mut expected = Vec::with_capacity(n);
        let mut digests = Vec::with_capacity(n);
        for s in &plan.subjects {
            let span = store
                .store
                .stored_span(&s.digest)
                .map_err(|e| Error::sequence(&s.digest, format!("{e:#}")))?;
            if let (Some(want), Some(got)) = (s.stored_len, span.len)
                && want != got
            {
                return Err(Error::sequence(
                    &s.digest,
                    format!("pack index says {got} bytes, expected {want}"),
                ));
            }
            let file = match span.len {
                Some(_) => *packs.entry(span.path).or_insert_with_key(|p| {
                    files.push(p.clone());
                    files.len() as u32 - 1
                }),
                None => {
                    files.push(span.path);
                    files.len() as u32 - 1
                }
            };
            spans.push(Span {
                file,
                offset: span.offset,
                len: span.len,
            });
            expected.push(s.stored_len.or(span.len));
            digests.push(s.digest.clone());
        }
        let shared = (0..files.len()).map(|_| OnceLock::new()).collect();
        Ok(MappedFiles {
            files,
            shared,
            spans,
            expected,
            digests,
            cache: Mutex::new(MapCache {
                maps: HashMap::new(),
                order: std::collections::VecDeque::new(),
                limit: Self::DEFAULT_CACHE,
            }),
        })
    }

    fn path(&self, oid: u32) -> &PathBuf {
        &self.files[self.spans[oid as usize].file as usize]
    }

    /// An I/O error naming the file; for a pack span also the offset and the digest.
    fn io_error(&self, oid: u32, e: std::io::Error) -> Error {
        let span = &self.spans[oid as usize];
        match span.len {
            None => Error::io(self.path(oid), e),
            Some(len) => Error::sequence(
                &self.digests[oid as usize],
                format!(
                    "{len} bytes at offset {} of {}: {e}",
                    span.offset,
                    self.path(oid).display()
                ),
            ),
        }
    }

    /// The shared open handle of a subject's pack.
    fn pack_file(&self, oid: u32) -> Result<&File, Error> {
        let cell = &self.shared[self.spans[oid as usize].file as usize];
        if let Some(f) = cell.get() {
            return Ok(f);
        }
        let p = self.path(oid);
        let f = File::open(p).map_err(|e| Error::io(p, e))?;
        // a racing thread may have opened it too; either handle reads the same sealed pack
        let _ = cell.set(f);
        Ok(cell.get().expect("set above"))
    }

    fn map(&self, oid: u32) -> Result<Arc<Mmap>, Error> {
        let mut c = self.cache.lock().unwrap();
        if let Some(m) = c.maps.get(&oid) {
            return Ok(m.clone());
        }
        let span = &self.spans[oid as usize];
        let p = self.path(oid);
        // SAFETY: the store is read-only while a search runs; a concurrent writer would be a
        // data race the caller must avoid (same contract as every other reader of the store).
        let m = match span.len {
            // only the record, never the whole pack; memmap2 aligns the offset to a page
            // itself, so the mapping derefs to exactly the record
            Some(len) => {
                let f = self.pack_file(oid)?;
                unsafe {
                    MmapOptions::new()
                        .offset(span.offset)
                        .len(len as usize)
                        .map(f)
                }
                .map_err(|e| self.io_error(oid, e))?
            }
            None => {
                let f = File::open(p).map_err(|e| Error::io(p, e))?;
                unsafe { Mmap::map(&f) }.map_err(|e| Error::io(p, e))?
            }
        };
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
            Some(0) => Ok(SubjectBytes::Owned(Vec::new())),
            Some(n) if n < Self::MMAP_MIN => {
                let span = &self.spans[oid as usize];
                let mut v = vec![0u8; n as usize];
                if span.in_pack() {
                    // one positioned read on the pack's shared handle, no open/close
                    self.pack_file(oid)?
                        .read_exact_at(&mut v, span.offset)
                        .map_err(|e| self.io_error(oid, e))?;
                    return Ok(SubjectBytes::Owned(v));
                }
                let p = self.path(oid);
                let mut f = File::open(p).map_err(|e| Error::io(p, e))?;
                // exactly the indexed length, with no size query (a `statx` per file was a
                // tenth of the search); a shorter file is an error, extra bytes are ignored
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

    /// A reader of a subject's stored record, for `read_at`: its `.seq` file, opened, or its
    /// pack's shared handle.
    pub fn reader(&self, oid: u32) -> Result<SubjectReader<'_>, Error> {
        let span = &self.spans[oid as usize];
        let file = if span.in_pack() {
            ReaderFile::Shared(self.pack_file(oid)?)
        } else {
            let p = self.path(oid);
            ReaderFile::Owned(File::open(p).map_err(|e| Error::io(p, e))?)
        };
        Ok(SubjectReader {
            file,
            base: span.offset,
        })
    }

    /// Bytes `off..off + out.len()` of a subject's stored record (`r` from `reader`), with a
    /// positioned read (no seek, so a shared pack handle is safe across threads). The search
    /// reads a long subject this way, piece by piece, rather than through its mapping: on
    /// GPFS every 4 KB page of a mapping costs a page fault (about 5 µs, more with many
    /// threads), while a read is one copy.
    pub fn read_at(
        &self,
        oid: u32,
        r: &SubjectReader<'_>,
        off: u64,
        out: &mut [u8],
    ) -> Result<(), Error> {
        r.file()
            .read_exact_at(out, r.base + off)
            .map_err(|e| self.io_error(oid, e))
    }
}

/// The access path for a store: `MappedFiles` over the store's own records, or, for a remote
/// store, over the cache after the subjects are pulled into it ([`OpenedStore::pull`]).
pub fn open(store: &mut OpenedStore, plan: &SubjectPlan) -> Result<MappedFiles, Error> {
    let digests: Vec<String> = plan.subjects.iter().map(|s| s.digest.clone()).collect();
    store.pull(&digests)?;
    MappedFiles::new(store, plan)
}

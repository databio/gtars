//! Where subjects' stored bytes come from: the store's own sequence files (a local store) or
//! the cache a remote store's selection was pulled into ([`OpenedStore::pull`]).
//!
//! [`OpenedStore::pull`]: super::OpenedStore::pull

use crate::error::Error;
use gtars_refget::store::{PrefetchItem, PrefetchOptions, PrefetchOrder, ReadonlyRefgetStore};

/// Stored bytes (packed, raw or zstd, as the store keeps them) of one subject.
pub type StoredBytes = Vec<u8>;

/// Stored bytes from [`SeqFetch::read_ahead`]; dropping them frees their share of its memory.
pub type HeldBytes = Box<dyn AsRef<[u8]> + Send>;

/// Stored bytes of many subjects, in the order asked for.
pub type HeldStream<'a> = Box<dyn Iterator<Item = Result<HeldBytes, Error>> + Send + 'a>;

/// Where a subject's stored bytes come from. The writer, `verify`, the N-count scan and the
/// reference FASTA see only this trait; test stubs implement it over `Vec<u8>`s.
pub trait SeqFetch: Sync {
    /// One digest's stored bytes.
    fn read(&self, digest: &str) -> Result<StoredBytes, Error>;

    /// The stored bytes of `digests`, in order. This reads them one at a time when asked;
    /// [`StoreFiles`] reads ahead on its own threads, so reading overlaps the caller's work.
    fn read_ahead(&self, digests: &[String]) -> Result<HeldStream<'_>, Error> {
        let digests = digests.to_vec();
        Ok(Box::new(digests.into_iter().map(move |d| {
            self.read(&d).map(|v| Box::new(v) as HeldBytes)
        })))
    }
}

/// A prefetched sequence, as [`HeldBytes`].
struct Held(PrefetchItem);

impl AsRef<[u8]> for Held {
    fn as_ref(&self) -> &[u8] {
        self.0.bytes()
    }
}

/// The one real [`SeqFetch`]: each sequence's file in the store (local) or in the cache
/// (remote). A remote sequence that was not pulled yet is pulled first, on its own.
pub struct StoreFiles<'a> {
    store: &'a ReadonlyRefgetStore,
    remote: bool,
}

impl<'a> StoreFiles<'a> {
    /// The files of an opened store. The sequence index must be loaded.
    pub fn new(store: &'a super::OpenedStore) -> Self {
        StoreFiles {
            store: &store.store,
            remote: store.is_remote(),
        }
    }
}

impl SeqFetch for StoreFiles<'_> {
    fn read(&self, digest: &str) -> Result<StoredBytes, Error> {
        let path = |s: &ReadonlyRefgetStore| {
            s.sequence_file_path(digest)
                .ok_or_else(|| Error::sequence(digest, "the store has no file for it"))
        };
        let mut p = path(self.store)?;
        if self.remote && !p.is_file() {
            super::pull(self.store, &[digest.to_string()])?;
            p = path(self.store)?;
        }
        std::fs::read(&p).map_err(|e| Error::io(&p, e))
    }

    /// gtars-refget's prefetch: 8 reader threads holding at most 1 GiB, in request order.
    /// A remote sequence missing from the cache is downloaded and added to it.
    fn read_ahead(&self, digests: &[String]) -> Result<HeldStream<'_>, Error> {
        let opts = PrefetchOptions {
            order: PrefetchOrder::Input,
            write_cache: self.remote,
            ..Default::default()
        };
        let run = self.store.prefetch(digests, opts).map_err(Error::Store)?;
        Ok(Box::new(run.map(|r| {
            r.map(|i| Box::new(Held(i)) as HeldBytes)
                .map_err(Error::Store)
        })))
    }
}

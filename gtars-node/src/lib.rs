use std::io::Read;
use std::sync::Mutex;

use napi::bindgen_prelude::*;
use napi::threadsafe_function::{
    ErrorStrategy, ThreadsafeFunction, ThreadsafeFunctionCallMode,
};
use napi::JsFunction;
use napi_derive::napi;

use gtars_refget::store::{FastaImportOptions, RefgetStore as RustRefgetStore};

fn to_napi_err(e: anyhow::Error) -> napi::Error {
    napi::Error::from_reason(format!("{:#}", e))
}

fn lock_err(e: std::sync::PoisonError<impl std::fmt::Debug>) -> napi::Error {
    napi::Error::from_reason(format!("Lock poisoned: {:?}", e))
}

// -- Metadata return types (plain JS objects via #[napi(object)]) --

#[napi(object, js_name = "SequenceMetadata")]
pub struct SequenceMetadataJs {
    pub name: String,
    pub length: u32,
    pub sha512t24u: String,
    pub md5: String,
}

#[napi(object, js_name = "CollectionMetadata")]
pub struct CollectionMetadataJs {
    pub digest: String,
    pub n_sequences: u32,
    pub names_digest: String,
    pub sequences_digest: String,
    pub lengths_digest: String,
}

#[napi(object, js_name = "StoreStats")]
pub struct StoreStatsJs {
    pub n_sequences: u32,
    pub n_sequences_loaded: u32,
    pub n_collections: u32,
    pub n_collections_loaded: u32,
    pub storage_mode: String,
}

// -- Main wrapper class --

#[napi]
pub struct RefgetStore {
    inner: Mutex<RustRefgetStore>,
}

#[napi]
impl RefgetStore {
    // -- Factory constructors --

    #[napi(factory)]
    pub fn open_local(path: String) -> Result<Self> {
        let store = RustRefgetStore::open_local(&path).map_err(to_napi_err)?;
        Ok(Self {
            inner: Mutex::new(store),
        })
    }

    #[napi(factory)]
    pub fn open_remote(cache_path: String, url: String) -> Result<Self> {
        let store =
            RustRefgetStore::open_remote(&cache_path, &url).map_err(to_napi_err)?;
        Ok(Self {
            inner: Mutex::new(store),
        })
    }

    #[napi(factory)]
    pub fn in_memory() -> Self {
        Self {
            inner: Mutex::new(RustRefgetStore::in_memory()),
        }
    }

    // -- Sequence retrieval (auto-loads stubs from disk/remote) --

    #[napi]
    pub fn get_sequence(&self, digest: String) -> Result<String> {
        let mut store = self.inner.lock().map_err(lock_err)?;
        store.load_sequence(&digest).map_err(to_napi_err)?;
        let record = store.get_sequence(&digest).map_err(to_napi_err)?;
        record
            .decode()
            .ok_or_else(|| napi::Error::from_reason("Failed to decode sequence"))
    }

    #[napi]
    pub fn get_substring(&self, digest: String, start: u32, end: u32) -> Result<String> {
        let mut store = self.inner.lock().map_err(lock_err)?;
        store.load_sequence(&digest).map_err(to_napi_err)?;
        store
            .get_substring(&digest, start as usize, end as usize)
            .map_err(to_napi_err)
    }

    #[napi]
    pub fn get_sequence_by_name(
        &self,
        collection_digest: String,
        name: String,
    ) -> Result<String> {
        let mut store = self.inner.lock().map_err(lock_err)?;
        // Load the collection first to populate the name lookup
        store.load_collection(&collection_digest).map_err(to_napi_err)?;
        let record = store
            .get_sequence_by_name(&collection_digest, &name)
            .map_err(to_napi_err)?;
        // Get the digest so we can load the sequence data
        let seq_digest = record.metadata().sha512t24u.clone();
        store.load_sequence(&seq_digest).map_err(to_napi_err)?;
        let record = store
            .get_sequence_by_name(&collection_digest, &name)
            .map_err(to_napi_err)?;
        record
            .decode()
            .ok_or_else(|| napi::Error::from_reason("Failed to decode sequence"))
    }

    // -- Listing / metadata --

    #[napi]
    pub fn list_collections(
        &self,
        page: Option<u32>,
        page_size: Option<u32>,
    ) -> Result<Vec<CollectionMetadataJs>> {
        let store = self.inner.lock().map_err(lock_err)?;
        let paged = store
            .list_collections(
                page.unwrap_or(0) as usize,
                page_size.map_or(usize::MAX, |s| s as usize),
                &[],
            )
            .map_err(to_napi_err)?;
        Ok(paged
            .results
            .into_iter()
            .map(|m| CollectionMetadataJs {
                digest: m.digest,
                n_sequences: m.n_sequences as u32,
                names_digest: m.names_digest,
                sequences_digest: m.sequences_digest,
                lengths_digest: m.lengths_digest,
            })
            .collect())
    }

    #[napi]
    pub fn list_sequences(&self) -> Result<Vec<SequenceMetadataJs>> {
        let store = self.inner.lock().map_err(lock_err)?;
        Ok(store
            .list_sequences()
            .into_iter()
            .map(|m| SequenceMetadataJs {
                name: m.name,
                length: m.length as u32,
                sha512t24u: m.sha512t24u,
                md5: m.md5,
            })
            .collect())
    }

    #[napi]
    pub fn get_collection_metadata(
        &self,
        digest: String,
    ) -> Result<Option<CollectionMetadataJs>> {
        let store = self.inner.lock().map_err(lock_err)?;
        Ok(store.get_collection_metadata(&digest).map(|m| {
            CollectionMetadataJs {
                digest: m.digest.clone(),
                n_sequences: m.n_sequences as u32,
                names_digest: m.names_digest.clone(),
                sequences_digest: m.sequences_digest.clone(),
                lengths_digest: m.lengths_digest.clone(),
            }
        }))
    }

    #[napi]
    pub fn get_sequence_metadata(
        &self,
        digest: String,
    ) -> Result<Option<SequenceMetadataJs>> {
        let store = self.inner.lock().map_err(lock_err)?;
        Ok(store.get_sequence_metadata(&digest).map(|m| {
            SequenceMetadataJs {
                name: m.name.clone(),
                length: m.length as u32,
                sha512t24u: m.sha512t24u.clone(),
                md5: m.md5.clone(),
            }
        }))
    }

    #[napi]
    pub fn stats(&self) -> Result<StoreStatsJs> {
        let store = self.inner.lock().map_err(lock_err)?;
        let s = store.stats();
        Ok(StoreStatsJs {
            n_sequences: s.n_sequences as u32,
            n_sequences_loaded: s.n_sequences_loaded as u32,
            n_collections: s.n_collections as u32,
            n_collections_loaded: s.n_collections_loaded as u32,
            storage_mode: s.storage_mode,
        })
    }

    // -- Cache management --

    #[napi]
    pub fn ensure_decoded(&self, digest: String) -> Result<()> {
        let mut store = self.inner.lock().map_err(lock_err)?;
        store.ensure_decoded(&digest).map_err(to_napi_err)
    }

    #[napi]
    pub fn clear_decoded_cache(&self) -> Result<()> {
        let mut store = self.inner.lock().map_err(lock_err)?;
        store.clear_decoded_cache();
        Ok(())
    }

    // -- Store building --

    #[napi]
    pub fn add_fasta(&self, fasta_path: String) -> Result<()> {
        let mut store = self.inner.lock().map_err(lock_err)?;
        store
            .add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new())
            .map_err(to_napi_err)?;
        Ok(())
    }

    #[napi]
    pub fn write(&self, path: String) -> Result<()> {
        let store = self.inner.lock().map_err(lock_err)?;
        store.write_store_to_dir(&path, None).map_err(to_napi_err)
    }

    // -- Seqcol comparison --

    /// Compare two sequence collections by digest.
    /// Returns the comparison result as a JSON string (call JSON.parse() on the JS side).
    #[napi]
    pub fn compare(
        &self,
        digest_a: String,
        digest_b: String,
    ) -> Result<String> {
        let mut store = self.inner.lock().map_err(lock_err)?;
        let comparison = store.compare(&digest_a, &digest_b).map_err(to_napi_err)?;
        serde_json::to_string(&comparison)
            .map_err(|e| napi::Error::from_reason(format!("Serialization error: {}", e)))
    }

    // -- Streaming --

    /// Stream a (sub)sequence to JavaScript callbacks.
    ///
    /// Spawns a background thread that reads the sequence in chunks and
    /// invokes `on_chunk(Buffer)` for each chunk, `on_error(message)` on
    /// failure, or `on_end()` on successful completion. Exactly one of
    /// `on_end` or `on_error` is invoked as the terminal callback. The
    /// JS-side `streamSequence` wrapper constructs a `Readable` from these
    /// callbacks.
    #[napi]
    pub fn stream_sequence_internal(
        &self,
        digest: String,
        start: Option<u32>,
        end: Option<u32>,
        on_chunk: JsFunction,
        on_error: JsFunction,
        on_end: JsFunction,
    ) -> Result<()> {
        // Build the reader under the mutex, then drop the lock so the worker
        // thread never contends with JS callers.
        //
        // Note: we deliberately do NOT call `load_sequence` here. Pulling the
        // entire sequence into memory before streaming defeats the O(1) memory
        // guarantee of `stream_sequence` — peak memory would grow to roughly
        // 2x the encoded sequence for chromosome-sized reads. `stream_sequence`
        // already falls back to the on-disk / remote byte source when the
        // record is a `Stub`, so no pre-loading is required.
        let reader: Box<dyn Read + Send> = {
            let store = self.inner.lock().map_err(lock_err)?;
            store
                .stream_sequence(&digest, start.map(|v| v as u64), end.map(|v| v as u64))
                .map_err(to_napi_err)?
        };

        // Build threadsafe functions for the three JS callbacks.
        let tsfn_chunk: ThreadsafeFunction<Vec<u8>, ErrorStrategy::Fatal> = on_chunk
            .create_threadsafe_function(0, |ctx| {
                let buf: Vec<u8> = ctx.value;
                ctx.env.create_buffer_with_data(buf).map(|b| vec![b.into_raw()])
            })?;
        let tsfn_error: ThreadsafeFunction<String, ErrorStrategy::Fatal> = on_error
            .create_threadsafe_function(0, |ctx| {
                let s: String = ctx.value;
                ctx.env.create_string(&s).map(|v| vec![v])
            })?;
        let tsfn_end: ThreadsafeFunction<(), ErrorStrategy::Fatal> = on_end
            .create_threadsafe_function(0, |_ctx| Ok(Vec::<napi::JsUndefined>::new()))?;

        std::thread::spawn(move || {
            let mut reader = reader;
            let mut buf = [0u8; 65536];
            loop {
                match reader.read(&mut buf) {
                    Ok(0) => {
                        tsfn_end.call((), ThreadsafeFunctionCallMode::NonBlocking);
                        break;
                    }
                    Ok(n) => {
                        let chunk = buf[..n].to_vec();
                        // Blocking so the tsfn queue applies backpressure.
                        tsfn_chunk.call(chunk, ThreadsafeFunctionCallMode::Blocking);
                    }
                    Err(e) => {
                        tsfn_error.call(
                            format!("stream decode error: {}", e),
                            ThreadsafeFunctionCallMode::NonBlocking,
                        );
                        break;
                    }
                }
            }
        });

        Ok(())
    }
}

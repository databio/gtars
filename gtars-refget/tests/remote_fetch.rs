//! Remote RefgetStores over a local HTTP fixture: parallel fetch, the shared
//! cache, resume, verification, metadata batching and freshness. (Ranged
//! downloads of big files are tested in the crate, where the piece size can
//! be made tiny.)

mod support;

use std::collections::HashMap;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::time::Duration;

use gtars_refget::store::{
    ExportOptions, ExportSelection, FastaImportOptions, FetchConfig, FhrMetadata, ItemSource,
    LabelTemplate, PrefetchOptions, PrefetchOrder, RefgetStore, RemoteMissing, StorageMode,
    SyncStrategy,
};
use support::http_fixture::{FaultKind, Fixture};
use tempfile::{TempDir, tempdir};

fn cfg() -> FetchConfig {
    FetchConfig {
        jobs: 4,
        ..FetchConfig::default()
    }
}

fn seq(alphabet: &[u8], n: usize, salt: usize) -> String {
    (0..n)
        .map(|i| alphabet[(i * 31 + salt * 7 + i / 3) % alphabet.len()] as char)
        .collect()
}

/// Sequences of every alphabet and several sizes.
fn fasta_text() -> String {
    let mut out = String::new();
    for (i, n) in [24usize, 28, 32, 92].iter().enumerate() {
        out.push_str(&format!(">acgt{}\n{}\n", i, seq(b"ACGT", *n, i)));
    }
    out.push_str(&format!(">n3bit\n{}\n", seq(b"ACGTN", 40, 1)));
    out.push_str(&format!(">iupac\n{}\n", seq(b"ACGTRYSWKM", 33, 2)));
    out.push_str(&format!(">prot\n{}\n", seq(b"MKVLAEGHST", 29, 3)));
    out
}

struct Built {
    _dir: TempDir,
    store: PathBuf,
    collections: Vec<String>,
}

fn build(mode: StorageMode, fastas: &[String]) -> Built {
    let dir = tempdir().unwrap();
    let store_path = dir.path().join("store");
    let mut store = RefgetStore::on_disk(&store_path).unwrap();
    store.set_encoding_mode(mode).unwrap();
    let mut collections = Vec::new();
    for (i, text) in fastas.iter().enumerate() {
        let fa = dir.path().join(format!("in{}.fa", i));
        fs::write(&fa, text).unwrap();
        let (meta, _) = store
            .add_sequence_collection_from_fasta(&fa, FastaImportOptions::new())
            .unwrap();
        collections.push(meta.digest);
    }
    store.write().unwrap();
    Built {
        _dir: dir,
        store: store_path,
        collections,
    }
}

fn digests(store_path: &Path) -> Vec<String> {
    let mut s = RefgetStore::open_local(store_path).unwrap();
    let mut d: Vec<String> = s
        .list_sequences()
        .unwrap()
        .into_iter()
        .map(|m| m.sha512t24u)
        .collect();
    d.sort();
    d
}

fn local_bytes(store_path: &Path) -> HashMap<String, Vec<u8>> {
    let mut s = RefgetStore::open_local(store_path).unwrap();
    let ds = digests(store_path);
    s.prefetch(&ds, PrefetchOptions::default())
        .unwrap()
        .map(|i| {
            let i = i.unwrap();
            (i.digest.clone(), i.bytes().to_vec())
        })
        .collect()
}

fn open(cache: &Path, url: &str, c: &FetchConfig) -> RefgetStore {
    let mut s = RefgetStore::open_remote(cache, url, c).unwrap();
    s.set_quiet(true);
    s
}

fn seq_requests(fx: &Fixture) -> usize {
    fx.requests_matching(".seq").len()
}

fn files_under(dir: &Path) -> usize {
    walk(dir).len()
}

fn walk(dir: &Path) -> Vec<PathBuf> {
    let mut out = Vec::new();
    let mut stack = vec![dir.to_path_buf()];
    while let Some(d) = stack.pop() {
        let Ok(entries) = fs::read_dir(&d) else {
            continue;
        };
        for e in entries.flatten() {
            if e.path().is_dir() {
                stack.push(e.path());
            } else {
                out.push(e.path());
            }
        }
    }
    out
}

fn export(store: &mut RefgetStore, sel: &ExportSelection, label: &str, out: &Path) -> String {
    let opts = ExportOptions {
        label: LabelTemplate::parse(label).unwrap(),
        ..Default::default()
    };
    store.export_fasta(sel, out, &opts).unwrap();
    fs::read_to_string(out).unwrap()
}

// ============================================================================
// Byte identity
// ============================================================================

#[test]
fn remote_bytes_match_local_for_every_mode() {
    for mode in [StorageMode::Encoded, StorageMode::Raw, StorageMode::Zstd] {
        let b = build(mode, &[fasta_text()]);
        let fx = Fixture::serve(&b.store, &["s"]);
        let expected = local_bytes(&b.store);
        let ds = digests(&b.store);
        let cache = tempdir().unwrap();

        for order in [PrefetchOrder::Any, PrefetchOrder::Input] {
            let mut remote = open(cache.path(), &fx.url("s"), &cfg());
            let got: Vec<_> = remote
                .prefetch(
                    &ds,
                    PrefetchOptions {
                        order,
                        ..Default::default()
                    },
                )
                .unwrap()
                .map(|i| i.unwrap())
                .collect();
            assert_eq!(got.len(), ds.len());
            for (pos, item) in got.iter().enumerate() {
                assert_eq!(
                    item.bytes(),
                    expected[&item.digest].as_slice(),
                    "{mode:?} {order:?} {}",
                    item.digest
                );
                if order == PrefetchOrder::Input {
                    assert_eq!(item.index, pos);
                    assert_eq!(item.digest, ds[pos]);
                }
            }
        }

        // load_sequence goes through prefetch with cache writes on.
        let mut remote = open(cache.path(), &fx.url("s"), &cfg());
        let mut local = RefgetStore::open_local(&b.store).unwrap();
        for d in &ds {
            remote.load_sequence(d).unwrap();
            let r = remote.get_sequence(d).unwrap();
            let l = local.get_sequence(d).unwrap();
            assert_eq!(r.sequence(), l.sequence(), "{mode:?} {d}");
        }

        // Exports match byte for byte.
        let out = tempdir().unwrap();
        for (sel, label) in [
            (ExportSelection::All, "{seq}"),
            (ExportSelection::All, "{col}/{name}"),
        ] {
            let mut remote = open(tempdir().unwrap().path(), &fx.url("s"), &cfg());
            let mut local = RefgetStore::open_local(&b.store).unwrap();
            local.set_quiet(true);
            let r = export(&mut remote, &sel, label, &out.path().join("r.fa"));
            let l = export(&mut local, &sel, label, &out.path().join("l.fa"));
            assert_eq!(r, l, "{mode:?} {label}");
            assert!(r.starts_with('>'));
        }
    }
}

#[test]
fn cold_prefetch_makes_the_planned_requests_concurrently() {
    let b = build(StorageMode::Encoded, &[fasta_text()]);
    let fx = Fixture::serve(&b.store, &["s"]);
    let expected = local_bytes(&b.store);
    let ds = digests(&b.store);
    let planned = expected.values().filter(|v| !v.is_empty()).count();
    fx.add_fault(".seq", FaultKind::Delay(Duration::from_millis(40)), 1000);
    let cache = tempdir().unwrap();
    let mut remote = open(cache.path(), &fx.url("s"), &cfg());
    let counters = remote
        .cache_sequences(&ds, PrefetchOptions::default())
        .unwrap();
    assert_eq!(seq_requests(&fx), planned);
    assert!(
        fx.peak_in_flight() > 1 && fx.peak_in_flight() <= 4,
        "peak {}",
        fx.peak_in_flight()
    );
    assert_eq!(
        counters.sequences_fetched as usize,
        expected.values().filter(|v| !v.is_empty()).count()
    );
    assert_eq!(
        remote.fetch_stats().requests,
        remote.fetch_stats().requests.max(planned as u64)
    );
}

// ============================================================================
// Retries, 403, ranges
// ============================================================================

#[test]
fn prefetch_retries_then_succeeds() {
    let b = build(StorageMode::Encoded, &[fasta_text()]);
    let fx = Fixture::serve(&b.store, &["s"]);
    let ds = digests(&b.store);
    let cache = tempdir().unwrap();
    let mut remote = open(cache.path(), &fx.url("s"), &cfg());
    remote.load_sequence_index().unwrap();
    fx.add_fault(".seq", FaultKind::Status(503), 2);
    let run = remote.prefetch(&ds, PrefetchOptions::default()).unwrap();
    let mut run = run;
    for item in run.by_ref() {
        item.unwrap();
    }
    assert_eq!(run.counters().retries, 2);
}

#[test]
fn forbidden_means_missing() {
    let b = build(StorageMode::Encoded, &[fasta_text()]);
    // Like an S3 bucket that denies listing: a missing key answers 403.
    let fx = Fixture::with_options(&b.store, &["s"], true, true);
    let cache = tempdir().unwrap();
    let mut remote = open(cache.path(), &fx.url("s"), &cfg());
    // No FHR sidecars on the server: each counts as not_found.
    let r = remote.pull_fhr(None, SyncStrategy::KeepOurs).unwrap();
    assert_eq!(r.not_found, 1);
    // A missing .seq is a clear "not found or access denied (HTTP 403)".
    let ds = digests(&b.store);
    remote.load_sequence_index().unwrap();
    let victim = ds
        .iter()
        .find(|d| remote.sequence_file_path(d).is_some())
        .unwrap()
        .clone();
    let rel = remote.sequence_relative_path(&victim).unwrap();
    fs::remove_file(b.store.join(&rel)).unwrap();
    let err = remote.load_sequence(&victim).unwrap_err();
    assert!(
        format!("{err:#}").contains("not found or access denied (HTTP 403)"),
        "{err:#}"
    );
    assert_eq!(err.downcast_ref::<RemoteMissing>().unwrap().status, 403);

    // Missing rgstore.json.
    let err =
        RefgetStore::open_remote(cache.path(), format!("{}nope", fx.url("s")), &cfg()).unwrap_err();
    assert!(err.to_string().contains("rgstore.json: not found"), "{err}");
}

#[test]
fn server_ignoring_range_serves_whole_files_but_not_partial_reads() {
    let b = build(StorageMode::Raw, &[fasta_text()]);
    let fx = Fixture::with_options(&b.store, &["s"], false, false);
    let expected = local_bytes(&b.store);
    let ds = digests(&b.store);
    let cache = tempdir().unwrap();
    let mut remote = open(cache.path(), &fx.url("s"), &cfg());
    // Small files come down whole, with no Range header.
    for item in remote.prefetch(&ds, PrefetchOptions::default()).unwrap() {
        let item = item.unwrap();
        assert_eq!(item.bytes(), expected[&item.digest].as_slice());
    }
    let reqs = fx.requests_matching(".seq");
    assert!(
        !reqs.is_empty() && reqs.iter().all(|r| r.range.is_none()),
        "{reqs:?}"
    );
    // A partial read needs Range, so it is refused.
    let d = ds.iter().find(|d| expected[*d].len() > 10).unwrap();
    let err = remote.get_substring(d, 2, 5).unwrap_err().to_string();
    assert!(err.contains("ignored the Range header"), "{err}");
}

// ============================================================================
// Verification and the cache
// ============================================================================

#[test]
fn flipped_byte_fails_verification_and_caches_nothing() {
    let b = build(StorageMode::Encoded, &[fasta_text()]);
    let fx = Fixture::serve(&b.store, &["s"]);
    let cache = tempdir().unwrap();
    let mut remote = open(cache.path(), &fx.url("s"), &cfg());
    let d = digests(&b.store)
        .into_iter()
        .find(|d| local_bytes(&b.store)[d].len() > 4)
        .unwrap();
    fx.add_fault(".seq", FaultKind::FlipByte, 100);
    let err = format!(
        "{:#}",
        remote
            .cache_sequences(std::slice::from_ref(&d), PrefetchOptions::default())
            .unwrap_err()
    );
    assert!(err.contains(".seq") && err.contains("hash"), "{err}");
    assert_eq!(files_under(&cache.path().join("sequences")), 0);

    // verify: false delivers the (bad) bytes without hashing.
    let mut remote = open(
        cache.path(),
        &fx.url("s"),
        &FetchConfig {
            verify: false,
            ..cfg()
        },
    );
    let item = remote
        .prefetch(std::slice::from_ref(&d), PrefetchOptions::default())
        .unwrap()
        .next()
        .unwrap()
        .unwrap();
    assert_ne!(item.bytes(), local_bytes(&b.store)[&d].as_slice());
}

#[test]
fn cache_hits_and_second_origin_share_sequences() {
    let b = build(StorageMode::Encoded, &[fasta_text()]);
    let fx = Fixture::serve(&b.store, &["a", "b"]);
    let ds = digests(&b.store);
    let n_nonempty = local_bytes(&b.store)
        .values()
        .filter(|v| !v.is_empty())
        .count() as u64;
    let cache = tempdir().unwrap();

    let mut a = open(cache.path(), &fx.url("a"), &cfg());
    a.cache_sequences(&ds, PrefetchOptions::default()).unwrap();
    fx.clear_log();
    let second = a.cache_sequences(&ds, PrefetchOptions::default()).unwrap();
    assert_eq!(seq_requests(&fx), 0);
    assert_eq!(second.cache_hits, n_nonempty);

    let mut b2 = open(cache.path(), &fx.url("b"), &cfg());
    let counters = b2.cache_sequences(&ds, PrefetchOptions::default()).unwrap();
    assert_eq!(seq_requests(&fx), 0, "the second origin reuses every .seq");
    assert_eq!(counters.cache_hits, n_nonempty);
    let b_meta: Vec<_> = fx
        .requests()
        .into_iter()
        .filter(|r| r.prefix == "b")
        .collect();
    assert!(
        b_meta.iter().any(|r| r.path == "rgstore.json")
            && b_meta.iter().any(|r| r.path == "sequences.rgsi")
    );
}

#[test]
fn without_write_cache_nothing_is_written_and_bad_cache_files_are_replaced() {
    let b = build(StorageMode::Encoded, &[fasta_text()]);
    let fx = Fixture::serve(&b.store, &["s"]);
    let ds = digests(&b.store);
    let expected = local_bytes(&b.store);
    let cache = tempdir().unwrap();
    let mut remote = open(cache.path(), &fx.url("s"), &cfg());
    for item in remote.prefetch(&ds, PrefetchOptions::default()).unwrap() {
        item.unwrap();
    }
    let out = cache.path().join("o.fa");
    export(&mut remote, &ExportSelection::All, "{seq}", &out);
    assert_eq!(files_under(&cache.path().join("sequences")), 0);
    assert_eq!(files_under(&cache.path().join("tmp")), 0);

    // A wrong-length cache file is deleted and fetched again.
    let d = ds.iter().find(|d| expected[*d].len() > 3).unwrap();
    let path = remote.sequence_file_path(d).unwrap();
    fs::create_dir_all(path.parent().unwrap()).unwrap();
    fs::write(&path, b"xx").unwrap();
    remote
        .cache_sequences(std::slice::from_ref(d), PrefetchOptions::default())
        .unwrap();
    assert_eq!(fs::read(&path).unwrap(), expected[d]);
}

// ============================================================================
// Resume, repeats, two writers
// ============================================================================

fn big_store() -> (Built, String, Vec<u8>) {
    let b = build(
        StorageMode::Raw,
        &[format!(">big\n{}\n", seq(b"ACGT", 23, 9))],
    );
    let d = digests(&b.store).remove(0);
    let bytes = local_bytes(&b.store)[&d].clone();
    (b, d, bytes)
}

#[test]
fn interrupted_cache_fill_resumes_with_only_missing_sequences() {
    let b = build(StorageMode::Encoded, &[fasta_text()]);
    let fx = Fixture::serve(&b.store, &["s"]);
    let ds = digests(&b.store);
    let expected = local_bytes(&b.store);
    let n_nonempty = expected.values().filter(|v| !v.is_empty()).count() as u64;
    // One request at a time, largest first: the smallest sequence goes last.
    let one = FetchConfig {
        jobs: 1,
        retries: 0,
        ..cfg()
    };
    let victim = ds
        .iter()
        .enumerate()
        .filter(|(_, d)| !expected[*d].is_empty())
        .min_by_key(|(i, d)| (expected[*d].len(), std::cmp::Reverse(*i)))
        .unwrap()
        .1
        .clone();
    let cache = tempdir().unwrap();
    let mut remote = open(cache.path(), &fx.url("s"), &one);
    remote.load_sequence_index().unwrap();
    let rel = remote.sequence_relative_path(&victim).unwrap();
    let rel = rel.to_string_lossy();
    fx.add_fault(&rel, FaultKind::Status(500), 1000);
    assert!(
        remote
            .cache_sequences(&ds, PrefetchOptions::default())
            .is_err()
    );

    fx.clear_faults();
    fx.clear_log();
    let mut remote = open(cache.path(), &fx.url("s"), &one);
    let counters = remote
        .cache_sequences(&ds, PrefetchOptions::default())
        .unwrap();
    // The rerun downloads only the sequence the first run did not finish.
    assert_eq!(counters.cache_hits, n_nonempty - 1, "{counters:?}");
    assert_eq!(counters.sequences_fetched, 1);
    let reqs = fx.requests_matching(".seq");
    assert!(reqs.len() == 1 && reqs[0].path.contains(&*rel), "{reqs:?}");
    assert_eq!(
        fs::read(remote.sequence_file_path(&victim).unwrap()).unwrap(),
        expected[&victim]
    );
}

#[test]
fn input_order_delivers_every_repeat() {
    let (b, d, bytes) = big_store();
    let fx = Fixture::serve(&b.store, &["s"]);
    let cache = tempdir().unwrap();
    let mut remote = open(cache.path(), &fx.url("s"), &cfg());
    let items: Vec<_> = remote
        .prefetch(
            &[d.clone(), d.clone()],
            PrefetchOptions {
                order: PrefetchOrder::Input,
                ..Default::default()
            },
        )
        .unwrap()
        .map(|i| i.unwrap())
        .collect();
    assert_eq!(
        items.iter().map(|i| i.index).collect::<Vec<_>>(),
        vec![0, 1]
    );
    assert!(items.iter().all(|i| i.bytes() == bytes.as_slice()));
}

#[test]
fn two_writers_fill_one_cache() {
    let (b, d, bytes) = big_store();
    let fx = Fixture::serve(&b.store, &["s"]);
    let cache = tempdir().unwrap();
    std::thread::scope(|s| {
        for _ in 0..2 {
            s.spawn(|| {
                let mut remote = open(cache.path(), &fx.url("s"), &cfg());
                remote
                    .cache_sequences(std::slice::from_ref(&d), PrefetchOptions::default())
                    .unwrap();
            });
        }
    });
    let remote = open(cache.path(), &fx.url("s"), &cfg());
    let _ = remote;
    let path = walk(&cache.path().join("sequences")).remove(0);
    assert_eq!(fs::read(path).unwrap(), bytes);
}

// ============================================================================
// Partial reads
// ============================================================================

#[test]
fn partial_reads_fetch_ranges_and_persist_nothing() {
    let b = build(
        StorageMode::Encoded,
        &[format!(">chr\n{}\n", seq(b"ACGT", 400, 4))],
    );
    let fx = Fixture::serve(&b.store, &["s"]);
    let d = digests(&b.store).remove(0);
    let local = RefgetStore::open_local(&b.store).unwrap();
    let cache = tempdir().unwrap();
    let mut remote = open(cache.path(), &fx.url("s"), &cfg());
    remote.load_sequence_index().unwrap();
    assert_eq!(
        remote.get_substring(&d, 10, 50).unwrap(),
        local.get_substring(&d, 10, 50).unwrap()
    );
    let mut streamed = String::new();
    remote
        .stream_sequence(&d, Some(5), Some(300))
        .unwrap()
        .read_to_string(&mut streamed)
        .unwrap();
    assert_eq!(streamed, local.get_substring(&d, 5, 300).unwrap());
    assert_eq!(files_under(&cache.path().join("sequences")), 0);

    // Many ranges at once promote the whole sequence into the cache.
    let ranges: Vec<(usize, usize)> = (0..20).map(|i| (i * 10, i * 10 + 7)).collect();
    let got = remote.get_substrings(&d, &ranges).unwrap();
    let want = local.get_substrings(&d, &ranges).unwrap();
    assert_eq!(got, want);
    assert_eq!(files_under(&cache.path().join("sequences")), 1);
    let items: Vec<_> = remote
        .prefetch(std::slice::from_ref(&d), PrefetchOptions::default())
        .unwrap()
        .map(|i| i.unwrap().source)
        .collect();
    assert_eq!(items, vec![ItemSource::Cache]);
}

// ============================================================================
// Metadata
// ============================================================================

fn many_collections(n: usize, fhr: bool) -> Built {
    let fastas: Vec<String> = (0..n)
        .map(|i| format!(">c{i}\n{}\n", seq(b"ACGT", 12 + i, i)))
        .collect();
    let b = build(StorageMode::Encoded, &fastas);
    if fhr {
        let mut s = RefgetStore::open_local(&b.store).unwrap();
        for d in &b.collections {
            s.set_fhr_metadata(
                d,
                FhrMetadata {
                    genome: Some(format!("g-{d}")),
                    ..Default::default()
                },
            )
            .unwrap();
        }
        s.write().unwrap();
    }
    b
}

#[test]
fn collection_and_fhr_metadata_download_concurrently() {
    let b = many_collections(20, true);
    let fx = Fixture::serve(&b.store, &["s"]);
    let cache = tempdir().unwrap();
    let mut remote = open(cache.path(), &fx.url("s"), &cfg());
    fx.add_fault(
        "collections/",
        FaultKind::Delay(Duration::from_millis(30)),
        1000,
    );
    fx.add_fault("fhr/", FaultKind::Delay(Duration::from_millis(30)), 1000);
    fx.clear_log();
    remote.load_all_collections().unwrap();
    assert_eq!(fx.requests_matching("collections/").len(), 20);
    assert!(fx.peak_in_flight() > 1, "peak {}", fx.peak_in_flight());
    fx.clear_log();
    let r = remote.pull_fhr(None, SyncStrategy::KeepOurs).unwrap();
    assert_eq!(r.pulled, 20);
    assert!(fx.peak_in_flight() > 1, "peak {}", fx.peak_in_flight());
    // Already cached: KeepOurs still loads each sidecar into memory.
    let mut again = open(cache.path(), &fx.url("s"), &cfg());
    let r = again.pull_fhr(None, SyncStrategy::KeepOurs).unwrap();
    assert_eq!(r.skipped, 20);
    assert!(again.get_fhr_metadata(&b.collections[0]).is_some());
    // load_collections batches a selection.
    let mut sel = open(tempdir().unwrap().path(), &fx.url("s"), &cfg());
    sel.load_collections(&b.collections[..5]).unwrap();
    assert!(
        b.collections[..5]
            .iter()
            .all(|d| sel.is_collection_loaded(d))
    );
}

#[test]
fn sequence_index_mismatch_is_caught() {
    let b = many_collections(3, false);
    let fx = Fixture::serve(&b.store, &["s"]);
    let cache = tempdir().unwrap();
    let mut remote = open(cache.path(), &fx.url("s"), &cfg());
    fx.clear_log();
    remote.load_sequence_index().unwrap();
    assert_eq!(fx.requests_matching("sequences.rgsi").len(), 1);
    assert_eq!(remote.list_sequences().unwrap().len(), 3);

    // The server's index no longer matches its manifest: one refetch, then an error.
    let idx = b.store.join("sequences.rgsi");
    let mut text = fs::read_to_string(&idx).unwrap();
    text.push_str("# changed\n");
    fs::write(&idx, text).unwrap();
    let cache2 = tempdir().unwrap();
    let mut remote = open(cache2.path(), &fx.url("s"), &cfg());
    fx.clear_log();
    let err = remote.load_sequence_index().unwrap_err().to_string();
    assert!(err.contains("does not match rgstore.json"), "{err}");
    assert_eq!(fx.requests_matching("sequences.rgsi").len(), 2);
}

#[test]
fn reopen_after_server_change_refreshes_metadata_but_not_sequences() {
    let b = many_collections(2, false);
    let fx = Fixture::serve(&b.store, &["s"]);
    let cache = tempdir().unwrap();
    let ds = digests(&b.store);
    {
        let mut remote = open(cache.path(), &fx.url("s"), &cfg());
        remote.load_all_collections().unwrap();
        remote
            .cache_sequences(&ds, PrefetchOptions::default())
            .unwrap();
    }
    // Add a collection on the server.
    {
        let mut w = RefgetStore::open_local(&b.store).unwrap();
        let fa = b.store.parent().unwrap().join("more.fa");
        fs::write(&fa, format!(">extra\n{}\n", seq(b"ACGT", 50, 77))).unwrap();
        w.add_sequence_collection_from_fasta(&fa, FastaImportOptions::new())
            .unwrap();
        w.write().unwrap();
    }
    fx.clear_log();
    let mut remote = open(cache.path(), &fx.url("s"), &cfg());
    remote.load_all_collections().unwrap();
    assert!(!fx.requests_matching("collections.rgci").is_empty());
    let mut paths: Vec<String> = fx
        .requests_matching("collections/")
        .into_iter()
        .map(|r| r.path)
        .collect();
    paths.sort();
    paths.dedup();
    assert_eq!(paths.len(), 3);
    remote
        .cache_sequences(&ds, PrefetchOptions::default())
        .unwrap();
    assert_eq!(seq_requests(&fx), 0);
}

#[test]
fn offline_reopen_uses_the_cache() {
    let b = many_collections(1, false);
    let fx = Fixture::serve(&b.store, &["s"]);
    let url = fx.url("s");
    let cache = tempdir().unwrap();
    let ds = digests(&b.store);
    {
        let mut remote = open(cache.path(), &url, &cfg());
        remote
            .cache_sequences(&ds, PrefetchOptions::default())
            .unwrap();
    }
    drop(fx);
    let mut remote = open(cache.path(), &url, &cfg());
    let items: Vec<_> = remote
        .prefetch(&ds, PrefetchOptions::default())
        .unwrap()
        .map(|i| i.unwrap().source)
        .collect();
    assert_eq!(items, vec![ItemSource::Cache]);
}

//! Multi-PROCESS concurrent-writer tests for RefgetStore.
//!
//! These deliberately spawn real OS processes rather than threads. The whole
//! point of the store lock is coordination BETWEEN processes — SLURM jobs,
//! separate `gtars refget build` invocations, Python callers — and a
//! thread-based test would exercise none of it: an in-process `Mutex` would
//! serialize the writers whether or not the file lock works, and the `/proc`
//! liveness check and steal-by-rename paths would never be reached.
//!
//! Each child re-executes THIS test binary with a marker env var set, which is
//! the standard way to get a real subprocess out of a cargo test without
//! shipping a separate helper binary. The binary runs with `harness = false`
//! (see Cargo.toml) so [`main`] can hand off to the child role immediately,
//! before any test-runner machinery starts.

use std::path::{Path, PathBuf};
use std::process::Command;

use gtars_refget::store::{FastaImportOptions, RefgetStore, lock_status};

/// Env var that turns a re-exec of this binary into a writer child.
const CHILD_STORE: &str = "GTARS_TEST_CHILD_STORE";
/// Env var naming the collection the child should add.
const CHILD_NAME: &str = "GTARS_TEST_CHILD_NAME";
/// Env var telling the child to die mid-commit instead of finishing.
const CHILD_ABORT_HOLDING_LOCK: &str = "GTARS_TEST_CHILD_ABORT_HOLDING_LOCK";
/// Env var turning a re-exec into a READER that loops until told to stop.
const CHILD_READ_UNTIL: &str = "GTARS_TEST_CHILD_READ_UNTIL";
/// Env var turning a re-exec into a REMOVER of one collection (with orphans).
const CHILD_REMOVE: &str = "GTARS_TEST_CHILD_REMOVE";

/// A distinct 32-base sequence per child, so every writer contributes rows no
/// other writer could have produced.
fn fasta_for(name: &str) -> String {
    let bases = ["ACGT", "AACC", "GGTT", "CCGG", "TTAA", "GCGC", "ATAT", "TGCA"];
    let idx = name.bytes().map(|b| b as usize).sum::<usize>() % bases.len();
    format!(">chr_{}\n{}\n", name, bases[idx].repeat(8))
}

/// Runs in the CHILD process: open the shared store, add one collection, commit.
fn run_as_child(store_dir: &str, name: &str) -> ! {
    let work = std::env::temp_dir().join(format!("gtars-child-{}-{}", name, std::process::id()));
    std::fs::create_dir_all(&work).unwrap();
    let fasta = work.join(format!("{}.fa", name));
    std::fs::write(&fasta, fasta_for(name)).unwrap();

    let mut store = RefgetStore::on_disk(store_dir).expect("child could not open store");
    store.set_quiet(true);
    store
        .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
        .expect("child could not import fasta");

    if std::env::var(CHILD_ABORT_HOLDING_LOCK).is_ok() {
        // Take the lock and die without releasing it, exactly as a SIGKILLed or
        // OOM-killed SLURM job would. `process::abort` skips every destructor,
        // so the lockfile is left on disk.
        store.lock_for_batch("deliberately-abandoned").unwrap();
        std::process::abort();
    }

    store.write().expect("child could not commit");
    std::process::exit(0);
}

/// Spawn a writer child against `store_dir`.
fn spawn_writer(store_dir: &Path, name: &str, abort: bool) -> std::process::Child {
    let exe = std::env::current_exe().expect("test binary path");
    let mut cmd = Command::new(exe);
    cmd.env(CHILD_STORE, store_dir).env(CHILD_NAME, name);
    if abort {
        cmd.env(CHILD_ABORT_HOLDING_LOCK, "1");
    }
    cmd.spawn().expect("failed to spawn writer child")
}

/// Runs in the CHILD process: hammer `open_local` until the sentinel file
/// appears, asserting the store is never observed in a torn state.
///
/// Exits non-zero on the FIRST bad read, so the parent's `status.success()`
/// check is a real assertion about every read the child performed.
fn run_as_reader(store_dir: &str, sentinel: &str, min_collections: usize) -> ! {
    let store_dir = Path::new(store_dir);
    let sentinel = Path::new(sentinel);
    let mut reads = 0usize;

    while !sentinel.exists() {
        reads += 1;
        // A commit publishes by rename, so every open must see a complete,
        // parseable store -- never a truncated index or a manifest describing
        // bytes that have not landed.
        let store = match RefgetStore::open_local(store_dir) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("reader saw a torn store after {} reads: {}", reads, e);
                std::process::exit(2);
            }
        };
        // The index must never SHRINK below what was already published: a
        // partially-written file that still happens to parse would show up here.
        let seen = store.stats().n_collections;
        if seen < min_collections {
            eprintln!(
                "reader saw {} collections, fewer than the {} already published",
                seen, min_collections
            );
            std::process::exit(3);
        }
    }

    if reads == 0 {
        eprintln!("reader never got a chance to read");
        std::process::exit(4);
    }
    println!("reader completed {} clean reads", reads);
    std::process::exit(0);
}

/// Runs in the CHILD process: remove one collection (and its orphan sequences)
/// from the shared store, then exit. Used to make the removal a genuinely
/// separate process from the writer holding the stale snapshot.
fn run_as_remover(store_dir: &str, digest: &str) -> ! {
    let mut store = RefgetStore::open_local(store_dir).expect("remover could not open store");
    store.set_quiet(true);
    let removed = store
        .remove_collection(digest, true)
        .expect("remover could not remove collection");
    assert!(removed, "the collection to remove was not in the store");
    drop(store);
    std::process::exit(0);
}

fn main() {
    // Child roles first: a re-exec must never fall through into the tests.
    if let (Ok(store), Ok(digest)) = (std::env::var(CHILD_STORE), std::env::var(CHILD_REMOVE)) {
        run_as_remover(&store, &digest);
    }
    if let (Ok(store), Ok(sentinel)) = (std::env::var(CHILD_STORE), std::env::var(CHILD_READ_UNTIL))
    {
        let min: usize = std::env::var(CHILD_NAME)
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0);
        run_as_reader(&store, &sentinel, min);
    }
    if let (Ok(store), Ok(name)) = (std::env::var(CHILD_STORE), std::env::var(CHILD_NAME)) {
        run_as_child(&store, &name);
    }

    let tests: &[(&str, fn())] = &[
        (
            "concurrent_processes_all_land_in_the_index",
            concurrent_processes_all_land_in_the_index,
        ),
        // Tier-1 (same host, dead pid) staleness needs /proc; elsewhere the
        // lock only expires via the 120s heartbeat threshold, which this test
        // deliberately does not wait for.
        #[cfg(target_os = "linux")]
        (
            "abandoned_lock_from_a_dead_process_is_broken",
            abandoned_lock_from_a_dead_process_is_broken,
        ),
        (
            "reader_never_observes_a_torn_store",
            reader_never_observes_a_torn_store,
        ),
        (
            "another_process_removal_is_not_resurrected",
            another_process_removal_is_not_resurrected,
        ),
    ];

    println!("\nrunning {} tests", tests.len());
    for (name, f) in tests {
        f();
        println!("test {} ... ok", name);
    }
    println!("\ntest result: ok. {} passed; 0 failed\n", tests.len());
}

fn seed_store(store_dir: &Path, work: &Path) -> String {
    let fasta = work.join("seed.fa");
    std::fs::write(&fasta, ">chr_seed\nAAAACCCCGGGGTTTT\n").unwrap();
    let mut store = RefgetStore::on_disk(store_dir).unwrap();
    store.set_quiet(true);
    let (meta, _) = store
        .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
        .unwrap();
    store.write().unwrap();
    meta.digest
}

fn collection_digests_on_disk(store_dir: &Path) -> Vec<String> {
    let rgci = std::fs::read_to_string(store_dir.join("collections.rgci")).unwrap();
    rgci.lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .filter_map(|l| l.split('\t').next().map(String::from))
        .collect()
}

/// N real processes build into ONE store directory simultaneously. Every
/// collection must appear in the final `collections.rgci`.
///
/// This is the 2026-07-23 incident as an executable test: four concurrent
/// `genome_init` jobs, one of which lost the race and was dropped from both
/// indexes without any writer reporting an error.
fn concurrent_processes_all_land_in_the_index() {
    let work = tempfile::tempdir().unwrap();
    let store_dir = work.path().join("store");
    let seed = seed_store(&store_dir, work.path());

    const WRITERS: usize = 4;
    let children: Vec<_> = (0..WRITERS)
        .map(|i| spawn_writer(&store_dir, &format!("w{}", i), false))
        .collect();

    for mut child in children {
        let status = child.wait().expect("child did not run");
        assert!(status.success(), "a writer process failed: {:?}", status);
    }

    let digests = collection_digests_on_disk(&store_dir);
    assert_eq!(
        digests.len(),
        WRITERS + 1,
        "expected the seed plus {} writers in collections.rgci, found {}: {:?}",
        WRITERS,
        digests.len(),
        digests
    );
    assert!(digests.contains(&seed), "the seed collection was dropped");

    // The store must still open, and the sequence index must carry one row per
    // writer's unique chromosome.
    let reopened = RefgetStore::open_local(&store_dir).unwrap();
    let names: std::collections::HashSet<String> = reopened
        .list_sequences()
        .iter()
        .map(|m| m.name.clone())
        .collect();
    for i in 0..WRITERS {
        let expected = format!("chr_w{}", i);
        assert!(
            names.contains(&expected),
            "sequence {} is missing from sequences.rgsi; have {:?}",
            expected,
            names
        );
    }

    // No lock or temp file may outlive the writers.
    assert!(lock_status(&store_dir).unwrap().is_none(), "lock left behind");
    assert!(
        transient_files(&store_dir).is_empty(),
        "transient files left behind: {:?}",
        transient_files(&store_dir)
    );
}

/// A writer killed while holding the lock must not wedge the store: the next
/// writer detects the dead pid on this host and breaks the lock immediately.
#[cfg(target_os = "linux")]
fn abandoned_lock_from_a_dead_process_is_broken() {
    let work = tempfile::tempdir().unwrap();
    let store_dir = work.path().join("store");
    let seed = seed_store(&store_dir, work.path());

    let mut killer = spawn_writer(&store_dir, "doomed", true);
    let status = killer.wait().expect("child did not run");
    assert!(!status.success(), "the child was supposed to abort");

    // The abandoned lockfile is still there, naming a pid that no longer exists.
    let stranded = lock_status(&store_dir)
        .unwrap()
        .expect("the aborted child should have left its lock behind");
    assert_eq!(stranded.operation, "deliberately-abandoned");

    // Tier-1 staleness (same host, dead pid) breaks it without waiting out the
    // 120s heartbeat threshold.
    let start = std::time::Instant::now();
    let mut store = RefgetStore::open_local(&store_dir).unwrap();
    let fasta = work.path().join("after.fa");
    std::fs::write(&fasta, ">chr_after\nGGGGTTTTAAAACCCC\n").unwrap();
    store
        .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
        .unwrap();
    store.write().unwrap();
    assert!(
        start.elapsed() < std::time::Duration::from_secs(30),
        "breaking a definitively-dead lock should be immediate, took {:?}",
        start.elapsed()
    );
    drop(store);

    // The store is intact and the stale lock is gone.
    let digests = collection_digests_on_disk(&store_dir);
    assert!(digests.contains(&seed), "the store lost data across the abort");
    assert!(lock_status(&store_dir).unwrap().is_none());
    let reopened = RefgetStore::open_local(&store_dir).unwrap();
    assert!(
        reopened
            .list_sequences()
            .iter()
            .any(|m| m.name == "chr_after"),
        "the recovering writer's commit did not land"
    );
}

fn transient_files(dir: &Path) -> Vec<PathBuf> {
    std::fs::read_dir(dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| {
            e.file_name()
                .to_str()
                .is_some_and(gtars_refget::store::is_transient_store_file)
        })
        .map(|e| e.path())
        .collect()
}

/// A reader process loops `open_local` while a writer commits repeatedly over a
/// non-trivial index. Because every whole-file write is published by
/// `rename(2)`, the reader must never observe a truncated index, a half-written
/// manifest, or a store that fails to parse.
///
/// Readers take NO lock — that is the design, not an oversight: reads must never
/// block on a shared filesystem. Atomic publish is what makes it safe.
fn reader_never_observes_a_torn_store() {
    let work = tempfile::tempdir().unwrap();
    let store_dir = work.path().join("store");
    let sentinel = work.path().join("STOP");

    // Seed with enough collections that an index write is not a single small
    // buffer -- a torn read needs a file big enough to catch mid-write.
    const SEEDED: usize = 40;
    let mut store = RefgetStore::on_disk(&store_dir).unwrap();
    store.set_quiet(true);
    for i in 0..SEEDED {
        let fasta = work.path().join(format!("seed{}.fa", i));
        let mut content = String::new();
        for j in 0..25 {
            content.push_str(&format!(">chr{}_{}\n{}\n", i, j, "ACGT".repeat(64)));
        }
        std::fs::write(&fasta, content).unwrap();
        store
            .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
            .unwrap();
    }
    store.write().unwrap();

    // Start the reader, then keep committing underneath it.
    let exe = std::env::current_exe().expect("test binary path");
    let mut reader = Command::new(exe)
        .env(CHILD_STORE, &store_dir)
        .env(CHILD_READ_UNTIL, &sentinel)
        .env(CHILD_NAME, SEEDED.to_string())
        .spawn()
        .expect("failed to spawn reader child");

    for i in 0..25 {
        let fasta = work.path().join(format!("live{}.fa", i));
        std::fs::write(&fasta, format!(">chr_live_{}\n{}\n", i, "GGTT".repeat(64))).unwrap();
        store
            .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
            .unwrap();
        store.write().unwrap();
    }

    std::fs::write(&sentinel, b"stop").unwrap();
    let status = reader.wait().expect("reader child did not run");
    assert!(
        status.success(),
        "the reader observed a torn store (exit {:?})",
        status.code()
    );
}

/// A removal performed by a DIFFERENT process must survive the next commit of a
/// writer that opened before it.
///
/// This is the 2026-07-28 finding as an executable test, and the ordering is the
/// whole point: the remover's commit completes before the writer's begins, so
/// the lock never even contends. Merge-at-commit still lost it, because the
/// writer wrote back the collection stub `open_local` had loaded for reading —
/// re-adding index rows for `.seq` files the remover had already unlinked.
fn another_process_removal_is_not_resurrected() {
    let work = tempfile::tempdir().unwrap();
    let store_dir = work.path().join("store");

    // Seed: a collection to keep, and a victim with its own sequences + alias.
    let keep_fasta = work.path().join("keep.fa");
    std::fs::write(&keep_fasta, ">chr_keep\nAAAACCCCGGGGTTTT\n").unwrap();
    let victim_fasta = work.path().join("victim.fa");
    std::fs::write(&victim_fasta, ">chr_v1\nGGGGTTTTAAAACCCC\n>chr_v2\nTTTTAAAACCCCGGGG\n").unwrap();

    let mut store = RefgetStore::on_disk(&store_dir).unwrap();
    store.set_quiet(true);
    let (keep, _) = store
        .add_sequence_collection_from_fasta(&keep_fasta, FastaImportOptions::new())
        .unwrap();
    let (victim, _) = store
        .add_sequence_collection_from_fasta(&victim_fasta, FastaImportOptions::new())
        .unwrap();
    store
        .add_collection_alias("refgenie", "victim", &victim.digest)
        .unwrap();
    store.write().unwrap();
    store.load_collection(&victim.digest).unwrap();
    let victim_sequences: Vec<String> = store
        .get_collection(&victim.digest)
        .unwrap()
        .sequences
        .iter()
        .map(|s| s.metadata().sha512t24u.clone())
        .collect();
    drop(store);

    // The writer opens NOW, while the victim is still in the index.
    let mut writer = RefgetStore::open_local(&store_dir).unwrap();
    writer.set_quiet(true);
    assert!(writer.get_collection_metadata(&victim.digest).is_some());

    // A separate process removes the victim and exits. Fully sequenced: the
    // removal is finished and committed before the writer touches anything.
    let exe = std::env::current_exe().expect("test binary path");
    let status = Command::new(exe)
        .env(CHILD_STORE, &store_dir)
        .env(CHILD_REMOVE, &victim.digest)
        .status()
        .expect("failed to spawn remover child");
    assert!(status.success(), "the remover process failed: {:?}", status);
    assert!(
        !collection_digests_on_disk(&store_dir).contains(&victim.digest),
        "precondition: the remover's commit must have landed"
    );

    // Only now does the writer commit its own, unrelated addition.
    let added_fasta = work.path().join("added.fa");
    std::fs::write(&added_fasta, ">chr_added\nACACACACGTGTGTGT\n").unwrap();
    let (added, _) = writer
        .add_sequence_collection_from_fasta(&added_fasta, FastaImportOptions::new())
        .unwrap();
    writer.write().unwrap();
    drop(writer);

    let digests = collection_digests_on_disk(&store_dir);
    assert!(
        !digests.contains(&victim.digest),
        "the other process's removal was undone: {} is back in collections.rgci",
        victim.digest
    );
    assert!(digests.contains(&keep.digest), "the untouched collection was dropped");
    assert!(digests.contains(&added.digest), "the writer's own addition was dropped");

    // Its sequence rows must be gone too -- a resurrected row points at a `.seq`
    // file the remover unlinked, which is corruption, not merely stale metadata.
    let rgsi = std::fs::read_to_string(store_dir.join("sequences.rgsi")).unwrap();
    for seq in &victim_sequences {
        assert!(
            !rgsi.contains(seq.as_str()),
            "orphan sequence {} came back into sequences.rgsi",
            seq
        );
    }

    // And the alias must not resolve to a collection that no longer exists.
    let reopened = RefgetStore::open_local(&store_dir).unwrap();
    assert!(
        reopened
            .get_collection_metadata_by_alias("refgenie", "victim")
            .is_none(),
        "the removed collection's alias came back"
    );

    assert!(lock_status(&store_dir).unwrap().is_none(), "lock left behind");
}

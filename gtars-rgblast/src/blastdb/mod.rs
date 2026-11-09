//! Writing (and reading back) BLAST version 5 nucleotide databases from a `SubjectPlan`.

pub mod alias;
pub mod asn;
pub mod date;
pub mod lmdb;
pub mod njs;
pub mod oidfiles;
pub mod output;
pub mod preflight;
pub mod reader;
pub mod sidecars;
pub mod verify;
pub mod writer;

use crate::convert::{BaseCounts, ConvertedBatch};
use crate::error::Error;
use crate::pipeline::run_ordered_with;
use crate::plan::{SkipReason, Subject, SubjectPlan};
use crate::progress::{ProgressEvent, Throttle};
use crate::select::nstats::{NCount, NstatsOptions};
use crate::store::decode;
use crate::store::fetch::{HeldBytes, SeqFetch};
pub use date::BlastDate;
use gtars_refget::store::{FetchCounters, StorageMode};
use output::{OutputSet, PartialOutput, with_ext};
use preflight::{PreflightInput, PreflightSummary, Prepared};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use std::sync::Mutex;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

/// How to write a database (the writer's part of `BuildOptions`, resolved).
#[derive(Debug, Clone)]
pub struct WriteOptions {
    pub out: PathBuf,
    pub title: String,
    pub max_file_size: u64,
    /// Converted bytes held in memory by in-flight batches (fetched bytes have their own
    /// budget, the store's `RemoteOptions::memory_bytes`).
    pub budget_bytes: u64,
    /// Conversion workers.
    pub threads: usize,
    /// Stored bytes per batch of consecutive OIDs (`--batch-mb`).
    pub batch_bytes: u64,
    pub force: bool,
    pub check_free_space: bool,
    pub date: BlastDate,
    pub progress_interval_ms: u64,
    /// For the preflight summary: where the store is, its label, local or remote.
    pub store: String,
    pub store_label: String,
    pub remote: bool,
    /// The selection line of the preflight summary (empty when unknown).
    pub selection: String,
    /// Where per-sequence N counts (a by-product of conversion) are cached.
    pub nstats: NstatsOptions,
    /// The manifest to write, with `volumes`, `files`, `totals` and `date` filled in by the
    /// writer; None writes no manifest (`view` then refuses the database).
    pub manifest: Option<sidecars::Manifest>,
}

#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
pub struct Timings {
    pub plan_s: f64,
    pub sequences_s: f64,
    pub finish_s: f64,
    /// Summed over workers.
    pub cpu_read_s: f64,
    /// Summed over workers.
    pub cpu_convert_s: f64,
    /// On the writing thread.
    pub write_s: f64,
}

/// What a remote build downloaded (gtars-refget's `FetchCounters` of the pull).
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct NetworkReport {
    /// Body bytes received, failed attempts included.
    pub bytes: u64,
    /// HTTP requests, retries included.
    pub requests: u64,
    pub retries: u64,
    /// Whole sequences that came over the network.
    pub sequences_fetched: u64,
    /// Sequences read from the shared cache instead.
    pub cache_hits: u64,
    pub cache_hit_bytes: u64,
    /// Wall time of the fetch.
    pub elapsed_s: f64,
    /// `bytes` over `elapsed_s`, in MB/s (10^6 bytes).
    pub mb_per_s: f64,
}

impl From<&FetchCounters> for NetworkReport {
    fn from(c: &FetchCounters) -> Self {
        NetworkReport {
            bytes: c.bytes_downloaded,
            requests: c.requests,
            retries: c.retries,
            sequences_fetched: c.sequences_fetched,
            cache_hits: c.cache_hits,
            cache_hit_bytes: c.cache_hit_bytes,
            elapsed_s: c.elapsed_secs,
            mb_per_s: c.rate_mb_s(),
        }
    }
}

impl NetworkReport {
    /// gtars-refget's summary line, e.g. `network: 1362000000 bytes in 346 requests (0
    /// retries), 27.6 MB/s; cache hits 0 (0.00 GB)`.
    pub fn line(&self) -> String {
        FetchCounters {
            requests: self.requests,
            retries: self.retries,
            bytes_downloaded: self.bytes,
            sequences_fetched: self.sequences_fetched,
            cache_hits: self.cache_hits,
            cache_hit_bytes: self.cache_hit_bytes,
            elapsed_secs: self.elapsed_s,
            ..FetchCounters::default()
        }
        .network_line()
    }
}

/// One taxon of a database: its taxid, species name and number of collections.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct TaxonCount {
    pub taxid: i32,
    pub species: Option<String>,
    pub collections: u64,
}

/// The taxa of a plan's collections (taxid > 0), ascending by taxid.
pub fn taxa_of(plan: &SubjectPlan) -> Vec<TaxonCount> {
    let mut m: std::collections::BTreeMap<i32, TaxonCount> = std::collections::BTreeMap::new();
    for c in &plan.collections {
        if let Some(t) = c.taxid.filter(|t| *t > 0) {
            let e = m.entry(t).or_insert_with(|| TaxonCount {
                taxid: t,
                species: c.species.clone(),
                collections: 0,
            });
            e.collections += 1;
            if e.species.is_none() {
                e.species = c.species.clone();
            }
        }
    }
    m.into_values().collect()
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BuildSummary {
    pub preflight: PreflightSummary,
    /// Taxa the deflines carry (taxid > 0): taxid, species, number of collections.
    pub taxa: Vec<TaxonCount>,
    /// Database files, sorted (what `.njs` lists, plus `.njs`).
    pub files: Vec<PathBuf>,
    pub sidecars: Vec<PathBuf>,
    pub bytes_written: u64,
    pub ambiguous_bases: u64,
    pub volumes: u64,
    pub warnings: Vec<String>,
    pub timings: Timings,
    pub network: Option<NetworkReport>,
    pub elapsed_s: f64,
}

/// Encode headers and compute the preflight summary (writes nothing).
pub fn plan_preflight(plan: &SubjectPlan, opts: &WriteOptions) -> (Prepared, PreflightSummary) {
    let prep = preflight::prepare(plan);
    let date = opts.date.blast();
    let s = preflight::compute(
        plan,
        &prep,
        &PreflightInput {
            store: opts.store.clone(),
            store_label: opts.store_label.clone(),
            remote: opts.remote,
            out: &opts.out,
            title: &opts.title,
            date: &date,
            max_file_size: opts.max_file_size,
        },
    );
    let mut s = s;
    s.selection = opts.selection.clone();
    (prep, s)
}

/// The warnings a plan produces, with their sidecar paths.
pub fn warnings(plan: &SubjectPlan, base: &std::path::Path) -> Vec<String> {
    use crate::progress::commas;
    let name = |ext: &str| crate::blastdb::writer::file_name(&with_ext(base, ext));
    let mut w = Vec::new();
    for (reason, what) in [
        (
            SkipReason::Oversize,
            "sequences of 2^31 bp or more, which BLAST can't store",
        ),
        (SkipReason::NonNucleotide, "non-nucleotide sequences"),
        (SkipReason::Empty, "empty sequences"),
    ] {
        let n = plan.n_skipped(reason);
        if n > 0 {
            w.push(format!(
                "skipped {} {what} (listed in {})",
                commas(n),
                name("skipped.tsv")
            ));
        }
    }
    if !plan.sanitized.is_empty() {
        w.push(format!(
            "replaced characters BLAST ids can't hold in {} names or aliases (listed in {})",
            commas(plan.sanitized.len() as u64),
            name("sanitized.tsv")
        ));
    }
    let n = plan.n_dash_leading_collections();
    if n > 0 {
        w.push(format!(
            "{} collection digests start with '-'; look their sequences up as gnl|<digest>|<name>",
            commas(n)
        ));
    }
    if !plan.ambiguous_aliases.is_empty() {
        w.push(format!(
            "{} collection aliases name more than one collection and are not used as ids",
            commas(plan.ambiguous_aliases.len() as u64)
        ));
    }
    let n = plan.n_orphans();
    if n > 0 {
        w.push(format!(
            "included {} sequences that are in no collection (ids SQ.<digest> only)",
            commas(n)
        ));
    }
    w
}

fn emit(progress: &mut dyn FnMut(&ProgressEvent) -> bool, e: ProgressEvent) -> Result<(), Error> {
    if progress(&e) {
        Ok(())
    } else {
        Err(Error::Cancelled)
    }
}

/// What the sequence phase measured.
pub(crate) struct SeqPhase {
    /// Base counts of each subject written, in order.
    pub counts: Vec<BaseCounts>,
    pub sequences_s: f64,
    pub cpu_read_s: f64,
    pub cpu_convert_s: f64,
    pub write_s: f64,
}

/// Stored bytes a subject is expected to take (its length when unknown).
fn stored_size(s: &Subject) -> u64 {
    s.stored_len.unwrap_or(s.length as u64)
}

/// Converted bytes a subject takes in memory (ncbi2na plus a little for ambiguities).
fn out_bytes(s: &Subject) -> u64 {
    s.length as u64 / 4 + 64
}

/// A worker's reused buffers: one converted subject.
#[derive(Default)]
struct Scratch {
    seq: Vec<u8>,
    amb: Vec<u8>,
}

/// How the sequence phase reads and batches.
pub(crate) struct SeqSource<'a> {
    pub mode: StorageMode,
    pub fetch: &'a dyn SeqFetch,
    /// Conversion workers.
    pub threads: usize,
    /// Converted bytes held by in-flight batches.
    pub budget: u64,
    /// Stored bytes per batch (`--batch-mb`).
    pub batch_bytes: u64,
}

/// Read, convert and append `subjects` to `w`, in order. The fetch reads the subjects' stored
/// bytes ahead in OID order; `run_ordered_with` hands each worker the next batch of
/// consecutive subjects with their bytes, the worker decodes it into one reused
/// [`ConvertedBatch`] arena, and the writer encodes each subject's header and keys from its
/// deflines as it writes it.
pub(crate) fn write_sequences(
    subjects: &[Subject],
    src: &SeqSource,
    w: &mut writer::DbWriter,
    progress_interval_ms: u64,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<SeqPhase, Error> {
    let t_seq = Instant::now();
    let fetch = src.fetch;
    let mode = src.mode;
    let workers = src.threads.max(1);
    let batches = crate::pipeline::plan_batches(
        subjects.iter().map(stored_size),
        src.batch_bytes.max(1),
        crate::pipeline::batch_subject_cap(subjects.len(), workers),
    );
    let pool: Mutex<Vec<ConvertedBatch>> = Mutex::new(Vec::new());
    let t_read = AtomicU64::new(0);
    let t_conv = AtomicU64::new(0);
    let digests: Vec<String> = subjects.iter().map(|s| s.digest.clone()).collect();
    let mut stream = fetch.read_ahead(&digests)?;
    // batch i with its subjects' stored bytes, taken from the stream in order
    let take = |i: usize| {
        let b = batches.get(i)?;
        let t = Instant::now();
        let held = b.oids.clone().map(|oid| {
            let d = &subjects[oid as usize].digest;
            match stream.next() {
                Some(r) => r.map_err(|e| e.with_digest(d)),
                None => Err(Error::sequence(d, "the read ended early")),
            }
        });
        let held = held.collect::<Result<Vec<_>, _>>();
        t_read.fetch_add(t.elapsed().as_micros() as u64, Ordering::Relaxed);
        Some(held.map(|h| (b, h)))
    };
    let cost = |(b, _): &(&crate::pipeline::OidBatch, Vec<HeldBytes>)| {
        subjects[b.oids.start as usize..b.oids.end as usize]
            .iter()
            .map(out_bytes)
            .sum()
    };
    let work = |sc: &mut Scratch, (b, held): (&crate::pipeline::OidBatch, Vec<HeldBytes>)| {
        let mut out = pool.lock().unwrap().pop().unwrap_or_default();
        out.reset(b.oids.start);
        for (s, stored) in subjects[b.oids.start as usize..b.oids.end as usize]
            .iter()
            .zip(&held)
        {
            let stored: &[u8] = (**stored).as_ref();
            if let Some(n) = s.stored_len
                && stored.len() as u64 != n
            {
                return Err(Error::sequence(
                    &s.digest,
                    format!("got {} stored bytes, expected {n}", stored.len()),
                ));
            }
            let t = Instant::now();
            let n = decode::to_blast(mode, s.alphabet, stored, s.length, &mut sc.seq, &mut sc.amb)
                .map_err(|e| e.with_digest(&s.digest))?;
            t_conv.fetch_add(t.elapsed().as_micros() as u64, Ordering::Relaxed);
            out.push(&sc.seq, &sc.amb, n)?;
        }
        Ok(out)
    };
    let total = subjects.len() as u64;
    let bases_total: u64 = subjects.iter().map(|s| s.length as u64).sum();
    let (mut done, mut bases_done) = (0u64, 0u64);
    let mut counts = Vec::with_capacity(subjects.len());
    let mut t_write = 0f64;
    let mut throttle = Throttle::new(progress_interval_ms);
    run_ordered_with(
        take,
        workers,
        src.budget,
        cost,
        Scratch::default,
        work,
        |_, batch: ConvertedBatch| {
            for i in 0..batch.len() {
                let s = &subjects[batch.first_oid as usize + i];
                let (seq, amb, n) = batch.subject(i);
                let t = Instant::now();
                w.add(s.length, seq, amb, &s.deflines)?;
                t_write += t.elapsed().as_secs_f64();
                counts.push(n);
                done += 1;
                bases_done += s.length as u64;
                if throttle.ready() {
                    emit(
                        progress,
                        ProgressEvent::Sequences {
                            done,
                            total,
                            bases_done,
                            bases_total,
                            elapsed_s: t_seq.elapsed().as_secs_f64(),
                        },
                    )?;
                }
            }
            pool.lock().unwrap().push(batch);
            Ok(())
        },
    )?;
    let sequences_s = t_seq.elapsed().as_secs_f64();
    emit(
        progress,
        ProgressEvent::Sequences {
            done,
            total,
            bases_done,
            bases_total,
            elapsed_s: sequences_s,
        },
    )?;
    Ok(SeqPhase {
        counts,
        sequences_s,
        cpu_read_s: t_read.into_inner() as f64 / 1e6,
        cpu_convert_s: t_conv.into_inner() as f64 / 1e6,
        write_s: t_write,
    })
}

/// Each plan collection's N fraction: Σ N / Σ length over the rows in this database, from the
/// per-OID counts. The plan's whole-collection value is used only when that sum can't be made
/// (a row was left out as oversize, or counts are missing). Rows come first so the column does
/// not change with what the N-count cache already holds.
pub fn collection_n_fractions(plan: &SubjectPlan, counts: &[BaseCounts]) -> Vec<Option<f64>> {
    let mut sums = vec![(0u64, 0u64); plan.collections.len()];
    for (oid, s) in plan.subjects.iter().enumerate() {
        let n = counts.get(oid).map(|c| c.n).unwrap_or(0);
        for o in &s.occurrences {
            let e = &mut sums[o.collection as usize];
            e.0 += n;
            e.1 += s.length as u64;
        }
    }
    plan.collections
        .iter()
        .zip(sums)
        .map(|(c, (n, len))| {
            (c.n_oversize == 0 && len > 0 && counts.len() == plan.subjects.len())
                .then(|| n as f64 / len as f64)
                .or(c.n_fraction)
        })
        .collect()
}

/// Append conversion by-products to the N-count cache; a failure is a warning, not an error.
pub(crate) fn cache_by_products(
    opts: &NstatsOptions,
    subjects: &[Subject],
    counts: &[BaseCounts],
    warnings: &mut Vec<String>,
) {
    if !opts.write_cache {
        return;
    }
    let Some(cache) = opts.cache() else {
        return;
    };
    let entries: Vec<(String, NCount)> = subjects
        .iter()
        .zip(counts)
        .filter(|(s, _)| crate::select::nstats::needs_read(s.alphabet))
        .map(|(s, c)| {
            (
                s.digest.clone(),
                NCount {
                    length: s.length as u64,
                    n: c.n,
                    ambiguous: c.ambiguous,
                },
            )
        })
        .collect();
    if let Err(e) = cache.append(&entries) {
        warnings.push(format!(
            "could not cache N counts in {}: {e}",
            cache.dir.display()
        ));
    }
}

/// The manifest's volume list from the writer's result.
pub(crate) fn manifest_volumes(
    names: &[(String, u32)],
    letters: &[u64],
) -> Vec<sidecars::ManifestVolume> {
    let mut first = 0;
    names
        .iter()
        .zip(letters)
        .map(|((n, k), l)| {
            let v = sidecars::ManifestVolume {
                stem: n.clone(),
                first_oid: first,
                n_oids: *k,
                letters: *l,
            };
            first += k;
            v
        })
        .collect()
}

/// Write the database for `plan`, reading sequence bytes from `fetch`. Preflight first (output
/// set, free space), then the sequences through `run_ordered`, then the LMDB, OID files,
/// sidecars, manifest and `.njs`. On any error or cancel, every file written is removed.
pub fn write_blastdb(
    plan: &SubjectPlan,
    fetch: &dyn SeqFetch,
    opts: &WriteOptions,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<BuildSummary, Error> {
    let t0 = Instant::now();
    if opts.max_file_size == 0 {
        return Err(Error::InvalidOptions(
            "--max-file-size must be greater than 0".into(),
        ));
    }
    if plan.subjects.is_empty() {
        return Err(Error::InvalidOptions(format!(
            "the selection has no nucleotide sequences to write ({} skipped)",
            plan.skipped.len()
        )));
    }
    let (prep, mut summary) = plan_preflight(plan, opts);
    let set = OutputSet::new(&opts.out);
    set.check(opts.force)?;
    preflight::check_space(&mut summary, opts.check_free_space)?;
    emit(progress, ProgressEvent::Planned(summary.clone()))?;
    if opts.force {
        set.remove_existing()?;
    }
    set.create_dir()?;

    let guard = PartialOutput::new();
    let mut w = writer::DbWriter::new(
        writer::WriterOptions {
            base: opts.out.clone(),
            title: opts.title.clone(),
            max_file_size: opts.max_file_size,
            date: opts.date.blast(),
        },
        &guard,
    );
    let phase = write_sequences(
        &plan.subjects,
        &SeqSource {
            mode: plan.store_mode,
            fetch,
            threads: opts.threads,
            budget: opts.budget_bytes,
            batch_bytes: opts.batch_bytes,
        },
        &mut w,
        opts.progress_interval_ms,
        progress,
    )?;

    // --- finishing ---
    let t_fin = Instant::now();
    let mut written = w.finish()?;
    let mut files = written.files.clone();
    emit(
        progress,
        ProgressEvent::Finishing {
            stage: "oid files".into(),
        },
    )?;
    // `.nos` reads the keys in OID order, so the OID files come before the LMDB sorts them
    let oid_files = oidfiles::write_oid_files(
        &opts.out,
        written.n_oids,
        &written.keys,
        &prep.taxids,
        &guard,
    )?;
    emit(
        progress,
        ProgressEvent::Finishing {
            stage: "lmdb".into(),
        },
    )?;
    let ndb = with_ext(&opts.out, "ndb");
    guard.record(&ndb);
    lmdb::write_ndb(&ndb, &written.volumes, &mut written.keys)?;
    drop(written.keys);
    files.push(ndb);
    files.extend(oid_files);
    emit(
        progress,
        ProgressEvent::Finishing {
            stage: "sidecars".into(),
        },
    )?;
    let mut extra_warnings = Vec::new();
    cache_by_products(
        &opts.nstats,
        &plan.subjects,
        &phase.counts,
        &mut extra_warnings,
    );
    let fractions = collection_n_fractions(plan, &phase.counts);
    let mut sidecars = sidecars::write_all(&opts.out, plan, Some(&fractions), &guard)?;
    sidecars.push(sidecars::write_oids(
        &opts.out,
        &sidecars::oid_rows(plan, &phase.counts),
        &guard,
    )?);
    let taxids = oidfiles::distinct_taxids(&prep.taxids);
    if let Some(m) = &opts.manifest {
        let mut m = m.clone();
        m.title = opts.title.clone();
        m.date = opts.date.blast();
        m.volumes = manifest_volumes(&written.volumes, &written.volume_letters);
        m.files = files
            .iter()
            .map(|f| writer::file_name(f))
            .filter(|n| !n.ends_with(".nal"))
            .collect();
        m.files.sort();
        m.totals = plan.totals;
        sidecars.push(sidecars::write_manifest(&opts.out, &m, Some(&guard))?);
    }
    emit(
        progress,
        ProgressEvent::Finishing {
            stage: "njs".into(),
        },
    )?;
    let njs_path = with_ext(&opts.out, "njs");
    guard.record(&njs_path);
    njs::write(&njs::NjsInput {
        base: &opts.out,
        title: &opts.title,
        letters: written.total_letters,
        sequences: written.n_oids as u64,
        last_updated: opts.date.iso(),
        volumes: written.volumes.len() as u64,
        files: &files,
        taxids: &taxids,
    })?;
    files.push(njs_path);
    files.sort();
    let finish_s = t_fin.elapsed().as_secs_f64();
    guard.commit();

    let mut bytes_written = 0u64;
    for f in files.iter().chain(sidecars.iter()) {
        bytes_written += std::fs::metadata(f).map_err(|e| Error::io(f, e))?.len();
    }
    let mut warnings = warnings(plan, &opts.out);
    warnings.extend(extra_warnings);
    let s = BuildSummary {
        preflight: summary,
        taxa: taxa_of(plan),
        files,
        sidecars,
        bytes_written,
        ambiguous_bases: phase.counts.iter().map(|c| c.ambiguous).sum(),
        volumes: written.volumes.len() as u64,
        warnings: warnings.clone(),
        timings: Timings {
            plan_s: 0.0,
            sequences_s: phase.sequences_s,
            finish_s,
            cpu_read_s: phase.cpu_read_s,
            cpu_convert_s: phase.cpu_convert_s,
            write_s: phase.write_s,
        },
        network: None,
        elapsed_s: t0.elapsed().as_secs_f64(),
    };
    // The database is complete; a cancel from here on no longer removes it.
    for wn in warnings {
        progress(&ProgressEvent::Warning(wn));
    }
    Ok(s)
}

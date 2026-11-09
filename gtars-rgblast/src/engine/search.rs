//! The library entry point: `SearchParams` in, `SearchResults` out.
//!
//! Orchestration mirrors BLAST+ 2.17.0 `app/blast/blastn_app.cpp` (the query batch loop with
//! `CBatchSizeMixer` from `blast_app_util.hpp`/`.cpp` and `CBlastInput::GetNextSeqBatch`),
//! then per batch `algo/blast/api/prelim_stage.cpp` / `traceback_stage.cpp`: query setup,
//! `BLAST_MainSetUp` (score block), lookup table, `BLAST_GapAlignSetUp` (parameters), the
//! subject loop of `BLAST_PreliminarySearchEngine` feeding the collector, then
//! `BLAST_ComputeTraceback`. Threads split the database as blastn's `-mt_mode 0` does, by work
//! items (runs of subjects by letters; a long subject alone) taken from an atomic counter;
//! see `build_work_items` for why neither the split nor the timing changes output.
//! `BLAST_GetSubjectTotals` (blast_setup.c) decides the database statistics: the STATS
//! override when given, else the real totals.
//!
//! Query batching matters for the output: every batch is its own `CLocalBlast` run with its
//! own concatenated query block, lookup table (type and width depend on the batch's total
//! length) and parameters, and `-task blastn`'s `s_NuclUngappedExtend` scores 4-base blocks
//! that straddle the sentinel before a query that is not first in its block. blastn's
//! `CBlastnAppArgs::GetQueryBatchSize()` returns 0 (no fixed size) unless the `BATCH_SIZE`
//! environment variable is set, so the adaptive mixer decides: target = max(1,000,000,
//! db_total_length / 3000) saved initial hits per batch; the first batch is target / 200
//! letters; a batch takes whole queries while the letters read so far are below the batch
//! size; after each batch the ratio (good_init_extends + 1) / batch_size is mixed 0.3 / 0.7
//! into the running ratio and the next size is target / ratio, clamped to
//! [100, chunk - 1000] (chunk = `SplitQuery_GetChunkSize`: 1,000,000 for blastn and
//! blastn-short, 5,000,000 for megablast); a clamp resets the running ratio. A single query
//! longer than that chunk is refused (`Error::InvalidOptions`): blastn would split it
//! (`split_query.c`), which is not ported.
//! `good_init_extends` is `Blast_UngappedStatsUpdate`'s sum of `init_hitlist->total` over the
//! word finder calls (one per subject chunk), read back through `CLocalBlast::GetNumExtensions`.

use super::hsp::{EditScript, HspResults, calc_length_and_gaps, evalue_compare_hsps};
use super::lookup::Lookup;
use super::options::MAX_DBSEQ_LEN;
use super::options::{BlastOptions, ContextStats, DustSpec, Parameters, Strand, Task};
use super::prelim::{Scratch, SearchContext, collect, finish_subject, merge_chunk, search_chunk};
use super::query::{QueryBlock, QuerySet};
use super::stat::ScoreBlk;
pub use super::subject::DbStatsOverride;
use super::subject::{ChunkBuf, StoreSubjects};
use super::traceback::compute_traceback;
use crate::error::Error;
use crate::progress::ProgressEvent;
use serde::{Deserialize, Serialize};
use std::ops::Range;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Mutex, mpsc};
use std::time::Instant;

/// Every accepted blastn flag, typed. `Default` = `-task megablast` defaults.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SearchParams {
    pub task: Task,
    /// `-evalue`: None = the task's default (10; 1000 for blastn-short).
    pub evalue: Option<f64>,
    pub word_size: Option<u32>,
    /// `-dust`: None = the task's default, Some(None) = `no`, Some(Some(spec)) = `yes`/explicit.
    pub dust: Option<Option<DustSpec>>,
    pub soft_masking: bool,
    pub strand: Strand,
    pub perc_identity: f64,
    pub max_target_seqs: u32,
    pub max_hsps: Option<u32>,
    /// `STATS_NSEQ` + `STATS_TOTLEN`: replaces both statistics terms (not the query batching,
    /// which reads the real total length).
    pub stats: Option<DbStatsOverride>,
    pub threads: usize,
    pub progress_interval_ms: u64,
    /// blastn's `BATCH_SIZE` environment variable: a fixed query batch size in letters.
    /// `None` or `Some(0)` = blastn's default, the adaptive `CBatchSizeMixer`.
    #[serde(default)]
    pub query_batch_size: Option<u32>,
    /// Subject bases per preliminary work item (consecutive subjects are grouped until this
    /// many letters or `MAX_BATCH_SUBJECTS`). `None` = `default_batch_letters`. Output does
    /// not depend on it.
    #[serde(default)]
    pub batch_letters: Option<u64>,
}

impl Default for SearchParams {
    fn default() -> Self {
        SearchParams {
            task: Task::Megablast,
            evalue: None,
            word_size: None,
            dust: None,
            soft_masking: true,
            strand: Strand::Both,
            perc_identity: 0.0,
            max_target_seqs: 500,
            max_hsps: None,
            stats: None,
            threads: 1,
            progress_interval_ms: 1000,
            query_batch_size: None,
            batch_letters: None,
        }
    }
}

/// Stack of every worker thread (the main thread's default on Linux; the ported loops are
/// iterative, this only keeps deep frames of the aligners off a 2 MiB default).
pub const WORKER_STACK: usize = 8 << 20;

/// `RGBLASTN_DEBUG_SHUFFLE=<seed>`: permute a work list with a seeded xorshift before the
/// threads take from it. Output must not change (the determinism tests set it); read per call
/// so a test can turn it on and off.
pub fn debug_shuffle<T>(items: &mut [T]) {
    let Some(seed) = std::env::var_os("RGBLASTN_DEBUG_SHUFFLE") else {
        return;
    };
    let mut x: u64 = seed
        .to_str()
        .and_then(|s| s.trim().parse::<u64>().ok())
        .unwrap_or(0x9E37_79B9_7F4A_7C15)
        | 1;
    // Fisher-Yates with xorshift64*
    for i in (1..items.len()).rev() {
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        let j = (x.wrapping_mul(0x2545_F491_4F6C_DD1D) % (i as u64 + 1)) as usize;
        items.swap(i, j);
    }
}

/// Most subjects in one work item.
pub const MAX_BATCH_SUBJECTS: u32 = 4096;
/// Bounds of the default letters per work item: `clamp(total / (16 * threads), 1 Mbp, 64 Mbp)`.
pub const MIN_BATCH_LETTERS: u64 = 1_000_000;
pub const MAX_BATCH_LETTERS: u64 = 64_000_000;

/// The default letters per work item: about sixteen items per thread, clamped to
/// [`MIN_BATCH_LETTERS`, `MAX_BATCH_LETTERS`].
pub fn default_batch_letters(total_length: i64, threads: usize) -> u64 {
    (total_length.max(0) as u64 / (16 * threads.max(1) as u64))
        .clamp(MIN_BATCH_LETTERS, MAX_BATCH_LETTERS)
}

/// The preliminary work items (plan step 25) for a database whose subjects (in OID order)
/// have `lengths`: consecutive subjects grouped while the group holds under `batch_letters`
/// letters and under `MAX_BATCH_SUBJECTS` subjects (a subject alone always fits), and a
/// subject longer than `MAX_DBSEQ_LEN` alone, searched chunk by chunk by one thread as
/// blastn's threads do. Longest first, so the long subjects do not finish last. Neither the
/// grouping nor the order changes the output: `Blast_HitListUpdate` is order-independent.
pub fn build_work_items(lengths: &[u32], batch_letters: u64) -> Vec<Range<u32>> {
    let n = lengths.len() as u32;
    let len = |oid: u32| lengths[oid as usize] as u64;
    let long = |oid: u32| len(oid) > MAX_DBSEQ_LEN as u64;
    let mut items = Vec::new();
    let mut oid = 0u32;
    while oid < n {
        let mut end = oid + 1;
        let mut letters = len(oid);
        while !long(oid)
            && end < n
            && end - oid < MAX_BATCH_SUBJECTS
            && !long(end)
            && letters + len(end) <= batch_letters
        {
            letters += len(end);
            end += 1;
        }
        items.push(oid..end);
        oid = end;
    }
    items.sort_by_cached_key(|r| std::cmp::Reverse(r.clone().map(len).sum::<u64>()));
    items
}

/// The preliminary results (`hspfilter_collector.c`'s writer, `prelim::collect`; its
/// `Blast_HitListUpdate` is order-independent, so the threads feed it as they finish) plus
/// the diagnostics the threads fold together.
struct Collector {
    results: HspResults,
    prelim_hitlist_size: i32,
    /// `Blast_DiagnosticsUpdate`: the threads' `good_init_extends` summed (wrapping `Int4`).
    good_init_extends: i32,
    subjects_with_hits: u64,
}

/// What the preliminary threads of one query batch share.
struct PrelimShared<'a> {
    ctx: &'a SearchContext<'a>,
    subjects: &'a StoreSubjects,
    items: Vec<Range<u32>>,
    next: AtomicUsize,
    cancel: AtomicBool,
    collector: Mutex<Collector>,
    max_subject_length: i32,
}

/// Messages from a preliminary worker to the main thread (progress and errors; results go
/// through the collector).
enum PrelimMsg {
    Done { subjects: u32, bases: u64 },
    Error(Error),
}

fn prelim_worker(sh: &PrelimShared, tx: &mpsc::Sender<PrelimMsg>) {
    let mut scratch = Scratch::new(sh.ctx, sh.max_subject_length);
    let mut buf = ChunkBuf::default();
    loop {
        if sh.cancel.load(Ordering::Relaxed) {
            break;
        }
        let i = sh.next.fetch_add(1, Ordering::Relaxed);
        if i >= sh.items.len() {
            break;
        }
        match run_item(sh, sh.items[i].clone(), &mut scratch, &mut buf) {
            Ok(msg) => {
                if tx.send(msg).is_err() {
                    break;
                }
            }
            Err(e) => {
                let _ = tx.send(PrelimMsg::Error(e));
                break;
            }
        }
    }
    let mut c = sh.collector.lock().unwrap();
    c.good_init_extends = c.good_init_extends.wrapping_add(scratch.good_init_extends);
}

/// One work item: each subject decoded chunk by chunk into this thread's buffer, each chunk
/// searched and merged in order (`s_BlastSearchEngineCore`), the lists handed to the
/// collector.
fn run_item(
    sh: &PrelimShared,
    oids: Range<u32>,
    scratch: &mut Scratch,
    buf: &mut ChunkBuf,
) -> Result<PrelimMsg, Error> {
    let ctx = sh.ctx;
    let mut lists = Vec::new();
    let mut bases = 0u64;
    for oid in oids.clone() {
        let mut combined = None;
        sh.subjects
            .for_each_chunk(oid, buf, |seq, offset, chunk_len| {
                if let Some(list) = search_chunk(ctx, seq, offset, chunk_len, oid as i32, scratch) {
                    merge_chunk(ctx, list, &mut combined, offset);
                }
            })?;
        lists.extend(finish_subject(ctx, combined));
        bases += sh.subjects.seq_len(oid as i32) as u64;
    }
    {
        let mut c = sh.collector.lock().unwrap();
        c.subjects_with_hits += lists.len() as u64;
        let size = c.prelim_hitlist_size;
        for l in lists {
            collect(&mut c.results, l, size);
        }
    }
    Ok(PrelimMsg::Done {
        subjects: oids.end - oids.start,
        bases,
    })
}

/// `BLAST_PreliminarySearchEngine` over the work items on `threads` threads. Returns the
/// preliminary results and `good_init_extends` for the mixer.
fn run_prelim(
    ctx: &SearchContext,
    subjects: &StoreSubjects,
    items: Vec<Range<u32>>,
    threads: usize,
    prelim_hitlist_size: i32,
    rp: &mut RunProgress,
) -> Result<(HspResults, i32), Error> {
    let num_seqs = subjects.num_seqs().max(0) as u32;
    let n_items = items.len();
    let sh = PrelimShared {
        ctx,
        subjects,
        items,
        next: AtomicUsize::new(0),
        cancel: AtomicBool::new(false),
        collector: Mutex::new(Collector {
            results: HspResults::new(ctx.qb.num_queries()),
            prelim_hitlist_size,
            good_init_extends: 0,
            subjects_with_hits: 0,
        }),
        max_subject_length: subjects.max_len(),
    };
    let workers = threads.clamp(1, n_items.max(1));
    let mut first_error: Option<Error> = None;
    let mut cancelled = false;
    let mut subjects_done = 0u32;
    let (tx, rx) = mpsc::channel::<PrelimMsg>();
    std::thread::scope(|s| {
        for t in 0..workers {
            let tx = tx.clone();
            let sh = &sh;
            std::thread::Builder::new()
                .name(format!("rgblastn-{t}"))
                .stack_size(WORKER_STACK)
                .spawn_scoped(s, move || prelim_worker(sh, &tx))
                .expect("spawn search worker");
        }
        drop(tx);
        loop {
            let wait = rp.interval.saturating_sub(rp.last_report.elapsed());
            let mut finished = false;
            match rx.recv_timeout(wait) {
                Ok(PrelimMsg::Done { subjects, bases }) => {
                    subjects_done += subjects;
                    rp.bases_done += bases;
                }
                Ok(PrelimMsg::Error(e)) => {
                    if first_error.is_none() {
                        first_error = Some(e);
                    }
                    sh.cancel.store(true, Ordering::Relaxed);
                }
                Err(mpsc::RecvTimeoutError::Timeout) => {}
                Err(mpsc::RecvTimeoutError::Disconnected) => finished = true,
            }
            if !cancelled
                && first_error.is_none()
                && (rp.last_report.elapsed() >= rp.interval || finished)
            {
                rp.last_report = Instant::now();
                let hits = sh.collector.lock().unwrap().subjects_with_hits;
                if !(rp.progress)(&ProgressEvent::SearchProgress {
                    subjects_done,
                    subjects_total: num_seqs,
                    bases_done: rp.bases_done,
                    hits,
                    elapsed_s: rp.t0.elapsed().as_secs_f64(),
                }) {
                    cancelled = true;
                    sh.cancel.store(true, Ordering::Relaxed);
                }
            }
            if finished {
                break;
            }
        }
    });
    if let Some(e) = first_error {
        return Err(e);
    }
    if cancelled {
        return Err(Error::Cancelled);
    }
    let c = sh.collector.into_inner().unwrap();
    rp.subjects_with_hits += c.subjects_with_hits;
    Ok((c.results, c.good_init_extends))
}

impl SearchParams {
    /// `CBlastOptionsFactory::CreateTask` + the CLI overrides of `blast_args.cpp`.
    pub fn to_options(&self) -> Result<BlastOptions, Error> {
        let mut o = BlastOptions::for_task(self.task);
        if let Some(e) = self.evalue {
            o.expect_value = e;
        }
        if let Some(w) = self.word_size {
            o.word_size = w as i32;
        }
        if let Some(d) = self.dust {
            o.dust = d;
        }
        o.mask_at_hash = self.soft_masking;
        o.strand = self.strand;
        o.percent_identity = self.perc_identity;
        o.hitlist_size = self.max_target_seqs as i32;
        o.max_hsps_per_subject = self.max_hsps.unwrap_or(0) as i32;
        o.validate()?;
        Ok(o)
    }
}

/// One final HSP, with blastn's printed coordinates (1-based; `s_start > s_end` on the minus
/// strand, the query always plus).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Hsp {
    pub score: i32,
    pub bit_score: f64,
    pub evalue: f64,
    pub num_ident: i32,
    pub align_length: i32,
    pub gaps: i32,
    pub gap_opens: i32,
    pub q_start: i32,
    pub q_end: i32,
    pub s_start: i32,
    pub s_end: i32,
    pub minus: bool,
    /// Engine coordinates (context-relative query offsets, 0-based half-open) and the
    /// edit script, for `qseq`/`sseq`/`btop` later.
    pub query_offset: i32,
    pub query_end: i32,
    pub subject_offset: i32,
    pub subject_end: i32,
    pub context: i32,
    pub edit_script: Option<EditScript>,
    /// Aligned letters in printed orientation, filled by `format::alignment::fill` when an
    /// output field needs them (`qseq`, `sseq`, `btop`, `positive`, `ppos`).
    pub qseq: Option<String>,
    pub sseq: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Hit {
    pub oid: u32,
    pub hsps: Vec<Hsp>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct QueryResult {
    pub id: String,
    pub defline: String,
    pub length: u32,
    /// In blastn output order.
    pub hits: Vec<Hit>,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct SearchSummary {
    pub queries: u32,
    pub subjects: u32,
    pub total_length: u64,
    pub db_length_used: i64,
    pub db_num_seqs_used: i32,
    pub hits: u64,
    pub hsps: u64,
    pub elapsed_s: f64,
    pub prelim_s: f64,
    pub traceback_s: f64,
    /// The lookup table(s) built: one description, or the distinct ones joined with " | "
    /// when the query batches chose differently.
    pub lookup: String,
    /// Queries per query batch, in order (`blastn_app.cpp`'s batch loop).
    #[serde(default)]
    pub batches: Vec<u32>,
    /// `CLocalBlast::GetNumExtensions()` per batch: the saved initial hits that fed the mixer.
    #[serde(default)]
    pub good_init_extends: Vec<i32>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SearchResults {
    pub queries: Vec<QueryResult>,
    pub stats: Vec<ContextStats>,
    pub summary: SearchSummary,
}

/// `SplitQuery_GetChunkSize` (local_blast.cpp) for the task's program: `eBlastn` (blastn and
/// blastn-short) 1,000,000; `eMegablast` 5,000,000. The mixer's maximum batch is 1,000
/// letters below it, and no single query may be longer (`search`).
pub fn query_chunk_size(task: Task) -> i32 {
    match task {
        Task::Blastn | Task::BlastnShort => 1_000_000,
        Task::Megablast => 5_000_000,
    }
}

/// Refuse a query longer than the task's chunk size: blastn splits it (`split_query.c`),
/// which is not ported, so its output could differ.
pub fn check_query_lengths(task: Task, queries: &QuerySet) -> Result<(), Error> {
    let limit = query_chunk_size(task) as usize;
    match queries.queries.iter().find(|q| q.seq.len() > limit) {
        Some(q) => Err(Error::InvalidOptions(format!(
            "query {} is {} bases, longer than the {limit} rgblastn searches in one piece; split it, or build a database with rgmakeblastdb and use blastn",
            q.id,
            q.seq.len()
        ))),
        None => Ok(()),
    }
}

/// `(Int4) x` for a double as gcc compiles it on x86-64 (`cvttsd2si` to a 32-bit register):
/// truncation toward zero, and the "integer indefinite" `INT_MIN` for NaN or any value out of
/// range. The mixer reaches that case when a batch saved no initial hits (target / ratio is
/// then in the billions), and the result clamps to the minimum batch, not the maximum.
fn c_int4(x: f64) -> i32 {
    if x.is_nan() || x >= 2_147_483_648.0 || x <= -2_147_483_649.0 {
        i32::MIN
    } else {
        x.trunc() as i32
    }
}

/// `CBatchSizeMixer` (blast_app_util.hpp / .cpp), member for member.
#[derive(Debug, Clone, PartialEq)]
pub struct BatchSizeMixer {
    /// `m_TargetHits`
    target_hits: i32,
    /// `m_Ratio` (-1 = none yet)
    ratio: f64,
    /// `m_BatchSize`
    batch_size: i32,
    /// `k_MaxBatchSize`
    max_batch_size: i32,
}

/// `k_MixIn`
const MIX_IN: f64 = 0.3;
/// `k_MinBatchSize`
const MIN_BATCH_SIZE: i32 = 100;
/// `k_MinTargetHits`
const MIN_TARGET_HITS: i32 = 1_000_000;

impl BatchSizeMixer {
    /// `CBatchSizeMixer(max_batch_size)`
    pub fn new(max_batch_size: i32) -> Self {
        BatchSizeMixer {
            target_hits: 1_000_000,
            ratio: -1.0,
            batch_size: 5000,
            max_batch_size,
        }
    }

    /// `SetTargetHits`
    pub fn set_target_hits(&mut self, target: i32) {
        self.target_hits = if target < MIN_TARGET_HITS {
            MIN_TARGET_HITS
        } else {
            target
        };
        self.batch_size = self.target_hits / 200;
    }

    /// `GetBatchSize(hits)`: `hits` < 0 only clamps (the first call).
    pub fn batch_size(&mut self, hits: i32) -> i32 {
        if hits >= 0 {
            let ratio = 1.0 * (hits.wrapping_add(1)) as f64 / self.batch_size as f64;
            self.ratio = if self.ratio < 0.0 {
                ratio
            } else {
                MIX_IN * ratio + (1.0 - MIX_IN) * self.ratio
            };
            self.batch_size = c_int4(1.0 * self.target_hits as f64 / self.ratio);
        }
        if self.batch_size > self.max_batch_size {
            self.batch_size = self.max_batch_size;
            self.ratio = -1.0;
        } else if self.batch_size < MIN_BATCH_SIZE {
            self.batch_size = MIN_BATCH_SIZE;
            self.ratio = -1.0;
        }
        self.batch_size
    }

    pub fn target_hits(&self) -> i32 {
        self.target_hits
    }

    pub fn ratio(&self) -> f64 {
        self.ratio
    }
}

/// The query batch sizes of a run, fixed (`BATCH_SIZE`) or from the mixer.
enum BatchSizing {
    Fixed(u32),
    Adaptive(BatchSizeMixer),
}

impl BatchSizing {
    /// `blastn_app.cpp`: `GetQueryBatchSize()` when set, else the mixer with
    /// `SetTargetHits(formatter.GetDbTotalLength() / 3000)` (the database's real total
    /// length, `CSeqDB::GetTotalLength`, not the STATS override) when that total is positive.
    fn new(params: &SearchParams, db_total_length: i64) -> (BatchSizing, u32) {
        match params.query_batch_size {
            Some(n) if n > 0 => (BatchSizing::Fixed(n), n),
            _ => {
                let mut mixer = BatchSizeMixer::new(query_chunk_size(params.task) - 1000);
                if db_total_length > 0 {
                    // Int8 total_len / 3000 passed to an Int4 parameter
                    mixer.set_target_hits((db_total_length / 3000) as i32);
                }
                let first = mixer.batch_size(-1) as u32;
                (BatchSizing::Adaptive(mixer), first)
            }
        }
    }

    /// `input.SetBatchSize(mixer.GetBatchSize(lcl_blast.GetNumExtensions()))` after a batch
    /// (a fixed size never changes).
    fn next(&mut self, good_init_extends: i32) -> u32 {
        match self {
            BatchSizing::Fixed(n) => *n,
            BatchSizing::Adaptive(m) => m.batch_size(good_init_extends) as u32,
        }
    }
}

/// `CBlastInput::GetNextSeqBatch`: the end of the batch starting at `start`, taking whole
/// queries while the letters read so far are below `batch_size` (so always at least one).
pub fn batch_end(lengths: &[usize], start: usize, batch_size: u32) -> usize {
    let mut end = start;
    let mut size_read = 0u64;
    while size_read < batch_size as u64 && end < lengths.len() {
        size_read += lengths[end] as u64;
        end += 1;
    }
    end
}

/// `BLAST_GetSubjectTotals`: (total_length, num_seqs) with the STATS override.
pub fn subject_totals(s: &StoreSubjects) -> (i64, i32) {
    let mut total_length = s.total_length_stats();
    if total_length <= 0 {
        total_length = s.total_length();
    }
    let mut num_seqs = -1;
    if total_length > 0 {
        num_seqs = s.num_seqs_stats();
        if num_seqs <= 0 {
            num_seqs = s.num_seqs();
        }
    }
    (total_length, num_seqs)
}

/// Query block, lookup table, score block and parameters for a search (the setup half).
pub struct Setup {
    pub options: BlastOptions,
    pub qb: QueryBlock,
    pub lookup: Lookup,
    pub sbp: ScoreBlk,
    pub params: Parameters,
    pub ctx: Vec<ContextStats>,
}

pub fn setup(
    params: &SearchParams,
    queries: &QuerySet,
    subjects: &StoreSubjects,
) -> Result<Setup, Error> {
    let options = params.to_options()?;
    let mut qb = QueryBlock::build(queries, &options);
    let letters: Vec<Option<&[u8]>> = qb
        .contexts
        .iter()
        .enumerate()
        .map(|(i, c)| {
            if c.is_valid {
                Some(qb.context_seq(i))
            } else {
                None
            }
        })
        .collect();
    let (sbp, valid) = ScoreBlk::new(
        options.reward,
        options.penalty,
        options.gap_open,
        options.gap_extend,
        options.gapped_calculation,
        &letters,
    )?;
    for (c, v) in qb.contexts.iter_mut().zip(valid) {
        if !v {
            c.is_valid = false;
        }
    }
    let lookup = Lookup::build(&mut qb, options.word_size, options.mask_at_hash)?;
    let mut ctx: Vec<ContextStats> = qb
        .contexts
        .iter()
        .map(|c| ContextStats {
            query_index: c.query_index,
            frame: c.frame,
            query_length: c.query_length,
            is_valid: c.is_valid,
            ..Default::default()
        })
        .collect();
    let (total_length, num_seqs) = subject_totals(subjects);
    let p = Parameters::new(
        &options,
        &sbp,
        &mut ctx,
        total_length,
        num_seqs,
        subjects.avg_len(),
    )?;
    Ok(Setup {
        options,
        qb,
        lookup,
        sbp,
        params: p,
        ctx,
    })
}

/// Progress bookkeeping shared by the batches of one run.
struct RunProgress<'a> {
    t0: Instant,
    interval: std::time::Duration,
    last_report: Instant,
    bases_done: u64,
    subjects_with_hits: u64,
    progress: &'a mut dyn FnMut(&ProgressEvent) -> bool,
}

/// One query batch's results (`CLocalBlast::Run` + the formatter's per-batch output).
struct BatchResults {
    queries: Vec<QueryResult>,
    stats: Vec<ContextStats>,
    prelim_s: f64,
    traceback_s: f64,
    hits: u64,
    hsps: u64,
    /// Empty for a batch without a valid context.
    lookup: String,
    db_length_used: Option<i64>,
    db_num_seqs_used: Option<i32>,
    /// `CLocalBlast::GetNumExtensions()`
    good_init_extends: i32,
}

/// The start of `BLAST_MainSetUp`'s message when no context of the block is valid
/// (`Blast_ScoreBlkKbpUngappedCalc` status 1; `engine::stat`).
const NO_VALID_CONTEXT: &str = "Could not calculate ungapped Karlin-Altschul parameters";

/// Run a search: blastn's query batch loop, each batch a full search of the database on
/// `params.threads` threads (`run_prelim`, `compute_traceback`), outputs concatenated in
/// query order.
pub fn search(
    params: &SearchParams,
    queries: &QuerySet,
    subjects: &StoreSubjects,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<SearchResults, Error> {
    let t0 = Instant::now();
    // validate the options once, before any batch (CBlastOptionsHandle::Validate)
    params.to_options()?;
    check_query_lengths(params.task, queries)?;
    let num_seqs = subjects.num_seqs();
    let total_bases = subjects.total_length() as u64;
    if !progress(&ProgressEvent::SearchStarted {
        queries: queries.queries.len() as u32,
        subjects: num_seqs as u32,
        bases: total_bases,
    }) {
        return Err(Error::Cancelled);
    }
    let mut rp = RunProgress {
        t0,
        interval: std::time::Duration::from_millis(params.progress_interval_ms.max(1)),
        last_report: Instant::now(),
        bases_done: 0,
        subjects_with_hits: 0,
        progress,
    };
    let nq = queries.queries.len();
    let lengths: Vec<usize> = queries.queries.iter().map(|q| q.seq.len()).collect();
    let (mut sizing, mut batch_size) = BatchSizing::new(params, subjects.total_length());
    let mut out = SearchResults {
        queries: Vec::with_capacity(nq),
        stats: Vec::with_capacity(nq * 2),
        summary: SearchSummary {
            queries: nq as u32,
            subjects: num_seqs as u32,
            total_length: total_bases,
            ..Default::default()
        },
    };
    let mut lookups: Vec<String> = Vec::new();
    let mut start = 0usize;
    while start < nq {
        let end = batch_end(&lengths, start, batch_size);
        let batch = QuerySet {
            queries: queries.queries[start..end].to_vec(),
        };
        let r = search_batch(params, &batch, start, subjects, &mut rp)?;
        out.queries.extend(r.queries);
        out.stats.extend(r.stats);
        out.summary.hits += r.hits;
        out.summary.hsps += r.hsps;
        out.summary.prelim_s += r.prelim_s;
        out.summary.traceback_s += r.traceback_s;
        if let (Some(l), Some(n)) = (r.db_length_used, r.db_num_seqs_used) {
            out.summary.db_length_used = l;
            out.summary.db_num_seqs_used = n;
        }
        out.summary.batches.push((end - start) as u32);
        out.summary.good_init_extends.push(r.good_init_extends);
        if !r.lookup.is_empty() && !lookups.contains(&r.lookup) {
            lookups.push(r.lookup);
        }
        batch_size = sizing.next(r.good_init_extends);
        start = end;
    }
    out.summary.lookup = lookups.join(" | ");
    out.summary.elapsed_s = t0.elapsed().as_secs_f64();
    let _ = (rp.progress)(&ProgressEvent::SearchFinished(out.summary.clone()));
    Ok(out)
}

/// One batch: `CLocalBlast::Run` on `queries` (query indices `qstart..`), prelim stage then
/// traceback, and the per-query results in blastn's order.
fn search_batch(
    params: &SearchParams,
    queries: &QuerySet,
    qstart: usize,
    subjects: &StoreSubjects,
    rp: &mut RunProgress,
) -> Result<BatchResults, Error> {
    let t0 = Instant::now();
    let st = match setup(params, queries, subjects) {
        Ok(st) => st,
        Err(Error::InvalidOptions(msg)) if msg.starts_with(NO_VALID_CONTEXT) => {
            // BLAST_MainSetUp's message is eBlastSevWarning, so CSetupFactory::CreateScoreBlock
            // does not throw: blastn prints it once per query of the batch
            // ("Warning: [blastn] <id>: ..."; the binary adds no prefix to SearchWarning), the
            // batch finds nothing, and CLocalBlast::GetNumExtensions hands the mixer 0
            for q in &queries.queries {
                let w = format!("{}: {msg}", q.id);
                if !(rp.progress)(&ProgressEvent::SearchWarning(w)) {
                    return Err(Error::Cancelled);
                }
            }
            let qb = QueryBlock::build(queries, &params.to_options()?);
            return Ok(BatchResults {
                queries: queries
                    .queries
                    .iter()
                    .map(|q| QueryResult {
                        id: q.id.clone(),
                        defline: q.defline.clone(),
                        length: q.seq.len() as u32,
                        hits: Vec::new(),
                    })
                    .collect(),
                stats: qb
                    .contexts
                    .iter()
                    .map(|c| ContextStats {
                        query_index: c.query_index + qstart as u32,
                        frame: c.frame,
                        query_length: c.query_length,
                        is_valid: false,
                        ..Default::default()
                    })
                    .collect(),
                prelim_s: 0.0,
                traceback_s: 0.0,
                hits: 0,
                hsps: 0,
                lookup: String::new(),
                db_length_used: None,
                db_num_seqs_used: None,
                good_init_extends: 0,
            });
        }
        Err(e) => return Err(e),
    };
    let eff: Vec<i64> = st.ctx.iter().map(|c| c.eff_searchsp).collect();
    let sctx = SearchContext {
        qb: &st.qb,
        lookup: &st.lookup,
        sbp: &st.sbp,
        params: &st.params,
        options: &st.options,
        eff_searchsp: eff,
    };
    let nq = queries.queries.len();
    let prelim_size =
        super::hsp::prelim_hitlist_size(st.options.hitlist_size, st.options.gapped_calculation);
    let threads = params.threads.max(1);
    let batch_letters = match params.batch_letters {
        Some(n) if n > 0 => n,
        _ => default_batch_letters(subjects.total_length(), threads),
    };
    let mut items = build_work_items(&subjects.lengths(), batch_letters);
    debug_shuffle(&mut items);
    let (prelim, good_init_extends) = run_prelim(&sctx, subjects, items, threads, prelim_size, rp)?;
    let prelim_s = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();
    let mut tb_progress =
        |done: u32, total: u32| (rp.progress)(&ProgressEvent::Traceback { done, total });
    let results = compute_traceback(
        &sctx,
        prelim,
        subjects,
        threads,
        subjects.max_len(),
        &mut tb_progress,
    )?;
    let traceback_s = t1.elapsed().as_secs_f64();
    let mut stats = st.ctx.clone();
    for c in &mut stats {
        c.query_index += qstart as u32;
    }
    let mut out = BatchResults {
        queries: Vec::with_capacity(nq),
        stats,
        prelim_s,
        traceback_s,
        hits: 0,
        hsps: 0,
        lookup: format!(
            "{:?} width {} step {}",
            st.lookup.lut_type, st.lookup.lut_word_length, st.lookup.scan_step
        ),
        db_length_used: Some(st.params.db_length),
        db_num_seqs_used: Some(st.params.db_num_seqs),
        good_init_extends,
    };
    for (qi, q) in queries.queries.iter().enumerate() {
        let qlen = q.seq.len() as i32;
        let mut hits = Vec::new();
        if let Some(hl) = &results.hitlists[qi] {
            for l in &hl.hsplists {
                // blast_seqalign.cpp sorts each list by e-value before building Seq-aligns
                let mut hsps: Vec<&super::hsp::Hsp> = l.hsps.iter().collect();
                hsps.sort_by(|a, b| evalue_compare_hsps(a, b));
                let mut outh = Vec::with_capacity(hsps.len());
                for h in hsps {
                    let (align_length, gaps, gap_opens) = calc_length_and_gaps(h);
                    let minus = h.query.frame < 0;
                    let (q_start, q_end, s_start, s_end) = if minus {
                        (
                            qlen - h.query.end + 1,
                            qlen - h.query.offset,
                            h.subject.end,
                            h.subject.offset + 1,
                        )
                    } else {
                        (
                            h.query.offset + 1,
                            h.query.end,
                            h.subject.offset + 1,
                            h.subject.end,
                        )
                    };
                    outh.push(Hsp {
                        score: h.score,
                        bit_score: h.bit_score,
                        evalue: h.evalue,
                        num_ident: h.num_ident,
                        align_length,
                        gaps,
                        gap_opens,
                        q_start,
                        q_end,
                        s_start,
                        s_end,
                        minus,
                        query_offset: h.query.offset,
                        query_end: h.query.end,
                        subject_offset: h.subject.offset,
                        subject_end: h.subject.end,
                        context: h.context,
                        edit_script: h.gap_info.clone(),
                        qseq: None,
                        sseq: None,
                    });
                }
                out.hsps += outh.len() as u64;
                out.hits += 1;
                hits.push(Hit {
                    oid: l.oid as u32,
                    hsps: outh,
                });
            }
        }
        out.queries.push(QueryResult {
            id: q.id.clone(),
            defline: q.defline.clone(),
            length: qlen as u32,
            hits,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chunk_size_per_program() {
        assert_eq!(query_chunk_size(Task::Blastn), 1_000_000);
        assert_eq!(query_chunk_size(Task::BlastnShort), 1_000_000);
        assert_eq!(query_chunk_size(Task::Megablast), 5_000_000);
    }

    #[test]
    fn long_query_is_refused() {
        let q = |id: &str, len: usize| super::super::query::Query {
            id: id.into(),
            defline: id.into(),
            seq: vec![0; len],
        };
        let set = QuerySet {
            queries: vec![q("a", 10), q("long", 1_000_001)],
        };
        let e = check_query_lengths(Task::Blastn, &set).unwrap_err();
        assert!(
            e.to_string()
                .starts_with("query long is 1000001 bases, longer than the 1000000"),
            "{e}"
        );
        assert!(check_query_lengths(Task::BlastnShort, &set).is_err());
        // megablast's chunk is 5,000,000
        check_query_lengths(Task::Megablast, &set).unwrap();
        let set = QuerySet {
            queries: vec![q("m", 5_000_000)],
        };
        check_query_lengths(Task::Megablast, &set).unwrap();
        let set = QuerySet {
            queries: vec![q("m", 5_000_001)],
        };
        assert!(check_query_lengths(Task::Megablast, &set).is_err());
    }

    #[test]
    fn int4_cast_is_cvttsd2si() {
        assert_eq!(c_int4(5147.9), 5147);
        assert_eq!(c_int4(-5147.9), -5147);
        assert_eq!(c_int4(2_147_483_647.5), i32::MAX);
        assert_eq!(c_int4(2_147_483_648.0), i32::MIN);
        assert_eq!(c_int4(5.3e9), i32::MIN);
        assert_eq!(c_int4(f64::NAN), i32::MIN);
        assert_eq!(c_int4(f64::INFINITY), i32::MIN);
        assert_eq!(c_int4(-2_147_483_648.0), i32::MIN);
        assert_eq!(c_int4(-3e9), i32::MIN);
    }

    #[test]
    fn mixer_first_batch_hg38_and_small_db() {
        // hg38: 3,088,286,401 / 3000 = 1,029,428 target hits -> 5,147 letters
        let mut m = BatchSizeMixer::new(1_000_000 - 1000);
        m.set_target_hits((3_088_286_401i64 / 3000) as i32);
        assert_eq!(m.target_hits(), 1_029_428);
        assert_eq!(m.batch_size(-1), 5147);
        assert_eq!(m.ratio(), -1.0);
        // a database under 3 Gbp keeps the 1,000,000 minimum -> 5,000 letters
        let mut m = BatchSizeMixer::new(5_000_000 - 1000);
        m.set_target_hits((248_956_422i64 / 3000) as i32);
        assert_eq!(m.target_hits(), 1_000_000);
        assert_eq!(m.batch_size(-1), 5000);
        // no SetTargetHits (total length 0): the constructor's 5,000
        let mut m = BatchSizeMixer::new(999_000);
        assert_eq!(m.batch_size(-1), 5000);
    }

    #[test]
    fn mixer_mixes_ratios_and_clamps() {
        let mut m = BatchSizeMixer::new(999_000);
        m.set_target_hits(1_029_428);
        assert_eq!(m.batch_size(-1), 5147);
        // first batch: ratio taken as is
        let r1 = 1.0 * (750_001 + 1) as f64 / 5147.0;
        let b1 = (1_029_428.0 / r1) as i32;
        assert_eq!(m.batch_size(750_001), b1);
        assert_eq!(m.ratio(), r1);
        // second batch: 0.3 / 0.7 mixing against the previous requested size
        let r2 = 0.3 * (1.0 * (20_000 + 1) as f64 / b1 as f64) + (1.0 - 0.3) * r1;
        let b2 = (1_029_428.0 / r2) as i32;
        assert_eq!(m.batch_size(20_000), b2);
        assert_eq!(m.ratio(), r2);
        // clamp to the maximum resets the ratio
        let mut m = BatchSizeMixer::new(999_000);
        m.set_target_hits(1_000_000);
        m.batch_size(-1);
        assert_eq!(m.batch_size(4), 999_000); // 1e6 / (5 / 5000) = 1e9
        assert_eq!(m.ratio(), -1.0);
        // clamp to the minimum resets the ratio
        let mut m = BatchSizeMixer::new(999_000);
        m.set_target_hits(1_000_000);
        m.batch_size(-1);
        assert_eq!(m.batch_size(2_000_000_000), 100); // 1e6 / 400,000 = 2
        assert_eq!(m.ratio(), -1.0);
        // after a reset the next ratio is taken as is again: (99 + 1) / 100 = 1 -> 1e6, which
        // the maximum clamp cuts to 999,000 (and resets); then 100 / 999,000 overflows -> 100
        assert_eq!(m.batch_size(99), 999_000);
        assert_eq!(m.ratio(), -1.0);
        assert_eq!(m.batch_size(99), 100);
    }

    #[test]
    fn mixer_zero_hits_overflows_to_the_minimum_batch() {
        // 1e6 / (1 / 5000) = 5e9 does not fit an Int4: cvttsd2si gives INT_MIN, which the
        // minimum clamp turns into 100 (not the maximum), and the ratio resets
        let mut m = BatchSizeMixer::new(999_000);
        m.set_target_hits(1_000_000);
        m.batch_size(-1);
        assert_eq!(m.batch_size(0), 100);
        assert_eq!(m.ratio(), -1.0);
        // one hit on the hg38 first batch: 1,029,428 / (2 / 5147) = 2.65e9, same overflow
        let mut m = BatchSizeMixer::new(999_000);
        m.set_target_hits(1_029_428);
        m.batch_size(-1);
        assert_eq!(m.batch_size(1), 100);
        // two hits: 1,029,428 / (3 / 5147) = 1.77e9 fits, then the maximum clamp
        let mut m = BatchSizeMixer::new(999_000);
        m.set_target_hits(1_029_428);
        m.batch_size(-1);
        assert_eq!(m.batch_size(2), 999_000);
    }

    #[test]
    fn batch_takes_whole_queries_until_the_size_is_reached() {
        // the hg38 query set against 5,147: q01..q06 (970 + 756 + 1089 + 330 + 1861 = 5006
        // < 5147, + 1170 = 6176)
        let l = [970, 756, 1089, 330, 1861, 1170, 1169, 628, 1490, 1458];
        assert_eq!(batch_end(&l, 0, 5147), 6);
        assert_eq!(batch_end(&l, 6, 5147), 10);
        assert_eq!(batch_end(&l, 6, 1000), 7);
        assert_eq!(batch_end(&l, 7, 100), 8);
        // always at least one query, even when it alone exceeds the size
        assert_eq!(batch_end(&l, 4, 100), 5);
        assert_eq!(batch_end(&l, 10, 100), 10);
        // exactly reaching the size stops
        assert_eq!(batch_end(&[100, 100, 100], 0, 200), 2);
    }

    #[test]
    fn sizing_fixed_or_adaptive() {
        let mut p = SearchParams {
            task: Task::Blastn,
            ..Default::default()
        };
        let (mut s, first) = BatchSizing::new(&p, 3_088_286_401);
        assert_eq!(first, 5147);
        assert_eq!(s.next(0), 100);
        p.query_batch_size = Some(100_000);
        let (mut s, first) = BatchSizing::new(&p, 3_088_286_401);
        assert_eq!(first, 100_000);
        assert_eq!(s.next(0), 100_000);
        p.query_batch_size = Some(0);
        let (_, first) = BatchSizing::new(&p, 0);
        assert_eq!(first, 5000);
    }
}

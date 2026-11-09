//! Parallel map with in-order delivery and a bound on bytes held in flight, and the build's
//! work items: batches of consecutive OIDs.

use crate::error::Error;
use std::collections::BTreeMap;
use std::ops::Range;
use std::sync::{Condvar, Mutex, mpsc};

/// Most subjects in one batch.
pub const BATCH_MAX_SUBJECTS: usize = 4096;
/// Default stored bytes in one batch (`--batch-mb 8`).
pub const DEFAULT_BATCH_BYTES: u64 = 8 << 20;

/// A run of consecutive OIDs and their stored bytes: one work item of the build.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OidBatch {
    pub oids: Range<u32>,
    pub bytes: u64,
}

/// Group subjects, given their stored sizes in OID order, into batches of consecutive OIDs.
/// A batch closes before a subject that would take it past `max_bytes` or `max_subjects`; a
/// subject larger than `max_bytes` is a batch of its own. OIDs are numbered from 0.
pub fn plan_batches(
    sizes: impl IntoIterator<Item = u64>,
    max_bytes: u64,
    max_subjects: usize,
) -> Vec<OidBatch> {
    let max_subjects = max_subjects.max(1) as u32;
    let mut out = Vec::new();
    let (mut start, mut bytes, mut oid) = (0u32, 0u64, 0u32);
    for size in sizes {
        let n = oid - start;
        if n > 0 && (n >= max_subjects || bytes + size > max_bytes) {
            out.push(OidBatch {
                oids: start..oid,
                bytes,
            });
            (start, bytes) = (oid, 0);
        }
        bytes += size;
        oid += 1;
    }
    if oid > start {
        out.push(OidBatch {
            oids: start..oid,
            bytes,
        });
    }
    out
}

/// The subject cap for `n` subjects on `workers` workers: [`BATCH_MAX_SUBJECTS`], lowered so
/// a small selection still gives each worker about four batches.
pub fn batch_subject_cap(n: usize, workers: usize) -> usize {
    n.div_ceil(4 * workers.max(1)).clamp(1, BATCH_MAX_SUBJECTS)
}

struct Gate {
    next_write: usize,
    inflight: u64,
    failed: bool,
}

/// Run `work` on every item using `threads` workers; call `sink` with results in item
/// order, on the calling thread. A worker only starts item i when it is the next one to
/// be written or when the summed `cost` of started-but-unwritten items stays within
/// `budget`. The first error (from `work` or `sink`) stops the workers and is returned.
pub fn run_ordered<I: Sync, O: Send>(
    items: &[I],
    threads: usize,
    budget: u64,
    cost: impl Fn(&I) -> u64 + Sync,
    work: impl Fn(&I) -> Result<O, Error> + Sync,
    sink: impl FnMut(usize, O) -> Result<(), Error>,
) -> Result<(), Error> {
    let take = |i: usize| items.get(i).map(Ok);
    run_ordered_with(
        take,
        threads,
        budget,
        |i| cost(i),
        || (),
        |_, i| work(i),
        sink,
    )
}

/// [`run_ordered`] with per-worker state: each worker calls `init` once and hands the state
/// to every `work` call it makes (a reused decode buffer, for example). Item i is `take(i)`,
/// called under a lock in index order, so it may read a stream; `None` ends the items.
pub fn run_ordered_with<I: Send, O: Send, S>(
    take: impl FnMut(usize) -> Option<Result<I, Error>> + Send,
    threads: usize,
    budget: u64,
    cost: impl Fn(&I) -> u64 + Sync,
    init: impl Fn() -> S + Sync,
    work: impl Fn(&mut S, I) -> Result<O, Error> + Sync,
    mut sink: impl FnMut(usize, O) -> Result<(), Error>,
) -> Result<(), Error> {
    // the next index to take, and `take` until the items end or one fails
    let source = Mutex::new((0usize, Some(take)));
    let gate = (
        Mutex::new(Gate {
            next_write: 0,
            inflight: 0,
            failed: false,
        }),
        Condvar::new(),
    );
    let (tx, rx) = mpsc::channel::<(usize, u64, Result<O, Error>)>();
    std::thread::scope(|s| {
        for _ in 0..threads.max(1) {
            let tx = tx.clone();
            let (source, gate, init, work, cost) = (&source, &gate, &init, &work, &cost);
            s.spawn(move || {
                let mut state = init();
                loop {
                    if gate.0.lock().unwrap().failed {
                        break;
                    }
                    let (i, item) = {
                        let mut src = source.lock().unwrap();
                        let i = src.0;
                        match src.1.as_mut().and_then(|take| take(i)) {
                            Some(Ok(item)) => {
                                src.0 += 1;
                                (i, item)
                            }
                            end => {
                                // the sink returns the error when it reaches item i
                                src.1 = None;
                                if let Some(Err(e)) = end {
                                    let _ = tx.send((i, 0, Err(e)));
                                }
                                break;
                            }
                        }
                    };
                    let c = cost(&item);
                    {
                        let mut g = gate.0.lock().unwrap();
                        while !(g.failed || i == g.next_write || g.inflight + c <= budget) {
                            g = gate.1.wait(g).unwrap();
                        }
                        if g.failed {
                            break;
                        }
                        g.inflight += c;
                    }
                    let r = work(&mut state, item);
                    if tx.send((i, c, r)).is_err() {
                        break;
                    }
                }
            });
        }
        drop(tx);
        let mut pending: BTreeMap<usize, (u64, Result<O, Error>)> = BTreeMap::new();
        let mut want = 0usize;
        let res: Result<(), Error> = (|| {
            for (i, c, r) in rx.iter() {
                pending.insert(i, (c, r));
                while let Some((c, r)) = pending.remove(&want) {
                    sink(want, r?)?;
                    let mut g = gate.0.lock().unwrap();
                    g.inflight -= c;
                    g.next_write = want + 1;
                    gate.1.notify_all();
                    drop(g);
                    want += 1;
                }
            }
            Ok(())
        })();
        if res.is_err() {
            gate.0.lock().unwrap().failed = true;
            gate.1.notify_all();
        }
        res
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};
    use std::time::Duration;

    fn delay(i: u64) -> Duration {
        // deterministic pseudo-random delays, 0..3 ms
        let x = i.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 60;
        Duration::from_micros(x * 200)
    }

    #[test]
    fn plan_batches_cuts_by_bytes_and_count() {
        let b = |v: &[u64], bytes, n| -> Vec<(u32, u32, u64)> {
            plan_batches(v.iter().copied(), bytes, n)
                .into_iter()
                .map(|b| (b.oids.start, b.oids.end, b.bytes))
                .collect()
        };
        assert!(b(&[], 10, 4).is_empty());
        assert_eq!(b(&[3, 3, 3, 3], 10, 4), vec![(0, 3, 9), (3, 4, 3)]);
        assert_eq!(b(&[1; 5], 100, 2), vec![(0, 2, 2), (2, 4, 2), (4, 5, 1)]);
        // a subject larger than the byte cap is a batch of its own
        assert_eq!(
            b(&[1, 50, 1], 10, 4),
            vec![(0, 1, 1), (1, 2, 50), (2, 3, 1)]
        );
        assert_eq!(b(&[2, 2], 1, 4), vec![(0, 1, 2), (1, 2, 2)]);
        // batches cover every OID once, in order
        let sizes: Vec<u64> = (0..10_000u64).map(|i| (i * 7919) % 5000).collect();
        let v = plan_batches(sizes.iter().copied(), 8 << 10, BATCH_MAX_SUBJECTS);
        let mut next = 0;
        for x in &v {
            assert_eq!(x.oids.start, next);
            assert!(x.oids.end > x.oids.start);
            let sum: u64 = sizes[x.oids.start as usize..x.oids.end as usize]
                .iter()
                .sum();
            assert_eq!(sum, x.bytes);
            assert!(x.bytes <= 8 << 10 || x.oids.len() == 1);
            next = x.oids.end;
        }
        assert_eq!(next, 10_000);
        assert_eq!(batch_subject_cap(0, 8), 1);
        assert_eq!(batch_subject_cap(30, 2), 4);
        assert_eq!(batch_subject_cap(384_354, 16), BATCH_MAX_SUBJECTS);
    }

    #[test]
    fn run_ordered_with_keeps_one_state_per_worker() {
        let items: Vec<u64> = (0..500).collect();
        let inits = AtomicU64::new(0);
        let mut got = Vec::new();
        run_ordered_with(
            |i| items.get(i).map(Ok),
            4,
            1 << 30,
            |_| 1,
            || {
                inits.fetch_add(1, Ordering::SeqCst);
                Vec::<u64>::new()
            },
            |buf, &i| {
                buf.push(i);
                std::thread::sleep(delay(i) / 10);
                Ok((i, buf.len()))
            },
            |idx, (i, _)| {
                assert_eq!(idx as u64, i);
                got.push(i);
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(got, items);
        assert!(inits.load(Ordering::SeqCst) <= 4);
    }

    #[test]
    fn run_ordered_with_returns_a_take_error_in_order() {
        let mut got = Vec::new();
        let r = run_ordered_with(
            |i| match i {
                0..5 => Some(Ok(i)),
                5 => Some(Err(Error::sequence("d", "unreadable"))),
                _ => panic!("took item {i} after the error"),
            },
            3,
            1 << 30,
            |_| 1,
            || (),
            |_, i| Ok(i),
            |_, i| {
                got.push(i);
                Ok(())
            },
        );
        assert!(matches!(r, Err(Error::Sequence { .. })));
        assert_eq!(got, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn run_ordered_preserves_order_under_random_delays() {
        let items: Vec<u64> = (0..200).collect();
        let mut got = Vec::new();
        run_ordered(
            &items,
            8,
            1 << 30,
            |_| 1,
            |&i| {
                std::thread::sleep(delay(i));
                Ok(i * 2)
            },
            |idx, o| {
                assert_eq!(o, idx as u64 * 2);
                got.push(o);
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(got, (0..200).map(|i| i * 2).collect::<Vec<_>>());
    }

    #[test]
    fn run_ordered_respects_budget() {
        let items: Vec<u64> = (0..100).collect();
        let inflight = AtomicU64::new(0);
        let peak = AtomicU64::new(0);
        run_ordered(
            &items,
            8,
            30,
            |_| 10,
            |&i| {
                let now = inflight.fetch_add(10, Ordering::SeqCst) + 10;
                peak.fetch_max(now, Ordering::SeqCst);
                std::thread::sleep(delay(i));
                Ok(i)
            },
            |_, _| {
                inflight.fetch_sub(10, Ordering::SeqCst);
                Ok(())
            },
        )
        .unwrap();
        // At most budget/cost items in flight, plus the one allowed because it is next.
        assert!(
            peak.load(Ordering::SeqCst) <= 40,
            "peak {}",
            peak.load(Ordering::SeqCst)
        );
    }

    #[test]
    fn run_ordered_stops_on_error() {
        let items: Vec<u64> = (0..10_000).collect();
        let started = AtomicU64::new(0);
        let mut sunk = 0;
        let r = run_ordered(
            &items,
            4,
            40,
            |_| 10,
            |&i| {
                started.fetch_add(1, Ordering::SeqCst);
                if i == 5 {
                    Err(Error::sequence("x", "boom"))
                } else {
                    Ok(i)
                }
            },
            |_, _| {
                sunk += 1;
                Ok(())
            },
        );
        assert!(matches!(r, Err(Error::Sequence { .. })));
        assert_eq!(sunk, 5);
        assert!(started.load(Ordering::SeqCst) < 100);

        // An error from the sink stops it too.
        let r = run_ordered(
            &items,
            4,
            40,
            |_| 10,
            |&i| Ok(i),
            |i, _| {
                if i == 3 {
                    Err(Error::Cancelled)
                } else {
                    Ok(())
                }
            },
        );
        assert!(matches!(r, Err(Error::Cancelled)));
    }
}

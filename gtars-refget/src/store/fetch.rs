//! Parallel HTTP fetching for remote RefgetStores.
//!
//! One [`Fetcher`] per remote store. It keeps a pool of keep-alive
//! connections (ureq, blocking, on plain threads), retries transient failures
//! with jittered exponential backoff, and counts everything it does in
//! [`FetchCounters`]. Metadata and small sequence files come down as one
//! whole-file GET each; a prefetch splits sequence files bigger than
//! `CHUNK_BYTES` (8 MiB) into ranged GETs that run in parallel.
//!
//! S3 serves 0.7-2 MB/s per connection, so throughput comes from many
//! connections at once, not from any single one.

use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Duration;

use serde::Serialize;

// ============================================================================
// Public configuration and counters
// ============================================================================

/// How a remote store downloads files.
#[derive(Debug, Clone)]
pub struct FetchConfig {
    /// Concurrent HTTP requests. Default 32.
    pub jobs: usize,
    /// Extra attempts per request after the first. Default 5.
    pub retries: u32,
    /// Hash every whole sequence fetched over the network and compare it with
    /// its stored digest. Default true.
    pub verify: bool,
}

impl Default for FetchConfig {
    fn default() -> Self {
        Self {
            jobs: 32,
            retries: 5,
            verify: true,
        }
    }
}

/// What a fetcher or a prefetch run did.
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct FetchCounters {
    /// HTTP requests sent, retries included.
    pub requests: u64,
    /// Requests that were retries of an earlier attempt.
    pub retries: u64,
    /// Body bytes received, failed attempts included.
    pub bytes_downloaded: u64,
    /// Whole sequences that came over the network.
    pub sequences_fetched: u64,
    /// Sequences served from the shared cache.
    pub cache_hits: u64,
    /// Bytes of those cache hits.
    pub cache_hit_bytes: u64,
    /// Sequences read from a local store's files.
    pub local_reads: u64,
    /// Bytes of those local reads.
    pub local_read_bytes: u64,
    /// Requests answered "not found" (HTTP 404, 410 or 403).
    pub missing: u64,
    /// Wall time covered by these counters.
    pub elapsed_secs: f64,
}

impl FetchCounters {
    /// Download rate in MB/s (10^6 bytes) over `elapsed_secs`.
    pub fn rate_mb_s(&self) -> f64 {
        if self.elapsed_secs > 0.0 {
            self.bytes_downloaded as f64 / 1e6 / self.elapsed_secs
        } else {
            0.0
        }
    }

    /// One summary line, e.g.
    /// `network: 1362000000 bytes in 346 requests (0 retries), 27.6 MB/s; cache hits 0 (0.00 GB)`.
    pub fn network_line(&self) -> String {
        format!(
            "network: {} bytes in {} requests ({} {}), {:.1} MB/s; cache hits {} ({:.2} GB)",
            self.bytes_downloaded,
            self.requests,
            self.retries,
            if self.retries == 1 {
                "retry"
            } else {
                "retries"
            },
            self.rate_mb_s(),
            self.cache_hits,
            self.cache_hit_bytes as f64 / 1e9
        )
    }
}

/// A required remote file is not there or may not be read: HTTP 404, 410 or
/// 403 (S3 answers 403 for a missing key when listing is denied). Carried
/// inside the `anyhow::Error`, so callers can tell "missing" from other
/// failures with `err.downcast_ref::<RemoteMissing>()`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RemoteMissing {
    pub url: String,
    pub status: u16,
}

impl std::fmt::Display for RemoteMissing {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let what = if self.status == 403 {
            "not found or access denied"
        } else {
            "not found"
        };
        write!(f, "{}: {} (HTTP {})", self.url, what, self.status)
    }
}

impl std::error::Error for RemoteMissing {}

/// The lock-free counters behind a [`FetchCounters`] snapshot. A fetcher
/// keeps one for its lifetime, a prefetch one for its run.
#[derive(Debug, Default)]
pub(crate) struct AtomicCounters {
    pub(crate) requests: AtomicU64,
    pub(crate) retries: AtomicU64,
    pub(crate) bytes_downloaded: AtomicU64,
    pub(crate) sequences_fetched: AtomicU64,
    pub(crate) cache_hits: AtomicU64,
    pub(crate) cache_hit_bytes: AtomicU64,
    pub(crate) local_reads: AtomicU64,
    pub(crate) local_read_bytes: AtomicU64,
    pub(crate) missing: AtomicU64,
}

/// Picks one counter of an [`AtomicCounters`].
pub(crate) type Field = fn(&AtomicCounters) -> &AtomicU64;

impl AtomicCounters {
    pub(crate) fn add(&self, field: Field, n: u64) {
        field(self).fetch_add(n, Ordering::Relaxed);
    }

    pub(crate) fn snapshot(&self, elapsed: Duration) -> FetchCounters {
        let g = |a: &AtomicU64| a.load(Ordering::Relaxed);
        FetchCounters {
            requests: g(&self.requests),
            retries: g(&self.retries),
            bytes_downloaded: g(&self.bytes_downloaded),
            sequences_fetched: g(&self.sequences_fetched),
            cache_hits: g(&self.cache_hits),
            cache_hit_bytes: g(&self.cache_hit_bytes),
            local_reads: g(&self.local_reads),
            local_read_bytes: g(&self.local_read_bytes),
            missing: g(&self.missing),
            elapsed_secs: elapsed.as_secs_f64(),
        }
    }
}

// ============================================================================
// Fetcher
// ============================================================================

#[cfg(feature = "http")]
pub(crate) use client::{Fetcher, Got};

#[cfg(feature = "http")]
mod client {
    use super::*;
    use crate::store::readonly::ReadonlyRefgetStore;
    use anyhow::{Result, anyhow};
    use std::io::Read;
    use std::sync::atomic::{AtomicBool, AtomicUsize};
    use std::sync::{Arc, Mutex};
    use std::time::Instant;

    /// Sequence files bigger than this download as parallel ranged GETs of this
    /// size; a remote stream reads this many bytes per request.
    const CHUNK_BYTES: u64 = 8 << 20;
    /// First retry wait, doubled per attempt up to [`MAX_BACKOFF`], with +/-50% jitter.
    pub(super) const BACKOFF: Duration = Duration::from_millis(500);
    /// Longest wait between two attempts of one request (before jitter).
    pub(super) const MAX_BACKOFF: Duration = Duration::from_secs(8);
    /// Longest `Retry-After` the fetcher honors.
    pub(super) const MAX_RETRY_AFTER: Duration = Duration::from_secs(30);
    /// TCP/TLS connect timeout.
    const CONNECT_TIMEOUT: Duration = Duration::from_secs(20);
    /// Timeout for the response headers and for each body read.
    const READ_TIMEOUT: Duration = Duration::from_secs(60);

    /// The answer to a GET: the body (the whole file, or the requested range),
    /// or "not there".
    #[derive(Debug)]
    pub(crate) enum Got {
        Found(Vec<u8>),
        Missing { status: u16 },
    }

    /// What one attempt produced.
    enum Step {
        Done(Got),
        Retry {
            why: String,
            retry_after: Option<Duration>,
        },
    }

    fn retry(why: String) -> Step {
        Step::Retry {
            why,
            retry_after: None,
        }
    }

    /// One logical GET.
    struct Req<'a> {
        url: &'a str,
        /// `[start, end)` byte range, sent as `Range: bytes=start-(end-1)`.
        range: Option<(u64, u64)>,
        /// Expected body length, when known.
        len: Option<u64>,
    }

    /// Downloads files from one remote store.
    pub(crate) struct Fetcher {
        /// Normalized store URL with a trailing `/`.
        base: String,
        agent: ureq::Agent,
        cfg: FetchConfig,
        /// Ranged-download piece size: [`CHUNK_BYTES`] outside tests.
        chunk: u64,
        stats: AtomicCounters,
        created: Instant,
        jitter: AtomicU64,
    }

    impl std::fmt::Debug for Fetcher {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            f.debug_struct("Fetcher")
                .field("base", &self.base)
                .field("cfg", &self.cfg)
                .field("counters", &self.counters())
                .finish()
        }
    }

    impl Fetcher {
        pub(crate) fn new(base: &str, cfg: FetchConfig) -> Result<Self> {
            let mut cfg = cfg;
            cfg.jobs = cfg.jobs.max(1);
            let conns = cfg.jobs * 2;
            let agent: ureq::Agent = ureq::Agent::config_builder()
                .http_status_as_error(false)
                .accept_encoding("identity")
                .max_idle_connections(conns)
                .max_idle_connections_per_host(conns)
                .timeout_connect(Some(CONNECT_TIMEOUT))
                .timeout_recv_response(Some(READ_TIMEOUT))
                .timeout_per_read(Some(READ_TIMEOUT))
                .build()
                .into();
            let seed = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_nanos() as u64)
                .unwrap_or(0x9E37_79B9_7F4A_7C15);
            Ok(Self {
                base: format!("{}/", base.trim_end_matches('/')),
                agent,
                cfg,
                chunk: CHUNK_BYTES,
                stats: AtomicCounters::default(),
                created: Instant::now(),
                jitter: AtomicU64::new(seed | 1),
            })
        }

        /// The same fetcher with a smaller ranged-download piece, so tests
        /// can split tiny files.
        #[cfg(test)]
        pub(crate) fn with_chunk(mut self, chunk: u64) -> Self {
            self.chunk = chunk;
            self
        }

        pub(crate) fn config(&self) -> &FetchConfig {
            &self.cfg
        }

        pub(crate) fn chunk(&self) -> u64 {
            self.chunk
        }

        pub(crate) fn stats(&self) -> &AtomicCounters {
            &self.stats
        }

        /// Cumulative counters since this fetcher was created.
        pub(crate) fn counters(&self) -> FetchCounters {
            self.stats.snapshot(self.created.elapsed())
        }

        /// Full URL of a store-relative path, after checking it for traversal.
        pub(crate) fn url(&self, rel: &str) -> Result<String> {
            ReadonlyRefgetStore::sanitize_relative_path(rel)?;
            Ok(format!("{}{}", self.base, rel.trim_start_matches('/')))
        }

        /// Count into this fetcher's counters and, when given, a run's.
        fn count(&self, sink: Option<&AtomicCounters>, field: Field, n: u64) {
            self.stats.add(field, n);
            if let Some(s) = sink {
                s.add(field, n);
            }
        }

        /// Wait before attempt number `attempt` (1-based count of retries).
        pub(super) fn backoff(&self, attempt: u32, retry_after: Option<Duration>) -> Duration {
            if let Some(ra) = retry_after {
                return ra.min(MAX_RETRY_AFTER);
            }
            let base = BACKOFF
                .saturating_mul(1u32 << (attempt - 1).min(16))
                .min(MAX_BACKOFF);
            // xorshift64* for +/-50% jitter, so 32 workers don't retry in lockstep.
            let mut x = self.jitter.load(Ordering::Relaxed);
            x ^= x >> 12;
            x ^= x << 25;
            x ^= x >> 27;
            self.jitter.store(x, Ordering::Relaxed);
            let r = (x.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 11) as f64 / (1u64 << 53) as f64;
            base.mul_f64(0.5 + r)
        }

        // ------------------------------------------------------------------------
        // Requests
        // ------------------------------------------------------------------------

        /// One logical GET, tried up to `attempts` times. Counts every attempt,
        /// retry and body byte in the fetcher's counters (and in `sink`).
        fn request(
            &self,
            req: &Req<'_>,
            sink: Option<&AtomicCounters>,
            attempts: u32,
        ) -> Result<Got> {
            let mut last = String::new();
            let mut retry_after = None;
            for attempt in 0..attempts {
                if attempt > 0 {
                    self.count(sink, |s| &s.retries, 1);
                    std::thread::sleep(self.backoff(attempt, retry_after));
                }
                self.count(sink, |s| &s.requests, 1);
                match self.attempt(req, sink)? {
                    Step::Done(got) => return Ok(got),
                    Step::Retry {
                        why,
                        retry_after: ra,
                    } => {
                        last = why;
                        retry_after = ra;
                    }
                }
            }
            Err(anyhow!(
                "{}: failed after {} attempt{}: {}",
                req.url,
                attempts,
                if attempts == 1 { "" } else { "s" },
                last
            ))
        }

        fn attempt(&self, req: &Req<'_>, sink: Option<&AtomicCounters>) -> Result<Step> {
            let mut rb = self.agent.get(req.url);
            if let Some((start, end)) = req.range {
                rb = rb.header("Range", format!("bytes={}-{}", start, end - 1));
            }
            let resp = match rb.call() {
                Ok(r) => r,
                Err(
                    e @ (ureq::Error::BadUri(_)
                    | ureq::Error::InvalidProxyUrl
                    | ureq::Error::Http(_)),
                ) => return Err(anyhow!("{}: {}", req.url, e)),
                Err(e) => return Ok(retry(e.to_string())),
            };
            let status = resp.status().as_u16();
            let header = |name: &str| {
                resp.headers()
                    .get(name)
                    .and_then(|v| v.to_str().ok())
                    .map(|s| s.trim().to_string())
            };
            match (status, req.range) {
                (403 | 404 | 410, _) => {
                    self.count(sink, |s| &s.missing, 1);
                    return Ok(Step::Done(Got::Missing { status }));
                }
                (408 | 429 | 500 | 502 | 503 | 504, _) => {
                    return Ok(Step::Retry {
                        why: format!("HTTP {}", status),
                        retry_after: header("Retry-After")
                            .and_then(|v| v.parse::<u64>().ok())
                            .map(Duration::from_secs),
                    });
                }
                (200, None) => {}
                (200, Some(_)) => {
                    return Err(anyhow!(
                        "{}: the server ignored the Range header (HTTP 200); \
                         remote stores need a server that supports range requests",
                        req.url
                    ));
                }
                (206, Some((start, end))) => {
                    let cr = header("Content-Range").and_then(|v| parse_content_range(&v));
                    if cr.is_none_or(|(s, e, _)| s != start || e + 1 != end) {
                        return Ok(retry(format!(
                            "Content-Range {:?} does not match requested bytes {}-{}",
                            header("Content-Range").unwrap_or_default(),
                            start,
                            end - 1
                        )));
                    }
                }
                _ => return Err(anyhow!("{}: HTTP {}", req.url, status)),
            }
            let content_length = header("Content-Length").and_then(|v| v.parse::<u64>().ok());
            let cap = req.len.or(content_length).unwrap_or(64 * 1024);
            let mut body = Vec::with_capacity(cap.min(1 << 30) as usize);
            let read = resp.into_body().into_reader().read_to_end(&mut body);
            let got = body.len() as u64;
            self.count(sink, |s| &s.bytes_downloaded, got);
            if let Err(e) = read {
                return Ok(retry(format!(
                    "body read failed after {} bytes: {}",
                    got, e
                )));
            }
            if let Some(n) = content_length
                && got != n
            {
                return Ok(retry(format!("got {} bytes, Content-Length {}", got, n)));
            }
            if let Some(n) = req.len
                && got != n
            {
                if req.range.is_none() && content_length.is_some() {
                    return Err(anyhow!(
                        "{}: remote file is {} bytes, expected {}",
                        req.url,
                        got,
                        n
                    ));
                }
                return Ok(retry(format!("got {} bytes, expected {}", got, n)));
            }
            Ok(Step::Done(Got::Found(body)))
        }

        // ------------------------------------------------------------------------
        // Crate operations
        // ------------------------------------------------------------------------

        /// GET a whole file (metadata).
        pub(crate) fn get(&self, rel: &str) -> Result<Got> {
            let url = self.url(rel)?;
            let req = Req {
                url: &url,
                range: None,
                len: None,
            };
            self.request(&req, None, 1 + self.cfg.retries)
        }

        /// One whole-file GET with no retries, for a quick check that has a
        /// fallback (the manifest of a store that is cached).
        pub(crate) fn get_quick(&self, rel: &str) -> Result<Got> {
            let url = self.url(rel)?;
            let req = Req {
                url: &url,
                range: None,
                len: None,
            };
            self.request(&req, None, 1)
        }

        /// GET a whole sequence file, `len` bytes when known. A missing file
        /// is an error carrying [`RemoteMissing`].
        pub(crate) fn get_file(
            &self,
            rel: &str,
            len: Option<u64>,
            sink: Option<&AtomicCounters>,
        ) -> Result<Vec<u8>> {
            self.get_bytes(rel, None, len, sink)
        }

        /// GET bytes `[start, end)` of a file. A server that ignores `Range`
        /// is refused; a missing file is an error carrying [`RemoteMissing`].
        pub(crate) fn get_range(&self, rel: &str, start: u64, end: u64) -> Result<Vec<u8>> {
            self.get_range_with(rel, start, end, None)
        }

        pub(crate) fn get_range_with(
            &self,
            rel: &str,
            start: u64,
            end: u64,
            sink: Option<&AtomicCounters>,
        ) -> Result<Vec<u8>> {
            if start >= end {
                return Ok(Vec::new());
            }
            self.get_bytes(rel, Some((start, end)), Some(end - start), sink)
        }

        fn get_bytes(
            &self,
            rel: &str,
            range: Option<(u64, u64)>,
            len: Option<u64>,
            sink: Option<&AtomicCounters>,
        ) -> Result<Vec<u8>> {
            let url = self.url(rel)?;
            let req = Req {
                url: &url,
                range,
                len,
            };
            match self.request(&req, sink, 1 + self.cfg.retries)? {
                Got::Found(b) => Ok(b),
                Got::Missing { status } => Err(RemoteMissing { url, status }.into()),
            }
        }

        /// Stream bytes `[start, end)` of a file, one ranged GET of `chunk`
        /// bytes at a time. The first piece comes down here, so a missing file
        /// or a server that ignores `Range` fails now.
        pub(crate) fn open_range(
            self: &Arc<Self>,
            rel: &str,
            start: u64,
            end: u64,
        ) -> Result<Box<dyn Read + Send>> {
            let mut reader = RangeReader {
                fetcher: Arc::clone(self),
                rel: rel.to_string(),
                pos: start,
                end,
                buf: Vec::new(),
                at: 0,
            };
            reader.fill()?;
            Ok(Box::new(reader))
        }

        /// GET many whole files in parallel. Results keep input order.
        pub(crate) fn get_many(&self, rels: &[String]) -> Result<Vec<Got>> {
            let results: Vec<Mutex<Option<Result<Got>>>> =
                rels.iter().map(|_| Mutex::new(None)).collect();
            let next = AtomicUsize::new(0);
            let failed = AtomicBool::new(false);
            std::thread::scope(|scope| {
                for _ in 0..self.cfg.jobs.min(rels.len()) {
                    scope.spawn(|| {
                        while !failed.load(Ordering::Relaxed) {
                            let i = next.fetch_add(1, Ordering::Relaxed);
                            if i >= rels.len() {
                                break;
                            }
                            let r = self.get(&rels[i]);
                            if r.is_err() {
                                failed.store(true, Ordering::Relaxed);
                            }
                            *results[i].lock().unwrap() = Some(r);
                        }
                    });
                }
            });
            let mut out = Vec::with_capacity(rels.len());
            for (i, slot) in results.into_iter().enumerate() {
                match slot.into_inner().unwrap() {
                    Some(r) => out.push(r?),
                    // Skipped after another request failed; that error comes first.
                    None => return Err(anyhow!("download of {} was cancelled", rels[i])),
                }
            }
            Ok(out)
        }
    }

    /// Parse `bytes START-END/TOTAL` (TOTAL may be `*`).
    pub(super) fn parse_content_range(v: &str) -> Option<(u64, u64, Option<u64>)> {
        let rest = v.strip_prefix("bytes")?.trim_start();
        let (range, total) = rest.split_once('/')?;
        let (s, e) = range.trim().split_once('-')?;
        let s: u64 = s.trim().parse().ok()?;
        let e: u64 = e.trim().parse().ok()?;
        if e < s {
            return None;
        }
        let total = match total.trim() {
            "*" => None,
            t => Some(t.parse().ok()?),
        };
        Some((s, e, total))
    }

    /// Reads `[pos, end)` of one remote file a piece at a time, each piece a
    /// retried ranged GET.
    struct RangeReader {
        fetcher: Arc<Fetcher>,
        rel: String,
        pos: u64,
        end: u64,
        buf: Vec<u8>,
        at: usize,
    }

    impl RangeReader {
        /// Fetch the next piece once the current one is used up.
        fn fill(&mut self) -> Result<()> {
            if self.at == self.buf.len() && self.pos < self.end {
                let to = self.end.min(self.pos + self.fetcher.chunk);
                self.buf = self.fetcher.get_range(&self.rel, self.pos, to)?;
                (self.pos, self.at) = (to, 0);
            }
            Ok(())
        }
    }

    impl Read for RangeReader {
        fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
            self.fill().map_err(std::io::Error::other)?;
            let n = out.len().min(self.buf.len() - self.at);
            out[..n].copy_from_slice(&self.buf[self.at..self.at + n]);
            self.at += n;
            Ok(n)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(feature = "http")]
    #[test]
    fn content_range_parses() {
        use super::client::parse_content_range;
        assert_eq!(
            parse_content_range("bytes 0-99/1000"),
            Some((0, 99, Some(1000)))
        );
        assert_eq!(parse_content_range("bytes 5-9/*"), Some((5, 9, None)));
        assert_eq!(parse_content_range("bytes 9-5/10"), None);
        assert_eq!(parse_content_range("items 0-1/2"), None);
    }

    #[cfg(feature = "http")]
    #[test]
    fn backoff_doubles_and_caps_with_jitter() {
        use super::client::{BACKOFF, MAX_BACKOFF, MAX_RETRY_AFTER};
        let f = Fetcher::new("http://127.0.0.1:1/x", FetchConfig::default()).unwrap();
        for attempt in 1..10 {
            let base = BACKOFF.saturating_mul(1 << (attempt - 1)).min(MAX_BACKOFF);
            let w = f.backoff(attempt, None);
            assert!(
                w >= base.mul_f64(0.5) && w <= base.mul_f64(1.5),
                "{attempt}: {w:?}"
            );
        }
        assert_eq!(f.backoff(1, Some(Duration::from_secs(90))), MAX_RETRY_AFTER);
    }

    #[test]
    fn missing_error_downcasts() {
        let e: anyhow::Error = RemoteMissing {
            url: "http://x/rgstore.json".into(),
            status: 403,
        }
        .into();
        let e = e.context("opening remote");
        let m = e.downcast_ref::<RemoteMissing>().unwrap();
        assert_eq!(m.status, 403);
        assert_eq!(
            m.to_string(),
            "http://x/rgstore.json: not found or access denied (HTTP 403)"
        );
    }
}

#[cfg(all(test, feature = "http"))]
mod http_tests {
    use super::*;
    use std::io::Read;
    use std::sync::Arc;
    use std::time::Instant;

    use crate::store::http_fixture::{FaultKind, Fixture};

    fn serve(files: &[(&str, &[u8])]) -> (tempfile::TempDir, Fixture) {
        let dir = tempfile::tempdir().unwrap();
        for (rel, body) in files {
            let p = dir.path().join(rel);
            std::fs::create_dir_all(p.parent().unwrap()).unwrap();
            std::fs::write(p, body).unwrap();
        }
        let fx = Fixture::serve(dir.path(), &["s"]);
        (dir, fx)
    }

    fn fetcher(fx: &Fixture) -> Fetcher {
        let cfg = FetchConfig {
            jobs: 4,
            ..FetchConfig::default()
        };
        Fetcher::new(&fx.url("s"), cfg).unwrap()
    }

    fn body(n: usize) -> Vec<u8> {
        (0..n).map(|i| (i * 7 % 251) as u8).collect()
    }

    #[test]
    fn ranged_reads_refuse_a_server_that_ignores_range() {
        let full = body(100);
        let (_d, fx) = serve(&[("f.seq", &full)]);
        fx.set_honor_range(false);
        let f = Arc::new(fetcher(&fx));
        for err in [
            f.get_range("f.seq", 10, 20).unwrap_err(),
            f.open_range("f.seq", 10, 20).err().unwrap(),
        ] {
            assert!(
                err.to_string().contains("ignored the Range header"),
                "{err}"
            );
        }
        // Whole-file GETs send no Range, so they still work.
        assert_eq!(f.get_file("f.seq", Some(100), None).unwrap(), full);
    }

    #[test]
    fn whole_files_are_one_get_each() {
        let full = body(100);
        let (_d, fx) = serve(&[
            ("meta.txt", &full),
            ("small.txt", b"abc"),
            ("empty.txt", b""),
        ]);
        let f = fetcher(&fx);
        assert!(matches!(f.get("meta.txt").unwrap(), Got::Found(b) if b == full));
        let reqs = fx.requests_matching("meta.txt");
        assert!(reqs.len() == 1 && reqs[0].range.is_none(), "{reqs:?}");
        assert!(matches!(f.get("small.txt").unwrap(), Got::Found(b) if b == b"abc"));
        assert!(matches!(f.get("empty.txt").unwrap(), Got::Found(b) if b.is_empty()));
        assert!(matches!(
            f.get("nope.txt").unwrap(),
            Got::Missing { status: 404 }
        ));
        // A wrong expected size is caught without retrying.
        let err = f
            .get_file("meta.txt", Some(99), None)
            .unwrap_err()
            .to_string();
        assert!(err.contains("100 bytes, expected 99"), "{err}");
        assert_eq!(fx.requests_matching("meta.txt").len(), 2);
    }

    #[test]
    fn ranges_return_exactly_the_requested_bytes() {
        let full = body(1000);
        let (_d, fx) = serve(&[("big.seq", &full)]);
        let f = Arc::new(fetcher(&fx));
        assert_eq!(f.get_range("big.seq", 100, 164).unwrap(), &full[100..164]);
        assert!(f.get_range("big.seq", 5, 5).unwrap().is_empty());
        let err = f.get_range("nope.seq", 0, 10).unwrap_err();
        assert_eq!(err.downcast_ref::<RemoteMissing>().unwrap().status, 404);
    }

    #[test]
    fn retries_transient_statuses_and_truncation() {
        let full = body(50);
        let (_d, fx) = serve(&[("a.seq", &full)]);
        let f = fetcher(&fx);
        fx.add_fault("a.seq", FaultKind::Status(503), 1);
        assert_eq!(f.get_file("a.seq", Some(50), None).unwrap(), full);
        assert_eq!(f.counters().retries, 1);
        fx.add_fault("a.seq", FaultKind::Truncate, 1);
        assert_eq!(f.get_file("a.seq", Some(50), None).unwrap(), full);
        assert_eq!(f.counters().retries, 2);
    }

    #[test]
    fn honors_retry_after() {
        let (_d, fx) = serve(&[("a.seq", b"xyz")]);
        let f = fetcher(&fx);
        fx.add_fault("a.seq", FaultKind::StatusRetryAfter(429, 1), 1);
        let t = Instant::now();
        assert_eq!(f.get_file("a.seq", Some(3), None).unwrap(), b"xyz");
        assert!(
            t.elapsed() >= Duration::from_millis(900),
            "{:?}",
            t.elapsed()
        );
    }

    #[test]
    fn permanent_failure_names_url_and_attempts() {
        let (_d, fx) = serve(&[("a.seq", b"xyz")]);
        let f = Fetcher::new(
            &fx.url("s"),
            FetchConfig {
                retries: 1,
                ..FetchConfig::default()
            },
        )
        .unwrap();
        fx.add_fault("a.seq", FaultKind::Status(500), 100);
        let err = f.get_file("a.seq", Some(3), None).unwrap_err().to_string();
        assert!(
            err.contains("/s/a.seq") && err.contains("2 attempts") && err.contains("HTTP 500"),
            "{err}"
        );
        assert_eq!(fx.requests_matching("a.seq").len(), 2);
    }

    #[test]
    fn forbidden_is_missing() {
        let (_d, fx) = serve(&[("x.seq", b"xyz")]);
        let f = fetcher(&fx);
        fx.add_fault("x.seq", FaultKind::Status(403), 2);
        assert!(matches!(
            f.get("x.seq").unwrap(),
            Got::Missing { status: 403 }
        ));
        let err = f.get_file("x.seq", Some(3), None).unwrap_err();
        assert_eq!(err.downcast_ref::<RemoteMissing>().unwrap().status, 403);
        assert!(
            err.to_string()
                .ends_with("x.seq: not found or access denied (HTTP 403)"),
            "{err}"
        );
        assert_eq!(f.counters().missing, 2);
    }

    #[test]
    fn get_many_keeps_order() {
        let (_d, fx) = serve(&[("a", b"1"), ("b", b"22"), ("c", b"333")]);
        let f = fetcher(&fx);
        let rels: Vec<String> = ["c", "zz", "a", "b"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let got = f.get_many(&rels).unwrap();
        let as_text: Vec<Option<Vec<u8>>> = got
            .into_iter()
            .map(|g| if let Got::Found(b) = g { Some(b) } else { None })
            .collect();
        assert_eq!(
            as_text,
            vec![
                Some(b"333".to_vec()),
                None,
                Some(b"1".to_vec()),
                Some(b"22".to_vec())
            ]
        );
    }

    #[test]
    fn range_reader_streams_in_pieces_and_retries_a_drop() {
        let full = body(200);
        let (_d, fx) = serve(&[("r.seq", &full)]);
        let f = Arc::new(fetcher(&fx).with_chunk(64));
        fx.add_fault("r.seq", FaultKind::Truncate, 1);
        let mut got = Vec::new();
        f.open_range("r.seq", 20, 180)
            .unwrap()
            .read_to_end(&mut got)
            .unwrap();
        assert_eq!(got, &full[20..180]);
        // 160 bytes in pieces of 64, plus the retried first piece.
        assert_eq!(fx.requests_matching("r.seq").len(), 4);
        assert_eq!(f.counters().retries, 1);
    }
}

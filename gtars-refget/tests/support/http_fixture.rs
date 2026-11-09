//! A small HTTP server for remote-store tests (standard library only).
//!
//! Serves one folder at `http://127.0.0.1:<port>/<prefix>/` for each prefix
//! (two prefixes give two "origins" over the same files). Supports byte
//! ranges (or ignoring them), 404 or 403 for a missing file, fault injection
//! and a request log.
#![allow(dead_code)]

use std::io::{Read, Write};
use std::net::{TcpListener, TcpStream};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Duration;

/// What a fault does to a matching request.
#[derive(Debug, Clone)]
pub enum FaultKind {
    /// Answer with this status and an empty body.
    Status(u16),
    /// Answer with this status and a `Retry-After` header (seconds).
    StatusRetryAfter(u16, u64),
    /// Send the full `Content-Length`, then close after half the body.
    Truncate,
    /// Flip one byte in the middle of the body.
    FlipByte,
    /// Wait before answering normally.
    Delay(Duration),
}

#[derive(Debug, Clone)]
struct Fault {
    path_contains: String,
    kind: FaultKind,
    remaining: usize,
}

/// One logged request.
#[derive(Debug, Clone)]
pub struct Logged {
    /// Path without the leading prefix, e.g. `sequences/ab/abc.seq`.
    pub path: String,
    pub prefix: String,
    pub range: Option<String>,
    pub status: u16,
}

struct State {
    root: PathBuf,
    prefixes: Vec<String>,
    honor_range: AtomicBool,
    /// Answer 403 for a missing file, as S3 does when listing is denied.
    listing_denied: bool,
    faults: Mutex<Vec<Fault>>,
    log: Mutex<Vec<Logged>>,
    in_flight: AtomicUsize,
    peak: AtomicUsize,
}

pub struct Fixture {
    pub port: u16,
    state: Arc<State>,
    stop: Arc<AtomicBool>,
    handle: Option<std::thread::JoinHandle<()>>,
}

impl Fixture {
    /// Serve `root` under each prefix, honoring ranges, 404 for a missing file.
    pub fn serve(root: &Path, prefixes: &[&str]) -> Self {
        Self::with_options(root, prefixes, true, false)
    }

    pub fn with_options(
        root: &Path,
        prefixes: &[&str],
        honor_range: bool,
        listing_denied: bool,
    ) -> Self {
        let listener = TcpListener::bind("127.0.0.1:0").expect("bind");
        let port = listener.local_addr().unwrap().port();
        let state = Arc::new(State {
            root: root.to_path_buf(),
            prefixes: prefixes
                .iter()
                .map(|p| p.trim_matches('/').to_string())
                .collect(),
            honor_range: AtomicBool::new(honor_range),
            listing_denied,
            faults: Mutex::new(Vec::new()),
            log: Mutex::new(Vec::new()),
            in_flight: AtomicUsize::new(0),
            peak: AtomicUsize::new(0),
        });
        let stop = Arc::new(AtomicBool::new(false));
        let handle = {
            let state = Arc::clone(&state);
            let stop = Arc::clone(&stop);
            std::thread::spawn(move || {
                for stream in listener.incoming() {
                    if stop.load(Ordering::Relaxed) {
                        break;
                    }
                    let Ok(stream) = stream else { continue };
                    let state = Arc::clone(&state);
                    std::thread::spawn(move || handle_conn(stream, &state));
                }
            })
        };
        Fixture {
            port,
            state,
            stop,
            handle: Some(handle),
        }
    }

    /// Base URL of one prefix, with a trailing `/`.
    pub fn url(&self, prefix: &str) -> String {
        format!(
            "http://127.0.0.1:{}/{}/",
            self.port,
            prefix.trim_matches('/')
        )
    }

    /// Apply `kind` to the next `hits` requests whose path contains `path_contains`.
    pub fn add_fault(&self, path_contains: &str, kind: FaultKind, hits: usize) {
        self.state.faults.lock().unwrap().push(Fault {
            path_contains: path_contains.to_string(),
            kind,
            remaining: hits,
        });
    }

    pub fn clear_faults(&self) {
        self.state.faults.lock().unwrap().clear();
    }

    pub fn requests(&self) -> Vec<Logged> {
        self.state.log.lock().unwrap().clone()
    }

    /// Requests whose path contains `s`.
    pub fn requests_matching(&self, s: &str) -> Vec<Logged> {
        self.requests()
            .into_iter()
            .filter(|r| r.path.contains(s))
            .collect()
    }

    pub fn clear_log(&self) {
        self.state.log.lock().unwrap().clear();
        self.state.peak.store(0, Ordering::Relaxed);
    }

    /// Most requests in flight at once since the last `clear_log`.
    pub fn peak_in_flight(&self) -> usize {
        self.state.peak.load(Ordering::Relaxed)
    }

    pub fn set_honor_range(&self, on: bool) {
        self.state.honor_range.store(on, Ordering::Relaxed);
    }
}

impl Drop for Fixture {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::Relaxed);
        let _ = TcpStream::connect(("127.0.0.1", self.port));
        if let Some(h) = self.handle.take() {
            let _ = h.join();
        }
    }
}

fn read_request(stream: &mut TcpStream) -> Option<String> {
    let mut buf = Vec::new();
    let mut chunk = [0u8; 4096];
    loop {
        let n = stream.read(&mut chunk).ok()?;
        if n == 0 {
            return None;
        }
        buf.extend_from_slice(&chunk[..n]);
        if buf.windows(4).any(|w| w == b"\r\n\r\n") {
            return Some(String::from_utf8_lossy(&buf).into_owned());
        }
        if buf.len() > 64 * 1024 {
            return None;
        }
    }
}

fn status_text(code: u16) -> &'static str {
    match code {
        200 => "OK",
        206 => "Partial Content",
        403 => "Forbidden",
        404 => "Not Found",
        416 => "Range Not Satisfiable",
        500 => "Internal Server Error",
        503 => "Service Unavailable",
        _ => "Status",
    }
}

fn handle_conn(mut stream: TcpStream, state: &State) {
    let _ = stream.set_read_timeout(Some(Duration::from_secs(10)));
    // One request per connection; every answer says `Connection: close`.
    let Some(req) = read_request(&mut stream) else {
        return;
    };
    let now = state.in_flight.fetch_add(1, Ordering::SeqCst) + 1;
    state.peak.fetch_max(now, Ordering::SeqCst);
    let status = respond(&mut stream, state, &req);
    let _ = status;
    let _ = stream.flush();
    state.in_flight.fetch_sub(1, Ordering::SeqCst);
}

fn respond(stream: &mut TcpStream, state: &State, req: &str) -> u16 {
    let mut lines = req.split("\r\n");
    let request_line = lines.next().unwrap_or("");
    let target = request_line
        .split_whitespace()
        .nth(1)
        .unwrap_or("/")
        .to_string();
    let range: Option<String> = lines.find_map(|l| {
        let (k, v) = l.split_once(':')?;
        k.trim()
            .eq_ignore_ascii_case("range")
            .then(|| v.trim().to_string())
    });
    let path = target
        .split('?')
        .next()
        .unwrap_or("")
        .trim_start_matches('/')
        .to_string();
    let (prefix, rel) = match path.split_once('/') {
        Some((p, r)) => (p.to_string(), r.to_string()),
        None => (path.clone(), String::new()),
    };

    let log = |status: u16| {
        state.log.lock().unwrap().push(Logged {
            path: rel.clone(),
            prefix: prefix.clone(),
            range: range.clone(),
            status,
        });
    };

    let fault = {
        let mut faults = state.faults.lock().unwrap();
        faults
            .iter_mut()
            .find(|f| f.remaining > 0 && path.contains(&f.path_contains))
            .map(|f| {
                f.remaining -= 1;
                f.kind.clone()
            })
    };
    let send = |stream: &mut TcpStream, code: u16, headers: &[(String, String)], body: &[u8]| {
        let mut head = format!(
            "HTTP/1.1 {} {}\r\nConnection: close\r\nContent-Length: {}\r\n",
            code,
            status_text(code),
            body.len()
        );
        for (k, v) in headers {
            head.push_str(&format!("{}: {}\r\n", k, v));
        }
        head.push_str("\r\n");
        let _ = stream.write_all(head.as_bytes());
        let _ = stream.write_all(body);
    };
    match &fault {
        Some(FaultKind::Status(code)) => {
            log(*code);
            send(stream, *code, &[], b"");
            return *code;
        }
        Some(FaultKind::StatusRetryAfter(code, secs)) => {
            log(*code);
            send(
                stream,
                *code,
                &[("Retry-After".into(), secs.to_string())],
                b"",
            );
            return *code;
        }
        Some(FaultKind::Delay(d)) => std::thread::sleep(*d),
        _ => {}
    }

    let missing_code = if state.listing_denied { 403 } else { 404 };
    let file = state.root.join(&rel);
    let body = if !state.prefixes.contains(&prefix) || rel.is_empty() || rel.contains("..") {
        None
    } else {
        std::fs::read(&file).ok()
    };
    let Some(mut body) = body else {
        log(missing_code);
        send(stream, missing_code, &[], b"");
        return missing_code;
    };
    let len = body.len() as u64;
    let mut code = 200;
    let mut headers = vec![("Accept-Ranges".to_string(), "bytes".to_string())];
    if let (Some(r), true) = (&range, state.honor_range.load(Ordering::Relaxed)) {
        let spec = r.trim_start_matches("bytes=");
        let (s, e) = spec.split_once('-').unwrap_or((spec, ""));
        let s: u64 = s.trim().parse().unwrap_or(0);
        let e: u64 = e.trim().parse().unwrap_or(u64::MAX);
        if s >= len {
            log(416);
            send(
                stream,
                416,
                &[("Content-Range".into(), format!("bytes */{}", len))],
                b"",
            );
            return 416;
        }
        let e = e.min(len - 1);
        body = body[s as usize..=e as usize].to_vec();
        headers.push(("Content-Range".into(), format!("bytes {}-{}/{}", s, e, len)));
        code = 206;
    }
    match fault {
        Some(FaultKind::Truncate) => {
            let mut head = format!(
                "HTTP/1.1 {} {}\r\nConnection: close\r\nContent-Length: {}\r\n",
                code,
                status_text(code),
                body.len()
            );
            for (k, v) in &headers {
                head.push_str(&format!("{}: {}\r\n", k, v));
            }
            head.push_str("\r\n");
            log(code);
            let _ = stream.write_all(head.as_bytes());
            let _ = stream.write_all(&body[..body.len() / 2]);
            let _ = stream.shutdown(std::net::Shutdown::Both);
            return code;
        }
        Some(FaultKind::FlipByte) if !body.is_empty() => {
            let mid = body.len() / 2;
            body[mid] ^= 0x01;
        }
        _ => {}
    }
    log(code);
    send(stream, code, &headers, &body);
    code
}

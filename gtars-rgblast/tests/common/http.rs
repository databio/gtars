//! A tiny static HTTP server with `Range` support that answers 403 for missing keys, like
//! an S3 bucket that denies listing.

use std::io::{BufRead, BufReader, Read, Write};
use std::net::{TcpListener, TcpStream};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

pub struct Server {
    /// `http://127.0.0.1:<port>`
    pub url: String,
    pub requests: Arc<AtomicU64>,
    pub ranged: Arc<AtomicU64>,
}

fn handle(
    mut s: TcpStream,
    root: &Path,
    requests: &AtomicU64,
    ranged: &AtomicU64,
) -> std::io::Result<()> {
    let mut r = BufReader::new(s.try_clone()?);
    let mut line = String::new();
    r.read_line(&mut line)?;
    let path = line.split_whitespace().nth(1).unwrap_or("/").to_string();
    let mut range: Option<(u64, u64)> = None;
    loop {
        let mut h = String::new();
        if r.read_line(&mut h)? == 0 || h == "\r\n" || h == "\n" {
            break;
        }
        let lower = h.to_ascii_lowercase();
        if let Some(v) = lower.strip_prefix("range: bytes=") {
            let (a, b) = v.trim().split_once('-').unwrap();
            range = Some((a.parse().unwrap(), b.parse().unwrap()));
        }
    }
    requests.fetch_add(1, Ordering::SeqCst);
    let rel = path.trim_start_matches('/');
    let file = root.join(rel);
    let body = if rel.contains("..") || !file.is_file() {
        None
    } else {
        std::fs::read(&file).ok()
    };
    match body {
        None => {
            let msg = b"<Error><Code>AccessDenied</Code></Error>";
            write!(
                s,
                "HTTP/1.1 403 Forbidden\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                msg.len()
            )?;
            s.write_all(msg)?;
        }
        Some(b) => match range {
            Some((a, e)) => {
                ranged.fetch_add(1, Ordering::SeqCst);
                let e = e.min(b.len() as u64 - 1);
                let part = &b[a as usize..=e as usize];
                write!(
                    s,
                    "HTTP/1.1 206 Partial Content\r\nContent-Length: {}\r\nContent-Range: bytes {a}-{e}/{}\r\nConnection: close\r\n\r\n",
                    part.len(),
                    b.len()
                )?;
                s.write_all(part)?;
            }
            None => {
                write!(
                    s,
                    "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                    b.len()
                )?;
                s.write_all(&b)?;
            }
        },
    }
    s.flush()?;
    // drain anything left so the client sees a clean close
    let _ = s.shutdown(std::net::Shutdown::Write);
    let mut sink = [0u8; 64];
    let _ = r.read(&mut sink);
    Ok(())
}

/// Serve `root` on a free local port (the server lives until the test process exits).
pub fn serve(root: &Path) -> Server {
    let l = TcpListener::bind("127.0.0.1:0").unwrap();
    let url = format!("http://{}", l.local_addr().unwrap());
    let requests = Arc::new(AtomicU64::new(0));
    let ranged = Arc::new(AtomicU64::new(0));
    let root: PathBuf = root.to_path_buf();
    let (rq, rg) = (requests.clone(), ranged.clone());
    std::thread::spawn(move || {
        for s in l.incoming().flatten() {
            let (root, rq, rg) = (root.clone(), rq.clone(), rg.clone());
            std::thread::spawn(move || {
                let _ = handle(s, &root, &rq, &rg);
            });
        }
    });
    Server {
        url,
        requests,
        ranged,
    }
}

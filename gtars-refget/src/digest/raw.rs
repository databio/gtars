//! Raw-byte digests and a fixed-size, heap-free digest text type.
//!
//! The store keeps digests as raw bytes (24 for sha512t24u, 16 for MD5) and
//! turns them back into text only on output. [`Digest32`] is that text: both a
//! sha512t24u (base64url of 24 bytes) and an MD5 (hex of 16 bytes) are exactly
//! 32 ASCII chars, so one inline `[u8; 32]` holds either with no allocation.

use serde::{Serialize, Serializer};
use std::fmt;
use std::ops::Deref;

/// Raw sha512t24u digest: the first 24 bytes of a SHA-512.
pub type RawSha = [u8; 24];
/// Raw MD5 digest.
pub type RawMd5 = [u8; 16];

/// A 32-char ASCII digest held inline. Derefs to `&str`.
#[repr(transparent)]
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Digest32([u8; 32]);

impl Digest32 {
    pub fn as_str(&self) -> &str {
        // SAFETY: only constructed from base64url or hex output, both ASCII.
        unsafe { std::str::from_utf8_unchecked(&self.0) }
    }
}

impl Deref for Digest32 {
    type Target = str;
    fn deref(&self) -> &str {
        self.as_str()
    }
}

impl AsRef<str> for Digest32 {
    fn as_ref(&self) -> &str {
        self.as_str()
    }
}

impl AsRef<[u8]> for Digest32 {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl fmt::Display for Digest32 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl fmt::Debug for Digest32 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Debug::fmt(self.as_str(), f)
    }
}

impl Serialize for Digest32 {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.serialize_str(self.as_str())
    }
}

impl PartialEq<str> for Digest32 {
    fn eq(&self, other: &str) -> bool {
        self.as_str() == other
    }
}

impl PartialEq<&str> for Digest32 {
    fn eq(&self, other: &&str) -> bool {
        self.as_str() == *other
    }
}

impl PartialEq<String> for Digest32 {
    fn eq(&self, other: &String) -> bool {
        self.as_str() == other.as_str()
    }
}

impl From<Digest32> for String {
    fn from(d: Digest32) -> String {
        d.as_str().to_string()
    }
}

const B64URL: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_";
const HEX: &[u8; 16] = b"0123456789abcdef";

fn b64url_value(c: u8) -> Option<u8> {
    match c {
        b'A'..=b'Z' => Some(c - b'A'),
        b'a'..=b'z' => Some(c - b'a' + 26),
        b'0'..=b'9' => Some(c - b'0' + 52),
        b'-' => Some(62),
        b'_' => Some(63),
        _ => None,
    }
}

fn hex_value(c: u8) -> Option<u8> {
    match c {
        b'0'..=b'9' => Some(c - b'0'),
        b'a'..=b'f' => Some(c - b'a' + 10),
        _ => None,
    }
}

/// Strip a GA4GH `SQ.` prefix, if present.
fn strip_sq(s: &[u8]) -> &[u8] {
    if s.len() > 3 && s[..3].eq_ignore_ascii_case(b"SQ.") { &s[3..] } else { s }
}

/// Decode a 32-char base64url sha512t24u (optionally `SQ.`-prefixed) to raw bytes.
pub fn sha_to_raw(s: impl AsRef<[u8]>) -> Option<RawSha> {
    let s = strip_sq(s.as_ref());
    if s.len() != 32 {
        return None;
    }
    let mut out = [0u8; 24];
    for (i, chunk) in s.chunks_exact(4).enumerate() {
        let a = b64url_value(chunk[0])? as u32;
        let b = b64url_value(chunk[1])? as u32;
        let c = b64url_value(chunk[2])? as u32;
        let d = b64url_value(chunk[3])? as u32;
        let v = (a << 18) | (b << 12) | (c << 6) | d;
        out[i * 3] = (v >> 16) as u8;
        out[i * 3 + 1] = (v >> 8) as u8;
        out[i * 3 + 2] = v as u8;
    }
    Some(out)
}

/// Decode a 32-char lowercase-hex MD5 to raw bytes.
pub fn md5_to_raw(s: impl AsRef<[u8]>) -> Option<RawMd5> {
    let s = s.as_ref();
    if s.len() != 32 {
        return None;
    }
    let mut out = [0u8; 16];
    for i in 0..16 {
        out[i] = (hex_value(s[2 * i])? << 4) | hex_value(s[2 * i + 1])?;
    }
    Some(out)
}

/// Encode raw sha512t24u bytes as 32-char base64url text.
pub fn raw_to_sha(raw: &RawSha) -> Digest32 {
    let mut out = [0u8; 32];
    for (i, chunk) in raw.chunks_exact(3).enumerate() {
        let v = ((chunk[0] as u32) << 16) | ((chunk[1] as u32) << 8) | chunk[2] as u32;
        out[i * 4] = B64URL[(v >> 18) as usize & 63];
        out[i * 4 + 1] = B64URL[(v >> 12) as usize & 63];
        out[i * 4 + 2] = B64URL[(v >> 6) as usize & 63];
        out[i * 4 + 3] = B64URL[v as usize & 63];
    }
    Digest32(out)
}

/// Encode raw MD5 bytes as 32-char lowercase hex text.
pub fn raw_to_md5(raw: &RawMd5) -> Digest32 {
    let mut out = [0u8; 32];
    for (i, b) in raw.iter().enumerate() {
        out[2 * i] = HEX[(b >> 4) as usize];
        out[2 * i + 1] = HEX[(b & 15) as usize];
    }
    Digest32(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::digest::algorithms::{md5, sha512t24u};

    #[test]
    fn sha_round_trip_matches_algorithm_output() {
        for input in ["", "ACGT", "the quick brown fox", "NNNNNNNNNNNN"] {
            let text = sha512t24u(input);
            let raw = sha_to_raw(&text).expect("valid sha");
            assert_eq!(raw_to_sha(&raw).as_str(), text);
            let prefixed = format!("SQ.{}", text);
            assert_eq!(sha_to_raw(&prefixed), Some(raw));
        }
    }

    #[test]
    fn md5_round_trip_matches_algorithm_output() {
        for input in ["", "ACGT", "the quick brown fox"] {
            let text = md5(input);
            let raw = md5_to_raw(&text).expect("valid md5");
            assert_eq!(raw_to_md5(&raw).as_str(), text);
        }
    }

    #[test]
    fn malformed_digests_are_rejected() {
        assert!(sha_to_raw("short").is_none());
        assert!(sha_to_raw("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!").is_none());
        assert!(md5_to_raw("zz000000000000000000000000000000").is_none());
        assert!(md5_to_raw("ABCDEF00000000000000000000000000").is_none());
    }
}

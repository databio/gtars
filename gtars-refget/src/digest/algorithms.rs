//! Core digest algorithms - WASM-safe, no filesystem dependencies.

use md5::Md5;
use serde_json::Value;
use sha2::{Digest, Sha512};

/// Processes a given string to compute its GA4GH sha512t24u digest.
///
/// This function processes a given string to compute its GA4GH sha512t24u digest. The input string
/// is processed in chunks of 800 bytes, and the digest is computed incrementally. The final digest
/// is a 24-byte string encoded using base64url encoding. You can provide either a string slice or
/// a byte slice as input.
///
/// # Arguments
///
/// * `input` - The input string to be processed, as a string slice or byte slice.
///
/// # Returns
///
/// A string SHA-512 digest of the input string.
pub fn sha512t24u<T: AsRef<[u8]>>(input: T) -> String {
    let mut sha512_hasher_box = Box::new(Sha512::new());
    for chunk in input.as_ref().chunks(1024) {
        sha512_hasher_box.as_mut().update(chunk);
    }
    base64_url::encode(&sha512_hasher_box.as_mut().finalize_reset()[0..24])
}

/// Process a string to compute its md5 digest
///
/// # Arguments
///
/// * `input` - The input string to be processed, as a string slice or byte slice.
///
/// # Returns
///
/// A string MD5 digest of the input string.
pub fn md5<T: AsRef<[u8]>>(input: T) -> String {
    let mut hasher = Md5::new();
    for chunk in input.as_ref().chunks(1024) {
        hasher.update(chunk);
    }
    format!("{:x}", hasher.finalize())
}

/// Incremental hasher producing the same sha512t24u / md5 digests as
/// [`sha512t24u`] / [`md5`], fed in arbitrarily-sized chunks.
///
/// WASM-safe (no filesystem dependencies). Used by `store::verify` to hash
/// streamed sequence bytes without buffering the whole sequence in memory.
pub struct SequenceHasher {
    sha512: Sha512,
    md5: Option<Md5>,
    len: u64,
}

impl SequenceHasher {
    /// Create a new hasher. When `with_md5` is `false`, [`Self::finalize`]
    /// always returns `None` for the md5 digest and the md5 state is never
    /// updated (saving the work when the caller does not need it).
    pub fn new(with_md5: bool) -> Self {
        Self {
            sha512: Sha512::new(),
            md5: if with_md5 { Some(Md5::new()) } else { None },
            len: 0,
        }
    }

    /// Feed the next chunk of bytes. Chunk size is arbitrary; the result is
    /// independent of how the input is chunked.
    pub fn update(&mut self, bytes: &[u8]) {
        self.sha512.update(bytes);
        if let Some(md5) = self.md5.as_mut() {
            md5.update(bytes);
        }
        self.len += bytes.len() as u64;
    }

    /// Finalize and return `(sha512t24u, md5, total_bytes_fed)`. `md5` is
    /// `None` when this hasher was created with `with_md5: false`.
    pub fn finalize(self) -> (String, Option<String>, u64) {
        let sha512t24u = base64_url::encode(&self.sha512.finalize()[0..24]);
        let md5 = self.md5.map(|hasher| format!("{:x}", hasher.finalize()));
        (sha512t24u, md5, self.len)
    }
}

/// Apply RFC-8785 JSON Canonicalization Scheme (JCS) to a JSON value
///
/// This function canonicalizes a JSON value according to RFC-8785 by:
/// 1. Removing whitespace
/// 2. Sorting object keys lexicographically
/// 3. Using consistent number formatting
/// 4. Ensuring UTF-8 encoding
///
/// # Arguments
/// * `value` - The JSON value to canonicalize
///
/// # Returns
/// A canonicalized JSON string
pub fn canonicalize_json(value: &Value) -> String {
    match value {
        Value::Null => "null".to_string(),
        Value::Bool(b) => b.to_string(),
        Value::Number(n) => {
            // RFC-8785 requires specific number formatting
            if let Some(i) = n.as_i64() {
                i.to_string()
            } else if let Some(u) = n.as_u64() {
                u.to_string()
            } else if let Some(f) = n.as_f64() {
                // Format floating point numbers without unnecessary trailing zeros
                if f.fract() == 0.0 {
                    format!("{:.0}", f)
                } else {
                    // Remove trailing zeros after decimal point
                    let formatted = format!("{}", f);
                    formatted
                        .trim_end_matches('0')
                        .trim_end_matches('.')
                        .to_string()
                }
            } else {
                n.to_string()
            }
        }
        Value::String(s) => {
            // Escape string according to JSON rules
            serde_json::to_string(s).unwrap()
        }
        Value::Array(arr) => {
            let elements: Vec<String> = arr.iter().map(canonicalize_json).collect();
            format!("[{}]", elements.join(","))
        }
        Value::Object(obj) => {
            // Sort keys lexicographically (RFC-8785 requirement)
            let mut sorted_keys: Vec<&String> = obj.keys().collect();
            sorted_keys.sort();

            let pairs: Vec<String> = sorted_keys
                .iter()
                .map(|key| {
                    let key_str = serde_json::to_string(key).unwrap();
                    let value_str = canonicalize_json(&obj[*key]);
                    format!("{}:{}", key_str, value_str)
                })
                .collect();

            format!("{{{}}}", pairs.join(","))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn digests_sha512t24u_bytes() {
        let digest = sha512t24u(b"hello world");
        assert_eq!(digest, "MJ7MSJwS1utMxA9QyQLytNDtd-5RGnx6");
    }

    #[test]
    fn digests_sha512t24u_str() {
        let digest = sha512t24u("hello world");
        assert_eq!(digest, "MJ7MSJwS1utMxA9QyQLytNDtd-5RGnx6");
    }

    #[test]
    fn digests_md5() {
        let digest = md5("hello world");
        assert_eq!(digest, "5eb63bbbe01eeed093cb22bb8f5acdc3");
    }

    #[test]
    fn sequence_hasher_matches_oneshot_regardless_of_chunking() {
        let data: Vec<u8> = (0..5000u32).map(|i| (i % 251) as u8).collect();
        let expected_sha = sha512t24u(&data);
        let expected_md5 = md5(&data);

        let mut hasher = SequenceHasher::new(true);
        let mut offset = 0usize;
        let mut chunk_len = 1usize;
        while offset < data.len() {
            let end = (offset + chunk_len).min(data.len());
            hasher.update(&data[offset..end]);
            offset = end;
            // Pseudo-random, varying chunk sizes (deterministic, no RNG dependency).
            chunk_len = chunk_len * 7 % 613 + 1;
        }
        let (sha512t24u_out, md5_out, len) = hasher.finalize();

        assert_eq!(sha512t24u_out, expected_sha);
        assert_eq!(md5_out, Some(expected_md5));
        assert_eq!(len, data.len() as u64);
    }

    #[test]
    fn sequence_hasher_without_md5_returns_none() {
        let mut hasher = SequenceHasher::new(false);
        hasher.update(b"ACGTACGT");
        let (_sha, md5_out, len) = hasher.finalize();
        assert_eq!(md5_out, None);
        assert_eq!(len, 8);
    }
}

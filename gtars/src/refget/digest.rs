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
}

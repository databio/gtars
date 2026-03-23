/// Size of the fixed-length byte key used for digest lookups.
/// Must be large enough to hold any GA4GH digest with prefix (e.g. "SQ." = 35 bytes,
/// "SQC." = 36 bytes).
pub const DIGEST_KEY_LEN: usize = 48;

/// Fixed-length byte array used as a HashMap key for digest lookups.
pub type DigestKey = [u8; DIGEST_KEY_LEN];

/// Trait for converting types to a fixed-length byte key suitable for hashing.
///
/// Uses a blanket implementation for any type implementing `AsRef<[u8]>`,
/// which covers `str`, `String`, `&[u8]`, `Vec<u8>`, etc. -- no need for
/// separate `str`/`String` implementations.
///
/// # Example
///
/// use crate::refget::hashkeyable::HashKeyable;
/// let key: DigestKey = "example".to_key();
///
pub trait HashKeyable {
    fn to_key(&self) -> DigestKey;
}

impl<T: AsRef<[u8]>> HashKeyable for T {
    fn to_key(&self) -> DigestKey {
        let mut key = [0u8; DIGEST_KEY_LEN];
        let bytes = self.as_ref();
        let len = std::cmp::min(bytes.len(), DIGEST_KEY_LEN);
        key[..len].copy_from_slice(&bytes[..len]);
        key
    }
}

/// Convert a [`DigestKey`] back to a digest string.
///
/// This is the inverse of [`HashKeyable::to_key()`]: it finds the first null byte
/// (or uses all bytes if none) and converts the prefix to a UTF-8 string.
pub(crate) fn key_to_digest_string(key: &DigestKey) -> String {
    let len = key.iter().position(|&b| b == 0).unwrap_or(DIGEST_KEY_LEN);
    String::from_utf8_lossy(&key[..len]).to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ga4gh_digest_roundtrip() {
        let digest = "SQ.KH8TQq_rHNXuRAfzbjTTnRL_k0ZMznaM";
        assert_eq!(digest.len(), 35);
        let key: DigestKey = digest.to_key();
        let recovered = key_to_digest_string(&key);
        assert_eq!(recovered, digest, "digest was truncated during roundtrip");
    }
}

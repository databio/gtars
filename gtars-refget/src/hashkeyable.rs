/// Trait for converting types to a 32-byte key suitable for hashing.
///
/// Uses a blanket implementation for any type implementing `AsRef<[u8]>`,
/// which covers `str`, `String`, `&[u8]`, `Vec<u8>`, etc. -- no need for
/// separate `str`/`String` implementations.
///
/// # Example
///
/// use crate::refget::hashkeyable::HashKeyable;
/// let key: [u8; 32] = "example".to_key();
///
pub trait HashKeyable {
    fn to_key(&self) -> [u8; 32];
}

impl<T: AsRef<[u8]>> HashKeyable for T {
    fn to_key(&self) -> [u8; 32] {
        let mut key = [0u8; 32];
        let bytes = self.as_ref();
        let len = std::cmp::min(bytes.len(), 32);
        key[..len].copy_from_slice(&bytes[..len]);
        key
    }
}

/// Convert a `[u8; 32]` key back to a digest string.
///
/// This is the inverse of [`HashKeyable::to_key()`]: it finds the first null byte
/// (or uses all 32 bytes if none) and converts the prefix to a UTF-8 string.
pub(crate) fn key_to_digest_string(key: &[u8; 32]) -> String {
    let len = key.iter().position(|&b| b == 0).unwrap_or(32);
    String::from_utf8_lossy(&key[..len]).to_string()
}

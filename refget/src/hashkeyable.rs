/// Trait for converting types to a 32-byte key suitable for hashing.
///
/// This trait provides to_key(), which converts the object into a fixed-size byte array
/// that can be used as a key in hashmaps or for hashing purposes.
/// The trait is implemented here for `str` and `String`. This allows me to use String
/// for easy display and human readability, but using the fixed-size byte array for hashing,
/// which is more efficient and avoids issues with variable-length strings.
///
/// # Example
///
/// use crate::refget::hashkeyable::HashKeyable;
/// let key: [u8; 32] = "example".to_key();
///
pub trait HashKeyable {
    fn to_key(&self) -> [u8; 32];
}

// impl HashKeyable for str {
//     fn to_key(&self) -> [u8; 32] {
//         let mut key = [0u8; 32];
//         let bytes = self.as_bytes();
//         let len = std::cmp::min(bytes.len(), 32);
//         key[..len].copy_from_slice(&bytes[..len]);
//         key
//     }
// }

// impl HashKeyable for String {
//     fn to_key(&self) -> [u8; 32] {
//         self.as_str().to_key()
//     }
// }
// Implementing HashKeyable for str and String is not necessary since we can use a generic implementation for any type that implements AsRef<[u8]>.

impl<T: AsRef<[u8]>> HashKeyable for T {
    fn to_key(&self) -> [u8; 32] {
        let mut key = [0u8; 32];
        let bytes = self.as_ref();
        let len = std::cmp::min(bytes.len(), 32);
        key[..len].copy_from_slice(&bytes[..len]);
        key
    }
}

//! Auto-decompressing writer for push-based decompression.
//!
//! This module provides `AutoDecompressWriter`, a wrapper that auto-detects gzip
//! compression from magic bytes and transparently decompresses data written to it.
//!
//! # Design Note
//!
//! This is the push-based (write) counterpart to `get_dynamic_reader` in gtars-core.
//! While `get_dynamic_reader` is pull-based (you read decompressed data out),
//! `AutoDecompressWriter` is push-based (you write compressed data in).
//!
//! This utility lives in the digest module to avoid filesystem feature gating,
//! since it works with any `Write` impl and doesn't require file paths.
//! Eventually it should move to a gtars-io module when that exists.
//!
//! # Example
//! ```
//! use std::io::Write;
//! use gtars_refget::digest::auto_decompress::AutoDecompressWriter;
//!
//! let mut output = Vec::new();
//! let mut writer = AutoDecompressWriter::new(&mut output);
//!
//! // Write data (plain or gzipped - auto-detected from magic bytes)
//! writer.write_all(b"Hello, world!").unwrap();
//! writer.finish().unwrap();
//!
//! assert_eq!(&output, b"Hello, world!");
//! ```

use flate2::write::GzDecoder;
use std::io::{self, Write};

/// Internal state for format detection
enum WriterState<W: Write> {
    /// Haven't detected format yet
    Detecting(W),
    /// Plain data (uncompressed) - pass through directly
    Plain(W),
    /// Gzipped data - decompress via GzDecoder
    Gzipped(GzDecoder<W>),
    /// Consumed (after finish())
    Consumed,
}

/// A writer wrapper that auto-detects gzip compression and decompresses transparently.
///
/// On the first write, examines magic bytes to detect gzip format:
/// - If bytes start with 0x1f 0x8b (gzip magic), wraps inner writer in GzDecoder
/// - Otherwise, passes data through directly to inner writer
///
/// This is useful for streaming scenarios where you receive data in chunks
/// and don't know the format ahead of time.
///
/// # Type Parameters
/// * `W` - The underlying writer that receives decompressed data
pub struct AutoDecompressWriter<W: Write> {
    state: WriterState<W>,
}

impl<W: Write> AutoDecompressWriter<W> {
    /// Create a new auto-decompressing writer wrapping the given writer.
    ///
    /// Format detection happens on the first write.
    pub fn new(inner: W) -> Self {
        Self {
            state: WriterState::Detecting(inner),
        }
    }

    /// Finish writing and return the inner writer.
    ///
    /// For gzipped data, this flushes the decompressor and verifies the gzip footer.
    /// Must be called after all data has been written.
    ///
    /// # Returns
    /// The inner writer, or an error if decompression failed.
    pub fn finish(self) -> io::Result<W> {
        match self.state {
            WriterState::Detecting(inner) => Ok(inner),
            WriterState::Plain(inner) => Ok(inner),
            WriterState::Gzipped(decoder) => decoder.finish(),
            WriterState::Consumed => Err(io::Error::new(
                io::ErrorKind::Other,
                "Writer already consumed",
            )),
        }
    }

    /// Check if gzip compression was detected.
    ///
    /// Returns `None` if no data has been written yet (format not yet detected).
    pub fn is_gzipped(&self) -> Option<bool> {
        match &self.state {
            WriterState::Detecting(_) => None,
            WriterState::Plain(_) => Some(false),
            WriterState::Gzipped(_) => Some(true),
            WriterState::Consumed => None,
        }
    }

    /// Get a reference to the inner writer (for plain mode only).
    ///
    /// Returns `None` if gzipped or not yet detected.
    pub fn get_ref(&self) -> Option<&W> {
        match &self.state {
            WriterState::Plain(inner) => Some(inner),
            WriterState::Detecting(inner) => Some(inner),
            WriterState::Gzipped(decoder) => Some(decoder.get_ref()),
            WriterState::Consumed => None,
        }
    }
}

impl<W: Write> Write for AutoDecompressWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        // Detect format on first non-empty write
        if let WriterState::Detecting(_) = &self.state {
            if !buf.is_empty() {
                // Take ownership of inner writer temporarily
                let inner = match std::mem::replace(&mut self.state, WriterState::Consumed) {
                    WriterState::Detecting(w) => w,
                    _ => unreachable!(),
                };

                // Check for gzip magic bytes (0x1f, 0x8b)
                let is_gzipped = buf.len() >= 2 && buf[0] == 0x1f && buf[1] == 0x8b;

                if is_gzipped {
                    self.state = WriterState::Gzipped(GzDecoder::new(inner));
                } else {
                    self.state = WriterState::Plain(inner);
                }
            }
        }

        // Write to appropriate target
        match &mut self.state {
            WriterState::Detecting(_) => Ok(0), // Empty buffer, nothing to write
            WriterState::Plain(inner) => inner.write(buf),
            WriterState::Gzipped(decoder) => decoder.write(buf),
            WriterState::Consumed => Err(io::Error::new(
                io::ErrorKind::Other,
                "Writer already consumed",
            )),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match &mut self.state {
            WriterState::Detecting(inner) => inner.flush(),
            WriterState::Plain(inner) => inner.flush(),
            WriterState::Gzipped(decoder) => decoder.flush(),
            WriterState::Consumed => Ok(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;

    #[test]
    fn test_plain_passthrough() {
        let mut output = Vec::new();
        {
            let mut writer = AutoDecompressWriter::new(&mut output);
            writer.write_all(b"Hello, world!").unwrap();
            writer.finish().unwrap();
        }
        assert_eq!(&output, b"Hello, world!");
    }

    #[test]
    fn test_gzip_decompression() {
        // Compress some data
        let original = b"Hello, gzip world!";
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(original).unwrap();
        let compressed = encoder.finish().unwrap();

        // Decompress via AutoDecompressWriter
        let mut output = Vec::new();
        {
            let mut writer = AutoDecompressWriter::new(&mut output);
            writer.write_all(&compressed).unwrap();
            writer.finish().unwrap();
        }
        assert_eq!(&output, original);
    }

    #[test]
    fn test_gzip_chunked() {
        // Compress some data
        let original = b"Hello, this is a longer message for chunked testing!";
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(original).unwrap();
        let compressed = encoder.finish().unwrap();

        // Decompress in small chunks
        let mut output = Vec::new();
        {
            let mut writer = AutoDecompressWriter::new(&mut output);
            for chunk in compressed.chunks(5) {
                writer.write_all(chunk).unwrap();
            }
            writer.finish().unwrap();
        }
        assert_eq!(&output, original);
    }

    #[test]
    fn test_is_gzipped_detection() {
        // Plain data
        let mut output = Vec::new();
        let mut writer = AutoDecompressWriter::new(&mut output);
        assert_eq!(writer.is_gzipped(), None); // Not yet detected
        writer.write_all(b"plain").unwrap();
        assert_eq!(writer.is_gzipped(), Some(false));
        drop(writer);

        // Gzipped data
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(b"data").unwrap();
        let compressed = encoder.finish().unwrap();

        let mut output = Vec::new();
        let mut writer = AutoDecompressWriter::new(&mut output);
        writer.write_all(&compressed).unwrap();
        assert_eq!(writer.is_gzipped(), Some(true));
    }

    #[test]
    fn test_empty_write() {
        let mut output = Vec::new();
        let mut writer = AutoDecompressWriter::new(&mut output);
        writer.write_all(b"").unwrap(); // Empty write
        assert_eq!(writer.is_gzipped(), None); // Still not detected
        writer.write_all(b"data").unwrap();
        assert_eq!(writer.is_gzipped(), Some(false));
        let _ = writer.finish();
        assert_eq!(&output, b"data");
    }

    #[test]
    fn test_single_byte_non_gzip() {
        // Single byte that's not gzip magic
        let mut output = Vec::new();
        {
            let mut writer = AutoDecompressWriter::new(&mut output);
            writer.write_all(b"X").unwrap();
            writer.finish().unwrap();
        }
        assert_eq!(&output, b"X");
    }

    #[test]
    fn test_get_ref() {
        let mut output = Vec::new();
        let writer = AutoDecompressWriter::new(&mut output);
        assert!(writer.get_ref().is_some());
    }
}

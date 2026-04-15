//! Streaming decoder for bit-packed encoded sequences - WASM-safe.
//!
//! Unlike [`decode_substring_from_bytes`](super::encoder::decode_substring_from_bytes),
//! which requires the entire encoded byte slice to be buffered in memory,
//! [`StreamingDecoder`] consumes an [`std::io::Read`] source incrementally and
//! emits decoded ASCII bases via its own `Read` implementation. Only a tiny
//! bit buffer (`u64`) of internal state is used, so memory is bounded regardless
//! of sequence length.
//!
//! The decoder relies on the supplied [`Alphabet`]'s `decoding_array`, so `N`
//! and other ambiguity codes are handled identically to the in-memory path.
//!
//! # Bit ordering
//!
//! Matches the encoder: MSB-first within each byte (see
//! [`encode_sequence`](super::encoder::encode_sequence)).
//!
//! # Non-byte-aligned windows
//!
//! For alphabets where `bits_per_symbol` does not divide 8 (2-bit, 3-bit, 5-bit),
//! a requested `[start, end)` base window typically does not land on byte
//! boundaries. Callers should:
//!
//! 1. Feed the decoder the exact byte range
//!    `[floor(start * bps / 8), ceil(end * bps / 8))`.
//! 2. Pass `leading_skip_bits = (start * bps) mod 8` — the number of *bits* at
//!    the very front of the byte window that precede the requested base window
//!    (0..=7). The decoder consumes and discards them before emitting.
//! 3. Pass `bases_to_emit = end - start`.
//!
//! For the Raw / ASCII case (`bits_per_symbol == 8`), byte and base boundaries
//! coincide; pass `leading_skip = 0`.

use std::io::{self, Read};

use super::alphabet::Alphabet;

/// Streaming decoder for bit-packed encoded sequences.
///
/// Wraps any `impl Read` source and implements `Read` itself, emitting
/// decoded ASCII bases. Internal state is O(bits_per_symbol) — just a
/// 64-bit bit buffer plus a handful of counters.
pub struct StreamingDecoder<R: Read> {
    inner: R,
    alphabet: &'static Alphabet,
    bits_per_symbol: u8,
    bases_remaining: u64,
    /// Number of bits to discard at the front of the byte stream (0..=7).
    /// Consumed on the first `read` call.
    leading_skip_bits: u8,
    bit_buffer: u64,
    bit_buffer_len: u8,
}

impl<R: Read> StreamingDecoder<R> {
    /// Construct a new streaming decoder.
    ///
    /// - `inner` - Byte source providing the bit-packed encoded sequence.
    /// - `alphabet` - Alphabet whose `bits_per_symbol` and `decoding_array` drive decoding.
    /// - `leading_skip_bits` - Number of bits at the very start of the byte
    ///   stream to discard before the first emitted symbol. Must be `< 8`.
    ///   Should be `0` when `bits_per_symbol == 8`.
    /// - `bases_to_emit` - Total number of bases the decoder will emit via `Read`.
    pub fn new(
        inner: R,
        alphabet: &'static Alphabet,
        leading_skip_bits: u8,
        bases_to_emit: u64,
    ) -> Self {
        let bits_per_symbol = alphabet.bits_per_symbol as u8;
        debug_assert!(
            bits_per_symbol > 0 && bits_per_symbol <= 8,
            "StreamingDecoder supports 1..=8 bits per symbol, got {}",
            bits_per_symbol
        );
        debug_assert!(
            leading_skip_bits < 8,
            "leading_skip_bits must be < 8, got {}",
            leading_skip_bits
        );
        StreamingDecoder {
            inner,
            alphabet,
            bits_per_symbol,
            bases_remaining: bases_to_emit,
            leading_skip_bits,
            bit_buffer: 0,
            bit_buffer_len: 0,
        }
    }

    /// Ensure `bit_buffer` contains at least `min_bits` valid bits, reading
    /// more bytes from `inner` if needed. Returns `UnexpectedEof` if the
    /// source is exhausted before the requirement is met.
    fn refill(&mut self, min_bits: u8) -> io::Result<()> {
        while self.bit_buffer_len < min_bits {
            let mut byte = [0u8; 1];
            match self.inner.read(&mut byte) {
                Ok(0) => {
                    return Err(io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        "StreamingDecoder: source ended before all bases were decoded",
                    ));
                }
                Ok(_) => {
                    // Shift in 8 bits (MSB-first: new byte's bits go in the low end
                    // and become the "next" symbol's bits after older symbols are consumed
                    // from the top).
                    self.bit_buffer = (self.bit_buffer << 8) | (byte[0] as u64);
                    self.bit_buffer_len += 8;
                }
                Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
                Err(e) => return Err(e),
            }
        }
        Ok(())
    }

    /// Pull the next decoded base from the stream, or `Ok(None)` if the
    /// configured base budget has been exhausted.
    fn next_symbol(&mut self) -> io::Result<Option<u8>> {
        if self.bases_remaining == 0 {
            return Ok(None);
        }
        self.refill(self.bits_per_symbol)?;

        let shift = self.bit_buffer_len - self.bits_per_symbol;
        let mask: u64 = if self.bits_per_symbol == 64 {
            u64::MAX
        } else {
            (1u64 << self.bits_per_symbol) - 1
        };
        let code = ((self.bit_buffer >> shift) & mask) as u8;
        self.bit_buffer_len -= self.bits_per_symbol;
        // Mask off consumed bits so they don't leak into later symbols.
        let keep_mask: u64 = if self.bit_buffer_len == 0 {
            0
        } else {
            (1u64 << self.bit_buffer_len) - 1
        };
        self.bit_buffer &= keep_mask;

        self.bases_remaining -= 1;
        Ok(Some(self.alphabet.decoding_array[code as usize]))
    }
}

impl<R: Read> Read for StreamingDecoder<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        // Discard any leading bits from the start of the byte stream (the
        // bits that preceded the requested base window within the first byte).
        if self.leading_skip_bits > 0 {
            self.refill(self.leading_skip_bits)?;
            // Consume (discard) the top `leading_skip_bits` of the buffer.
            self.bit_buffer_len -= self.leading_skip_bits;
            let keep_mask: u64 = if self.bit_buffer_len == 0 {
                0
            } else {
                (1u64 << self.bit_buffer_len) - 1
            };
            self.bit_buffer &= keep_mask;
            self.leading_skip_bits = 0;
        }

        let mut written = 0;
        while written < buf.len() {
            match self.next_symbol()? {
                Some(base) => {
                    buf[written] = base;
                    written += 1;
                }
                None => break,
            }
        }
        Ok(written)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::digest::alphabet::{
        ASCII_ALPHABET, DNA_2BIT_ALPHABET, DNA_3BIT_ALPHABET, DNA_IUPAC_ALPHABET, PROTEIN_ALPHABET,
    };
    use crate::digest::encoder::{decode_substring_from_bytes, encode_sequence};
    use std::io::Cursor;

    /// Compute the tight byte window and leading-skip for a `[start, end)`
    /// base range, matching what `RefgetStore::stream_sequence` will do.
    fn byte_window(
        start: usize,
        end: usize,
        bits_per_symbol: usize,
    ) -> (usize, usize, u8) {
        let start_bit = start * bits_per_symbol;
        let end_bit = end * bits_per_symbol;
        let byte_start = start_bit / 8;
        let byte_end = end_bit.div_ceil(8);
        let leading_skip_bits = (start_bit - byte_start * 8) as u8;
        (byte_start, byte_end, leading_skip_bits)
    }

    fn run_range_test(
        sequence: &[u8],
        alphabet: &'static Alphabet,
        start: usize,
        end: usize,
    ) {
        let encoded = encode_sequence(sequence, alphabet);
        let (byte_start, byte_end, leading_skip) =
            byte_window(start, end, alphabet.bits_per_symbol);
        let slice = &encoded[byte_start..byte_end.min(encoded.len())];

        let expected = decode_substring_from_bytes(&encoded, start, end, alphabet);

        let mut decoder = StreamingDecoder::new(
            Cursor::new(slice.to_vec()),
            alphabet,
            leading_skip,
            (end - start) as u64,
        );
        let mut out = Vec::new();
        decoder.read_to_end(&mut out).expect("streaming read failed");
        assert_eq!(
            out, expected,
            "alphabet bps={} range [{}, {}) sequence={:?}",
            alphabet.bits_per_symbol,
            start,
            end,
            std::str::from_utf8(sequence).unwrap_or("<non-utf8>")
        );
    }

    // Test inputs per alphabet. Each sequence is chosen to exercise the decoder
    // across the full symbol set where possible.
    fn fixtures() -> Vec<(&'static Alphabet, &'static [u8])> {
        vec![
            (&DNA_2BIT_ALPHABET, b"ACGTACGTACGTACGT" as &[u8]),
            (&DNA_3BIT_ALPHABET, b"ACGTNRYXACGTNRYX" as &[u8]),
            (&DNA_IUPAC_ALPHABET, b"ACGTRYMKSWBDHVN-" as &[u8]),
            (&PROTEIN_ALPHABET, b"ACDEFGHIKLMNPQRSTVWY*X-" as &[u8]),
            (&ASCII_ALPHABET, b"Hello, World! 1234" as &[u8]),
        ]
    }

    #[test]
    fn test_full_sequence() {
        for (alphabet, seq) in fixtures() {
            run_range_test(seq, alphabet, 0, seq.len());
        }
    }

    #[test]
    fn test_aligned_start_aligned_end() {
        // For 2-bit, 4 symbols per byte — [4, 12) is byte-aligned.
        run_range_test(
            b"ACGTACGTACGTACGT",
            &DNA_2BIT_ALPHABET,
            4,
            12,
        );
    }

    #[test]
    fn test_unaligned_start() {
        run_range_test(
            b"ACGTACGTACGTACGT",
            &DNA_2BIT_ALPHABET,
            1,
            16,
        );
    }

    #[test]
    fn test_unaligned_end() {
        run_range_test(
            b"ACGTACGTACGTACGT",
            &DNA_2BIT_ALPHABET,
            0,
            15,
        );
    }

    #[test]
    fn test_unaligned_both_3bit() {
        run_range_test(
            b"ACGTNRYXACGTNRYX",
            &DNA_3BIT_ALPHABET,
            1,
            7,
        );
    }

    #[test]
    fn test_zero_length() {
        // Zero-length windows should emit nothing regardless of alphabet or
        // leading_skip_bits — the bit-buffer never needs refilling.
        for (alphabet, seq) in fixtures() {
            let encoded = encode_sequence(seq, alphabet);
            // A zero-length window at offset 5 — the decoder should read
            // zero bytes and produce an empty output.
            let mut decoder = StreamingDecoder::new(
                Cursor::new(encoded.clone()),
                alphabet,
                0,
                0,
            );
            let mut out = Vec::new();
            decoder.read_to_end(&mut out).expect("zero-length read");
            assert!(out.is_empty());
        }
    }

    #[test]
    fn test_single_base() {
        for (alphabet, seq) in fixtures() {
            for start in 0..seq.len() {
                run_range_test(seq, alphabet, start, start + 1);
            }
        }
    }

    #[test]
    fn test_start_at_end_minus_one() {
        for (alphabet, seq) in fixtures() {
            run_range_test(seq, alphabet, seq.len() - 1, seq.len());
        }
    }

    #[test]
    fn test_read_small_buf() {
        let sequence = b"ACGTACGTACGTACGT";
        let alphabet = &DNA_2BIT_ALPHABET;
        let encoded = encode_sequence(sequence, alphabet);
        let (byte_start, byte_end, leading_skip) =
            byte_window(1, 15, alphabet.bits_per_symbol);
        let slice = &encoded[byte_start..byte_end];

        let mut decoder = StreamingDecoder::new(
            Cursor::new(slice.to_vec()),
            alphabet,
            leading_skip,
            (15 - 1) as u64,
        );

        // Read into a 3-byte buffer repeatedly.
        let mut small_buf = [0u8; 3];
        let mut collected = Vec::new();
        loop {
            let n = decoder.read(&mut small_buf).unwrap();
            if n == 0 {
                break;
            }
            collected.extend_from_slice(&small_buf[..n]);
        }

        let expected = decode_substring_from_bytes(&encoded, 1, 15, alphabet);
        assert_eq!(collected, expected);
    }

    #[test]
    fn test_protein_5bit_unaligned() {
        let sequence = b"ACDEFGHIKLMNPQRSTVWY*X-";
        let alphabet = &PROTEIN_ALPHABET;
        // Test starts that are non-byte-aligned at 5 bits per symbol.
        for start in 0..sequence.len() {
            for end in (start + 1)..=sequence.len() {
                run_range_test(sequence, alphabet, start, end);
            }
        }
    }

    #[test]
    fn test_ascii_passthrough() {
        let sequence = b"The quick brown fox";
        let alphabet = &ASCII_ALPHABET;
        run_range_test(sequence, alphabet, 0, sequence.len());
        run_range_test(sequence, alphabet, 4, 9);
    }

    #[test]
    fn test_unexpected_eof() {
        // Provide a truncated source and ask for more bases than fit.
        let sequence = b"ACGTACGT";
        let alphabet = &DNA_2BIT_ALPHABET;
        let encoded = encode_sequence(sequence, alphabet);
        // Only feed 1 byte (4 bases worth) but ask for 8 bases.
        let mut decoder = StreamingDecoder::new(
            Cursor::new(encoded[..1].to_vec()),
            alphabet,
            0,
            8,
        );
        let mut buf = [0u8; 16];
        let err = decoder.read(&mut buf).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    }
}

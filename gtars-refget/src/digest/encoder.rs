//! Sequence encoding and decoding - WASM-safe, no filesystem dependencies.

use super::alphabet;
use super::alphabet::Alphabet;

/// A struct for encoding biological sequences into a compact, bit-packed format.
///
/// The `SequenceEncoder` is designed to efficiently encode sequences using a specified
/// alphabet type (e.g., DNA, protein, or ASCII). It uses a bit buffer to pack symbols
/// into bytes, minimizing memory usage while maintaining fast encoding performance.
///
/// # Fields
///
/// * `encoding_array` - A lookup table that maps input bytes to their encoded values.
/// * `bits_per_symbol` - The number of bits used to represent each symbol in the sequence.
/// * `encoded_sequence` - A vector that stores the resulting bit-packed sequence.
/// * `bit_pos` - The current bit position in the encoded sequence.
/// * `buffer` - An internal bit buffer used for packing symbols into bytes.
/// * `buffer_bits` - The number of bits currently stored in the buffer.
///
/// # Usage
///
/// 1. Create a new `SequenceEncoder` using the `new` method.
/// 2. Add sequence chunks with the `update` method.
/// 3. Call the `finalize` method to retrieve the encoded sequence as a `Vec<u8>`.
///
pub struct SequenceEncoder {
    alphabet: &'static alphabet::Alphabet, // alphabet used by this encoder
    encoded_sequence: Vec<u8>,             // resulting encoded sequence
    bit_pos: usize,                        // internal bit position in encoded sequence
    buffer: u64,                           // internal bit buffer
    buffer_bits: usize,                    // number of bits currently in buffer
}

impl SequenceEncoder {
    pub fn new(alphabet_type: alphabet::AlphabetType, length: usize) -> Self {
        let alphabet = alphabet::lookup_alphabet(&alphabet_type);
        let bits_per_symbol = alphabet.bits_per_symbol;
        let estimated_bytes = (length * bits_per_symbol).div_ceil(8);

        SequenceEncoder {
            alphabet,
            encoded_sequence: Vec::with_capacity(estimated_bytes),
            bit_pos: 0,
            buffer: 0,
            buffer_bits: 0,
        }
    }

    pub fn update(&mut self, sequence: &[u8]) {
        for &byte in sequence {
            let code = self.alphabet.encoding_array[byte as usize] as u64;
            self.buffer = (self.buffer << self.alphabet.bits_per_symbol) | code;
            self.buffer_bits += self.alphabet.bits_per_symbol;

            while self.buffer_bits >= 8 {
                self.buffer_bits -= 8;
                let out_byte = (self.buffer >> self.buffer_bits) as u8;
                self.encoded_sequence.push(out_byte);
                self.bit_pos += 8;
                self.buffer &= (1 << self.buffer_bits) - 1; // Mask to keep remaining bits
            }
        }
    }

    pub fn finalize(mut self) -> Vec<u8> {
        if self.buffer_bits > 0 {
            let out_byte = (self.buffer << (8 - self.buffer_bits)) as u8;
            self.encoded_sequence.push(out_byte);
            self.bit_pos += self.buffer_bits;
        }
        self.encoded_sequence
    }
}

/// Encodes a sequence using the specified encoding array.
///
/// This function takes a sequence of bytes and encodes it using the provided
/// encoding array. The encoding is done by converting each byte in the sequence
/// to its corresponding code in the encoding array, and then packing those
/// codes into a byte array. The number of bits used to represent each symbol
/// is specified by the `bits_per_symbol` parameter. The function returns
/// a bit-packed vector of bytes containing the encoded sequence.
///
/// **Bit Ordering: MSB-first (Most Significant Bit first)**
///
/// Example with 2-bit DNA encoding (A=00, C=01, G=10, T=11):
/// - Sequence "ACGT" → byte 0x1B (00011011)
/// - A (00) in bits 7-6, C (01) in bits 5-4, G (10) in bits 3-2, T (11) in bits 1-0
///
/// # Arguments
///
/// * `sequence` - The sequence to encode
/// * `alphabet` - The alphabet defining the encoding
///
/// # Returns
///
/// A vector containing the encoded (bit-packed) sequence
pub fn encode_sequence<T: AsRef<[u8]>>(sequence: T, alphabet: &Alphabet) -> Vec<u8> {
    let sequence = sequence.as_ref();
    let total_bits = sequence.len() * alphabet.bits_per_symbol;
    let mut bytes = vec![0u8; total_bits.div_ceil(8)];

    let mut bit_index = 0;
    for &byte in sequence {
        let code = alphabet.encoding_array[byte as usize];

        for i in (0..alphabet.bits_per_symbol).rev() {
            let bit = (code >> i) & 1;
            let byte_index = bit_index / 8;
            let bit_offset = 7 - (bit_index % 8); // MSB-first in each byte
            bytes[byte_index] |= bit << bit_offset;
            bit_index += 1;
        }
    }
    bytes
}

/// Returns the `[byte_start, byte_end)` range in the encoded file that contains
/// all bits for bases `[start, end)`.
///
/// `byte_start` is floor-aligned; `byte_end` is ceiling-aligned (exclusive).
/// Use this to compute an HTTP `Range` header (or a partial read) before fetching
/// only the bytes that cover a query window, then pass `byte_start` as the
/// `byte_offset` to [`decode_substring_from_bytes_at_offset`].
pub fn byte_range_for_bases(start: usize, end: usize, bits_per_symbol: usize) -> (usize, usize) {
    let byte_start = (start * bits_per_symbol) / 8;
    let byte_end = (end * bits_per_symbol + 7) / 8;
    (byte_start, byte_end)
}

/// Decodes a substring from a bit-packed encoded sequence.
///
/// Decodes a substring from a bit-packed sequence of bytes.
/// Decodes via the decoding array, which maps encoded
/// values back to original symbols.
///
/// # Arguments
///
/// * `encoded_bytes` - A slice of bytes containing the bit-packed encoded sequence.
/// * `start` - The start index of the substring (inclusive), in terms of symbols.
/// * `end` - The end index of the substring (exclusive), in terms of symbols.
/// * `bits_per_symbol` - The number of bits used to represent each symbol in the sequence.
/// * `decoding_array` - A lookup table that maps encoded values back to their original symbols.
///
/// # Returns
///
/// A `Vec<u8>` containing the decoded substring.
pub fn decode_substring_from_bytes(
    encoded_bytes: &[u8],
    start: usize,
    end: usize,
    alphabet: &Alphabet,
) -> Vec<u8> {
    decode_substring_from_bytes_at_offset(encoded_bytes, 0, start, end, alphabet)
}

/// Naive, bit-by-bit reference decoder. Kept private for differential testing.
///
/// This is the original implementation; the fast paths above must produce
/// byte-identical output to this function for all inputs.
#[cfg_attr(not(test), allow(dead_code))]
fn decode_substring_naive(
    encoded_bytes: &[u8],
    start: usize,
    end: usize,
    alphabet: &Alphabet,
) -> Vec<u8> {
    decode_substring_naive_at_offset(encoded_bytes, 0, start, end, alphabet)
}

/// Naive, bit-by-bit reference decoder honoring a `byte_offset`. Kept private
/// for differential testing against the fast accumulator path.
#[cfg_attr(not(test), allow(dead_code))]
fn decode_substring_naive_at_offset(
    bytes: &[u8],
    byte_offset: usize,
    start: usize,
    end: usize,
    alphabet: &Alphabet,
) -> Vec<u8> {
    let mut decoded = Vec::with_capacity(end - start);
    for i in start..end {
        let bit_offset = i * alphabet.bits_per_symbol;
        let mut code = 0u8;

        for j in 0..alphabet.bits_per_symbol {
            let bit_pos = bit_offset + j;
            let byte_index = bit_pos / 8;
            let bit_in_byte = 7 - (bit_pos % 8); // MSB0

            let bit = if byte_index >= byte_offset && byte_index - byte_offset < bytes.len() {
                (bytes[byte_index - byte_offset] >> bit_in_byte) & 1
            } else {
                0
            };

            code = (code << 1) | bit;
        }
        decoded.push(alphabet.decoding_array[code as usize]);
    }
    decoded
}

/// Decodes bases `[start, end)` from a partial encoded buffer.
///
/// `bytes` contains encoded data starting at `byte_offset` in the full encoded
/// file (i.e. `bytes[0]` corresponds to byte `byte_offset` of the complete
/// encoding). The function subtracts `byte_offset` from every computed byte
/// index so the bit arithmetic stays correct against the partial slice.
///
/// Use [`byte_range_for_bases`] with `alphabet.bits_per_symbol` to compute the
/// `(byte_offset, _)` range before fetching the partial buffer.
pub fn decode_substring_from_bytes_at_offset(
    bytes: &[u8],
    byte_offset: usize,
    start: usize,
    end: usize,
    alphabet: &Alphabet,
) -> Vec<u8> {
    if end <= start {
        return Vec::new();
    }
    let bits = alphabet.bits_per_symbol;
    let decoding_array = alphabet.decoding_array;

    // Specialized, byte-aligned 2-bit fast path: exactly 4 bases per byte.
    if bits == 2 {
        return decode_2bit(bytes, byte_offset, start, end, decoding_array);
    }
    if bits == 3 {
        return decode_3bit(bytes, byte_offset, start, end, decoding_array);
    }

    // General rolling-accumulator path for any bit width (covers anything else).
    // Produces byte-identical output to the naive decoder, including zero-fill
    // semantics past the end / before byte_offset.
    decode_rolling(bytes, byte_offset, start, end, bits, decoding_array)
}

/// Fetch the encoded byte for absolute byte index `abs_idx`, honoring the
/// `byte_offset` window and zero-filling anything outside `[byte_offset, byte_offset+len)`.
#[inline(always)]
fn fetch_byte(bytes: &[u8], byte_offset: usize, abs_idx: usize) -> u8 {
    if abs_idx >= byte_offset {
        let rel = abs_idx - byte_offset;
        if rel < bytes.len() {
            return bytes[rel];
        }
    }
    0
}

/// Specialized 2-bit decoder. 2 bits per symbol is byte-aligned (4 bases/byte),
/// so the aligned middle is a straight table lookup; partial head/tail bases are
/// handled with the rolling accumulator.
fn decode_2bit(
    bytes: &[u8],
    byte_offset: usize,
    start: usize,
    end: usize,
    decoding_array: &[u8; 256],
) -> Vec<u8> {
    // Precompute byte -> 4 decoded symbols (MSB0: base 0 is bits 7-6, etc).
    let mut table = [[0u8; 4]; 256];
    for (b, row) in table.iter_mut().enumerate() {
        for (k, slot) in row.iter_mut().enumerate() {
            let code = (b >> (6 - 2 * k)) & 0b11;
            *slot = decoding_array[code];
        }
    }

    let mut decoded = Vec::with_capacity(end - start);

    let mut i = start;
    // Head: align `i` up to a multiple of 4 (the next byte boundary).
    let aligned_start = i.div_ceil(4) * 4;
    let head_end = aligned_start.min(end);
    while i < head_end {
        let code = read_code(bytes, byte_offset, i * 2, 2);
        decoded.push(decoding_array[code as usize]);
        i += 1;
    }

    // Middle: whole bytes, 4 bases each, via the table.
    let aligned_end = end - (end % 4);
    while i < aligned_end {
        let abs_byte = (i * 2) / 8; // == i/4
        let byte = fetch_byte(bytes, byte_offset, abs_byte);
        decoded.extend_from_slice(&table[byte as usize]);
        i += 4;
    }

    // Tail: leftover bases.
    while i < end {
        let code = read_code(bytes, byte_offset, i * 2, 2);
        decoded.push(decoding_array[code as usize]);
        i += 1;
    }

    decoded
}

/// Specialized 3-bit decoder. 3 bits/symbol is unaligned, but lcm(3, 8) = 24
/// bits = 3 bytes = exactly 8 symbols, so the aligned middle (base indices that
/// are multiples of 8) decodes 8 symbols per 3 bytes with no inner refill loop.
/// Head/tail bases and any byte not fully in-bounds fall back to the scalar path.
fn decode_3bit(
    bytes: &[u8],
    byte_offset: usize,
    start: usize,
    end: usize,
    decoding_array: &[u8; 256],
) -> Vec<u8> {
    let mut decoded = Vec::with_capacity(end - start);

    let mut i = start;

    // Head: advance to the next base index that is a multiple of 8 (aligned to a
    // 3-byte boundary), using the scalar reader.
    let aligned_start = i.div_ceil(8) * 8;
    let head_end = aligned_start.min(end);
    while i < head_end {
        let code = read_code(bytes, byte_offset, i * 3, 3);
        decoded.push(decoding_array[code as usize]);
        i += 1;
    }

    // Middle: blocks of 8 symbols from 3 bytes. Only run while all 3 bytes are
    // fully in-bounds (abs byte index in [byte_offset, byte_offset+len)).
    let aligned_end = end - (end % 8);
    while i < aligned_end {
        let abs_byte = (i * 3) / 8; // == i/8*3
        // Need abs_byte, abs_byte+1, abs_byte+2 all in window.
        if abs_byte < byte_offset || abs_byte + 3 > byte_offset + bytes.len() {
            break;
        }
        let rel = abs_byte - byte_offset;
        let w = ((bytes[rel] as u32) << 16)
            | ((bytes[rel + 1] as u32) << 8)
            | (bytes[rel + 2] as u32);
        // Extract 8 codes of 3 bits, MSB first (bits 23..21, 20..18, ...).
        decoded.push(decoding_array[((w >> 21) & 0b111) as usize]);
        decoded.push(decoding_array[((w >> 18) & 0b111) as usize]);
        decoded.push(decoding_array[((w >> 15) & 0b111) as usize]);
        decoded.push(decoding_array[((w >> 12) & 0b111) as usize]);
        decoded.push(decoding_array[((w >> 9) & 0b111) as usize]);
        decoded.push(decoding_array[((w >> 6) & 0b111) as usize]);
        decoded.push(decoding_array[((w >> 3) & 0b111) as usize]);
        decoded.push(decoding_array[(w & 0b111) as usize]);
        i += 8;
    }

    // Tail (and any aligned blocks skipped because they touched the buffer edge):
    // scalar reader with zero-fill.
    while i < end {
        let code = read_code(bytes, byte_offset, i * 3, 3);
        decoded.push(decoding_array[code as usize]);
        i += 1;
    }

    decoded
}

/// Read a single `bits`-wide code starting at absolute bit position `bit_offset`
/// (MSB0), zero-filling bits outside the buffer window.
#[inline(always)]
fn read_code(bytes: &[u8], byte_offset: usize, bit_offset: usize, bits: usize) -> u8 {
    let mut code = 0u8;
    for j in 0..bits {
        let bit_pos = bit_offset + j;
        let abs_byte = bit_pos / 8;
        let bit_in_byte = 7 - (bit_pos % 8);
        let byte = fetch_byte(bytes, byte_offset, abs_byte);
        let bit = (byte >> bit_in_byte) & 1;
        code = (code << 1) | bit;
    }
    code
}

/// General rolling bit-accumulator decoder for arbitrary bit widths.
///
/// Maintains a `u64` window fed one byte at a time (MSB-first). Symbols are
/// extracted from the top of the window with shifts/masks, avoiding the
/// per-bit `/8` and `%8` of the naive path. Bits past the buffer window are
/// zero (we feed 0 bytes), matching the naive zero-fill semantics exactly.
fn decode_rolling(
    bytes: &[u8],
    byte_offset: usize,
    start: usize,
    end: usize,
    bits: usize,
    decoding_array: &[u8; 256],
) -> Vec<u8> {
    let n = end - start;
    let mut decoded = Vec::with_capacity(n);

    let first_bit = start * bits;
    // Absolute byte index we will load next into the accumulator.
    let mut next_byte = first_bit / 8;
    // How many bits at the front of that first byte we must discard.
    let skip_bits = first_bit % 8;

    let mut acc: u64 = 0;
    let mut acc_bits: usize = 0;

    // Prime the accumulator and drop the leading skip bits.
    while acc_bits < skip_bits {
        acc = (acc << 8) | fetch_byte(bytes, byte_offset, next_byte) as u64;
        next_byte += 1;
        acc_bits += 8;
    }
    acc_bits -= skip_bits;
    // Mask off the discarded high bits so only `acc_bits` low bits are valid.
    acc &= (1u64 << acc_bits) - 1;

    let mask = (1u64 << bits) - 1;

    // The bytes available directly from the slice are at absolute indices
    // `[byte_offset, byte_offset + bytes.len())`. While `next_byte` is in that
    // window we can index the slice directly (no zero-fill branch). We need at
    // least one full byte of headroom before refilling, so guard accordingly.
    //
    // Process the bulk where every refill byte is guaranteed in-bounds, then
    // fall back to the branchy path for the tail that touches the zero-fill
    // region (or the very edges).
    let slice_end_abs = byte_offset + bytes.len(); // exclusive, absolute

    // How many output symbols can be produced while only consuming bytes that
    // are strictly in-bounds? Each symbol needs `bits` bits; we may refill up to
    // (bits) extra bits worth (<= 1 byte) per symbol. Stay conservative: stop the
    // fast loop a few symbols before the slice runs out so the slow path handles
    // the boundary precisely.
    let mut produced = 0usize;
    if next_byte >= byte_offset && next_byte < slice_end_abs {
        // Largest number of symbols we can emit while next_byte stays < slice_end_abs.
        // After emitting k symbols we've consumed at most ceil((skip_bits=0 here is
        // already handled)+k*bits)/8 bytes beyond the primed state. Compute the bit
        // budget available from in-bounds bytes and derive a safe symbol count.
        let inbounds_bits = (slice_end_abs - next_byte) * 8 + acc_bits;
        // Reserve nothing extra: a symbol is safe iff after consuming it we never
        // had to refill from out-of-bounds. We can emit floor(inbounds_bits / bits)
        // symbols safely.
        let fast_syms = (inbounds_bits / bits).min(n);

        while produced < fast_syms {
            while acc_bits < bits {
                // SAFETY/CORRECTNESS: next_byte < slice_end_abs here because the
                // fast_syms budget guarantees enough in-bounds bits remain.
                let rel = next_byte - byte_offset;
                acc = (acc << 8) | bytes[rel] as u64;
                next_byte += 1;
                acc_bits += 8;
            }
            acc_bits -= bits;
            let code = ((acc >> acc_bits) & mask) as usize;
            decoded.push(decoding_array[code]);
            produced += 1;
        }
    }

    // Slow tail: edges and the zero-fill region beyond the buffer.
    for _ in produced..n {
        while acc_bits < bits {
            acc = (acc << 8) | fetch_byte(bytes, byte_offset, next_byte) as u64;
            next_byte += 1;
            acc_bits += 8;
        }
        acc_bits -= bits;
        let code = ((acc >> acc_bits) & mask) as usize;
        decoded.push(decoding_array[code]);
    }

    decoded
}

pub fn decode_string_from_bytes(
    encoded_bytes: &[u8],
    seq_len: usize,
    alphabet: &Alphabet,
) -> Vec<u8> {
    decode_substring_from_bytes(encoded_bytes, 0, seq_len, alphabet)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dna_2bit_encoding() {
        let alphabet = &alphabet::DNA_2BIT_ALPHABET;
        let sequence = b"ACGT";
        let encoded = encode_sequence(sequence, alphabet);
        let ans = [0b10, 0b01, 0b11, 0b00];
        let packed: Vec<u8> = ans
            .chunks(8 / alphabet.bits_per_symbol) // Number of symbols that fit in a byte
            .map(|chunk| {
                chunk
                    .iter()
                    .fold(0, |acc, &b| (acc << alphabet.bits_per_symbol) | b)
            })
            .collect();
        assert_eq!(encoded, packed);

        let decoded: Vec<u8> = decode_substring_from_bytes(&encoded, 0, sequence.len(), alphabet);
        assert_eq!(decoded, sequence);
    }

    #[test]
    fn test_dna_iupac_encoding() {
        let sequence = b"ACGTRYMK";
        let alphabet = &alphabet::DNA_IUPAC_ALPHABET;
        let encoded = encode_sequence(sequence, alphabet);
        let ans = [
            0b0001, 0b0010, 0b0100, 0b1000, 0b0101, 0b1010, 0b0011, 0b0111,
        ];
        let packed: Vec<u8> = ans
            .chunks(8 / alphabet.bits_per_symbol) // Number of symbols that fit in a byte
            .map(|chunk| {
                chunk
                    .iter()
                    .fold(0, |acc, &b| (acc << alphabet.bits_per_symbol) | b)
            })
            .collect();
        assert_eq!(encoded, packed);
        let decoded: Vec<u8> = decode_substring_from_bytes(&encoded, 0, sequence.len(), alphabet);
        assert_eq!(decoded, sequence);
    }

    #[test]
    fn test_protein_encoding() {
        let sequence = b"ACDEFGHIKLMNPQRSTVWY*X-";
        let alphabet = &alphabet::PROTEIN_ALPHABET;
        let encoded = encode_sequence(sequence, alphabet);
        // Don't want to re-implement bit-packing here for 5-bit symbols, so just check the length.
        assert_eq!(
            encoded.len(),
            (sequence.len() * alphabet.bits_per_symbol).div_ceil(8)
        );
        let decoded: Vec<u8> = decode_substring_from_bytes(&encoded, 0, sequence.len(), alphabet);
        assert_eq!(decoded, sequence);
    }

    #[test]
    fn test_ascii_encoding() {
        let sequence = b"Hello, World!";
        let alphabet = &alphabet::ASCII_ALPHABET;
        let encoded = encode_sequence(sequence, alphabet);
        let decoded = decode_substring_from_bytes(&encoded, 0, sequence.len(), alphabet);
        assert_eq!(decoded, sequence);
    }

    #[test]
    fn test_dna_3bit_encoding() {
        let sequence = b"ACGTNRYX"; // 8 chars * 3 bits/char = 24 bits.
        let alphabet = &alphabet::DNA_3BIT_ALPHABET;
        let encoded = encode_sequence(sequence, alphabet);
        // let ans =  vec![0b000, 0b001, 0b010, 0b011, 0b100, 0b101, 0b110, 0b111];
        let packed = vec![0b00000101, 0b00111001, 0b01110111]; //manually bit-packed above 3-bits
        assert_eq!(encoded, packed);
        let decoded: Vec<u8> = decode_substring_from_bytes(&encoded, 0, sequence.len(), alphabet);
        assert_eq!(decoded, sequence);
    }

    #[test]
    fn test_byte_range_for_bases() {
        // Dna2bit: 2 bits/base, 4 bases/byte
        assert_eq!(byte_range_for_bases(0, 4, 2), (0, 1));
        assert_eq!(byte_range_for_bases(1, 5, 2), (0, 2)); // start=1 is at bit 2, still byte 0
        assert_eq!(byte_range_for_bases(4, 8, 2), (1, 2));
        // Dna3bit: 3 bits/base
        assert_eq!(byte_range_for_bases(0, 8, 3), (0, 3));
        assert_eq!(byte_range_for_bases(3, 6, 3), (1, 3)); // bit 9 to bit 18; bytes 1..3
    }

    fn check_offset_roundtrip(sequence: &[u8], alphabet: &Alphabet, start: usize, end: usize) {
        let encoded = encode_sequence(sequence, alphabet);
        let (byte_start, byte_end) = byte_range_for_bases(start, end, alphabet.bits_per_symbol);
        let partial = &encoded[byte_start..byte_end];

        let from_partial =
            decode_substring_from_bytes_at_offset(partial, byte_start, start, end, alphabet);
        let from_full = decode_substring_from_bytes(&encoded, start, end, alphabet);

        assert_eq!(from_partial, from_full, "partial decode must match full decode");
        assert_eq!(from_partial, &sequence[start..end], "decode must match source");
    }

    #[test]
    fn test_decode_at_offset_dna_2bit() {
        let sequence = b"ACGTACGTACGTACGT";
        let alphabet = &alphabet::DNA_2BIT_ALPHABET;
        check_offset_roundtrip(sequence, alphabet, 0, 4);
        check_offset_roundtrip(sequence, alphabet, 1, 5);
        check_offset_roundtrip(sequence, alphabet, 4, 8);
        check_offset_roundtrip(sequence, alphabet, 5, 13);
        check_offset_roundtrip(sequence, alphabet, 3, 16);
    }

    #[test]
    fn test_decode_at_offset_dna_3bit() {
        let sequence = b"ACGTNRYXACGTNRYX";
        let alphabet = &alphabet::DNA_3BIT_ALPHABET;
        check_offset_roundtrip(sequence, alphabet, 0, 8);
        check_offset_roundtrip(sequence, alphabet, 3, 6);
        check_offset_roundtrip(sequence, alphabet, 5, 11);
        check_offset_roundtrip(sequence, alphabet, 2, 16);
    }

    /// Tiny deterministic xorshift64 RNG so the differential test is reproducible
    /// without pulling in an external rng crate.
    struct Rng(u64);
    impl Rng {
        fn new(seed: u64) -> Self {
            Rng(seed | 1)
        }
        fn next_u64(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }
        fn below(&mut self, bound: usize) -> usize {
            (self.next_u64() % bound as u64) as usize
        }
    }

    /// Differential test: the fast decoders must be byte-identical to the naive
    /// reference for both 2-bit (byte-aligned) and 3-bit (unaligned) alphabets,
    /// across many random (start, end) and byte_offset combinations, INCLUDING
    /// ranges that run at and beyond the buffer's base count (zero-fill region).
    #[test]
    fn test_fast_decoder_matches_naive_differential() {
        let mut rng = Rng::new(0xDEADBEEF_CAFEF00D);

        for &alphabet in &[
            &alphabet::DNA_2BIT_ALPHABET,
            &alphabet::DNA_3BIT_ALPHABET,
            &alphabet::DNA_IUPAC_ALPHABET, // 4-bit, also unaligned in spots
            &alphabet::PROTEIN_ALPHABET,   // 5-bit
        ] {
            let bits = alphabet.bits_per_symbol;

            // Try several buffer sizes, including tiny and empty-ish ones.
            for &n_bases in &[0usize, 1, 3, 4, 5, 7, 8, 17, 64, 257, 1000] {
                // Build a random encoded buffer of n_bases symbols by emitting
                // random raw bytes and encoding them. We need the raw bytes to be
                // valid symbols for the alphabet so encode/decode is meaningful,
                // but for the differential test we only compare fast vs naive, so
                // we can also just fabricate random encoded bytes directly.
                let n_bytes = (n_bases * bits).div_ceil(8);
                let mut encoded = vec![0u8; n_bytes];
                for b in encoded.iter_mut() {
                    *b = rng.next_u64() as u8;
                }

                // Max base count representable; allow start/end to exceed it to
                // exercise zero-fill.
                let max_base = n_bases + 16;

                for _ in 0..200 {
                    let a = rng.below(max_base + 1);
                    let b = rng.below(max_base + 1);
                    let (start, end) = if a <= b { (a, b) } else { (b, a) };

                    // Full-buffer variant.
                    let fast = decode_substring_from_bytes(&encoded, start, end, alphabet);
                    let naive = decode_substring_naive(&encoded, start, end, alphabet);
                    assert_eq!(
                        fast, naive,
                        "full decode mismatch: bits={bits} n_bases={n_bases} start={start} end={end}"
                    );

                    // _at_offset variant: choose a byte_offset and pass the
                    // corresponding partial slice. byte_offset can be anywhere
                    // in [0, n_bytes]; the slice is the remaining bytes.
                    let byte_offset = if n_bytes == 0 {
                        0
                    } else {
                        rng.below(n_bytes + 1)
                    };
                    // Also test partial slices that don't reach the end.
                    let slice_end = if byte_offset >= n_bytes {
                        byte_offset
                    } else {
                        byte_offset + rng.below(n_bytes - byte_offset + 1)
                    };
                    let partial = &encoded[byte_offset.min(n_bytes)..slice_end.min(n_bytes)];

                    let fast_off = decode_substring_from_bytes_at_offset(
                        partial, byte_offset, start, end, alphabet,
                    );
                    let naive_off = decode_substring_naive_at_offset(
                        partial, byte_offset, start, end, alphabet,
                    );
                    assert_eq!(
                        fast_off, naive_off,
                        "at_offset decode mismatch: bits={bits} n_bases={n_bases} \
                         byte_offset={byte_offset} slice_len={} start={start} end={end}",
                        partial.len()
                    );
                }

                // decode_string_from_bytes must also match the naive whole-seq decode.
                let fast_full = decode_string_from_bytes(&encoded, n_bases, alphabet);
                let naive_full = decode_substring_naive(&encoded, 0, n_bases, alphabet);
                assert_eq!(
                    fast_full, naive_full,
                    "decode_string mismatch: bits={bits} n_bases={n_bases}"
                );
            }
        }
    }
}

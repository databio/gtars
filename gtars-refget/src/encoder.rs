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
    let mut decoded = Vec::with_capacity(end - start);
    for i in start..end {
        let bit_offset = i * alphabet.bits_per_symbol;
        let mut code = 0u8;

        for j in 0..alphabet.bits_per_symbol {
            let bit_pos = bit_offset + j;
            let byte_index = bit_pos / 8;
            let bit_in_byte = 7 - (bit_pos % 8); // MSB0

            let bit = if byte_index < encoded_bytes.len() {
                (encoded_bytes[byte_index] >> bit_in_byte) & 1
            } else {
                0
            };

            code = (code << 1) | bit;
        }
        decoded.push(alphabet.decoding_array[code as usize]);
    }
    decoded
}

pub fn decode_string_from_bytes(
    encoded_bytes: &[u8],
    seq_len: usize,
    alphabet: &Alphabet,
) -> Vec<u8> {
    let mut decoded = Vec::with_capacity(seq_len);
    for i in 0..seq_len {
        let bit_offset = i * alphabet.bits_per_symbol;
        let mut code = 0u8;

        for j in 0..alphabet.bits_per_symbol {
            let bit_pos = bit_offset + j;
            let byte_index = bit_pos / 8;
            let bit_in_byte = 7 - (bit_pos % 8); // MSB0

            let bit = if byte_index < encoded_bytes.len() {
                (encoded_bytes[byte_index] >> bit_in_byte) & 1
            } else {
                0
            };

            code = (code << 1) | bit;
        }
        decoded.push(alphabet.decoding_array[code as usize]);
    }
    decoded
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
}

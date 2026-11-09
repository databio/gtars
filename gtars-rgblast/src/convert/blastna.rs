//! gtars packed record -> blastna letters, straight from the packed codes.
//!
//! The traceback stage reads subjects in blastna (A0 C1 G2 T3 R4 Y5 M6 K7 W8 S9 B10 D11 H12
//! V13 N14). Decoding through ASCII (`decode_substring_from_bytes`, then a per-letter table
//! lookup pushed one byte at a time) was 40% of a single-threaded hg38 megablast search; this
//! writes the letters directly with tables built from gtars' alphabet definitions at compile
//! time, several codes per lookup. The code layout (first code in the high bits of a byte,
//! 3-bit codes spanning bytes) is the one `convert::conv2/conv3/conv4` read.

use crate::engine::encoding::ASCII_TO_BLASTNA;
use gtars_refget::digest::{
    Alphabet, AlphabetType, DNA_2BIT_ALPHABET, DNA_3BIT_ALPHABET, DNA_IUPAC_ALPHABET,
};

/// gtars code -> blastna for the first `n` codes of an alphabet.
const fn code_table(a: &Alphabet, n: usize) -> [u8; 16] {
    let mut t = [0u8; 16];
    let mut code = 0;
    while code < n {
        let letter = a.decoding_array[code].to_ascii_uppercase();
        let v = ASCII_TO_BLASTNA[letter as usize];
        assert!(v != 255, "gtars alphabet symbol has no blastna value");
        t[code] = v;
        code += 1;
    }
    t
}

const B2: [u8; 16] = code_table(&DNA_2BIT_ALPHABET, 4);
const B3: [u8; 16] = code_table(&DNA_3BIT_ALPHABET, 8);
const B4: [u8; 16] = code_table(&DNA_IUPAC_ALPHABET, 16);

/// One 2-bit byte -> four letters.
const LUT2: [[u8; 4]; 256] = {
    let mut t = [[0u8; 4]; 256];
    let mut b = 0;
    while b < 256 {
        t[b] = [B2[b >> 6], B2[(b >> 4) & 3], B2[(b >> 2) & 3], B2[b & 3]];
        b += 1;
    }
    t
};

/// Twelve bits of 3-bit codes -> four letters.
const LUT3: [[u8; 4]; 4096] = {
    let mut t = [[0u8; 4]; 4096];
    let mut w = 0;
    while w < 4096 {
        t[w] = [B3[w >> 9], B3[(w >> 6) & 7], B3[(w >> 3) & 7], B3[w & 7]];
        w += 1;
    }
    t
};

/// One 4-bit byte -> two letters.
const LUT4: [[u8; 2]; 256] = {
    let mut t = [[0u8; 2]; 256];
    let mut b = 0;
    while b < 256 {
        t[b] = [B4[b >> 4], B4[b & 15]];
        b += 1;
    }
    t
};

#[inline(always)]
fn code3_at(packed: &[u8], p: usize) -> usize {
    let bit = 3 * p;
    let i = bit / 8;
    let hi = packed[i] as u32;
    let lo = *packed.get(i + 1).unwrap_or(&0) as u32;
    ((((hi << 8) | lo) >> (13 - (bit % 8))) & 7) as usize
}

/// Append letters `lo..hi` of a packed record to `out` in blastna. The caller has checked
/// that `packed` holds `len` codes of `alphabet` and that `lo <= hi <= len`.
pub fn packed_to_blastna(
    alphabet: AlphabetType,
    packed: &[u8],
    lo: usize,
    hi: usize,
    out: &mut Vec<u8>,
) {
    out.reserve(hi - lo + 1);
    let mut p = lo;
    match alphabet {
        AlphabetType::Dna2bit => {
            while p < hi && !p.is_multiple_of(4) {
                out.push(B2[((packed[p / 4] >> (6 - 2 * (p % 4))) & 3) as usize]);
                p += 1;
            }
            while p + 4 <= hi {
                out.extend_from_slice(&LUT2[packed[p / 4] as usize]);
                p += 4;
            }
            while p < hi {
                out.push(B2[((packed[p / 4] >> (6 - 2 * (p % 4))) & 3) as usize]);
                p += 1;
            }
        }
        AlphabetType::Dna3bit => {
            while p < hi && !p.is_multiple_of(8) {
                out.push(B3[code3_at(packed, p)]);
                p += 1;
            }
            while p + 8 <= hi {
                let b = &packed[3 * (p / 8)..3 * (p / 8) + 3];
                let w = (b[0] as usize) << 16 | (b[1] as usize) << 8 | b[2] as usize;
                out.extend_from_slice(&LUT3[w >> 12]);
                out.extend_from_slice(&LUT3[w & 0xFFF]);
                p += 8;
            }
            while p < hi {
                out.push(B3[code3_at(packed, p)]);
                p += 1;
            }
        }
        _ => {
            if p < hi && !p.is_multiple_of(2) {
                out.push(B4[(packed[p / 2] & 15) as usize]);
                p += 1;
            }
            while p + 2 <= hi {
                out.extend_from_slice(&LUT4[packed[p / 2] as usize]);
                p += 2;
            }
            if p < hi {
                out.push(B4[(packed[p / 2] >> 4) as usize]);
            }
        }
    }
}

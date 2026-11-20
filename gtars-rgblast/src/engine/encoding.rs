//! Sequence encodings used by the engine (BLAST+ `algo/blast/core/blast_encoding.c`,
//! `blast_util.c::GetReverseNuclSequence`).
//!
//! blastna: A0 C1 G2 T3 R4 Y5 M6 K7 W8 S9 B10 D11 H12 V13 N14, 15 = gap/sentinel.
//! ncbi2na: 4 bases per byte, first base in the high bits, A0 C1 G2 T3.

pub const NUCL_SENTINEL: u8 = 15;
pub const NUCL_MASK: u8 = 14;
pub const FENCE_SENTRY: u8 = 201;
pub const COMPRESSION_RATIO: i32 = 4;

/// `IUPACNA_TO_BLASTNA` for ASCII letters, both cases; U/u -> T. 255 = not a nucleotide letter.
pub const ASCII_TO_BLASTNA: [u8; 256] = {
    let mut t = [255u8; 256];
    let pairs: [(u8, u8); 16] = [
        (b'A', 0),
        (b'C', 1),
        (b'G', 2),
        (b'T', 3),
        (b'R', 4),
        (b'Y', 5),
        (b'M', 6),
        (b'K', 7),
        (b'W', 8),
        (b'S', 9),
        (b'B', 10),
        (b'D', 11),
        (b'H', 12),
        (b'V', 13),
        (b'N', 14),
        (b'U', 3),
    ];
    let mut i = 0;
    while i < 16 {
        let (c, v) = pairs[i];
        t[c as usize] = v;
        t[(c + 32) as usize] = v;
        i += 1;
    }
    t[b'-' as usize] = 15;
    t
};

pub const BLASTNA_TO_IUPAC: [u8; 16] = *b"ACGTRYMKWSBDHVN-";

/// Complement of a blastna code (via ncbi4na, as `GetReverseNuclSequence` does it).
pub const BLASTNA_COMPLEMENT: [u8; 16] = {
    // blastna -> ncbi4na
    let to4: [u8; 16] = [1, 2, 4, 8, 5, 10, 3, 12, 9, 6, 14, 13, 11, 7, 15, 0];
    // ncbi4na complement (GetReverseNuclSequence's conversion_table)
    let comp4: [u8; 16] = [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15];
    // ncbi4na -> blastna
    let from4: [u8; 16] = [15, 0, 1, 6, 2, 4, 9, 13, 3, 8, 5, 12, 7, 11, 10, 14];
    let mut t = [0u8; 16];
    let mut i = 0;
    while i < 16 {
        t[i] = from4[comp4[to4[i] as usize] as usize];
        i += 1;
    }
    t
};

/// Reverse complement of a blastna slice (FENCE_SENTRY kept in place).
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        out.push(if b == FENCE_SENTRY {
            FENCE_SENTRY
        } else {
            BLASTNA_COMPLEMENT[(b & 15) as usize]
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complement_table() {
        assert_eq!(BLASTNA_COMPLEMENT[0], 3);
        assert_eq!(BLASTNA_COMPLEMENT[1], 2);
        assert_eq!(BLASTNA_COMPLEMENT[4], 5); // R <-> Y
        assert_eq!(BLASTNA_COMPLEMENT[6], 7); // M <-> K
        assert_eq!(BLASTNA_COMPLEMENT[8], 8); // W
        assert_eq!(BLASTNA_COMPLEMENT[9], 9); // S
        assert_eq!(BLASTNA_COMPLEMENT[10], 13); // B <-> V
        assert_eq!(BLASTNA_COMPLEMENT[11], 12); // D <-> H
        assert_eq!(BLASTNA_COMPLEMENT[14], 14);
        assert_eq!(BLASTNA_COMPLEMENT[15], 15);
        assert_eq!(ASCII_TO_BLASTNA[b'u' as usize], 3);
        assert_eq!(ASCII_TO_BLASTNA[b'x' as usize], 255);
        assert_eq!(reverse_complement(&[0, 0, 1, 2]), vec![1, 2, 3, 3]);
    }
}

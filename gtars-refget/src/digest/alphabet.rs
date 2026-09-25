//! Alphabet types and encoding - WASM-safe, no filesystem dependencies.
//!
//! Each bit-packed alphabet is declared once, as a list of `(symbol, code)`
//! pairs. Its encoding table, decoding table and membership set are all built
//! from that list at compile time, so the three can never disagree.

use serde::{Deserialize, Serialize};
use std::fmt::Display;
use std::str::FromStr;

/// Encoding-table value for a byte that is not a member of the alphabet.
///
/// Every real code is below `0x80` (at most 5 bits), so the encoder can detect
/// a non-member byte by OR-ing codes together and testing the high bit.
pub const INVALID_CODE: u8 = 0xFF;

/// Represents an alphabet with its encoding and decoding arrays.
pub struct Alphabet {
    // Maps the alphabet type to its encoding and decoding arrays
    pub alphabet_type: AlphabetType,
    /// Maps an input byte to its code, or [`INVALID_CODE`] for non-members.
    pub encoding_array: &'static [u8; 256],
    pub decoding_array: &'static [u8; 256],
    /// `membership[b]` is true when byte `b` is a symbol of this alphabet.
    pub membership: &'static [bool; 256],
    pub bits_per_symbol: usize,
}

impl Alphabet {
    /// Returns true when `b` is a symbol of this alphabet (so it round-trips).
    pub fn contains(&self, b: u8) -> bool {
        self.membership[b as usize]
    }
}

/// Bit in the guesser accumulator for the Ascii alphabet (always set).
const ASCII_BIT: u8 = 1 << 4;
/// Initial guesser accumulator: every alphabet is still possible.
const ALL_ALPHABETS: u8 = 0x1F;
/// How many bytes the guesser folds before it checks for an early exit.
const GUESS_CHUNK: usize = 4096;

/// A struct to guess alphabet types based on the sequence content.
///
/// Selection rule: the result is the first alphabet in [`ALPHABET_ORDER`]
/// (Dna2bit, Dna3bit, DnaIupac, Protein, Ascii) that holds every byte seen.
/// The alphabets are not strictly nested (for example `X` is in Dna3bit and
/// Protein but not in DnaIupac), so the guesser intersects membership sets
/// instead of assuming each alphabet is a superset of the one before it.
///
/// Input must already be uppercased, as every ingest path does. Alphabets hold
/// uppercase symbols only, so a lowercase byte selects Ascii, which round-trips
/// it exactly.
///
/// This struct handles a sequence as a stream, to preserve memory when dealing
/// with large sequences, and can run alongside digesting.
pub struct AlphabetGuesser {
    /// Bit `i` is set while alphabet `ALPHABET_ORDER[i]` can still hold every
    /// byte seen so far.
    acc: u8,
}

impl AlphabetGuesser {
    /// Creates a new AlphabetGuesser with every alphabet still possible.
    #[allow(clippy::new_without_default)]
    pub fn new() -> Self {
        AlphabetGuesser { acc: ALL_ALPHABETS }
    }

    /// Folds `sequence` (already uppercased) into the guess.
    pub fn update(&mut self, sequence: &[u8]) {
        if self.acc == ASCII_BIT {
            return;
        }

        for chunk in sequence.chunks(GUESS_CHUNK) {
            let mut acc = self.acc;
            for &byte in chunk {
                acc &= ALPHABET_MASK[byte as usize];
            }
            self.acc = acc;
            if acc == ASCII_BIT {
                break;
            }
        }
    }

    /// The first alphabet in [`ALPHABET_ORDER`] that holds every byte seen.
    /// An empty sequence gives Dna2bit.
    pub fn guess(&self) -> AlphabetType {
        ALPHABET_ORDER[self.acc.trailing_zeros() as usize]
    }
}

impl Default for AlphabetGuesser {
    fn default() -> Self {
        Self::new()
    }
}

/// Represents the type of alphabet used for sequence encoding.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum AlphabetType {
    /// 2-bit DNA encoding (A, C, G, T only)
    Dna2bit,
    /// 3-bit DNA encoding (A, C, G, T, N, R, Y, X)
    Dna3bit,
    /// IUPAC DNA/RNA encoding: A C G T U plus ambiguity codes R Y S W K M B D H V N (4 bits/symbol)
    DnaIupac,
    /// Protein encoding (5 bits/symbol): the 20 standard amino acids, U
    /// (selenocysteine), O (pyrrolysine), the ambiguity codes B, Z, J, X, stop
    /// `*` and gaps `-` `.`
    Protein,
    /// ASCII encoding (for general text)
    Ascii,
    /// Unknown alphabet type
    Unknown,
}

impl Display for AlphabetType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AlphabetType::Dna2bit => write!(f, "dna2bit"),
            AlphabetType::Dna3bit => write!(f, "dna3bit"),
            AlphabetType::DnaIupac => write!(f, "dnaio"),
            AlphabetType::Protein => write!(f, "protein"),
            AlphabetType::Ascii => write!(f, "ASCII"),
            AlphabetType::Unknown => write!(f, "Unknown"),
        }
    }
}

impl FromStr for AlphabetType {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "dna2bit" => Ok(AlphabetType::Dna2bit),
            "dna3bit" => Ok(AlphabetType::Dna3bit),
            "dnaio" => Ok(AlphabetType::DnaIupac),
            "protein" => Ok(AlphabetType::Protein),
            "ascii" => Ok(AlphabetType::Ascii),
            "unknown" => Ok(AlphabetType::Unknown),
            _ => Err(()),
        }
    }
}

impl AlphabetType {
    /// Returns the number of bits used per symbol for this alphabet type
    pub fn bits_per_symbol(&self) -> usize {
        match self {
            AlphabetType::Dna2bit => 2,
            AlphabetType::Dna3bit => 3,
            AlphabetType::DnaIupac => 4,
            AlphabetType::Protein => 5,
            AlphabetType::Ascii => 8,
            AlphabetType::Unknown => 8, // Default to 8 bits for unknown
        }
    }
}

// ============================================================================
// Symbol lists: the single source of truth for each bit-packed alphabet.
//
// Each list is `(uppercase symbol, code)`. The encode table, decode table and
// membership set are all derived from it. Lists hold uppercase symbols only;
// ingest uppercases before guessing and encoding, and decoding yields
// uppercase.
//
// The layouts are frozen: encoded files already on disk use these codes.
// Never renumber an existing symbol. Only a code that no input ever produced
// may be given to a new symbol.
// ============================================================================

/// UCSC 2bit layout: T=00, C=01, A=10, G=11.
const DNA_2BIT_SYMBOLS: &[(u8, u8)] = &[(b'T', 0b00), (b'C', 0b01), (b'A', 0b10), (b'G', 0b11)];

/// 3-bit DNA: A C G T plus N, R, Y and X.
const DNA_3BIT_SYMBOLS: &[(u8, u8)] = &[
    (b'A', 0b000),
    (b'C', 0b001),
    (b'G', 0b010),
    (b'T', 0b011),
    (b'N', 0b100),
    (b'R', 0b101),
    (b'Y', 0b110),
    (b'X', 0b111),
];

/// IUPAC DNA/RNA symbols mapped to a 4-bit code (16 possible values).
///
/// **This is NOT an IUPAC bitmask.** The codes do not follow A=1, C=2, G=4,
/// T=8 union semantics -- if they did, K (G or T) would be `0b1100`, not
/// `0b0111`; B (C, G, or T) would be `0b1110`, not `0b1100`; N (any base)
/// would be `0b1111`, not `0b0000`; and so on. Treat this as an arbitrary
/// one-to-one mapping between 16 symbols and 16 codes.
///
/// Full code table:
/// `N=0000 A=0001 C=0010 M=0011 G=0100 R=0101 S=0110 K=0111 T=1000 W=1001
/// Y=1010 U=1011 B=1100 D=1101 H=1110 V=1111`.
const DNA_IUPAC_SYMBOLS: &[(u8, u8)] = &[
    (b'N', 0b0000), // any base
    (b'A', 0b0001), // adenine
    (b'C', 0b0010), // cytosine
    (b'M', 0b0011), // A or C
    (b'G', 0b0100), // guanine
    (b'R', 0b0101), // A or G
    (b'S', 0b0110), // G or C
    (b'K', 0b0111), // G or T
    (b'T', 0b1000), // thymine
    (b'W', 0b1001), // A or T
    (b'Y', 0b1010), // C or T
    (b'U', 0b1011), // RNA uracil; own code so it does not collapse to T
    (b'B', 0b1100), // not A
    (b'D', 0b1101), // not C
    (b'H', 0b1110), // not G
    (b'V', 0b1111), // not T
];

/// Protein: the 20 standard amino acids, stop, unknown, two gap symbols, and
/// the IUPAC-IUB extras U (selenocysteine), B (Asx), Z (Glx), O (pyrrolysine)
/// and J (Xle). The extras use codes that were unassigned before, so older
/// protein payloads decode unchanged. 0b11101-0b11111 stay unassigned.
const PROTEIN_SYMBOLS: &[(u8, u8)] = &[
    (b'A', 0b00000), // Alanine
    (b'C', 0b00001), // Cysteine
    (b'D', 0b00010), // Aspartic acid
    (b'E', 0b00011), // Glutamic acid
    (b'F', 0b00100), // Phenylalanine
    (b'G', 0b00101), // Glycine
    (b'H', 0b00110), // Histidine
    (b'I', 0b00111), // Isoleucine
    (b'K', 0b01000), // Lysine
    (b'L', 0b01001), // Leucine
    (b'M', 0b01010), // Methionine
    (b'N', 0b01011), // Asparagine
    (b'P', 0b01100), // Proline
    (b'Q', 0b01101), // Glutamine
    (b'R', 0b01110), // Arginine
    (b'S', 0b01111), // Serine
    (b'T', 0b10000), // Threonine
    (b'V', 0b10001), // Valine
    (b'W', 0b10010), // Tryptophan
    (b'Y', 0b10011), // Tyrosine
    (b'*', 0b10100), // Stop codon
    (b'X', 0b10101), // Unknown
    (b'-', 0b10110), // Gap
    (b'.', 0b10111), // Gap
    (b'U', 0b11000), // Selenocysteine
    (b'B', 0b11001), // Asx (D or N)
    (b'Z', 0b11010), // Glx (E or Q)
    (b'O', 0b11011), // Pyrrolysine
    (b'J', 0b11100), // Xle (L or I)
];

/// Encoding table from a symbol list. Non-members map to [`INVALID_CODE`].
const fn build_encoding(symbols: &[(u8, u8)]) -> [u8; 256] {
    let mut arr = [INVALID_CODE; 256];
    let mut i = 0;
    while i < symbols.len() {
        arr[symbols[i].0 as usize] = symbols[i].1;
        i += 1;
    }
    arr
}

/// Decoding table from a symbol list. Codes no symbol uses decode to
/// `default`; the encoder never writes them.
const fn build_decoding(symbols: &[(u8, u8)], default: u8) -> [u8; 256] {
    let mut arr = [default; 256];
    let mut i = 0;
    while i < symbols.len() {
        arr[symbols[i].1 as usize] = symbols[i].0;
        i += 1;
    }
    arr
}

/// Membership set from a symbol list. Never looks at code values, so a
/// symbol whose code is 0 (Protein `A`, DnaIupac `N`) is still a member.
const fn build_membership(symbols: &[(u8, u8)]) -> [bool; 256] {
    let mut arr = [false; 256];
    let mut i = 0;
    while i < symbols.len() {
        arr[symbols[i].0 as usize] = true;
        i += 1;
    }
    arr
}

/// Compile-time sanity check of a symbol list: no duplicate symbols, no
/// duplicate codes, every code fits in `bits` bits, no lowercase symbols.
const fn check_symbols(symbols: &[(u8, u8)], bits: u32) {
    let mut i = 0;
    while i < symbols.len() {
        let (sym, code) = symbols[i];
        assert!(code < (1u8 << bits), "code does not fit in bits_per_symbol");
        assert!(
            !sym.is_ascii_lowercase(),
            "symbol lists hold uppercase only"
        );
        let mut j = i + 1;
        while j < symbols.len() {
            assert!(symbols[j].0 != sym, "duplicate symbol");
            assert!(symbols[j].1 != code, "duplicate code");
            j += 1;
        }
        i += 1;
    }
}

const _: () = {
    check_symbols(DNA_2BIT_SYMBOLS, 2);
    check_symbols(DNA_3BIT_SYMBOLS, 3);
    check_symbols(DNA_IUPAC_SYMBOLS, 4);
    check_symbols(PROTEIN_SYMBOLS, 5);
};

const DNA_2BIT_ENCODING_ARRAY: [u8; 256] = build_encoding(DNA_2BIT_SYMBOLS);
const DNA_2BIT_DECODING_ARRAY: [u8; 256] = build_decoding(DNA_2BIT_SYMBOLS, b'N');
const DNA_2BIT_MEMBERSHIP: [bool; 256] = build_membership(DNA_2BIT_SYMBOLS);

const DNA_3BIT_ENCODING_ARRAY: [u8; 256] = build_encoding(DNA_3BIT_SYMBOLS);
const DNA_3BIT_DECODING_ARRAY: [u8; 256] = build_decoding(DNA_3BIT_SYMBOLS, b'X');
const DNA_3BIT_MEMBERSHIP: [bool; 256] = build_membership(DNA_3BIT_SYMBOLS);

const DNA_IUPAC_ENCODING_ARRAY: [u8; 256] = build_encoding(DNA_IUPAC_SYMBOLS);
const DNA_IUPAC_DECODING_ARRAY: [u8; 256] = build_decoding(DNA_IUPAC_SYMBOLS, b'N');
const DNA_IUPAC_MEMBERSHIP: [bool; 256] = build_membership(DNA_IUPAC_SYMBOLS);

const PROTEIN_ENCODING_ARRAY: [u8; 256] = build_encoding(PROTEIN_SYMBOLS);
const PROTEIN_DECODING_ARRAY: [u8; 256] = build_decoding(PROTEIN_SYMBOLS, b'X');
const PROTEIN_MEMBERSHIP: [bool; 256] = build_membership(PROTEIN_SYMBOLS);

const fn const_u8_array() -> [u8; 256] {
    let mut arr = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        arr[i] = i as u8;
        i += 1;
    }
    arr
}

// Simple 8-bit passthrough for ASCII: every byte is a member.
const ASCII_ENCODING_ARRAY: [u8; 256] = const_u8_array();
const ASCII_MEMBERSHIP: [bool; 256] = [true; 256];

pub const DNA_3BIT_ALPHABET: Alphabet = Alphabet {
    alphabet_type: AlphabetType::Dna3bit,
    bits_per_symbol: 3,
    encoding_array: &DNA_3BIT_ENCODING_ARRAY,
    decoding_array: &DNA_3BIT_DECODING_ARRAY,
    membership: &DNA_3BIT_MEMBERSHIP,
};

pub const DNA_2BIT_ALPHABET: Alphabet = Alphabet {
    alphabet_type: AlphabetType::Dna2bit,
    bits_per_symbol: 2,
    encoding_array: &DNA_2BIT_ENCODING_ARRAY,
    decoding_array: &DNA_2BIT_DECODING_ARRAY,
    membership: &DNA_2BIT_MEMBERSHIP,
};

pub const DNA_IUPAC_ALPHABET: Alphabet = Alphabet {
    alphabet_type: AlphabetType::DnaIupac,
    bits_per_symbol: 4,
    encoding_array: &DNA_IUPAC_ENCODING_ARRAY,
    decoding_array: &DNA_IUPAC_DECODING_ARRAY,
    membership: &DNA_IUPAC_MEMBERSHIP,
};

pub const PROTEIN_ALPHABET: Alphabet = Alphabet {
    alphabet_type: AlphabetType::Protein,
    bits_per_symbol: 5,
    encoding_array: &PROTEIN_ENCODING_ARRAY,
    decoding_array: &PROTEIN_DECODING_ARRAY,
    membership: &PROTEIN_MEMBERSHIP,
};

pub const ASCII_ALPHABET: Alphabet = Alphabet {
    alphabet_type: AlphabetType::Ascii,
    bits_per_symbol: 8,
    encoding_array: &ASCII_ENCODING_ARRAY,
    decoding_array: &ASCII_ENCODING_ARRAY,
    membership: &ASCII_MEMBERSHIP,
};

/// Look up the alphabet for a given alphabet type.
pub fn lookup_alphabet(alphabet_type: &AlphabetType) -> &'static Alphabet {
    match alphabet_type {
        AlphabetType::Dna2bit => &DNA_2BIT_ALPHABET,
        AlphabetType::Dna3bit => &DNA_3BIT_ALPHABET,
        AlphabetType::DnaIupac => &DNA_IUPAC_ALPHABET,
        AlphabetType::Protein => &PROTEIN_ALPHABET,
        AlphabetType::Ascii => &ASCII_ALPHABET,
        AlphabetType::Unknown => &ASCII_ALPHABET, // Default to ASCII for unknown
    }
}

/// The order in which alphabets are tried, smallest first. The guesser picks
/// the first one that holds every byte of the sequence.
pub const ALPHABET_ORDER: [AlphabetType; 5] = [
    AlphabetType::Dna2bit,
    AlphabetType::Dna3bit,
    AlphabetType::DnaIupac,
    AlphabetType::Protein,
    AlphabetType::Ascii,
];

/// For each byte, bit `i` is set when `ALPHABET_ORDER[i]` holds that byte.
/// Bit 4 (Ascii) is always set.
const ALPHABET_MASK: [u8; 256] = {
    let memberships: [&[bool; 256]; 5] = [
        &DNA_2BIT_MEMBERSHIP,
        &DNA_3BIT_MEMBERSHIP,
        &DNA_IUPAC_MEMBERSHIP,
        &PROTEIN_MEMBERSHIP,
        &ASCII_MEMBERSHIP,
    ];
    let mut mask = [0u8; 256];
    let mut b = 0;
    while b < 256 {
        let mut i = 0;
        while i < 5 {
            if memberships[i][b] {
                mask[b] |= 1 << i;
            }
            i += 1;
        }
        b += 1;
    }
    mask
};

/// Guesses the smallest alphabet that can represent every byte of `sequence`.
///
/// Returns the first alphabet in [`ALPHABET_ORDER`] that holds every byte.
/// `sequence` must already be uppercased; see [`AlphabetGuesser`].
pub fn guess_alphabet(sequence: &[u8]) -> AlphabetType {
    let mut guesser = AlphabetGuesser::new();
    guesser.update(sequence);
    guesser.guess()
}

#[cfg(test)]
mod tests {
    use super::{
        ALPHABET_ORDER, AlphabetGuesser, AlphabetType, DNA_IUPAC_DECODING_ARRAY,
        DNA_IUPAC_ENCODING_ARRAY, INVALID_CODE, guess_alphabet, lookup_alphabet,
    };

    #[test]
    fn test_dna_iupac_tables_are_inverse() {
        use std::collections::HashSet;

        for &b in b"ACGTURYSWKMBDHVN" {
            let code = DNA_IUPAC_ENCODING_ARRAY[b as usize];
            assert_eq!(
                DNA_IUPAC_DECODING_ARRAY[code as usize], b,
                "decode(encode({})) did not round-trip",
                b as char
            );
        }

        // Alphabets hold uppercase only: ingest uppercases before encoding,
        // so lowercase bytes are not members and must not get a code.
        for &b in b"acgturyswkmbdhvn" {
            assert_eq!(
                DNA_IUPAC_ENCODING_ARRAY[b as usize], INVALID_CODE,
                "lowercase {} must not be a DnaIupac member",
                b as char
            );
        }

        let codes: HashSet<u8> = b"ACGTURYSWKMBDHVN"
            .iter()
            .map(|&b| DNA_IUPAC_ENCODING_ARRAY[b as usize])
            .collect();
        assert_eq!(
            codes.len(),
            16,
            "the 16 IUPAC codes must be pairwise distinct"
        );
        assert!(codes.iter().all(|&c| c < 16));

        // Pin the frozen layout: guards against "fixing" the codes into real
        // IUPAC bitmasks, which would break existing on-disk stores.
        const EXPECTED_CODES: [(u8, u8); 16] = [
            (b'N', 0b0000),
            (b'A', 0b0001),
            (b'C', 0b0010),
            (b'M', 0b0011),
            (b'G', 0b0100),
            (b'R', 0b0101),
            (b'S', 0b0110),
            (b'K', 0b0111),
            (b'T', 0b1000),
            (b'W', 0b1001),
            (b'Y', 0b1010),
            (b'U', 0b1011),
            (b'B', 0b1100),
            (b'D', 0b1101),
            (b'H', 0b1110),
            (b'V', 0b1111),
        ];
        for &(symbol, expected_code) in EXPECTED_CODES.iter() {
            assert_eq!(
                DNA_IUPAC_ENCODING_ARRAY[symbol as usize], expected_code,
                "frozen layout changed for symbol {}",
                symbol as char
            );
        }
    }

    #[test]
    fn test_guess_alphabet() {
        assert_eq!(guess_alphabet(b"ACGT"), AlphabetType::Dna2bit);
        assert_eq!(guess_alphabet(b"ACGTNRY"), AlphabetType::Dna3bit);
        assert_eq!(guess_alphabet(b"ACGTRYMK"), AlphabetType::DnaIupac);
        // For protein test, use a sequence with characters that are only in the protein alphabet
        assert_eq!(guess_alphabet(b"EFILPQ"), AlphabetType::Protein);
        assert_eq!(guess_alphabet(b"Hello, World!"), AlphabetType::Ascii);

        assert_eq!(guess_alphabet(b"ACTGEG"), AlphabetType::Protein);

        // Test cases where the original guess_alphabet would fail but guess_alphabet succeeds
        // ACTGM: Contains 'M' which is IUPAC DNA, but original would stop at 3-bit
        assert_eq!(guess_alphabet(b"ACTGM"), AlphabetType::DnaIupac);

        // ACGTSKWV: Contains IUPAC codes that would force upgrade to IUPAC
        assert_eq!(guess_alphabet(b"ACGTSKWV"), AlphabetType::DnaIupac);

        // ACGTE: Contains 'E' which is protein-only, should jump directly to Protein
        assert_eq!(guess_alphabet(b"ACGTE"), AlphabetType::Protein);

        // ACGT*: Contains stop codon, should be protein
        assert_eq!(guess_alphabet(b"ACGT*"), AlphabetType::Protein);

        // ACGT-: Contains gap character, should be protein
        assert_eq!(guess_alphabet(b"ACGT-"), AlphabetType::Protein);

        // Lowercase is not in any bit-packed alphabet (input must be
        // uppercased first), so it falls to Ascii, which round-trips exactly.
        assert_eq!(guess_alphabet(b"actgEFIL"), AlphabetType::Ascii);

        // Test with non-standard characters that should force ASCII
        assert_eq!(guess_alphabet(b"ACGT123"), AlphabetType::Ascii);
        assert_eq!(guess_alphabet(b"ACGT@#$"), AlphabetType::Ascii);

        // Edge cases with single characters
        assert_eq!(guess_alphabet(b"A"), AlphabetType::Dna2bit);
        assert_eq!(guess_alphabet(b"N"), AlphabetType::Dna3bit);
        assert_eq!(guess_alphabet(b"M"), AlphabetType::DnaIupac);
        assert_eq!(guess_alphabet(b"E"), AlphabetType::Protein);
        assert_eq!(guess_alphabet(b"1"), AlphabetType::Ascii);
    }

    #[test]
    fn test_tables_are_inverse_for_every_alphabet() {
        for alphabet_type in ALPHABET_ORDER {
            let alphabet = lookup_alphabet(&alphabet_type);
            for b in 0u8..=255 {
                let code = alphabet.encoding_array[b as usize];
                if alphabet.membership[b as usize] {
                    assert!(alphabet.contains(b));
                    assert_eq!(
                        alphabet.decoding_array[code as usize], b,
                        "{alphabet_type:?}: decode(encode(0x{b:02x})) did not round-trip"
                    );
                } else if alphabet_type != AlphabetType::Ascii {
                    assert!(!alphabet.contains(b));
                    assert_eq!(
                        code, INVALID_CODE,
                        "{alphabet_type:?}: non-member 0x{b:02x} must encode to INVALID_CODE"
                    );
                }
            }
        }
        // Ascii holds every byte.
        let ascii = lookup_alphabet(&AlphabetType::Ascii);
        assert!((0u8..=255).all(|b| ascii.contains(b)));
    }

    #[test]
    fn test_guess_protein_extras_and_u_b() {
        // Protein-only letters (E, F, L) plus U: the chosen alphabet must hold U.
        assert_eq!(guess_alphabet(b"MEFLU"), AlphabetType::Protein);
        // M, G, C, R and U are all IUPAC nucleotide letters, so this is DnaIupac.
        assert_eq!(guess_alphabet(b"MGCRU"), AlphabetType::DnaIupac);
        assert_eq!(guess_alphabet(b"ACGTU"), AlphabetType::DnaIupac);
        assert_eq!(guess_alphabet(b"U"), AlphabetType::DnaIupac);
        assert_eq!(guess_alphabet(b"B"), AlphabetType::DnaIupac);
        assert_eq!(guess_alphabet(b"MEFBZOJ"), AlphabetType::Protein);
        assert_eq!(guess_alphabet(b"ACGTNX"), AlphabetType::Dna3bit);
        assert_eq!(guess_alphabet(b"actg"), AlphabetType::Ascii);
        assert_eq!(guess_alphabet(b""), AlphabetType::Dna2bit);
    }

    #[test]
    fn test_guess_all_iupac_nucleotides_roundtrips() {
        use crate::digest::encoder::{decode_string_from_bytes, encode_sequence};

        let seq = b"MCVTYHNGTGYC";
        let alphabet_type = guess_alphabet(seq);
        assert_eq!(alphabet_type, AlphabetType::DnaIupac);
        let alphabet = lookup_alphabet(&alphabet_type);
        let encoded = encode_sequence(seq, alphabet).unwrap();
        assert_eq!(decode_string_from_bytes(&encoded, seq.len(), alphabet), seq);
    }

    #[test]
    fn test_guess_single_byte_matches_naive_reference() {
        for b in 0u8..=255 {
            let expected = ALPHABET_ORDER
                .into_iter()
                .find(|t| lookup_alphabet(t).membership[b as usize])
                .unwrap();
            assert_eq!(
                guess_alphabet(&[b]),
                expected,
                "single byte 0x{b:02x} picked the wrong alphabet"
            );
        }

        // Every pair of bytes: the first alphabet that holds both.
        for a in 0u8..=255 {
            for b in 0u8..=255 {
                let expected = ALPHABET_ORDER
                    .into_iter()
                    .find(|t| {
                        let m = lookup_alphabet(t).membership;
                        m[a as usize] && m[b as usize]
                    })
                    .unwrap();
                assert_eq!(
                    guess_alphabet(&[a, b]),
                    expected,
                    "pair 0x{a:02x} 0x{b:02x} picked the wrong alphabet"
                );
            }
        }
    }

    #[test]
    fn test_alphabet_guesser_matches_guess_alphabet() {
        let test_cases: Vec<&[u8]> = vec![
            b"ACGT",
            b"ACGTNRY",
            b"ACGTRYMK",
            b"EFILPQ",
            b"Hello, World!",
            b"ACTGEG",
            b"ACTGM",
            b"ACGTSKWV",
            b"ACGTE",
            b"ACGT*",
            b"ACGT-",
            b"actgEFIL",
            b"ACGT123",
            b"ACGT@#$",
            b"A",
            b"N",
            b"M",
            b"E",
            b"1",
        ];

        for test_case in test_cases {
            let mut guesser = AlphabetGuesser::new();
            guesser.update(test_case);
            let guesser_result = guesser.guess();
            let function_result = guess_alphabet(test_case);

            assert_eq!(
                guesser_result,
                function_result,
                "AlphabetGuesser and guess_alphabet disagree on sequence: {:?}",
                std::str::from_utf8(test_case).unwrap_or("(invalid UTF-8)")
            );
        }
    }

    #[test]
    fn test_alphabet_guesser_streaming() {
        // Test that the guesser works correctly when fed data in chunks
        let mut guesser = AlphabetGuesser::new();

        // Start with DNA 2-bit
        guesser.update(b"ACGT");
        assert_eq!(guesser.guess(), AlphabetType::Dna2bit);

        // Add some 3-bit characters
        guesser.update(b"NRY");
        assert_eq!(guesser.guess(), AlphabetType::Dna3bit);

        // Add IUPAC characters
        guesser.update(b"MKS");
        assert_eq!(guesser.guess(), AlphabetType::DnaIupac);

        // Add protein characters
        guesser.update(b"EFIL");
        assert_eq!(guesser.guess(), AlphabetType::Protein);

        // Add ASCII characters
        guesser.update(b"123");
        assert_eq!(guesser.guess(), AlphabetType::Ascii);

        // Verify it matches the full sequence result
        let full_sequence = b"ACGTNRYMSKEFILP123";
        assert_eq!(guesser.guess(), guess_alphabet(full_sequence));
    }
}

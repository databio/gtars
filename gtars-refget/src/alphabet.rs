use serde::{Deserialize, Serialize};
use std::fmt::Display;
use std::str::FromStr;

/// Represents an alphabet with its encoding and decoding arrays.
pub struct Alphabet {
    // Maps the alphabet type to its encoding and decoding arrays
    pub alphabet_type: AlphabetType,
    pub encoding_array: &'static [u8; 256],
    pub decoding_array: &'static [u8; 256],
    pub bits_per_symbol: usize,
}

/// A struct to guess alphabet types based on the sequence content.
///
/// This struct is meant to handle a sequence as a stream, to preserve memory
/// when dealing with large sequences.
/// It is ameable to processing along with digesting the sequence.
pub struct AlphabetGuesser {
    alphabet_type: AlphabetType,
}

impl AlphabetGuesser {
    /// Creates a new AlphabetGuesser with the initial alphabet type set to Dna2bit.
    #[allow(clippy::new_without_default)]
    pub fn new() -> Self {
        AlphabetGuesser {
            alphabet_type: AlphabetType::Dna2bit,
        }
    }

    pub fn update(&mut self, sequence: &[u8]) {
        if self.alphabet_type == AlphabetType::Ascii {
            return;
        }

        for &byte in sequence {
            let byte_upper = byte.to_ascii_uppercase();
            let char_required_alphabet = get_minimum_alphabet_for_char(byte_upper);

            // Upgrade to the more general alphabet if needed
            if is_more_general_alphabet(char_required_alphabet, self.alphabet_type) {
                self.alphabet_type = char_required_alphabet;
            }

            // Early exit if we've reached ASCII
            if self.alphabet_type == AlphabetType::Ascii {
                break;
            }
        }
    }

    pub fn guess(&self) -> AlphabetType {
        self.alphabet_type
    }
}

/// Represents the type of alphabet used for sequence encoding.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum AlphabetType {
    /// 2-bit DNA encoding (A, C, G, T only)
    Dna2bit,
    /// 3-bit DNA encoding (A, C, G, T, N, R, Y, and others as X)
    Dna3bit,
    /// IUPAC DNA encoding (includes ambiguity codes)
    DnaIupac,
    /// Protein sequence encoding
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

/// A lookup table that maps ASCII characters representing DNA bases
/// (T, C, A, G, and their lowercase counterparts) to their 2-bit encoding.
/// This constant is used for efficient conversion of DNA sequences into a
/// more compact 2-bit representation.
const DNA_2BIT_ENCODING_ARRAY: [u8; 256] = {
    let mut arr = [0u8; 256];
    // UCSC 2bit encoding: T=00, C=01, A=10, G=11
    arr[b'T' as usize] = 0b00;
    arr[b't' as usize] = 0b00;
    arr[b'C' as usize] = 0b01;
    arr[b'c' as usize] = 0b01;
    arr[b'A' as usize] = 0b10;
    arr[b'a' as usize] = 0b10;
    arr[b'G' as usize] = 0b11;
    arr[b'g' as usize] = 0b11;
    arr
};

/// Reverse lookup: 2-bit value → representative DNA base
const DNA_2BIT_DECODING_ARRAY: [u8; 256] = {
    let mut arr = [b'N'; 256]; // Default to 'N' (unknown)
    arr[0b00] = b'T';
    arr[0b01] = b'C';
    arr[0b10] = b'A';
    arr[0b11] = b'G';
    arr
};

/// A lookup table that maps ASCII characters representing DNA and ambiguous
/// bases to their 3-bit encoding. This builds on the 2-bit encoding, adding
/// support for ambiguous bases (N, R, Y), and then encodes all other
/// characters as X.
const DNA_3BIT_ENCODING_ARRAY: [u8; 256] = {
    let mut arr = [0b111; 256]; // Default to 'X' (0b111)
    // Define the 8 possible values in 3-bit encoding
    arr[b'A' as usize] = 0b000;
    arr[b'a' as usize] = 0b000; // A
    arr[b'C' as usize] = 0b001;
    arr[b'c' as usize] = 0b001; // C
    arr[b'G' as usize] = 0b010;
    arr[b'g' as usize] = 0b010; // G
    arr[b'T' as usize] = 0b011;
    arr[b't' as usize] = 0b011; // T
    arr[b'N' as usize] = 0b100;
    arr[b'n' as usize] = 0b100; // N
    arr[b'R' as usize] = 0b101;
    arr[b'r' as usize] = 0b101; // R
    arr[b'Y' as usize] = 0b110;
    arr[b'y' as usize] = 0b110; // Y
    // All other characters will be encoded as 'X' (0b111)
    arr
};

const DNA_3BIT_DECODING_ARRAY: [u8; 256] = {
    let mut arr = [b'X'; 256]; // Default to 'X'
    arr[0b000] = b'A'; // A
    arr[0b001] = b'C'; // C
    arr[0b010] = b'G'; // G
    arr[0b011] = b'T'; // T
    arr[0b100] = b'N'; // N
    arr[0b101] = b'R'; // R
    arr[0b110] = b'Y'; // Y
    arr
};

/// A lookup table that maps ASCII characters representing DNA bases.
/// This encoding uses 4 bits to represent each base, allowing for
/// the representation of 16 different values.
const DNA_IUPAC_ENCODING_ARRAY: [u8; 256] = {
    let mut arr = [0u8; 256];
    arr[b'A' as usize] = 0b0001; // A
    arr[b'C' as usize] = 0b0010; // C
    arr[b'G' as usize] = 0b0100; // G
    arr[b'T' as usize] = 0b1000; // T
    arr[b'U' as usize] = 0b1000; // U (common in RNA)
    arr[b'R' as usize] = 0b0101; // A or G
    arr[b'Y' as usize] = 0b1010; // C or T
    arr[b'S' as usize] = 0b0110; // G or C
    arr[b'W' as usize] = 0b1001; // A or T
    arr[b'K' as usize] = 0b0111; // G or T
    arr[b'M' as usize] = 0b0011; // A or C
    arr[b'B' as usize] = 0b1100; // C or G or T
    arr[b'D' as usize] = 0b1101; // A or G or T
    arr[b'H' as usize] = 0b1110; // A or C or T
    arr[b'V' as usize] = 0b1111; // A or C or G
    arr[b'N' as usize] = 0b0000; // Any base
    // Add lowercase variants
    arr[b'a' as usize] = 0b0001;
    arr[b'c' as usize] = 0b0010;
    arr[b'g' as usize] = 0b0100;
    arr[b't' as usize] = 0b1000;
    arr[b'u' as usize] = 0b1000;
    arr[b'r' as usize] = 0b0101;
    arr[b'y' as usize] = 0b1010;
    arr[b's' as usize] = 0b0110;
    arr[b'w' as usize] = 0b1001;
    arr[b'k' as usize] = 0b0111;
    arr[b'm' as usize] = 0b0011;
    arr[b'b' as usize] = 0b1100;
    arr[b'd' as usize] = 0b1101;
    arr[b'h' as usize] = 0b1110;
    arr[b'v' as usize] = 0b1111;
    arr[b'n' as usize] = 0b0000;
    arr
};

/// Maps 4-bit IUPAC bit patterns (from 0b0000 to 0b1111) back to representative DNA characters.
/// Values outside the valid 4-bit range default to 'N'.
const DNA_IUPAC_DECODING_ARRAY: [u8; 256] = {
    let mut arr = [b'N'; 256]; // Default to 'N' for all
    arr[0b0000] = b'N'; // Any
    arr[0b0001] = b'A';
    arr[0b0010] = b'C';
    arr[0b0011] = b'M'; // A or C
    arr[0b0100] = b'G';
    arr[0b0101] = b'R'; // A or G
    arr[0b0110] = b'S'; // G or C
    arr[0b0111] = b'K'; // G or T
    arr[0b1000] = b'T';
    arr[0b1001] = b'W'; // A or T
    arr[0b1010] = b'Y'; // C or T
    arr[0b1011] = b'D'; // A or G or T
    arr[0b1100] = b'B'; // C or G or T
    arr[0b1101] = b'H'; // A or C or T
    arr[0b1110] = b'V'; // A or C or G
    arr[0b1111] = b'V'; // A or C or G (or use 'N' if you'd rather treat 0b1111 as invalid)
    arr
};

const PROTEIN_ENCODING_ARRAY: [u8; 256] = {
    let mut arr = [0u8; 256];
    // Standard amino acids (20) + special characters
    arr[b'A' as usize] = 0b00000;
    arr[b'a' as usize] = 0b00000; // Alanine
    arr[b'C' as usize] = 0b00001;
    arr[b'c' as usize] = 0b00001; // Cysteine
    arr[b'D' as usize] = 0b00010;
    arr[b'd' as usize] = 0b00010; // Aspartic acid
    arr[b'E' as usize] = 0b00011;
    arr[b'e' as usize] = 0b00011; // Glutamic acid
    arr[b'F' as usize] = 0b00100;
    arr[b'f' as usize] = 0b00100; // Phenylalanine
    arr[b'G' as usize] = 0b00101;
    arr[b'g' as usize] = 0b00101; // Glycine
    arr[b'H' as usize] = 0b00110;
    arr[b'h' as usize] = 0b00110; // Histidine
    arr[b'I' as usize] = 0b00111;
    arr[b'i' as usize] = 0b00111; // Isoleucine
    arr[b'K' as usize] = 0b01000;
    arr[b'k' as usize] = 0b01000; // Lysine
    arr[b'L' as usize] = 0b01001;
    arr[b'l' as usize] = 0b01001; // Leucine
    arr[b'M' as usize] = 0b01010;
    arr[b'm' as usize] = 0b01010; // Methionine
    arr[b'N' as usize] = 0b01011;
    arr[b'n' as usize] = 0b01011; // Asparagine
    arr[b'P' as usize] = 0b01100;
    arr[b'p' as usize] = 0b01100; // Proline
    arr[b'Q' as usize] = 0b01101;
    arr[b'q' as usize] = 0b01101; // Glutamine
    arr[b'R' as usize] = 0b01110;
    arr[b'r' as usize] = 0b01110; // Arginine
    arr[b'S' as usize] = 0b01111;
    arr[b's' as usize] = 0b01111; // Serine
    arr[b'T' as usize] = 0b10000;
    arr[b't' as usize] = 0b10000; // Threonine
    arr[b'V' as usize] = 0b10001;
    arr[b'v' as usize] = 0b10001; // Valine
    arr[b'W' as usize] = 0b10010;
    arr[b'w' as usize] = 0b10010; // Tryptophan
    arr[b'Y' as usize] = 0b10011;
    arr[b'y' as usize] = 0b10011; // Tyrosine
    arr[b'*' as usize] = 0b10100; // Stop codon
    arr[b'X' as usize] = 0b10101;
    arr[b'x' as usize] = 0b10101; // Unknown
    arr[b'-' as usize] = 0b10110; // Gap
    arr[b'.' as usize] = 0b10111; // Gap
    arr
};

const PROTEIN_DECODING_ARRAY: [u8; 256] = {
    let mut arr = [b'X'; 256]; // Default to 'X'
    arr[0b00000] = b'A'; // Alanine
    arr[0b00001] = b'C'; // Cysteine
    arr[0b00010] = b'D'; // Aspartic acid
    arr[0b00011] = b'E'; // Glutamic acid
    arr[0b00100] = b'F'; // Phenylalanine
    arr[0b00101] = b'G'; // Glycine
    arr[0b00110] = b'H'; // Histidine
    arr[0b00111] = b'I'; // Isoleucine
    arr[0b01000] = b'K'; // Lysine
    arr[0b01001] = b'L'; // Leucine
    arr[0b01010] = b'M'; // Methionine
    arr[0b01011] = b'N'; // Asparagine
    arr[0b01100] = b'P'; // Proline
    arr[0b01101] = b'Q'; // Glutamine
    arr[0b01110] = b'R'; // Arginine
    arr[0b01111] = b'S'; // Serine
    arr[0b10000] = b'T'; // Threonine
    arr[0b10001] = b'V'; // Valine
    arr[0b10010] = b'W'; // Tryptophan
    arr[0b10011] = b'Y'; // Tyrosine
    arr[0b10100] = b'*'; // Stop codon
    arr[0b10101] = b'X'; // Unknown
    arr[0b10110] = b'-'; // Gap
    arr[0b10111] = b'.'; // Gap
    // All other values default to 'X'
    arr
};

const fn const_u8_array() -> [u8; 256] {
    let mut arr = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        arr[i] = i as u8;
        i += 1;
    }
    arr
}

// Simple 8-bit passthrough for ASCII
const ASCII_ENCODING_ARRAY: [u8; 256] = const_u8_array();

pub const DNA_3BIT_ALPHABET: Alphabet = Alphabet {
    alphabet_type: AlphabetType::Dna3bit,
    bits_per_symbol: 3,
    encoding_array: &DNA_3BIT_ENCODING_ARRAY,
    decoding_array: &DNA_3BIT_DECODING_ARRAY,
};

pub const DNA_2BIT_ALPHABET: Alphabet = Alphabet {
    alphabet_type: AlphabetType::Dna2bit,
    bits_per_symbol: 2,
    encoding_array: &DNA_2BIT_ENCODING_ARRAY,
    decoding_array: &DNA_2BIT_DECODING_ARRAY,
};

pub const DNA_IUPAC_ALPHABET: Alphabet = Alphabet {
    alphabet_type: AlphabetType::DnaIupac,
    bits_per_symbol: 4,
    encoding_array: &DNA_IUPAC_ENCODING_ARRAY,
    decoding_array: &DNA_IUPAC_DECODING_ARRAY,
};

pub const PROTEIN_ALPHABET: Alphabet = Alphabet {
    alphabet_type: AlphabetType::Protein,
    bits_per_symbol: 5,
    encoding_array: &PROTEIN_ENCODING_ARRAY,
    decoding_array: &PROTEIN_DECODING_ARRAY,
};

pub const ASCII_ALPHABET: Alphabet = Alphabet {
    alphabet_type: AlphabetType::Ascii,
    bits_per_symbol: 8,
    encoding_array: &ASCII_ENCODING_ARRAY,
    decoding_array: &ASCII_ENCODING_ARRAY,
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

/// Guesses the most appropriate alphabet type for a given sequence.
///
/// Tries to determine if the sequence is DNA (2-bit or IUPAC), protein, or general ASCII.
///
/// # Arguments
///
/// * `sequence` - The sequence to analyze
///
/// # Returns
///
/// The most specific alphabet type that matches the sequence
pub fn guess_alphabet_fast(sequence: &[u8]) -> AlphabetType {
    // Start with the smallest alphabet
    let mut smallest_alphabet = AlphabetType::Dna2bit;

    for &byte in sequence {
        let byte_upper = byte.to_ascii_uppercase();

        // Check if the character fits in the current alphabet
        match smallest_alphabet {
            AlphabetType::Dna2bit => {
                // Check if valid for 2-bit DNA
                if !matches!(byte_upper, b'A' | b'C' | b'G' | b'T') {
                    // Upgrade to 3-bit DNA
                    smallest_alphabet = AlphabetType::Dna3bit;
                }
            }
            AlphabetType::Dna3bit => {
                // Check if valid for 3-bit DNA
                if !matches!(byte_upper, b'A' | b'C' | b'G' | b'T' | b'N' | b'R' | b'Y') {
                    // Upgrade to IUPAC DNA
                    smallest_alphabet = AlphabetType::DnaIupac;
                }
            }
            AlphabetType::DnaIupac => {
                // Check if valid for IUPAC DNA
                if DNA_IUPAC_ENCODING_ARRAY[byte as usize] == 0 && byte_upper != b'N' {
                    // Upgrade to Protein
                    smallest_alphabet = AlphabetType::Protein;
                }
            }
            AlphabetType::Protein => {
                // Check if valid for Protein
                if PROTEIN_ENCODING_ARRAY[byte as usize] == 0 && byte != b'-' && byte != b'*' {
                    // Upgrade to ASCII
                    smallest_alphabet = AlphabetType::Ascii;
                    // No need to check further characters
                    break;
                }
            }
            _ => break, // Already at ASCII, no need to check further
        }
    }

    smallest_alphabet
}

/// Guesses the most appropriate alphabet type for a given sequence (accurate version).
///
/// This version properly checks each character against all alphabet types to find
/// the minimum alphabet that can represent the entire sequence. It's slower than
/// the basic version but more accurate.
///
/// # Arguments
///
/// * `sequence` - The sequence to analyze
///
/// # Returns
///
/// The most specific alphabet type that can represent the entire sequence
pub fn guess_alphabet(sequence: &[u8]) -> AlphabetType {
    let mut required_alphabet = AlphabetType::Dna2bit;

    for &byte in sequence {
        let byte_upper = byte.to_ascii_uppercase();
        let char_required_alphabet = get_minimum_alphabet_for_char(byte_upper);

        // Upgrade to the more general alphabet if needed
        if is_more_general_alphabet(char_required_alphabet, required_alphabet) {
            required_alphabet = char_required_alphabet;
        }

        // Early exit if we've reached ASCII
        if required_alphabet == AlphabetType::Ascii {
            break;
        }
    }

    required_alphabet
}

/// Determines the minimum alphabet required to represent a single character.
fn get_minimum_alphabet_for_char(byte: u8) -> AlphabetType {
    // Check 2-bit DNA first (most restrictive)
    if matches!(byte, b'A' | b'C' | b'G' | b'T') {
        return AlphabetType::Dna2bit;
    }

    // Check 3-bit DNA
    if matches!(byte, b'N' | b'R' | b'Y') {
        return AlphabetType::Dna3bit;
    }

    // Check IUPAC DNA
    if DNA_IUPAC_ENCODING_ARRAY[byte as usize] != 0 || byte == b'N' {
        return AlphabetType::DnaIupac;
    }

    // Check Protein
    if PROTEIN_ENCODING_ARRAY[byte as usize] != 0 || byte == b'-' || byte == b'*' {
        return AlphabetType::Protein;
    }

    // Default to ASCII
    AlphabetType::Ascii
}

/// Checks if alphabet A is more general (can represent more characters) than alphabet B.
fn is_more_general_alphabet(a: AlphabetType, b: AlphabetType) -> bool {
    let alphabet_hierarchy = [
        AlphabetType::Dna2bit,
        AlphabetType::Dna3bit,
        AlphabetType::DnaIupac,
        AlphabetType::Protein,
        AlphabetType::Ascii,
    ];

    let pos_a = alphabet_hierarchy.iter().position(|&x| x == a).unwrap_or(4);
    let pos_b = alphabet_hierarchy.iter().position(|&x| x == b).unwrap_or(4);

    pos_a > pos_b
}

#[cfg(test)]
mod tests {
    use super::{AlphabetGuesser, AlphabetType, guess_alphabet, guess_alphabet_fast};

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

        // Mixed case with protein characters
        assert_eq!(guess_alphabet(b"actgEFIL"), AlphabetType::Protein);

        // Test with non-standard characters that should force ASCII
        assert_eq!(guess_alphabet(b"ACGT123"), AlphabetType::Ascii);
        assert_eq!(guess_alphabet(b"ACGT@#$"), AlphabetType::Ascii);

        // Edge cases with single characters
        assert_eq!(guess_alphabet(b"A"), AlphabetType::Dna2bit);
        assert_eq!(guess_alphabet(b"N"), AlphabetType::Dna3bit);
        assert_eq!(guess_alphabet(b"M"), AlphabetType::DnaIupac);
        assert_eq!(guess_alphabet(b"E"), AlphabetType::Protein);
        assert_eq!(guess_alphabet(b"1"), AlphabetType::Ascii);

        // Test that demonstrates the difference between original and accurate versions
        // The guess_alphabet_fast would return Dna3bit for "ACTGM" but, the better one returns DnaIupac
        assert_eq!(guess_alphabet_fast(b"ACTGM"), AlphabetType::Dna3bit); // Original is wrong
        assert_eq!(guess_alphabet(b"ACTGM"), AlphabetType::DnaIupac); // Accurate is correct
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

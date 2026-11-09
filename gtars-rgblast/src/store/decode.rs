//! Stored record -> BLAST bytes or ASCII, for every gtars storage mode.

use crate::convert;
use crate::error::Error;
use gtars_refget::digest::{AlphabetType, decode_substring_from_bytes, lookup_alphabet};
use gtars_refget::store::StorageMode;
use std::borrow::Cow;

/// Letters `lo..hi` of a stored record in blastna (A0 C1 G2 T3 R4 Y5 M6 K7 W8 S9 B10 D11 H12
/// V13 N14; `U` reads as `T`), followed by one sentinel (15), as `seqsrc_seqdb.cpp` hands the
/// traceback stage a subject (`eBlastEncodingNucleotide`): the window the traceback stage and
/// `format::alignment` read. `out` is cleared first
/// (callers reuse one buffer across windows). Encoded and raw records decode only the window;
/// a zstd record is decompressed whole first.
pub fn to_blastna_window(
    mode: StorageMode,
    alphabet: AlphabetType,
    stored: &[u8],
    len: u32,
    lo: usize,
    hi: usize,
    out: &mut Vec<u8>,
) -> Result<(), Error> {
    debug_assert!(lo <= hi && hi <= len as usize);
    out.clear();
    out.reserve(hi - lo + 1);
    to_blastna_into(mode, alphabet, stored, len, lo, hi, out)?;
    out.push(15u8);
    Ok(())
}

/// Letters `lo..hi` appended to `out` in blastna: packed records straight from their codes
/// (`convert::blastna`), ASCII records through the letter table.
fn to_blastna_into(
    mode: StorageMode,
    alphabet: AlphabetType,
    stored: &[u8],
    len: u32,
    lo: usize,
    hi: usize,
    out: &mut Vec<u8>,
) -> Result<(), Error> {
    match mode {
        StorageMode::Encoded => {
            let a = lookup_alphabet(&alphabet);
            let need = (len as u64 * a.bits_per_symbol as u64).div_ceil(8);
            if stored.len() as u64 != need {
                return Err(Error::sequence(
                    "",
                    format!(
                        "stored record is {} bytes, expected {need} for {len} bp of {alphabet}",
                        stored.len()
                    ),
                ));
            }
            if a.bits_per_symbol == 2 || a.bits_per_symbol == 3 || a.bits_per_symbol == 4 {
                convert::blastna::packed_to_blastna(alphabet, stored, lo, hi, out);
            } else {
                let ascii = decode_substring_from_bytes(stored, lo, hi, a);
                push_blastna(out, &ascii, lo)?;
            }
        }
        StorageMode::Raw => {
            check_raw_len(stored.len(), len)?;
            push_blastna(out, &stored[lo..hi].to_ascii_uppercase(), lo)?;
        }
        StorageMode::Zstd => {
            let ascii = unzstd(stored, len)?;
            push_blastna(out, &ascii[lo..hi].to_ascii_uppercase(), lo)?;
        }
    }
    Ok(())
}

/// ASCII letters -> blastna appended to `out`; `at` is the position of `ascii[0]` for the
/// error message.
fn push_blastna(out: &mut Vec<u8>, ascii: &[u8], at: usize) -> Result<(), Error> {
    for (i, &b) in ascii.iter().enumerate() {
        let v = crate::engine::encoding::ASCII_TO_BLASTNA[b as usize];
        if v == 255 {
            return Err(Error::sequence(
                "",
                format!(
                    "byte {:?} (0x{b:02x}) at position {} is not a nucleotide letter",
                    b as char,
                    at + i
                ),
            ));
        }
        out.push(v);
    }
    Ok(())
}

/// Stored record -> BLAST ncbi2na bytes + ambiguity table; returns the ambiguous and N base
/// counts. `Sequence` errors carry an empty digest; callers add it with `Error::with_digest`.
pub fn to_blast(
    mode: StorageMode,
    alphabet: AlphabetType,
    stored: &[u8],
    len: u32,
    seq: &mut Vec<u8>,
    amb: &mut Vec<u8>,
) -> Result<convert::BaseCounts, Error> {
    let (pieces, src) = to_blast_pieces(mode, alphabet, stored, len)?;
    pieces.whole(&src, seq, amb)
}

/// Stored record -> a piecewise conversion to ncbi2na (`convert::Pieces`) and the bytes it
/// reads: the packed codes, or the ASCII letters of a raw or zstd record (zstd decompressed
/// whole).
pub fn to_blast_pieces(
    mode: StorageMode,
    alphabet: AlphabetType,
    stored: &[u8],
    len: u32,
) -> Result<(convert::Pieces, Cow<'_, [u8]>), Error> {
    let (alphabet, src) = match mode {
        StorageMode::Encoded => (Some(alphabet), Cow::Borrowed(stored)),
        StorageMode::Raw => (None, Cow::Borrowed(stored)),
        StorageMode::Zstd => (None, Cow::Owned(unzstd(stored, len)?)),
    };
    Ok((convert::Pieces::new(alphabet, src.len(), len)?, src))
}

/// Stored record -> uppercase ASCII (verify slow path, reference FASTA).
pub fn to_ascii(
    mode: StorageMode,
    alphabet: AlphabetType,
    stored: &[u8],
    len: u32,
) -> Result<Vec<u8>, Error> {
    match mode {
        StorageMode::Encoded => {
            let a = lookup_alphabet(&alphabet);
            let need = (len as u64 * a.bits_per_symbol as u64).div_ceil(8);
            if stored.len() as u64 != need {
                return Err(Error::sequence(
                    "",
                    format!(
                        "stored record is {} bytes, expected {need} for {len} bp of {alphabet}",
                        stored.len()
                    ),
                ));
            }
            Ok(decode_substring_from_bytes(stored, 0, len as usize, a))
        }
        StorageMode::Raw => {
            check_raw_len(stored.len(), len)?;
            Ok(stored.to_ascii_uppercase())
        }
        StorageMode::Zstd => Ok(unzstd(stored, len)?.to_ascii_uppercase()),
    }
}

fn check_raw_len(got: usize, len: u32) -> Result<(), Error> {
    if got as u64 != len as u64 {
        return Err(Error::sequence(
            "",
            format!("stored record is {got} bytes, expected {len}"),
        ));
    }
    Ok(())
}

/// The whole ASCII record compressed with zstd (gtars `StorageMode::Zstd`).
fn unzstd(stored: &[u8], len: u32) -> Result<Vec<u8>, Error> {
    let ascii = zstd::decode_all(stored)
        .map_err(|e| Error::sequence("", format!("zstd decode failed: {e}")))?;
    check_raw_len(ascii.len(), len)?;
    Ok(ascii)
}

#[cfg(test)]
mod tests {
    use super::*;
    use gtars_refget::digest::encode_sequence;

    fn random_iupac(seed: u64, n: usize, letters: &[u8]) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                letters[((x >> 33) as usize) % letters.len()]
            })
            .collect()
    }

    /// The direct packed -> blastna decoder equals the ASCII route for every alphabet, whole
    /// records and windows at every start and end phase.
    #[test]
    fn packed_windows_match_the_ascii_route() {
        for (alphabet, letters) in [
            (AlphabetType::Dna2bit, &b"ACGT"[..]),
            (AlphabetType::Dna3bit, &b"ACGTNRYU"[..]),
            (AlphabetType::DnaIupac, &b"NACMGRSKTWYUBDHV"[..]),
        ] {
            for n in [0usize, 1, 5, 7, 8, 9, 23, 24, 25, 1000, 1003] {
                let ascii = random_iupac(n as u64 + 11, n, letters);
                let packed = encode_sequence(&ascii, lookup_alphabet(&alphabet)).unwrap();
                let mut want = vec![15u8];
                push_blastna(&mut want, &ascii, 0).unwrap();
                want.push(15);
                let mut out = Vec::new();
                for lo in 0..n.min(13) {
                    for hi in [lo, (lo + 1).min(n), (lo + 7).min(n), (lo + 9).min(n), n] {
                        to_blastna_window(
                            StorageMode::Encoded,
                            alphabet,
                            &packed,
                            n as u32,
                            lo,
                            hi,
                            &mut out,
                        )
                        .unwrap();
                        let mut w = want[1 + lo..1 + hi].to_vec();
                        w.push(15);
                        assert_eq!(out, w, "{alphabet} n={n} {lo}..{hi}");
                    }
                }
            }
        }
    }
}

//! VRS allele normalization.
//!
//! Port of bioutils `normalize()` for fully-justified (EXPAND mode) normalization.
//!
//! The reference sequence is accessed through the [`RefSeq`] trait, which exposes
//! only single-base reads (`base_at`) plus a range copy (`extend_range`). This lets
//! normalization run either over a decoded `&[u8]` (zero-cost) or directly over a
//! 2-bit-encoded buffer that decodes bases on the fly ([`EncodedSeq`]), without the
//! algorithm knowing the difference.

use gtars_refget::digest::alphabet::lookup_alphabet;
use gtars_refget::store::ReadonlyRefgetStore;

/// Random-access, single-base view over a reference sequence.
///
/// Implemented for decoded byte slices (`[u8]`) and for the bit-packed
/// [`EncodedSeq`]. Normalization only ever reads individual bases (during rolling)
/// and copies the small rolled-context ranges, so this minimal surface is enough.
pub trait RefSeq {
    /// Number of bases in the sequence.
    fn len(&self) -> usize;
    /// True if the sequence is empty.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// Decoded base at position `i` (0-based). Caller guarantees `i < len()`.
    fn base_at(&self, i: usize) -> u8;
    /// Append decoded bases for `start..end` onto `out`.
    fn extend_range(&self, start: usize, end: usize, out: &mut Vec<u8>);
}

impl RefSeq for [u8] {
    #[inline]
    fn len(&self) -> usize {
        <[u8]>::len(self)
    }
    #[inline]
    fn base_at(&self, i: usize) -> u8 {
        self[i]
    }
    #[inline]
    fn extend_range(&self, start: usize, end: usize, out: &mut Vec<u8>) {
        out.extend_from_slice(&self[start..end]);
    }
}

/// A bit-packed (e.g. 2-bit DNA) reference view that decodes single bases on the fly.
///
/// Holds a borrowed slice of the encoded bytes plus just enough of the alphabet
/// (`bits_per_symbol` + `decoding_array`) to decode any base by position. Decoding
/// mirrors `gtars_refget`'s `decode_substring_from_bytes` (MSB-first within each byte).
pub struct EncodedSeq<'a> {
    /// The bit-packed encoded bytes.
    pub bytes: &'a [u8],
    /// Number of symbols (bases) represented.
    pub length: usize,
    /// Bits used per symbol (2 for DNA_2BIT).
    pub bits_per_symbol: usize,
    /// Maps an encoded code back to its original symbol byte.
    pub decoding_array: &'a [u8; 256],
}

impl RefSeq for EncodedSeq<'_> {
    #[inline]
    fn len(&self) -> usize {
        self.length
    }
    #[inline]
    fn base_at(&self, i: usize) -> u8 {
        let bit_offset = i * self.bits_per_symbol;
        let mut code = 0u8;
        for j in 0..self.bits_per_symbol {
            let bit_pos = bit_offset + j;
            let byte_index = bit_pos / 8;
            let bit_in_byte = 7 - (bit_pos % 8); // MSB0, matches encoder
            let bit = if byte_index < self.bytes.len() {
                (self.bytes[byte_index] >> bit_in_byte) & 1
            } else {
                0
            };
            code = (code << 1) | bit;
        }
        self.decoding_array[code as usize]
    }
    #[inline]
    fn extend_range(&self, start: usize, end: usize, out: &mut Vec<u8>) {
        out.reserve(end - start);
        for i in start..end {
            out.push(self.base_at(i));
        }
    }
}

/// A reference view that is either already-decoded bytes (Raw-mode store) or a
/// bit-packed [`EncodedSeq`] (Encoded-mode store). Lets a caller hand `normalize`
/// a `RefSeq` without knowing the store's storage mode.
pub enum RefView<'a> {
    /// Already-decoded raw bytes (e.g. a Raw-mode store, or an mmap'd decoded file).
    Decoded(&'a [u8]),
    /// Bit-packed bytes, decoded per base on the fly.
    Encoded(EncodedSeq<'a>),
}

impl RefSeq for RefView<'_> {
    #[inline]
    fn len(&self) -> usize {
        match self {
            RefView::Decoded(s) => s.len(),
            RefView::Encoded(e) => e.len(),
        }
    }
    #[inline]
    fn base_at(&self, i: usize) -> u8 {
        match self {
            RefView::Decoded(s) => s[i],
            RefView::Encoded(e) => e.base_at(i),
        }
    }
    #[inline]
    fn extend_range(&self, start: usize, end: usize, out: &mut Vec<u8>) {
        match self {
            RefView::Decoded(s) => out.extend_from_slice(&s[start..end]),
            RefView::Encoded(e) => e.extend_range(start, end, out),
        }
    }
}

/// Errors building a [`RefView`] from a refget store.
#[derive(Debug, thiserror::Error)]
pub enum RefViewError {
    /// `get_sequence` failed for the digest (not in the store / lookup error).
    #[error("sequence not found in store")]
    NotFound,
    /// The record exists but its bytes are not resident (not loaded).
    #[error("sequence not resident in store")]
    NotResident,
}

/// Build a [`RefView`] over a resident sequence's bytes.
///
/// Whether the bytes are already-decoded (`len == length`) or 2-bit-packed
/// (shorter) is decided by length — robust to storage mode and to the legacy
/// mmap'd decoded-cache variant. The record must be resident (Full); load it
/// first for on-disk stores.
pub fn ref_view_for<'a>(
    store: &'a ReadonlyRefgetStore,
    raw_digest: &str,
) -> Result<RefView<'a>, RefViewError> {
    let rec = store
        .get_sequence(raw_digest)
        .map_err(|_| RefViewError::NotFound)?;
    let meta = rec.metadata();
    let bytes = rec.sequence().ok_or(RefViewError::NotResident)?;
    Ok(if bytes.len() == meta.length {
        RefView::Decoded(bytes)
    } else {
        let alphabet = lookup_alphabet(&meta.alphabet);
        RefView::Encoded(EncodedSeq {
            bytes,
            length: meta.length,
            bits_per_symbol: alphabet.bits_per_symbol,
            decoding_array: alphabet.decoding_array,
        })
    })
}

/// Result of normalizing an allele against a reference sequence.
#[derive(Debug, Clone, PartialEq)]
pub struct NormalizedAllele {
    pub start: u64,
    pub end: u64,
    pub allele: Vec<u8>,
}

/// Trim common prefix from a set of alleles (including the reference).
///
/// Returns `(trimmed_count, trimmed_alleles)`.
fn trim_left(alleles: &[&[u8]]) -> (usize, Vec<Vec<u8>>) {
    if alleles.is_empty() {
        return (0, Vec::new());
    }
    let min_len = alleles.iter().map(|a| a.len()).min().unwrap_or(0);
    let mut trimmed = 0;
    while trimmed < min_len {
        let ch = alleles[0][trimmed];
        if alleles.iter().all(|a| a[trimmed] == ch) {
            trimmed += 1;
        } else {
            break;
        }
    }
    let result = alleles.iter().map(|a| a[trimmed..].to_vec()).collect();
    (trimmed, result)
}

/// Trim common suffix from a set of alleles (including the reference).
///
/// Returns `(trimmed_count, trimmed_alleles)`.
fn trim_right(alleles: &[&[u8]]) -> (usize, Vec<Vec<u8>>) {
    if alleles.is_empty() {
        return (0, Vec::new());
    }
    let min_len = alleles.iter().map(|a| a.len()).min().unwrap_or(0);
    let mut trimmed = 0;
    while trimmed < min_len {
        let ch = alleles[0][alleles[0].len() - 1 - trimmed];
        if alleles.iter().all(|a| a[a.len() - 1 - trimmed] == ch) {
            trimmed += 1;
        } else {
            break;
        }
    }
    let result = alleles
        .iter()
        .map(|a| a[..a.len() - trimmed].to_vec())
        .collect();
    (trimmed, result)
}

/// Roll left: find how far alleles can be circularly shifted left in the reference.
///
/// `sequence` is the full reference, `ref_pos` is the current left boundary (0-based),
/// `bound` is the leftmost allowed position.
fn roll_left<R: RefSeq + ?Sized>(sequence: &R, alleles: &[&[u8]], ref_pos: usize, bound: usize) -> usize {
    // Collect non-empty alleles with their index and length
    let mut non_empty: Vec<(usize, usize)> = Vec::new();
    for (i, a) in alleles.iter().enumerate() {
        if !a.is_empty() {
            non_empty.push((i, a.len()));
        }
    }

    if non_empty.is_empty() || ref_pos <= bound {
        return 0;
    }

    let max_d = ref_pos - bound;
    let mut d = 0;
    while d < max_d {
        let seq_pos = ref_pos - 1 - d;
        let base = sequence.base_at(seq_pos);
        let mismatched = non_empty.iter().any(|&(i, len)| {
            let allele = alleles[i];
            // Circular index from the end of the allele
            let idx = if (d + 1) % len == 0 {
                0
            } else {
                len - ((d + 1) % len)
            };
            allele[idx] != base
        });
        if mismatched {
            break;
        }
        d += 1;
    }
    d
}

/// Roll right: find how far alleles can be circularly shifted right in the reference.
///
/// `sequence` is the full reference, `ref_pos` is the current right boundary (0-based),
/// `bound` is the rightmost allowed position.
fn roll_right<R: RefSeq + ?Sized>(sequence: &R, alleles: &[&[u8]], ref_pos: usize, bound: usize) -> usize {
    let mut non_empty: Vec<(usize, usize)> = Vec::new();
    for (i, a) in alleles.iter().enumerate() {
        if !a.is_empty() {
            non_empty.push((i, a.len()));
        }
    }

    if non_empty.is_empty() || ref_pos >= bound {
        return 0;
    }

    let max_d = bound - ref_pos;
    let mut d = 0;
    while d < max_d {
        let seq_pos = ref_pos + d;
        let base = sequence.base_at(seq_pos);
        let mismatched = non_empty.iter().any(|&(i, len)| {
            let allele = alleles[i];
            allele[d % len] != base
        });
        if mismatched {
            break;
        }
        d += 1;
    }
    d
}

/// Errors that can occur during allele normalization.
#[derive(Debug, Clone, PartialEq)]
pub enum NormalizeError {
    /// Start position exceeds sequence length or overflows usize.
    StartOutOfBounds { start: u64, seq_len: usize },
    /// Reference allele extends past the end of the sequence.
    RefAllelePastEnd { start: usize, ref_len: usize, seq_len: usize },
}

impl std::fmt::Display for NormalizeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NormalizeError::StartOutOfBounds { start, seq_len } => {
                write!(f, "start position {} exceeds sequence length {}", start, seq_len)
            }
            NormalizeError::RefAllelePastEnd { start, ref_len, seq_len } => {
                write!(
                    f,
                    "ref allele (start={}, len={}) extends past sequence length {}",
                    start, ref_len, seq_len
                )
            }
        }
    }
}

impl std::error::Error for NormalizeError {}

/// Normalize an allele against a reference sequence using fully-justified (EXPAND) mode.
///
/// This is the VRS normalization algorithm:
/// 1. Trim common prefix/suffix between ref and alt alleles
/// 2. Expand by rolling left and right through repeat regions
///
/// # Arguments
/// * `sequence` - Full reference chromosome sequence as bytes
/// * `start` - 0-based interbase start position of the variant
/// * `ref_allele` - Reference allele bytes
/// * `alt_allele` - Alternate allele bytes
///
/// # Returns
/// `NormalizedAllele` with adjusted start, end, and alternate allele sequence.
pub fn normalize(
    sequence: &[u8],
    start: u64,
    ref_allele: &[u8],
    alt_allele: &[u8],
) -> Result<NormalizedAllele, NormalizeError> {
    normalize_ref(sequence, start, ref_allele, alt_allele)
}

/// Generic over any [`RefSeq`] reference view (decoded `&[u8]` or [`EncodedSeq`]).
///
/// This is the real implementation; [`normalize`] is a `&[u8]` convenience wrapper.
/// Producing byte-identical results for the decoded and encoded views is what lets
/// us keep the reference 2-bit-encoded in memory without changing any VRS digest.
pub fn normalize_ref<R: RefSeq + ?Sized>(
    sequence: &R,
    start: u64,
    ref_allele: &[u8],
    alt_allele: &[u8],
) -> Result<NormalizedAllele, NormalizeError> {
    let seq_len = sequence.len();
    let s_usize = usize::try_from(start).map_err(|_| NormalizeError::StartOutOfBounds {
        start,
        seq_len,
    })?;
    let mut s = s_usize;
    let mut e = s.checked_add(ref_allele.len()).ok_or(NormalizeError::StartOutOfBounds {
        start,
        seq_len,
    })?;
    if e > seq_len {
        return Err(NormalizeError::RefAllelePastEnd {
            start: s,
            ref_len: ref_allele.len(),
            seq_len,
        });
    }

    // Step 1: Trim common prefix
    let alleles_for_trim: Vec<&[u8]> = vec![ref_allele, alt_allele];
    let (left_trimmed, trimmed_alleles) = trim_left(&alleles_for_trim);
    s += left_trimmed;

    // Step 2: Trim common suffix
    let trimmed_refs: Vec<&[u8]> = trimmed_alleles.iter().map(|a| a.as_slice()).collect();
    let (right_trimmed, trimmed_alleles2) = trim_right(&trimmed_refs);
    e -= right_trimmed;

    let ref_trimmed = &trimmed_alleles2[0];
    let alt_trimmed = &trimmed_alleles2[1];

    // Step 3: Expand (fully-justified normalization)
    // Roll left from start, roll right from end
    let bound_left = 0;
    let bound_right = seq_len;

    let alleles_for_roll: Vec<&[u8]> = vec![ref_trimmed.as_slice(), alt_trimmed.as_slice()];
    let left_roll = roll_left(sequence, &alleles_for_roll, s, bound_left);
    let right_roll = roll_right(sequence, &alleles_for_roll, e, bound_right);

    let new_start = s - left_roll;
    let new_end = e + right_roll;

    // Rebuild the alt allele with the expanded context
    // Left context: sequence[new_start..s], then alt_trimmed, then right context: sequence[e..new_end]
    let mut new_alt = Vec::with_capacity((s - new_start) + alt_trimmed.len() + (new_end - e));
    sequence.extend_range(new_start, s, &mut new_alt);
    new_alt.extend_from_slice(alt_trimmed);
    sequence.extend_range(e, new_end, &mut new_alt);

    Ok(NormalizedAllele {
        start: new_start as u64,
        end: new_end as u64,
        allele: new_alt,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trim_left_basic() {
        let alleles: Vec<&[u8]> = vec![b"ATCG", b"ATGG"];
        let (trimmed, result) = trim_left(&alleles);
        assert_eq!(trimmed, 2);
        assert_eq!(result, vec![b"CG".to_vec(), b"GG".to_vec()]);
    }

    #[test]
    fn test_trim_right_basic() {
        let alleles: Vec<&[u8]> = vec![b"ATCG", b"AGCG"];
        let (trimmed, result) = trim_right(&alleles);
        assert_eq!(trimmed, 2);
        assert_eq!(result, vec![b"AT".to_vec(), b"AG".to_vec()]);
    }

    #[test]
    fn test_trim_empty_allele() {
        // Insertion: ref is empty after trimming
        let alleles: Vec<&[u8]> = vec![b"A", b"AT"];
        let (left, result) = trim_left(&alleles);
        assert_eq!(left, 1);
        assert_eq!(result, vec![b"".to_vec(), b"T".to_vec()]);
    }

    #[test]
    fn test_normalize_snv() {
        // Simple SNV: no trimming or rolling needed
        let seq = b"ACGTACGT";
        let result = normalize(seq, 2, b"G", b"T").unwrap();
        assert_eq!(result.start, 2);
        assert_eq!(result.end, 3);
        assert_eq!(result.allele, b"T");
    }

    #[test]
    fn test_normalize_insertion_in_repeat() {
        // Insert "A" into a run of 4 A's: TAAAAG
        // VCF-style: pos=2 (1-based), ref=A, alt=AA → interbase start=1, ref=A, alt=AA
        let seq = b"TAAAAG";
        //          012345
        let result = normalize(seq, 1, b"A", b"AA").unwrap();
        // After trim: start=2, end=2, ref="", alt="A"
        // Roll left from 2: seq[1]='A' matches → 1 step, seq[0]='T' doesn't → stop
        // left_roll=1, new_start=1
        // Roll right from 2: seq[2]='A', seq[3]='A', seq[4]='A' all match, seq[5]='G' doesn't
        // right_roll=3, new_end=5
        assert_eq!(result.start, 1);
        assert_eq!(result.end, 5);
        // new alt = seq[1..2] + "A" + seq[2..5] = "A" + "A" + "AAA" = "AAAAA"
        assert_eq!(result.allele, b"AAAAA");
    }

    #[test]
    fn test_normalize_deletion() {
        // Delete one A from a run of 4 A's: TAAAAG
        let seq = b"TAAAAG";
        //          012345
        let result = normalize(seq, 1, b"AA", b"A").unwrap();
        // After trim: ref="A", alt="" (one A trimmed from left)
        // start=2, end=3
        // Roll left: seq[1]='A' matches → 1 step, seq[0]='T' doesn't → stop
        // left_roll=1, new_start=1
        // Roll right from 3: seq[3]='A', seq[4]='A' match, seq[5]='G' doesn't
        // right_roll=2, new_end=5
        assert_eq!(result.start, 1);
        assert_eq!(result.end, 5);
        // new alt = seq[1..2] + "" + seq[3..5] = "A" + "AA" = "AAA"
        assert_eq!(result.allele, b"AAA");
    }

    #[test]
    fn test_normalize_identity() {
        // ref == alt: no change needed
        let seq = b"ACGTACGT";
        let result = normalize(seq, 2, b"GT", b"GT").unwrap();
        // Fully trimmed away
        assert_eq!(result.start, 4);
        assert_eq!(result.end, 4);
        assert_eq!(result.allele, b"");
    }

    #[test]
    fn test_normalize_out_of_bounds() {
        let seq = b"ACGT";
        let result = normalize(seq, 10, b"G", b"T");
        assert!(result.is_err());
    }
}

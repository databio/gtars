//! VRS allele normalization.
//!
//! Port of bioutils `normalize()` for fully-justified (EXPAND mode) normalization.
//! Operates on `&[u8]` byte slices for zero-allocation sequence access.

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
fn roll_left(sequence: &[u8], alleles: &[&[u8]], ref_pos: usize, bound: usize) -> usize {
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
        let mismatched = non_empty.iter().any(|&(i, len)| {
            let allele = alleles[i];
            // Circular index from the end of the allele
            let idx = if (d + 1) % len == 0 {
                0
            } else {
                len - ((d + 1) % len)
            };
            allele[idx] != sequence[seq_pos]
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
fn roll_right(sequence: &[u8], alleles: &[&[u8]], ref_pos: usize, bound: usize) -> usize {
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
        let mismatched = non_empty.iter().any(|&(i, len)| {
            let allele = alleles[i];
            allele[d % len] != sequence[seq_pos]
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
    let s_usize = usize::try_from(start).map_err(|_| NormalizeError::StartOutOfBounds {
        start,
        seq_len: sequence.len(),
    })?;
    let mut s = s_usize;
    let mut e = s.checked_add(ref_allele.len()).ok_or(NormalizeError::StartOutOfBounds {
        start,
        seq_len: sequence.len(),
    })?;
    if e > sequence.len() {
        return Err(NormalizeError::RefAllelePastEnd {
            start: s,
            ref_len: ref_allele.len(),
            seq_len: sequence.len(),
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
    let bound_right = sequence.len();

    let alleles_for_roll: Vec<&[u8]> = vec![ref_trimmed.as_slice(), alt_trimmed.as_slice()];
    let left_roll = roll_left(sequence, &alleles_for_roll, s, bound_left);
    let right_roll = roll_right(sequence, &alleles_for_roll, e, bound_right);

    let new_start = s - left_roll;
    let new_end = e + right_roll;

    // Rebuild the alt allele with the expanded context
    // Left context: sequence[new_start..s], then alt_trimmed, then right context: sequence[e..new_end]
    let mut new_alt = Vec::with_capacity((s - new_start) + alt_trimmed.len() + (new_end - e));
    new_alt.extend_from_slice(&sequence[new_start..s]);
    new_alt.extend_from_slice(alt_trimmed);
    new_alt.extend_from_slice(&sequence[e..new_end]);

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

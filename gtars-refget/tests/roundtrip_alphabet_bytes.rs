//! Byte-level round-trip property test across every alphabet the ingest
//! pipeline can choose.
//!
//! The rule under test: **decoded output must equal the (uppercased)
//! original input**, for every byte value, in every alphabet context the
//! guesser can produce, through every decode path (whole-string, partial
//! from a full buffer, partial from a byte-windowed buffer, and the
//! streaming decoder at every window).
//!
//! Earlier tests (`encoder.rs`, `streaming_decoder.rs`) only ever compared
//! one decoder's output to *another* decoder's output -- e.g. the streaming
//! decoder against `decode_substring_from_bytes`. Both read the same
//! `decoding_array`, so if that table is wrong, both decoders are wrong in
//! the same way and agree with each other. Comparing decoders to each other
//! can never catch a bad table entry; only comparing against the original
//! input can. That is what every check in this file does.

use gtars_refget::digest::{
    Alphabet, AlphabetGuesser, SequenceEncoder, StreamingDecoder, byte_range_for_bases,
    decode_string_from_bytes, decode_substring_from_bytes, decode_substring_from_bytes_at_offset,
    encode_sequence, lookup_alphabet,
};
use std::io::{Cursor, Read};

/// What the ingest pipeline turns a raw byte into before it digests, guesses
/// or encodes (`fasta.rs` uppercases every sequence line).
fn expected(b: u8) -> u8 {
    b.to_ascii_uppercase()
}

/// Run the same alphabet-guessing pipeline the store import path uses:
/// `AlphabetGuesser::new()` + `update` + `guess()`, then `lookup_alphabet`.
fn pipeline_alphabet(seq: &[u8]) -> &'static Alphabet {
    let mut guesser = AlphabetGuesser::new();
    guesser.update(seq);
    lookup_alphabet(&guesser.guess())
}

/// Encode with both encoders the real pipeline uses: the streaming
/// `SequenceEncoder` (import path) and the bulk `encode_sequence`
/// (`set_encoding_mode` path). Assert they agree (this also catches drift
/// between the two encoders) and return the shared bytes.
fn encode_both(seq: &[u8], alphabet: &'static Alphabet) -> Vec<u8> {
    let mut enc = SequenceEncoder::new(alphabet.alphabet_type, seq.len());
    enc.update(seq);
    let streaming = enc.finalize();
    let bulk = encode_sequence(seq, alphabet);
    assert_eq!(
        streaming, bulk,
        "SequenceEncoder and encode_sequence disagree for {:?} under {:?}",
        String::from_utf8_lossy(seq),
        alphabet.alphabet_type
    );
    streaming
}

/// Decode `[start, end)` via `StreamingDecoder`, using the same byte-window
/// math `RefgetStore::stream_sequence` uses to compute the byte range and
/// leading-skip bits for a base range.
fn stream_decode(encoded: &[u8], alphabet: &'static Alphabet, start: usize, end: usize) -> Vec<u8> {
    let bps = alphabet.bits_per_symbol;
    let start_bit = start * bps;
    let end_bit = end * bps;
    let byte_start = start_bit / 8;
    let byte_end = end_bit.div_ceil(8);
    let leading_skip = (start_bit - byte_start * 8) as u8;
    let slice = &encoded[byte_start..byte_end.min(encoded.len())];

    let mut decoder = StreamingDecoder::new(
        Cursor::new(slice.to_vec()),
        alphabet,
        leading_skip,
        (end - start) as u64,
    );
    let mut out = Vec::new();
    decoder
        .read_to_end(&mut out)
        .expect("streaming decode failed");
    out
}

/// Round-trip `seq` (already uppercased) through every decode path listed in
/// the plan (a-f), collecting failures instead of panicking so one run
/// reports every broken byte/alphabet combination.
fn check_roundtrip(label: &str, seq: &[u8], failures: &mut Vec<String>) {
    let alphabet = pipeline_alphabet(seq);
    let encoded = encode_both(seq, alphabet);
    let len = seq.len();

    let mut record = |check: &str, expected: &[u8], actual: &[u8]| {
        failures.push(format!(
            "{label}: [{check}] alphabet={:?} seq={:?} raw_bytes={:?} -- expected {:?} got {:?}",
            alphabet.alphabet_type,
            String::from_utf8_lossy(seq),
            seq,
            String::from_utf8_lossy(expected),
            String::from_utf8_lossy(actual),
        ));
    };

    // a. decode_string_from_bytes must equal the whole input.
    let got = decode_string_from_bytes(&encoded, len, alphabet);
    if got != seq {
        record("a:decode_string_from_bytes", seq, &got);
    }

    // b. decode_substring_from_bytes(i, len) must equal seq[i..], every i.
    for i in 0..len {
        let got = decode_substring_from_bytes(&encoded, i, len, alphabet);
        if got != &seq[i..] {
            record("b:decode_substring_from_bytes", &seq[i..], &got);
        }
    }

    // c. decode_substring_from_bytes_at_offset (the partial-read path used
    // by `get_substring_from_disk`) must equal seq[i..i+1], every i.
    for i in 0..len {
        let (bs, be) = byte_range_for_bases(i, i + 1, alphabet.bits_per_symbol);
        let partial = &encoded[bs..be.min(encoded.len())];
        let got = decode_substring_from_bytes_at_offset(partial, bs, i, i + 1, alphabet);
        if got != &seq[i..i + 1] {
            record(
                "c:decode_substring_from_bytes_at_offset",
                &seq[i..i + 1],
                &got,
            );
        }
    }

    // d. Streaming decode of the whole sequence.
    let got = stream_decode(&encoded, alphabet, 0, len);
    if got != seq {
        record("d:stream_decode(full)", seq, &got);
    }

    // e. Streaming decode of every single-base window.
    for i in 0..len {
        let got = stream_decode(&encoded, alphabet, i, i + 1);
        if got != &seq[i..i + 1] {
            record("e:stream_decode(single)", &seq[i..i + 1], &got);
        }
    }

    // f. Streaming decode of every tail window.
    for i in 0..len {
        let got = stream_decode(&encoded, alphabet, i, len);
        if got != &seq[i..] {
            record("f:stream_decode(tail)", &seq[i..], &got);
        }
    }
}

/// Build the 9-position sweep of `base` with `extra` inserted at every
/// possible position (`0..=base.len()`), so a single extra byte lands at
/// every bit offset inside a byte for the 2-, 3-, 4- and 5-bit widths.
fn positional_variants(base: &[u8], extra: u8) -> Vec<Vec<u8>> {
    (0..=base.len())
        .map(|p| {
            let mut v = Vec::with_capacity(base.len() + 1);
            v.extend_from_slice(&base[..p]);
            v.push(extra);
            v.extend_from_slice(&base[p..]);
            v
        })
        .collect()
}

/// The test that would have caught every known bug: every byte value,
/// through every context the guesser can send it, through every decode path.
#[test]
fn every_byte_roundtrips_in_every_context() {
    let mut failures = Vec::new();

    for b in 0u8..=255 {
        let u = expected(b);

        // Bare: the byte alone, repeated, at lengths that hit the odd tail
        // for every bit width (2, 3, 4, 5, 8 bits/symbol).
        for &n in &[1usize, 2, 3, 9] {
            let seq = vec![u; n];
            check_roundtrip(&format!("bare x{n} byte=0x{b:02x}"), &seq, &mut failures);
        }

        // Nucleotide context.
        for seq in positional_variants(b"ACGTACGT", u) {
            check_roundtrip(&format!("nucleotide byte=0x{b:02x}"), &seq, &mut failures);
        }

        // Ambiguous-nucleotide context, pushed to at least Dna3bit.
        for seq in positional_variants(b"NRYACGT", u) {
            check_roundtrip(&format!("ambiguous(NRY) byte=0x{b:02x}"), &seq, &mut failures);
        }

        // Ambiguous-nucleotide context, pushed to at least DnaIupac; the
        // base string itself covers every IUPAC letter.
        for seq in positional_variants(b"ACGTRYSWKMBDHVN", u) {
            check_roundtrip(
                &format!("ambiguous(IUPAC) byte=0x{b:02x}"),
                &seq,
                &mut failures,
            );
        }

        // Protein context: E, F, I, L, P, Q are not IUPAC letters, so this
        // forces at least Protein.
        for seq in positional_variants(b"MEFILPQ", u) {
            check_roundtrip(&format!("protein byte=0x{b:02x}"), &seq, &mut failures);
        }

        // RNA context (U present).
        for seq in positional_variants(b"ACGUACGU", u) {
            check_roundtrip(&format!("rna byte=0x{b:02x}"), &seq, &mut failures);
        }
    }

    assert!(
        failures.is_empty(),
        "{} round-trip failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// The alphabet-aware guard: whatever alphabet the guesser picks for a given
/// input must actually be able to represent every byte in that input. This
/// is phrased independently of how alphabet selection works, so it does not
/// need updating when the selection rules change -- it just names the
/// broken table entry directly.
#[test]
fn guesser_never_picks_an_alphabet_that_cannot_hold_the_input() {
    let mut failures = Vec::new();

    let mut check = |seq: &[u8]| {
        let alphabet = pipeline_alphabet(seq);
        for &x in seq {
            let code = alphabet.encoding_array[x as usize];
            let back = alphabet.decoding_array[code as usize];
            if back != x {
                failures.push(format!(
                    "{:?} cannot represent '{}' (0x{:02x}): encodes to {} which decodes to '{}' (seq={:?})",
                    alphabet.alphabet_type,
                    x as char,
                    x,
                    code,
                    back as char,
                    String::from_utf8_lossy(seq),
                ));
            }
        }
    };

    for b in 0u8..=255 {
        let u = expected(b);

        for &n in &[1usize, 2, 3, 9] {
            check(&vec![u; n]);
        }
        for seq in positional_variants(b"MEFILPQ", u) {
            check(&seq);
        }
    }

    assert!(
        failures.is_empty(),
        "{} alphabet-selection failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// Pins the "decoding always yields uppercase" contract: for every non-ASCII
/// alphabet, encoding a lowercase letter that is a real member of that
/// alphabet must decode back to the uppercase form. This matters for callers
/// that skip the pipeline's own uppercasing step (e.g. `set_encoding_mode`
/// applied to hand-built records). ASCII is excluded: it passes bytes
/// through unchanged by design, so it must NOT uppercase.
#[test]
fn lowercase_input_to_non_ascii_alphabets_decodes_uppercase() {
    use gtars_refget::digest::AlphabetType;

    let mut failures = Vec::new();
    let alphabets = [
        lookup_alphabet(&AlphabetType::Dna2bit),
        lookup_alphabet(&AlphabetType::Dna3bit),
        lookup_alphabet(&AlphabetType::DnaIupac),
        lookup_alphabet(&AlphabetType::Protein),
    ];

    for alphabet in alphabets {
        for upper in b'A'..=b'Z' {
            // Only test letters that are real members of this alphabet.
            let code = alphabet.encoding_array[upper as usize];
            if alphabet.decoding_array[code as usize] != upper {
                continue;
            }
            let lower = upper.to_ascii_lowercase();
            let encoded = encode_sequence([lower], alphabet);
            let decoded = decode_string_from_bytes(&encoded, 1, alphabet);
            if decoded != vec![upper] {
                failures.push(format!(
                    "{:?}: lowercase '{}' decoded to {:?}, expected '{}'",
                    alphabet.alphabet_type, lower as char, decoded, upper as char
                ));
            }
        }
    }

    assert!(
        failures.is_empty(),
        "{} lowercase-decoding failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

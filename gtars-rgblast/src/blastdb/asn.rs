//! BER encoding of Blast-def-line-set, as NCBI's serializer writes it (indefinite lengths
//! for every constructed value).
//!
//! ```text
//! Blast-def-line ::= SEQUENCE {
//!   title [0] VisibleString OPTIONAL, seqid [1] SEQUENCE OF Seq-id,
//!   taxid [2] INTEGER OPTIONAL, ... }
//! Seq-id ::= CHOICE { local [0] Object-id, ..., general [10] Dbtag, ... }
//! Object-id ::= CHOICE { id [0] INTEGER, str [1] VisibleString }
//! Dbtag ::= SEQUENCE { db [0] VisibleString, tag [1] Object-id }
//! ```

use crate::plan::{Defline, SeqId, numeric_tag};

const END: [u8; 2] = [0, 0];

fn len(out: &mut Vec<u8>, n: usize) {
    if n < 0x80 {
        out.push(n as u8);
    } else {
        let b = (n as u64).to_be_bytes();
        let skip = b.iter().take_while(|&&x| x == 0).count();
        out.push(0x80 | (8 - skip) as u8);
        out.extend_from_slice(&b[skip..]);
    }
}

fn vstr(out: &mut Vec<u8>, s: &str) {
    out.push(0x1A);
    len(out, s.len());
    out.extend_from_slice(s.as_bytes());
}

fn int(out: &mut Vec<u8>, v: i64) {
    let b = v.to_be_bytes();
    let mut i = 0;
    while i < 7 && ((b[i] == 0 && b[i + 1] & 0x80 == 0) || (b[i] == 0xFF && b[i + 1] & 0x80 != 0)) {
        i += 1;
    }
    out.push(0x02);
    len(out, 8 - i);
    out.extend_from_slice(&b[i..]);
}

/// Open an indefinite-length constructed value with the given identifier byte.
fn open(out: &mut Vec<u8>, tag: u8) {
    out.push(tag);
    out.push(0x80);
}

fn objid_str(out: &mut Vec<u8>, s: &str) {
    open(out, 0xA1);
    vstr(out, s);
    out.extend_from_slice(&END);
}

/// An Object-id: makeblastdb stores a canonical all-digit string (see [`numeric_tag`]) as
/// `id` INTEGER, for local ids and general-id tags alike (checked against makeblastdb 2.17:
/// `lcl|123456` and `gnl|db|22` are integers; `lcl|0099`, `gnl|db|0` and `gnl|db|007` strings).
fn objid(out: &mut Vec<u8>, s: &str) {
    match numeric_tag(s) {
        Some(v) => {
            open(out, 0xA0);
            int(out, v);
            out.extend_from_slice(&END);
        }
        None => objid_str(out, s),
    }
}

fn seqid(out: &mut Vec<u8>, id: &SeqId) {
    match id {
        SeqId::Local(s) => {
            open(out, 0xA0);
            objid(out, s);
            out.extend_from_slice(&END);
        }
        SeqId::General { db, tag } => {
            open(out, 0xAA);
            open(out, 0x30);
            open(out, 0xA0);
            vstr(out, db);
            out.extend_from_slice(&END);
            open(out, 0xA1);
            objid(out, tag);
            out.extend_from_slice(&END);
            out.extend_from_slice(&END);
            out.extend_from_slice(&END);
        }
    }
}

/// Encode a Blast-def-line-set into `out` (cleared first).
pub fn encode(deflines: &[Defline], out: &mut Vec<u8>) {
    out.clear();
    open(out, 0x30);
    for d in deflines {
        open(out, 0x30);
        open(out, 0xA0);
        vstr(out, &d.title);
        out.extend_from_slice(&END);
        open(out, 0xA1);
        open(out, 0x30);
        for id in &d.ids {
            seqid(out, id);
        }
        out.extend_from_slice(&END);
        out.extend_from_slice(&END);
        if let Some(t) = d.taxid {
            open(out, 0xA2);
            int(out, t as i64);
            out.extend_from_slice(&END);
        }
        out.extend_from_slice(&END);
    }
    out.extend_from_slice(&END);
}

/// Bytes of a BER length field for `n`.
fn len_len(n: usize) -> usize {
    if n < 0x80 {
        1
    } else {
        1 + 8 - (n as u64).leading_zeros() as usize / 8
    }
}

fn vstr_len(s: &str) -> usize {
    1 + len_len(s.len()) + s.len()
}

fn int_len(v: i64) -> usize {
    let b = v.to_be_bytes();
    let mut i = 0;
    while i < 7 && ((b[i] == 0 && b[i + 1] & 0x80 == 0) || (b[i] == 0xFF && b[i + 1] & 0x80 != 0)) {
        i += 1;
    }
    1 + len_len(8 - i) + 8 - i
}

/// `open` + content + `END`.
const WRAP: usize = 4;

fn objid_len(s: &str) -> usize {
    WRAP + match numeric_tag(s) {
        Some(v) => int_len(v),
        None => vstr_len(s),
    }
}

fn seqid_len(id: &SeqId) -> usize {
    match id {
        SeqId::Local(s) => WRAP + objid_len(s),
        SeqId::General { db, tag } => WRAP + WRAP + (WRAP + vstr_len(db)) + (WRAP + objid_len(tag)),
    }
}

/// The length of what [`encode`] writes for `deflines`, without writing it.
pub fn encoded_len(deflines: &[Defline]) -> u32 {
    let mut n = WRAP;
    for d in deflines {
        n += WRAP;
        n += WRAP + vstr_len(&d.title);
        n += WRAP + WRAP + d.ids.iter().map(seqid_len).sum::<usize>();
        if let Some(t) = d.taxid {
            n += WRAP + int_len(t as i64);
        }
    }
    n as u32
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check(d: &[Defline]) {
        let mut out = Vec::new();
        encode(d, &mut out);
        assert_eq!(encoded_len(d) as usize, out.len(), "{d:?}");
    }

    #[test]
    fn encoded_len_equals_encode_len() {
        let g = |db: &str, tag: &str| SeqId::General {
            db: db.into(),
            tag: tag.into(),
        };
        check(&[]);
        for title_len in [0, 1, 127, 128, 255, 256, 65_535, 65_536, 70_000] {
            for taxid in [
                None,
                Some(0),
                Some(127),
                Some(128),
                Some(-1),
                Some(-129),
                Some(9606),
                Some(i32::MAX),
                Some(i32::MIN),
            ] {
                check(&[Defline {
                    ids: vec![SeqId::Local("SQ.abc".into()), g("db", "chr1")],
                    title: "t".repeat(title_len),
                    taxid,
                }]);
            }
        }
        let tags = [
            "0",
            "1",
            "22",
            "00",
            "007",
            "2147483647",
            "2147483648",
            "+1",
            "-1",
            "x",
        ];
        let ids: Vec<SeqId> = tags
            .iter()
            .flat_map(|t| [SeqId::Local(t.to_string()), g(&"d".repeat(200), t)])
            .collect();
        check(&[
            Defline {
                ids: ids.clone(),
                title: String::new(),
                taxid: None,
            },
            Defline {
                ids: vec![g("a", &"n".repeat(300))],
                title: "two".into(),
                taxid: Some(10090),
            },
        ]);
        // a set whose total crosses every length-field size
        let many: Vec<Defline> = (0..2_000)
            .map(|i| Defline {
                ids: vec![g("db", &format!("name{i}"))],
                title: format!("title {i}"),
                taxid: Some(i),
            })
            .collect();
        check(&many);
    }
}

//! Deflines and sequence ids, and the lookup keys BLAST derives from them.

use serde::{Deserialize, Serialize};

/// A BLAST sequence id. Only these two kinds are ever written, so the first id of
/// every defline (what blastn's `sallseqid` prints) is under our control.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SeqId {
    /// `lcl|<string>`
    Local(String),
    /// `gnl|<db>|<tag>`; an all-digit canonical tag is stored as an integer (see [`numeric_tag`]).
    General { db: String, tag: String },
}

/// One Blast-def-line.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Defline {
    pub ids: Vec<SeqId>,
    pub title: String,
    /// None omits the `[2] taxid` field (only the `SQ.` defline of a sequence whose
    /// collections disagree on the taxon); makeblastdb itself always writes one (0 = unknown).
    pub taxid: Option<i32>,
}

/// One id string that `blastdbcmd -entry` accepts for an OID.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct LookupForm {
    pub id: String,
    /// False for a `db:tag` form whose db starts with '-': it resolves inside BLAST
    /// but looks like a flag on a command line, so it is left out of `<db>.ids.tsv`.
    pub cli_safe: bool,
}

impl SeqId {
    /// The acc2oid key(s) WriteDB stores for this id (`GetSeqIdString(true)`), plus an
    /// upper-case copy when that differs. The flag marks keys that also go to the
    /// OID-to-seqids file (`.nos`): the original spelling only.
    pub fn lmdb_keys(&self, out: &mut Vec<(String, bool)>) {
        let k = self.lmdb_key();
        let u = k.to_uppercase();
        if u != k {
            out.push((u, false));
        }
        out.push((k, true));
    }

    /// The lengths of the keys [`lmdb_keys`](Self::lmdb_keys) gives, in the same order and
    /// with the same `.nos` flags, without building them.
    pub fn lmdb_key_lens(&self, mut f: impl FnMut(usize, bool)) {
        let (a, b) = match self {
            SeqId::Local(s) => (s.as_str(), None),
            SeqId::General { db, tag } => (db.as_str(), Some(tag.as_str())),
        };
        let (mut len, mut upper, mut differs) = (0, 0, false);
        for part in [Some(a), b.map(|_| ":"), b].into_iter().flatten() {
            len += part.len();
            if part.is_ascii() {
                upper += part.len();
                differs |= part.bytes().any(|c| c.is_ascii_lowercase());
            } else {
                for c in part.chars() {
                    let mut u = c.to_uppercase();
                    if u.len() != 1 || u.clone().next() != Some(c) {
                        differs = true;
                    }
                    upper += u.by_ref().map(char::len_utf8).sum::<usize>();
                }
            }
        }
        if differs {
            f(upper, false);
        }
        f(len, true);
    }

    /// Append the original-spelling acc2oid key ([`lmdb_key`](Self::lmdb_key)) to `out`.
    pub fn write_lmdb_key(&self, out: &mut Vec<u8>) {
        match self {
            SeqId::Local(s) => out.extend_from_slice(s.as_bytes()),
            SeqId::General { db, tag } => {
                out.extend_from_slice(db.as_bytes());
                out.push(b':');
                out.extend_from_slice(tag.as_bytes());
            }
        }
    }

    /// The original-spelling acc2oid key: the local string, or `db:tag`.
    pub fn lmdb_key(&self) -> String {
        match self {
            SeqId::Local(s) => s.clone(),
            SeqId::General { db, tag } => format!("{db}:{tag}"),
        }
    }

    /// FASTA spelling, for a FASTA that `makeblastdb -parse_seqids` turns into the same ids.
    pub fn fasta(&self) -> String {
        match self {
            SeqId::Local(s) => format!("lcl|{s}"),
            SeqId::General { db, tag } => format!("gnl|{db}|{tag}"),
        }
    }

    /// Every lookup string blastdbcmd accepts for this id.
    pub fn lookup_forms(&self, out: &mut Vec<LookupForm>) {
        match self {
            SeqId::Local(s) => out.push(LookupForm {
                id: s.clone(),
                cli_safe: true,
            }),
            SeqId::General { db, tag } => {
                out.push(LookupForm {
                    id: format!("gnl|{db}|{tag}"),
                    cli_safe: true,
                });
                out.push(LookupForm {
                    id: format!("{db}:{tag}"),
                    cli_safe: !db.starts_with('-'),
                });
            }
        }
    }

    /// The key makeblastdb's duplicate check compares (label, upper-cased for general ids).
    pub(crate) fn dup_key(&self) -> String {
        match self {
            SeqId::Local(s) => format!("lcl|{s}"),
            SeqId::General { .. } => self.lmdb_key().to_uppercase(),
        }
    }
}

/// A canonical positive decimal that fits an Int4, which NCBI's FASTA id parser turns into
/// an integer Object-id (for local ids and general-id tags alike): digits only, no leading
/// zero, from 1 to `i32::MAX`. Checked against makeblastdb 2.17: `1`, `22`, `2147483647` are
/// integers; `0`, `00`, `007`, `2147483648`, `+1`, `-1` stay strings.
pub fn numeric_tag(s: &str) -> Option<i64> {
    if s.is_empty() || s.len() > 10 || !s.bytes().all(|b| b.is_ascii_digit()) || s.starts_with('0')
    {
        return None;
    }
    let v: i64 = s.parse().ok()?;
    (v <= i32::MAX as i64).then_some(v)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numeric_tag_edges() {
        // (tag, stored as integer?, LMDB key for gnl|db|tag)
        let table: [(&str, Option<i64>); 11] = [
            ("0", None),
            ("1", Some(1)),
            ("22", Some(22)),
            ("00", None),
            ("007", None),
            ("0012", None),
            ("2147483647", Some(2_147_483_647)),
            ("2147483648", None),
            ("4294967296", None),
            ("+1", None),
            ("-1", None),
        ];
        for (tag, want) in table {
            assert_eq!(numeric_tag(tag), want, "{tag}");
            let id = SeqId::General {
                db: "asm".into(),
                tag: tag.into(),
            };
            assert_eq!(id.lmdb_key(), format!("asm:{tag}"));
        }
        assert_eq!(numeric_tag(""), None);
        assert_eq!(numeric_tag("12345678901"), None);
    }

    #[test]
    fn lmdb_key_lens_match_lmdb_keys() {
        let g = |db: &str, tag: &str| SeqId::General {
            db: db.into(),
            tag: tag.into(),
        };
        for id in [
            SeqId::Local("SQ.abC".into()),
            SeqId::Local("SQ.ABC".into()),
            SeqId::Local("123".into()),
            g("db", "chr1"),
            g("DB", "CHR1"),
            g("-mUDCQ", "dash36"),
            g("ß", "ﬁx"),
            g("É", "Ω"),
            SeqId::Local("straße".into()),
        ] {
            let mut keys = Vec::new();
            id.lmdb_keys(&mut keys);
            let mut lens = Vec::new();
            id.lmdb_key_lens(|l, save| lens.push((l, save)));
            let want: Vec<(usize, bool)> = keys.iter().map(|(k, s)| (k.len(), *s)).collect();
            assert_eq!(lens, want, "{id:?}");
            let mut b = Vec::new();
            id.write_lmdb_key(&mut b);
            assert_eq!(b, id.lmdb_key().into_bytes());
        }
    }

    #[test]
    fn lmdb_keys_uppercase_copy() {
        let mut k = Vec::new();
        SeqId::Local("SQ.abC".into()).lmdb_keys(&mut k);
        assert_eq!(
            k,
            vec![("SQ.ABC".to_string(), false), ("SQ.abC".to_string(), true)]
        );
        k.clear();
        SeqId::General {
            db: "ASM".into(),
            tag: "CHR1".into(),
        }
        .lmdb_keys(&mut k);
        assert_eq!(k, vec![("ASM:CHR1".to_string(), true)]);
    }

    #[test]
    fn fasta_spelling() {
        assert_eq!(SeqId::Local("SQ.x".into()).fasta(), "lcl|SQ.x");
        assert_eq!(
            SeqId::General {
                db: "hg38".into(),
                tag: "chr1".into()
            }
            .fasta(),
            "gnl|hg38|chr1"
        );
        let mut f = Vec::new();
        SeqId::General {
            db: "-abc".into(),
            tag: "chr1".into(),
        }
        .lookup_forms(&mut f);
        assert_eq!(
            f[0],
            LookupForm {
                id: "gnl|-abc|chr1".into(),
                cli_safe: true
            }
        );
        assert_eq!(
            f[1],
            LookupForm {
                id: "-abc:chr1".into(),
                cli_safe: false
            }
        );
    }
}

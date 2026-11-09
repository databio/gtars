//! Id spelling rules: sanitizing, repeated ids, uniqueness, and `db_tag`.

use super::defline::SeqId;
use std::collections::{HashMap, HashSet};

/// Replace the characters BLAST ids can't hold with `_`: space `>` `[` `]` `|` `"`,
/// `;` (blastn joins `sallseqid` with it), tab and ASCII control characters.
pub fn sanitize(s: &str) -> String {
    s.chars()
        .map(|c| {
            if matches!(c, ' ' | '>' | '[' | ']' | '|' | '"' | ';') || c.is_ascii_control() {
                '_'
            } else {
                c
            }
        })
        .collect()
}

/// Drop ids that repeat an earlier one in the same defline. General ids are compared
/// upper-cased (makeblastdb rejects repeated non-local ids that way), local ids exactly.
pub fn dedup_ids(ids: Vec<SeqId>) -> Vec<SeqId> {
    let mut seen = HashSet::new();
    ids.into_iter()
        .filter(|id| seen.insert(id.dup_key()))
        .collect()
}

/// The ids a collection contributes to its occurrence deflines.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CollectionTags {
    /// First store-unique alias string, else the digest.
    pub db_tag: String,
    /// Every store-unique alias string, in alias order (the first is `db_tag` when there is one).
    pub unique_aliases: Vec<String>,
    /// Alias strings of this collection left out of the ids because another collection uses them too.
    pub shared_aliases: Vec<String>,
}

/// Compute `db_tag` for every collection of a store from `(digest, sanitized alias strings in
/// alias order)`. An alias string is store-unique if, compared upper-cased, it appears in
/// exactly one collection's list and equals no collection digest. Repeats within one
/// collection's list (case-insensitive) count once.
pub fn db_tags(collections: &[(String, Vec<String>)]) -> HashMap<String, CollectionTags> {
    let digests: HashSet<String> = collections.iter().map(|(d, _)| d.to_uppercase()).collect();
    let mut owners: HashMap<String, usize> = HashMap::new();
    let lists: Vec<Vec<String>> = collections
        .iter()
        .map(|(_, aliases)| {
            let mut seen = HashSet::new();
            aliases
                .iter()
                .filter(|a| seen.insert(a.to_uppercase()))
                .cloned()
                .collect()
        })
        .collect();
    for list in &lists {
        for a in list {
            *owners.entry(a.to_uppercase()).or_default() += 1;
        }
    }
    collections
        .iter()
        .zip(lists)
        .map(|((digest, _), list)| {
            let (unique, shared): (Vec<String>, Vec<String>) = list.into_iter().partition(|a| {
                let u = a.to_uppercase();
                owners[&u] == 1 && !digests.contains(&u)
            });
            let db_tag = unique.first().cloned().unwrap_or_else(|| digest.clone());
            (
                digest.clone(),
                CollectionTags {
                    db_tag,
                    unique_aliases: unique,
                    shared_aliases: shared,
                },
            )
        })
        .collect()
}

/// Who an id belongs to, for the uniqueness check.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum Owner {
    /// (collection index, original name): general ids must each come from one such pair.
    Occurrence(u32, String),
    /// OID: local ids must each map to one OID.
    Oid(u32),
}

/// Every LMDB lookup key must name exactly one thing. `T` is extra information kept
/// with each id's first owner, for error messages.
pub struct UniqueIds<T = ()> {
    general: HashMap<String, (Owner, T)>,
    local: HashMap<String, (Owner, T)>,
}

impl<T> Default for UniqueIds<T> {
    fn default() -> Self {
        UniqueIds {
            general: HashMap::new(),
            local: HashMap::new(),
        }
    }
}

impl<T: Clone> UniqueIds<T> {
    /// Register an id. When another owner already holds it, returns that owner and its info.
    /// General ids are compared upper-cased (as makeblastdb's duplicate check does), local ids exactly.
    pub fn add(&mut self, id: &SeqId, owner: Owner, info: T) -> Result<(), (Owner, T)> {
        let (map, key) = match id {
            SeqId::General { .. } => (&mut self.general, id.lmdb_key().to_uppercase()),
            SeqId::Local(s) => (&mut self.local, s.clone()),
        };
        match map.get(&key) {
            Some((prev, info)) if *prev != owner => Err((prev.clone(), info.clone())),
            Some(_) => Ok(()),
            None => {
                map.insert(key, (owner, info));
                Ok(())
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn g(db: &str, tag: &str) -> SeqId {
        SeqId::General {
            db: db.into(),
            tag: tag.into(),
        }
    }

    #[test]
    fn sanitize_table() {
        let table = [
            ("chr1", "chr1"),
            ("ENA|KI270740|KI270740.1", "ENA_KI270740_KI270740.1"),
            ("chrX;alt", "chrX_alt"),
            ("a b", "a_b"),
            ("a>b[c]d\"e", "a_b_c_d_e"),
            ("tab\there", "tab_here"),
            ("ctl\u{1}\u{7f}x", "ctl__x"),
            ("HLA-A*01:01:01:01", "HLA-A*01:01:01:01"),
            ("-dash", "-dash"),
        ];
        for (a, b) in table {
            assert_eq!(sanitize(a), b, "{a:?}");
        }
    }

    #[test]
    fn uniqueness_rules() {
        let mut u: UniqueIds = UniqueIds::default();
        // general ids: case-insensitive, one (collection, name) each
        u.add(&g("X", "chrA"), Owner::Occurrence(0, "chrA".into()), ())
            .unwrap();
        u.add(&g("X", "chrA"), Owner::Occurrence(0, "chrA".into()), ())
            .unwrap();
        let e = u
            .add(&g("x", "CHRA"), Owner::Occurrence(0, "chra".into()), ())
            .unwrap_err();
        assert_eq!(e.0, Owner::Occurrence(0, "chrA".into()));
        assert!(
            u.add(&g("X", "chrA"), Owner::Occurrence(1, "chrA".into()), ())
                .is_err()
        );
        // local ids: exact, one OID each
        u.add(&SeqId::Local("SQ.abc".into()), Owner::Oid(0), ())
            .unwrap();
        u.add(&SeqId::Local("SQ.ABC".into()), Owner::Oid(1), ())
            .unwrap();
        assert_eq!(
            u.add(&SeqId::Local("SQ.abc".into()), Owner::Oid(2), ())
                .unwrap_err()
                .0,
            Owner::Oid(0)
        );
        // repeats within one defline are dropped, general ones case-insensitively
        let ids = vec![
            g("A", "c1"),
            g("a", "C1"),
            SeqId::Local("x".into()),
            SeqId::Local("x".into()),
            SeqId::Local("X".into()),
        ];
        assert_eq!(
            dedup_ids(ids),
            vec![
                g("A", "c1"),
                SeqId::Local("x".into()),
                SeqId::Local("X".into())
            ]
        );
    }

    #[test]
    fn db_tag_rule() {
        let cols = vec![
            (
                "DIGA".to_string(),
                vec!["hg38".to_string(), "HG38".to_string(), "grch38".to_string()],
            ),
            (
                "DIGB".to_string(),
                vec!["shared".to_string(), "mm10".to_string()],
            ),
            ("DIGC".to_string(), vec!["SHARED".to_string()]),
            ("DIGD".to_string(), vec!["diga".to_string()]),
            ("DIGE".to_string(), vec![]),
        ];
        let t = db_tags(&cols);
        assert_eq!(t["DIGA"].db_tag, "hg38");
        assert_eq!(t["DIGA"].unique_aliases, vec!["hg38", "grch38"]); // HG38 is a case repeat
        assert_eq!(t["DIGB"].db_tag, "mm10"); // "shared" is used by DIGC too (case-insensitive)
        assert_eq!(t["DIGB"].shared_aliases, vec!["shared"]);
        assert_eq!(t["DIGC"].db_tag, "DIGC");
        assert_eq!(t["DIGD"].db_tag, "DIGD"); // alias equal to a digest (case-insensitive)
        assert_eq!(t["DIGE"].db_tag, "DIGE");
    }
}

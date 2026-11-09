//! `SubjectPlan::build`: the id scheme and OID order (moved from the prototype's
//! `refget_items`).

use super::ids::{self, Owner, UniqueIds};
use super::{
    DbTotals, Defline, Occurrence, OnOversize, PlanOptions, PlannedCollection, ResolvedSelection,
    SanitizedId, SanitizedKind, SeqId, SkipReason, SkippedSequence, Subject, SubjectPlan, order,
    taxonomy, titles,
};
use crate::convert::nucleotide_bits;
use crate::error::{Error, IdOrigin};
use gtars_refget::store::{ReadonlyRefgetStore, StorageMode};
use std::collections::{BTreeSet, HashMap, HashSet};

struct Row {
    collection: u32,
    name: String,
    description: Option<String>,
}

/// Origin info kept with each registered id: (who, whether sanitizing changed it).
type Info = (IdOrigin, bool);

pub(super) fn build(
    store: &ReadonlyRefgetStore,
    sel: &ResolvedSelection,
    opts: &PlanOptions,
) -> Result<SubjectPlan, Error> {
    let mode = store.storage_mode();
    let calias = order::collection_aliases(store);
    let salias = order::sequence_aliases(store);

    // db_tag is computed once over every collection in the store, never just the selection.
    let all = store.list_collections(0, usize::MAX, &[])?.results;
    let lists: Vec<(String, Vec<String>)> = all
        .iter()
        .map(|m| {
            let s = calias
                .get(&m.digest)
                .map(|v| v.iter().map(|(_, a)| ids::sanitize(a)).collect())
                .unwrap_or_default();
            (m.digest.clone(), s)
        })
        .collect();
    let tags = ids::db_tags(&lists);

    let mut digests: Vec<String> = Vec::new();
    let mut seen = HashSet::new();
    for d in &sel.collections {
        if seen.insert(d.clone()) {
            digests.push(d.clone());
        }
    }
    order::sort_collections(&mut digests, &calias);

    let mut sanitized = Vec::new();
    let mut ambiguous: BTreeSet<String> = BTreeSet::new();
    let mut ambiguous_upper = HashSet::new();
    let mut collections = Vec::with_capacity(digests.len());
    for d in &digests {
        let t = tags
            .get(d)
            .ok_or_else(|| Error::Store(anyhow::anyhow!("collection {d} is not in the store")))?;
        let aliases = calias.get(d).cloned().unwrap_or_default();
        for (_, a) in &aliases {
            let s = ids::sanitize(a);
            if s != *a {
                sanitized.push(SanitizedId {
                    kind: SanitizedKind::CollectionAlias,
                    collection: Some(d.clone()),
                    original: a.clone(),
                    sanitized: s,
                });
            }
        }
        for a in &t.shared_aliases {
            if ambiguous_upper.insert(a.to_uppercase()) {
                ambiguous.insert(a.clone());
            }
        }
        let info = sel.info.get(d).cloned().unwrap_or_default();
        collections.push(PlannedCollection {
            digest: d.clone(),
            display: aliases
                .first()
                .map(|(_, a)| a.clone())
                .unwrap_or_else(|| d.clone()),
            aliases,
            db_tag: t.db_tag.clone(),
            taxid: info.taxid,
            species: info.species,
            assembly: info.assembly,
            n_fraction: info.n_fraction,
            n_oversize: 0,
            n_names: 0,
            n_sequences: 0,
            total_length: 0,
        });
    }

    // Walk collections in canonical order and their rows in .rgsi order; a sequence takes
    // the next OID at its first occurrence.
    let mut seq_order: Vec<String> = Vec::new();
    let mut occ: HashMap<String, Vec<Row>> = HashMap::new();
    for (ci, c) in collections.iter_mut().enumerate() {
        let coll = store.get_collection(&c.digest)?;
        for rec in &coll.sequences {
            let m = rec.metadata();
            c.n_names += 1;
            let e = occ.entry(m.sha512t24u.clone()).or_default();
            if e.is_empty() {
                seq_order.push(m.sha512t24u.clone());
            }
            e.push(Row {
                collection: ci as u32,
                name: m.name.clone(),
                description: m.description.clone(),
            });
        }
    }
    if sel.include_orphans {
        for m in store.list_sequences() {
            if !occ.contains_key(&m.sha512t24u) {
                seq_order.push(m.sha512t24u.clone());
                occ.insert(m.sha512t24u.clone(), Vec::new());
            }
        }
    }

    let mut subjects: Vec<Subject> = Vec::with_capacity(seq_order.len());
    let mut skipped = Vec::new();
    let mut unique: UniqueIds<Info> = UniqueIds::default();
    let mut seen_seq_alias = HashSet::new();
    for d in &seq_order {
        let rows = &occ[d];
        let meta = store.get_sequence_metadata(d).ok_or_else(|| {
            Error::Store(anyhow::anyhow!(
                "sequence {d} is in a collection but not in sequences.rgsi"
            ))
        })?;
        let mut cset: Vec<u32> = rows.iter().map(|r| r.collection).collect();
        cset.sort_unstable();
        cset.dedup();
        let n_coll = cset.len();
        let len = meta.length as u64;
        let bits = nucleotide_bits(meta.alphabet);
        let reason = if bits.is_none() {
            Some(SkipReason::NonNucleotide)
        } else if len == 0 {
            Some(SkipReason::Empty)
        } else if len > opts.max_subject_len.min(i32::MAX as u64) {
            Some(SkipReason::Oversize)
        } else {
            None
        };
        if let Some(reason) = reason {
            if reason == SkipReason::Oversize {
                for r in rows {
                    collections[r.collection as usize].n_oversize += 1;
                }
            }
            skipped.push(SkippedSequence {
                digest: d.clone(),
                length: len,
                alphabet: meta.alphabet,
                reason,
                n_collections: n_coll as u32,
            });
            continue;
        }
        let oid = subjects.len() as u32;
        let length = len as u32;

        // Defline 0: SQ.<digest>, then each sequence alias.
        let mut d0 = vec![SeqId::Local(format!("SQ.{d}"))];
        let mut d0_info = vec![(
            IdOrigin::Sequence {
                digest: d.clone(),
                alias: None,
            },
            false,
        )];
        for (_, a) in salias.get(d).into_iter().flatten() {
            let s = ids::sanitize(a);
            if s != *a && seen_seq_alias.insert(a.clone()) {
                sanitized.push(SanitizedId {
                    kind: SanitizedKind::SequenceAlias,
                    collection: None,
                    original: a.clone(),
                    sanitized: s.clone(),
                });
            }
            d0_info.push((
                IdOrigin::Sequence {
                    digest: d.clone(),
                    alias: Some(a.clone()),
                },
                s != *a,
            ));
            d0.push(SeqId::Local(s));
        }
        let mut deflines = Vec::with_capacity(rows.len() + 1);
        let mut kept = HashSet::new();
        let mut d0_ids = Vec::new();
        for (id, info) in d0.into_iter().zip(d0_info) {
            if !kept.insert(id.dup_key()) {
                continue;
            }
            register(&mut unique, &id, Owner::Oid(oid), info)?;
            d0_ids.push(id);
        }
        // title and taxid are filled in once the occurrences are known
        deflines.push(Defline {
            ids: d0_ids,
            title: String::new(),
            taxid: None,
        });

        // One defline per (collection, name) occurrence.
        let mut occurrences = Vec::with_capacity(rows.len());
        for r in rows {
            let c = &collections[r.collection as usize];
            let t = &tags[&c.digest];
            let id_name = ids::sanitize(&r.name);
            if id_name != r.name {
                sanitized.push(SanitizedId {
                    kind: SanitizedKind::Name,
                    collection: Some(c.digest.clone()),
                    original: r.name.clone(),
                    sanitized: id_name.clone(),
                });
            }
            let mut v = vec![SeqId::General {
                db: t.db_tag.clone(),
                tag: id_name.clone(),
            }];
            for a in &t.unique_aliases {
                if *a != t.db_tag {
                    v.push(SeqId::General {
                        db: a.clone(),
                        tag: id_name.clone(),
                    });
                }
            }
            if t.db_tag != c.digest {
                v.push(SeqId::General {
                    db: c.digest.clone(),
                    tag: id_name.clone(),
                });
            }
            let v = ids::dedup_ids(v);
            let info = (
                IdOrigin::Name {
                    collection: c.digest.clone(),
                    display: c.display.clone(),
                    name: r.name.clone(),
                },
                id_name != r.name,
            );
            for id in &v {
                register(
                    &mut unique,
                    id,
                    Owner::Occurrence(r.collection, r.name.clone()),
                    info.clone(),
                )?;
            }
            deflines.push(Defline {
                ids: v,
                title: titles::occurrence_title(&r.name, r.description.as_deref(), c),
                taxid: Some(taxonomy::collection_defline_taxid(c)),
            });
            occurrences.push(Occurrence {
                collection: r.collection,
                name: r.name.clone(),
                id_name,
            });
        }
        let stored_len = match mode {
            StorageMode::Encoded => Some((len * bits.unwrap()).div_ceil(8)),
            StorageMode::Raw => Some(len),
            StorageMode::Zstd => None,
        };
        let ctax: Vec<i32> = deflines[1..].iter().filter_map(|x| x.taxid).collect();
        deflines[0].taxid = taxonomy::sequence_defline_taxid(&ctax);
        let mut sub = Subject {
            digest: d.clone(),
            md5: meta.md5.clone(),
            alphabet: meta.alphabet,
            length,
            stored_len,
            deflines,
            occurrences,
        };
        sub.deflines[0].title = titles::sequence_title(&sub, &collections);
        subjects.push(sub);
    }

    let oversize: Vec<&SkippedSequence> = skipped
        .iter()
        .filter(|s| s.reason == SkipReason::Oversize)
        .collect();
    if opts.on_oversize == OnOversize::Fail && !oversize.is_empty() {
        return Err(Error::Oversize {
            count: oversize.len() as u64,
            examples: oversize
                .iter()
                .take(10)
                .map(|s| (s.digest.clone(), s.length))
                .collect(),
        });
    }

    for s in &subjects {
        for o in &s.occurrences {
            let c = &mut collections[o.collection as usize];
            c.n_sequences += 1;
            c.total_length += s.length as u64;
        }
    }

    Ok(SubjectPlan {
        store_mode: mode,
        totals: DbTotals::of(&subjects),
        collections,
        subjects,
        skipped,
        sanitized,
        ambiguous_aliases: ambiguous.into_iter().collect(),
    })
}

fn register(
    unique: &mut UniqueIds<Info>,
    id: &SeqId,
    owner: Owner,
    info: Info,
) -> Result<(), Error> {
    let changed = info.1;
    unique
        .add(id, owner, info.clone())
        .map_err(|(_, (first, first_changed))| Error::IdCollision {
            id: id.fasta(),
            first,
            second: info.0,
            sanitized: changed || first_changed,
        })
}

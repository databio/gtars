//! The canonical order of aliases, collections and sequences. All ordering flows from
//! [`alias_namespace_order`], which puts the readable `name` namespace first and the
//! accession namespaces (`insdc`, `refseq`) last.

use gtars_refget::store::ReadonlyRefgetStore;
use std::collections::HashMap;

/// digest -> [(namespace, alias)] in canonical alias order.
pub type ReverseAliases = HashMap<String, Vec<(String, String)>>;

/// Accession namespaces, ranked after every other namespace.
pub const ACCESSION_NAMESPACES: [&str; 2] = ["insdc", "refseq"];

/// The namespace ranking: `name` first; then the other namespaces in manifest order
/// (`rgstore.json` lists them; gtars keeps that order), except the accession namespaces;
/// then `insdc` and `refseq` last, in manifest order.
pub fn alias_namespace_order(manifest: &[String]) -> Vec<String> {
    let is_acc = |n: &String| ACCESSION_NAMESPACES.contains(&n.as_str());
    let mut out: Vec<String> = manifest.iter().filter(|n| *n == "name").cloned().collect();
    out.extend(
        manifest
            .iter()
            .filter(|n| *n != "name" && !is_acc(n))
            .cloned(),
    );
    out.extend(manifest.iter().filter(|n| is_acc(n)).cloned());
    out
}

/// Collection digest -> its aliases, namespaces in canonical order and aliases sorted by
/// byte order within a namespace. Built once per store with one pass per namespace; the
/// per-digest reverse scans of gtars are quadratic over a whole store.
pub fn collection_aliases(store: &ReadonlyRefgetStore) -> ReverseAliases {
    let mut out: ReverseAliases = HashMap::new();
    for ns in alias_namespace_order(store.available_alias_namespaces().collections) {
        let mut aliases = store.list_collection_aliases(&ns).unwrap_or_default();
        aliases.sort();
        for a in aliases {
            if let Some(m) = store.get_collection_metadata_by_alias(&ns, &a) {
                out.entry(m.digest.clone())
                    .or_default()
                    .push((ns.clone(), a));
            }
        }
    }
    out
}

/// Sequence digest -> its aliases, in the same order as [`collection_aliases`].
/// Needs the sequence index loaded.
pub fn sequence_aliases(store: &ReadonlyRefgetStore) -> ReverseAliases {
    let mut out: ReverseAliases = HashMap::new();
    for ns in alias_namespace_order(store.available_alias_namespaces().sequences) {
        let mut aliases = store.list_sequence_aliases(&ns).unwrap_or_default();
        aliases.sort();
        for a in aliases {
            if let Some(m) = store.get_sequence_metadata_by_alias(&ns, &a) {
                out.entry(m.sha512t24u.clone())
                    .or_default()
                    .push((ns.clone(), a));
            }
        }
    }
    out
}

/// Sort collections by (first alias in alias order, else the digest; then the digest).
pub fn sort_collections(digests: &mut [String], aliases: &ReverseAliases) {
    digests.sort_by_cached_key(|d| {
        let first = aliases
            .get(d)
            .and_then(|v| v.first())
            .map(|(_, a)| a.clone())
            .unwrap_or_else(|| d.clone());
        (first, d.clone())
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn namespace_ranking() {
        let v = |x: &[&str]| x.iter().map(|s| s.to_string()).collect::<Vec<_>>();
        // jungle's rgstore.json order
        assert_eq!(
            alias_namespace_order(&v(&["genome_assembly", "insdc", "name", "refseq"])),
            v(&["name", "genome_assembly", "insdc", "refseq"])
        );
        assert_eq!(
            alias_namespace_order(&v(&["refseq", "ucsc", "insdc", "ensembl"])),
            v(&["ucsc", "ensembl", "refseq", "insdc"])
        );
        assert_eq!(alias_namespace_order(&v(&["name"])), v(&["name"]));
        assert!(alias_namespace_order(&[]).is_empty());
    }

    #[test]
    fn collections_sort_by_first_alias_then_digest() {
        let mut a = ReverseAliases::new();
        a.insert("D2".into(), vec![("name".into(), "alpha".into())]);
        a.insert("D3".into(), vec![("name".into(), "B".into())]);
        let mut d = vec!["D1".to_string(), "D2".into(), "D3".into(), "A0".into()];
        sort_collections(&mut d, &a);
        // byte order: "A0" < "B" < "D1" < "alpha"
        assert_eq!(d, vec!["A0", "D3", "D1", "D2"]);
    }
}

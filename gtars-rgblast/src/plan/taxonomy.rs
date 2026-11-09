//! Taxonomy from FHR: parsing `taxon.uri`, the species and assembly a collection contributes,
//! `--taxid-map`, and the per-defline taxid rule.
//!
//! Only the taxid and the species name come from FHR; lineage and LCA come from an NCBI
//! taxdump in the report (plan 5).

use super::{CollectionInfo, Defline};
use crate::error::Error;
use gtars_refget::store::FhrMetadata;
use serde::{Deserialize, Serialize};
use std::path::Path;

/// One resolved `--taxid-map` line: a collection digest and the taxon it is given.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct TaxidMapEntry {
    pub collection: String,
    pub taxid: i32,
    pub species: Option<String>,
}

/// One raw `--taxid-map` line, before its selector is resolved.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TaxidMapLine {
    pub line: usize,
    pub selector: String,
    pub taxid: i32,
    pub species: Option<String>,
}

fn digits_to_taxid(s: &str) -> Option<i32> {
    if s.is_empty() || s.len() > 10 || !s.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let v: u64 = s.parse().ok()?;
    (1..=i32::MAX as u64).contains(&v).then_some(v as i32)
}

/// The NCBI taxid in an FHR `taxon.uri`, compared case-insensitively:
/// `http(s)://identifiers.org/taxonomy:<n>` (or `…/taxonomy/<n>`), `taxonomy:<n>`,
/// `ncbitaxon:<n>`, `http(s)://purl.obolibrary.org/obo/NCBITaxon_<n>`, an NCBI
/// `wwwtax.cgi?…id=<n>` URL, or a bare `<n>`. BLAST taxids are Int4, so 1 ≤ n ≤ 2147483647.
pub fn parse_taxon_uri(uri: &str) -> Option<i32> {
    let u = uri.trim().to_ascii_lowercase();
    if let Some(t) = digits_to_taxid(&u) {
        return Some(t);
    }
    let rest = u
        .strip_prefix("https://")
        .or_else(|| u.strip_prefix("http://"))
        .unwrap_or(&u);
    let rest = rest.strip_prefix("www.").unwrap_or(rest);
    for p in [
        "identifiers.org/taxonomy:",
        "identifiers.org/taxonomy/",
        "purl.obolibrary.org/obo/ncbitaxon_",
    ] {
        if let Some(n) = rest.strip_prefix(p) {
            return digits_to_taxid(n.trim_end_matches('/'));
        }
    }
    if rest.len() == u.len() {
        for p in ["taxonomy:", "ncbitaxon:"] {
            if let Some(n) = u.strip_prefix(p) {
                return digits_to_taxid(n);
            }
        }
    }
    if rest.starts_with("ncbi.nlm.nih.gov/") && rest.contains("wwwtax.cgi") {
        let q = rest.split_once('?')?.1;
        for kv in q.split(['&', ';']) {
            if let Some(v) = kv.strip_prefix("id=") {
                return digits_to_taxid(v);
            }
        }
    }
    None
}

/// What FHR gives a collection: taxid (None plus a warning when `taxon.uri` is not an NCBI
/// taxid), species (`taxon.name`, trimmed) and assembly (`accessionID.name`, else the first
/// `refseq` alias, else the first `insdc` alias). `aliases` are in canonical order.
pub fn info_from_fhr(
    display: &str,
    fhr: Option<&FhrMetadata>,
    aliases: &[(String, String)],
) -> (CollectionInfo, Option<String>) {
    let mut warning = None;
    let taxon = fhr.and_then(|f| f.taxon.as_ref());
    let taxid = match taxon.and_then(|t| t.uri.as_deref()) {
        Some(uri) => {
            let t = parse_taxon_uri(uri);
            if t.is_none() {
                warning = Some(format!(
                    "collection {display}: FHR taxon.uri '{uri}' is not an NCBI taxonomy id; taxid unknown"
                ));
            }
            t
        }
        None => None,
    };
    let species = taxon
        .and_then(|t| t.name.as_deref())
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(String::from);
    let first_in = |ns: &str| {
        aliases
            .iter()
            .find(|(n, _)| n == ns)
            .map(|(_, a)| a.clone())
    };
    let assembly = fhr
        .and_then(|f| f.accession_id.as_ref())
        .and_then(|a| a.name.as_deref())
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(String::from)
        .or_else(|| first_in("refseq"))
        .or_else(|| first_in("insdc"));
    (
        CollectionInfo {
            taxid,
            species,
            assembly,
            n_fraction: None,
        },
        warning,
    )
}

/// Read a `--taxid-map` file: `<collection selector><TAB><taxid>[<TAB><species name>]` per
/// line; blank lines and `#` comments are skipped.
pub fn read_taxid_map(path: &Path) -> Result<Vec<TaxidMapLine>, Error> {
    let text = std::fs::read_to_string(path).map_err(|e| Error::io(path, e))?;
    let mut out = Vec::new();
    for (i, line) in text.lines().enumerate() {
        let t = line.trim_end_matches('\r');
        if t.trim().is_empty() || t.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = t.split('\t').collect();
        let bad = |why: &str| {
            Error::InvalidOptions(format!(
                "--taxid-map {} line {}: {why} (expected <collection><TAB><taxid>[<TAB><species>])",
                path.display(),
                i + 1
            ))
        };
        if f.len() < 2 || f.len() > 3 || f[0].trim().is_empty() {
            return Err(bad("wrong number of fields"));
        }
        let taxid = digits_to_taxid(f[1].trim())
            .ok_or_else(|| bad(&format!("'{}' is not a taxid from 1 to 2147483647", f[1])))?;
        let species = f
            .get(2)
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .map(String::from);
        out.push(TaxidMapLine {
            line: i + 1,
            selector: f[0].trim().to_string(),
            taxid,
            species,
        });
    }
    Ok(out)
}

/// The taxid of a collection's defline: its collection's taxid, or 0 (makeblastdb's own value
/// for "unknown") when the collection has none.
pub fn collection_defline_taxid(c: &super::PlannedCollection) -> i32 {
    c.taxid.unwrap_or(0)
}

/// The taxid of an `SQ.` defline from its sequence's collection-defline taxids: the common
/// taxid when they all agree (0 included), else None (the field is omitted, so `staxids`
/// prints `9606;10090` and not `0;9606;10090`). A sequence in no collection gets 0.
pub fn sequence_defline_taxid(collection_taxids: &[i32]) -> Option<i32> {
    match collection_taxids.split_first() {
        None => Some(0),
        Some((first, rest)) => rest.iter().all(|t| t == first).then_some(*first),
    }
}

/// The sorted unique taxids of an OID's deflines that carry one (`.not` / `.nto` input).
pub fn taxid_set(deflines: &[Defline]) -> Vec<i32> {
    let mut v: Vec<i32> = deflines.iter().filter_map(|d| d.taxid).collect();
    v.sort_unstable();
    v.dedup();
    v
}

#[cfg(test)]
mod tests {
    use super::*;
    use gtars_refget::store::{FhrIdentifier, FhrTaxon};

    #[test]
    fn uri_forms() {
        let ok = [
            ("https://identifiers.org/taxonomy:9606", 9606),
            ("http://identifiers.org/taxonomy:10090", 10090),
            ("HTTPS://IDENTIFIERS.ORG/TAXONOMY:9606", 9606),
            ("https://identifiers.org/taxonomy/9606", 9606),
            ("taxonomy:32630", 32630),
            ("NCBITaxon:9606", 9606),
            ("ncbitaxon:9606", 9606),
            ("http://purl.obolibrary.org/obo/NCBITaxon_9606", 9606),
            ("https://purl.obolibrary.org/obo/NCBITaxon_10090", 10090),
            (
                "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=9606&lvl=3",
                9606,
            ),
            (
                "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=10090",
                10090,
            ),
            ("9606", 9606),
            (" 9606 ", 9606),
            ("2147483647", 2147483647),
            ("1", 1),
        ];
        for (u, t) in ok {
            assert_eq!(parse_taxon_uri(u), Some(t), "{u}");
        }
        for u in [
            "",
            "0",
            "2147483648",
            "-5",
            "Homo sapiens",
            "taxonomy:",
            "taxonomy:abc",
            "https://identifiers.org/taxonomy:",
            "https://example.org/taxonomy:9606",
            "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info",
            "ncbitaxon:0",
            "https://identifiers.org/taxonomy:99999999999",
        ] {
            assert_eq!(parse_taxon_uri(u), None, "{u}");
        }
    }

    fn fhr(name: Option<&str>, uri: Option<&str>, acc: Option<&str>) -> FhrMetadata {
        FhrMetadata {
            taxon: Some(FhrTaxon {
                name: name.map(String::from),
                uri: uri.map(String::from),
            }),
            accession_id: acc.map(|a| FhrIdentifier {
                name: Some(a.into()),
                url: None,
            }),
            ..Default::default()
        }
    }

    #[test]
    fn info_and_assembly_fallbacks() {
        let al = |v: &[(&str, &str)]| {
            v.iter()
                .map(|(a, b)| (a.to_string(), b.to_string()))
                .collect::<Vec<_>>()
        };
        let aliases = al(&[
            ("name", "hg38"),
            ("insdc", "GCA_1.1"),
            ("refseq", "GCF_1.1"),
            ("refseq", "GCF_2.1"),
        ]);
        let f = fhr(
            Some(" Homo sapiens "),
            Some("https://identifiers.org/taxonomy:9606"),
            Some("GCF_000001405.39"),
        );
        let (i, w) = info_from_fhr("hg38", Some(&f), &aliases);
        assert_eq!(i.taxid, Some(9606));
        assert_eq!(i.species.as_deref(), Some("Homo sapiens"));
        assert_eq!(i.assembly.as_deref(), Some("GCF_000001405.39"));
        assert!(w.is_none());
        // no accessionID: first refseq alias, then first insdc alias
        let f = fhr(Some("Homo sapiens"), Some("taxonomy:9606"), None);
        assert_eq!(
            info_from_fhr("hg38", Some(&f), &aliases)
                .0
                .assembly
                .as_deref(),
            Some("GCF_1.1")
        );
        assert_eq!(
            info_from_fhr("x", Some(&f), &al(&[("insdc", "GCA_9.1")]))
                .0
                .assembly
                .as_deref(),
            Some("GCA_9.1")
        );
        assert_eq!(
            info_from_fhr("x", None, &al(&[("name", "x")])).0.assembly,
            None
        );
        // no FHR at all
        let (i, w) = info_from_fhr("x", None, &[]);
        assert_eq!(i, CollectionInfo::default());
        assert!(w.is_none());
        // an unparseable uri warns and leaves the taxid unknown, the species stays
        let f = fhr(
            Some("Homo sapiens"),
            Some("https://example.org/human"),
            None,
        );
        let (i, w) = info_from_fhr("hg38", Some(&f), &[]);
        assert_eq!(
            (i.taxid, i.species.as_deref()),
            (None, Some("Homo sapiens"))
        );
        assert!(w.unwrap().contains("'https://example.org/human'"));
    }

    #[test]
    fn defline_taxid_rule() {
        assert_eq!(sequence_defline_taxid(&[9606, 9606]), Some(9606));
        assert_eq!(sequence_defline_taxid(&[9606, 10090]), None);
        assert_eq!(sequence_defline_taxid(&[0, 0]), Some(0));
        assert_eq!(sequence_defline_taxid(&[0, 9606]), None);
        assert_eq!(sequence_defline_taxid(&[]), Some(0));
        let d = |t: Option<i32>| Defline {
            ids: vec![],
            title: String::new(),
            taxid: t,
        };
        assert_eq!(
            taxid_set(&[d(None), d(Some(10090)), d(Some(9606)), d(Some(9606))]),
            vec![9606, 10090]
        );
        assert_eq!(taxid_set(&[d(Some(0))]), vec![0]);
    }

    #[test]
    fn taxid_map_lines() {
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("m.tsv");
        std::fs::write(
            &p,
            "# comment\n\nname:hg38\t9606\tHomo sapiens\nmm10\t10090\n",
        )
        .unwrap();
        let v = read_taxid_map(&p).unwrap();
        assert_eq!(v.len(), 2);
        assert_eq!(v[0].selector, "name:hg38");
        assert_eq!(
            (v[0].taxid, v[0].species.as_deref()),
            (9606, Some("Homo sapiens"))
        );
        assert_eq!(
            (v[1].line, v[1].taxid, v[1].species.clone()),
            (4, 10090, None)
        );
        for bad in ["x\t0\n", "x\tabc\n", "x\n", "x\t1\ty\tz\n"] {
            std::fs::write(&p, bad).unwrap();
            let e = read_taxid_map(&p).unwrap_err();
            assert!(matches!(e, Error::InvalidOptions(_)), "{bad:?}: {e}");
        }
    }
}

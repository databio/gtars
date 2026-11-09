//! Defline titles. `-outfmt 0` prints each defline as `<first id> <title>`, and `stitle` and
//! `salltitles` print titles without ids, so a title names the species and the assembly by
//! itself and never repeats the name that is already the id tag.

use super::{PlannedCollection, Subject};

/// Tabs become three spaces (titles end up in tab-separated tables).
fn clean(s: &str) -> String {
    s.replace('\t', "   ")
}

/// Title of a per-collection defline: `{species} {display}`, then ` ({assembly})` when the
/// assembly accession is known and differs from `display`, then ` {description}` when the row
/// has one. The species is left out when unknown; the name is never repeated.
pub fn occurrence_title(_name: &str, description: Option<&str>, c: &PlannedCollection) -> String {
    let mut s = String::new();
    if let Some(sp) = &c.species {
        s += sp;
        s.push(' ');
    }
    s += &c.display;
    if let Some(a) = &c.assembly
        && *a != c.display
    {
        s += &format!(" ({a})");
    }
    if let Some(d) = description
        && !d.is_empty()
    {
        s.push(' ');
        s += d;
    }
    clean(&s)
}

/// Title of defline 0: `{species list} {first name}; {n} assemblies; length:{L}; md5:{md5}`.
/// The species list holds the distinct known species in collection order (`, `-joined), else
/// `unknown species`; the first name is the row name in the first collection; n counts the
/// distinct selected collections holding the sequence (`1 assembly` when one). A sequence in
/// no collection: `no collection; length:{L}; md5:{md5}`.
pub fn sequence_title(subject: &Subject, collections: &[PlannedCollection]) -> String {
    let (length, md5) = (subject.length, &subject.md5);
    let Some(first) = subject.occurrences.first() else {
        return format!("no collection; length:{length}; md5:{md5}");
    };
    let mut colls: Vec<u32> = subject.occurrences.iter().map(|o| o.collection).collect();
    colls.sort_unstable();
    colls.dedup();
    let mut species: Vec<&str> = Vec::new();
    for &ci in &colls {
        if let Some(sp) = &collections[ci as usize].species
            && !species.contains(&sp.as_str())
        {
            species.push(sp);
        }
    }
    let sp = if species.is_empty() {
        "unknown species".to_string()
    } else {
        species.join(", ")
    };
    let n = colls.len();
    clean(&format!(
        "{sp} {}; {n} assembl{}; length:{length}; md5:{md5}",
        first.name,
        if n == 1 { "y" } else { "ies" }
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plan::Occurrence;
    use gtars_refget::digest::AlphabetType;

    fn coll(display: &str, species: Option<&str>, assembly: Option<&str>) -> PlannedCollection {
        PlannedCollection {
            digest: format!("D_{display}"),
            aliases: vec![],
            db_tag: display.into(),
            display: display.into(),
            taxid: None,
            species: species.map(String::from),
            assembly: assembly.map(String::from),
            n_fraction: None,
            n_oversize: 0,
            n_names: 0,
            n_sequences: 0,
            total_length: 0,
        }
    }

    #[test]
    fn occurrence_titles() {
        let hg = coll("hg38-p14-ucsc", Some("Homo sapiens"), None);
        assert_eq!(
            occurrence_title("chr1", None, &hg),
            "Homo sapiens hg38-p14-ucsc"
        );
        let grc = coll(
            "GRCh38.p13-fasta-genomic",
            Some("Homo sapiens"),
            Some("GCF_000001405.39"),
        );
        assert_eq!(
            occurrence_title("NC_000001.11", None, &grc),
            "Homo sapiens GRCh38.p13-fasta-genomic (GCF_000001405.39)"
        );
        // assembly equal to display is not repeated; unknown species is left out
        let same = coll("GCF_X.1", None, Some("GCF_X.1"));
        assert_eq!(occurrence_title("c1", Some(""), &same), "GCF_X.1");
        assert_eq!(
            occurrence_title("c1", Some("a\tcontig"), &same),
            "GCF_X.1 a   contig"
        );
    }

    #[test]
    fn sequence_titles() {
        let cs = vec![
            coll("hg38", Some("Homo sapiens"), None),
            coll("hg19", Some("Homo sapiens"), None),
            coll("mm10", Some("Mus musculus"), None),
            coll("x", None, None),
        ];
        let occ = |c: u32, n: &str| Occurrence {
            collection: c,
            name: n.into(),
            id_name: n.into(),
        };
        let mut s = Subject {
            digest: "d".into(),
            md5: "m5".into(),
            alphabet: AlphabetType::Dna2bit,
            length: 8,
            stored_len: None,
            deflines: vec![],
            occurrences: vec![
                occ(0, "chr1"),
                occ(1, "1"),
                occ(2, "decoy"),
                occ(2, "decoy2"),
            ],
        };
        assert_eq!(
            sequence_title(&s, &cs),
            "Homo sapiens, Mus musculus chr1; 3 assemblies; length:8; md5:m5"
        );
        s.occurrences = vec![occ(3, "u1")];
        assert_eq!(
            sequence_title(&s, &cs),
            "unknown species u1; 1 assembly; length:8; md5:m5"
        );
        s.occurrences.clear();
        assert_eq!(sequence_title(&s, &cs), "no collection; length:8; md5:m5");
    }
}

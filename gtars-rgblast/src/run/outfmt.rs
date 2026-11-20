//! `-outfmt` parsing and rewriting: `search` appends the columns `rgblast report` needs, and the
//! report maps `-outfmt 7` long column names back to their short names.

/// What `std` expands to.
pub const STD_FIELDS: [&str; 12] = [
    "qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart",
    "send", "evalue", "bitscore",
];

/// Columns `search` appends when missing: `sallseqid` (the fan-out source), `salltitles`
/// (original names), `staxids` (fallback taxa), `qlen` (coverage), `score` (exact ties).
pub const REPORT_FIELDS: [&str; 5] = ["qlen", "score", "sallseqid", "salltitles", "staxids"];

/// The fields origin always runs with (format 6: origin reads only the hit lines, and
/// rgblastn writes no other tabular format).
pub const ORIGIN_OUTFMT: &str = "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen score sallseqid salltitles staxids";

/// `-outfmt 7` `# Fields:` long names (blastn 2.17) and their short names.
pub const LONG_NAMES: &[(&str, &str)] = &[
    ("query id", "qseqid"),
    ("query gi", "qgi"),
    ("query acc.", "qacc"),
    ("query acc.ver", "qaccver"),
    ("query length", "qlen"),
    ("subject id", "sseqid"),
    ("subject ids", "sallseqid"),
    ("subject gi", "sgi"),
    ("subject gis", "sallgi"),
    ("subject acc.", "sacc"),
    ("subject acc.ver", "saccver"),
    ("subject accs.", "sallacc"),
    ("subject length", "slen"),
    ("q. start", "qstart"),
    ("q. end", "qend"),
    ("s. start", "sstart"),
    ("s. end", "send"),
    ("query seq", "qseq"),
    ("subject seq", "sseq"),
    ("evalue", "evalue"),
    ("bit score", "bitscore"),
    ("score", "score"),
    ("alignment length", "length"),
    ("% identity", "pident"),
    ("identical", "nident"),
    ("mismatches", "mismatch"),
    ("positives", "positive"),
    ("gap opens", "gapopen"),
    ("gaps", "gaps"),
    ("% positives", "ppos"),
    ("query/sbjct frames", "frames"),
    ("query frame", "qframe"),
    ("sbjct frame", "sframe"),
    ("BTOP", "btop"),
    ("subject tax id", "staxid"),
    ("subject sci name", "ssciname"),
    ("subject com names", "scomname"),
    ("subject blast name", "sblastname"),
    ("subject super kingdom", "sskingdom"),
    ("subject tax ids", "staxids"),
    ("subject sci names", "sscinames"),
    ("subject com names", "scomnames"),
    ("subject blast names", "sblastnames"),
    ("subject super kingdoms", "sskingdoms"),
    ("subject title", "stitle"),
    ("subject titles", "salltitles"),
    ("subject strand", "sstrand"),
    ("% query coverage per subject", "qcovs"),
    ("% query coverage per hsp", "qcovhsp"),
    ("% query coverage per uniq subject", "qcovus"),
];

/// The short name of a `# Fields:` long name (unknown names are returned unchanged).
pub fn short_name(long: &str) -> String {
    let l = long.trim();
    LONG_NAMES
        .iter()
        .find(|(k, _)| *k == l)
        .map(|(_, s)| s.to_string())
        .unwrap_or_else(|| l.to_string())
}

/// A field list in blastn syntax (`std` expanded), duplicates kept.
pub fn expand_fields(list: &str) -> Vec<String> {
    let mut v = Vec::new();
    for t in list.split_whitespace() {
        if t == "std" {
            v.extend(STD_FIELDS.iter().map(|s| s.to_string()));
        } else {
            v.push(t.to_string());
        }
    }
    v
}

/// What `rewrite` did to `-outfmt`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Outfmt {
    /// No `-outfmt`: blastn's pairwise default, passed unchanged.
    Absent,
    /// Tabular 6 or 7; `added` lists the appended fields.
    Tabular {
        format: u8,
        value: String,
        added: Vec<String>,
    },
    /// Another format, passed unchanged; `rgblast report` can't read it.
    Unreadable(String),
}

impl Outfmt {
    /// The stderr note for an unreadable format.
    pub fn note(&self) -> Option<String> {
        match self {
            Outfmt::Unreadable(v) => Some(format!(
                "rgblast: -outfmt {v} is not readable by 'rgblast report'; use 6 or 7"
            )),
            _ => None,
        }
    }
}

/// Rewrite the `-outfmt` value in `args` in place (see `rewrite_value`).
pub fn rewrite(args: &mut [String]) -> Outfmt {
    let Some(i) = args.iter().position(|a| a == "-outfmt") else {
        return Outfmt::Absent;
    };
    let Some(v) = args.get(i + 1).cloned() else {
        return Outfmt::Absent;
    };
    let r = rewrite_value(&v);
    if let Outfmt::Tabular { value, .. } = &r {
        args[i + 1] = value.clone();
    }
    r
}

/// For `6`/`7` (with or without a field list) append whichever `REPORT_FIELDS` are missing, so
/// the user's columns keep their positions; `6` alone becomes `6 std qlen score ...`.
pub fn rewrite_value(v: &str) -> Outfmt {
    let mut tok = v.split_whitespace();
    let Some(first) = tok.next() else {
        return Outfmt::Unreadable(v.to_string());
    };
    let format = match first {
        "6" => 6,
        "7" => 7,
        _ => return Outfmt::Unreadable(v.to_string()),
    };
    let rest: Vec<&str> = tok.collect();
    let mut out: Vec<String> = vec![first.to_string()];
    if rest.is_empty() {
        out.push("std".into());
    } else {
        out.extend(rest.iter().map(|s| s.to_string()));
    }
    let have = expand_fields(&out[1..].join(" "));
    let added: Vec<String> = REPORT_FIELDS
        .iter()
        .filter(|f| !have.iter().any(|h| h == *f))
        .map(|f| f.to_string())
        .collect();
    out.extend(added.iter().cloned());
    Outfmt::Tabular {
        format,
        value: out.join(" "),
        added,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn val(v: &str) -> String {
        match rewrite_value(v) {
            Outfmt::Tabular { value, .. } => value,
            o => panic!("{o:?}"),
        }
    }

    #[test]
    fn rewrites_tabular() {
        assert_eq!(val("6"), "6 std qlen score sallseqid salltitles staxids");
        assert_eq!(val("7"), "7 std qlen score sallseqid salltitles staxids");
        assert_eq!(
            val("6 qseqid sseqid bitscore staxids"),
            "6 qseqid sseqid bitscore staxids qlen score sallseqid salltitles"
        );
        // nothing to add, no duplicates
        assert_eq!(
            val("6 std qlen score sallseqid salltitles staxids"),
            "6 std qlen score sallseqid salltitles staxids"
        );
        assert_eq!(
            val("7 qlen std sallseqid"),
            "7 qlen std sallseqid score salltitles staxids"
        );
        assert_eq!(expand_fields("std qlen").len(), 13);
    }

    #[test]
    fn other_formats_pass() {
        for v in ["0", "5", "10", "10 qseqid", "11", "17"] {
            let o = rewrite_value(v);
            assert_eq!(o, Outfmt::Unreadable(v.into()));
            assert!(o.note().unwrap().contains("use 6 or 7"));
        }
        let mut a: Vec<String> = vec!["-query".into(), "q.fa".into()];
        assert_eq!(rewrite(&mut a), Outfmt::Absent);
        let mut a: Vec<String> = vec!["-outfmt".into(), "6".into(), "-evalue".into()];
        rewrite(&mut a);
        assert_eq!(a[1], "6 std qlen score sallseqid salltitles staxids");
    }

    #[test]
    fn long_names() {
        assert_eq!(short_name("subject ids"), "sallseqid");
        assert_eq!(short_name("% query coverage per hsp"), "qcovhsp");
        assert_eq!(short_name("query acc.ver"), "qaccver");
        assert_eq!(short_name("made up"), "made up");
    }
}

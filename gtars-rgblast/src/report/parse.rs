//! Reading blastn tabular output (`-outfmt 6` and `7`).

use crate::error::Error;
use crate::run::outfmt::{expand_fields, short_name};
use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

/// What `search` produces from `-outfmt 6`.
pub const DEFAULT_FIELDS: &str = "std qlen score sallseqid salltitles staxids";

pub const REQUIRED: [&str; 10] = [
    "pident",
    "length",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "sallseqid",
];

/// One HSP line. Printed numbers keep their raw strings for output.
#[derive(Debug, Clone, PartialEq)]
pub struct Hit {
    pub query: String,
    pub subject: Option<String>,
    pub pident: f64,
    pub pident_s: String,
    pub length: u64,
    pub qstart: u64,
    pub qend: u64,
    pub sstart: u64,
    pub send: u64,
    pub evalue: f64,
    pub evalue_s: String,
    pub bits: f64,
    pub bits_s: String,
    pub qlen: u64,
    pub score: Option<f64>,
    /// One id per defline (the first), in order.
    pub sallseqid: Vec<String>,
    pub salltitles: Option<String>,
    pub staxids: Vec<u32>,
}

impl Hit {
    /// Query coverage of this HSP, percent.
    pub fn qcov(&self) -> f64 {
        if self.qlen == 0 {
            return 0.0;
        }
        (self.qend.abs_diff(self.qstart) + 1) as f64 / self.qlen as f64 * 100.0
    }
}

#[derive(Debug, Clone, Default)]
pub struct HitsInput {
    pub hits: Vec<Hit>,
    /// Queries named by `# Query:` lines (first token), in order.
    pub queries: Vec<String>,
}

struct Columns {
    idx: HashMap<String, usize>,
    query: usize,
}

fn columns(fields: &[String], source: &str) -> Result<Columns, Error> {
    let mut idx = HashMap::new();
    for (i, f) in fields.iter().enumerate() {
        idx.entry(f.clone()).or_insert(i);
    }
    let mut missing: Vec<&str> = Vec::new();
    let query = idx.get("qaccver").or_else(|| idx.get("qseqid")).copied();
    if query.is_none() {
        missing.push("qseqid");
    }
    for r in REQUIRED {
        if !idx.contains_key(r) {
            missing.push(r);
        }
    }
    if !missing.is_empty() {
        return Err(Error::InvalidOptions(format!(
            "{source}: report needs columns: {} (add it to -outfmt, or produce hits with 'rgblast search')",
            missing.join(", ")
        )));
    }
    Ok(Columns {
        idx,
        query: query.unwrap(),
    })
}

/// Parse a hits stream. `fields` are the `-outfmt 6` columns (blastn syntax, `std` allowed);
/// `# Fields:` lines (`-outfmt 7`) replace them.
pub fn read_hits(r: impl BufRead, fields: &str, source: &str) -> Result<HitsInput, Error> {
    let mut cols: Option<Columns> = None;
    let mut given = Some(expand_fields(fields));
    let mut out = HitsInput::default();
    let bad = |n: usize, msg: String| Error::InvalidOptions(format!("{source}: line {n}: {msg}"));
    for (n, line) in r.lines().enumerate() {
        let n = n + 1;
        let line = line.map_err(|e| Error::io(source, e))?;
        if let Some(c) = line.strip_prefix('#') {
            let c = c.trim();
            if let Some(f) = c.strip_prefix("Fields:") {
                let f: Vec<String> = f.split(", ").map(short_name).collect();
                cols = Some(columns(&f, source)?);
                given = None;
            } else if let Some(q) = c.strip_prefix("Query:")
                && let Some(t) = q.split_whitespace().next()
            {
                out.queries.push(t.to_string());
            }
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        if cols.is_none() {
            cols = Some(columns(given.as_ref().unwrap(), source)?);
        }
        let c = cols.as_ref().unwrap();
        let f: Vec<&str> = line.split('\t').collect();
        let get = |name: &str| -> Option<&str> { c.idx.get(name).and_then(|&i| f.get(i).copied()) };
        let req = |name: &str| -> Result<&str, Error> {
            get(name).ok_or_else(|| bad(n, format!("no {name} column ({} fields)", f.len())))
        };
        let num_u = |name: &str| -> Result<u64, Error> {
            let v = req(name)?;
            v.trim()
                .parse()
                .map_err(|_| bad(n, format!("{name} '{v}' is not a number")))
        };
        let num_f = |name: &str| -> Result<(f64, String), Error> {
            let v = req(name)?.trim();
            let x: f64 = v
                .parse()
                .map_err(|_| bad(n, format!("{name} '{v}' is not a number")))?;
            Ok((x, v.to_string()))
        };
        let (pident, pident_s) = num_f("pident")?;
        let (evalue, evalue_s) = num_f("evalue")?;
        let (bits, bits_s) = num_f("bitscore")?;
        let score = match get("score") {
            Some(v) => Some(
                v.trim()
                    .parse::<f64>()
                    .map_err(|_| bad(n, format!("score '{v}' is not a number")))?,
            ),
            None => None,
        };
        let staxids = get("staxids")
            .map(|v| {
                v.split(';')
                    .filter_map(|t| t.trim().parse::<u32>().ok())
                    .filter(|t| *t > 0)
                    .collect()
            })
            .unwrap_or_default();
        out.hits.push(Hit {
            query: f
                .get(c.query)
                .map(|s| s.to_string())
                .ok_or_else(|| bad(n, "no query column".into()))?,
            subject: get("saccver").or_else(|| get("sseqid")).map(str::to_string),
            pident,
            pident_s,
            length: num_u("length")?,
            qstart: num_u("qstart")?,
            qend: num_u("qend")?,
            sstart: num_u("sstart")?,
            send: num_u("send")?,
            evalue,
            evalue_s,
            bits,
            bits_s,
            qlen: num_u("qlen")?,
            score,
            sallseqid: req("sallseqid")?
                .split(';')
                .filter(|s| !s.is_empty())
                .map(str::to_string)
                .collect(),
            salltitles: get("salltitles").map(str::to_string),
            staxids,
        });
    }
    Ok(out)
}

/// Read a hits file (`-` = stdin).
pub fn read_hits_file(path: &Path, fields: &str) -> Result<HitsInput, Error> {
    let name = path.display().to_string();
    if name == "-" {
        let stdin = std::io::stdin();
        return read_hits(stdin.lock(), fields, "<stdin>");
    }
    let f = std::fs::File::open(path).map_err(|e| Error::io(path, e))?;
    read_hits(std::io::BufReader::new(f), fields, &name)
}

/// Query names (first token of each defline) of a FASTA file, in order.
pub fn fasta_names(path: &Path) -> Result<Vec<String>, Error> {
    let t = std::fs::read_to_string(path).map_err(|e| Error::io(path, e))?;
    Ok(t.lines()
        .filter_map(|l| l.strip_prefix('>'))
        .map(|l| l.split_whitespace().next().unwrap_or("").to_string())
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    const LINE: &str = "q1\tSQ.abc\t100.000\t1000\t0\t0\t1\t1000\t2001\t3000\t0.0\t1847\t1000\t1000\tSQ.abc;gnl|hg|chrM\tt1<>t2\t9606";

    #[test]
    fn outfmt6_default_fields() {
        let h = read_hits(LINE.as_bytes(), DEFAULT_FIELDS, "x").unwrap();
        let h = &h.hits[0];
        assert_eq!(h.query, "q1");
        assert_eq!(h.sallseqid, vec!["SQ.abc", "gnl|hg|chrM"]);
        assert_eq!(
            (h.pident_s.as_str(), h.bits_s.as_str()),
            ("100.000", "1847")
        );
        assert_eq!(h.score, Some(1000.0));
        assert_eq!(h.staxids, vec![9606]);
        assert_eq!(h.qcov(), 100.0);
    }

    #[test]
    fn outfmt7_fields_and_queries() {
        let t = "# BLASTN 2.17.0+\n# Query: Q1 human chrM\n# Database: x\n# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, query length, subject ids\n# 1 hits found\nQ1\tSQ.a\t99.5\t200\t1\t0\t1\t200\t5\t204\t1e-90\t350\t400\tSQ.a\n# Query: Q2\n# 0 hits found\n";
        let h = read_hits(t.as_bytes(), DEFAULT_FIELDS, "x").unwrap();
        assert_eq!(h.queries, vec!["Q1", "Q2"]);
        assert_eq!(h.hits.len(), 1);
        assert_eq!(h.hits[0].qcov(), 50.0);
        assert_eq!(h.hits[0].score, None);
    }

    #[test]
    fn missing_column_is_named() {
        let e = read_hits(LINE.as_bytes(), "std sallseqid", "hits.tsv").unwrap_err();
        assert_eq!(
            e.to_string(),
            "hits.tsv: report needs columns: qlen (add it to -outfmt, or produce hits with 'rgblast search')"
        );
    }
}

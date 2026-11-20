//! Origin's query file: records renamed `>Q<index>` (the original defline kept after a space),
//! so stage 2 can take an exact subset, stdin works, and no-hit queries keep their order.

use crate::error::Error;
use std::collections::HashSet;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QueryRec {
    /// `Q1`, `Q2`, ...
    pub id: String,
    /// The original first token.
    pub name: String,
    /// The original defline (without `>`).
    pub defline: String,
    pub seq: String,
}

/// Read FASTA records from a reader.
pub fn parse_queries(r: impl BufRead, source: &str) -> Result<Vec<QueryRec>, Error> {
    let mut v: Vec<QueryRec> = Vec::new();
    for line in r.lines() {
        let line = line.map_err(|e| Error::io(source, e))?;
        let line = line.trim_end();
        if let Some(d) = line.strip_prefix('>') {
            let i = v.len() + 1;
            v.push(QueryRec {
                id: format!("Q{i}"),
                name: d.split_whitespace().next().unwrap_or("").to_string(),
                defline: d.to_string(),
                seq: String::new(),
            });
        } else if !line.is_empty() {
            match v.last_mut() {
                Some(q) => q.seq.push_str(line.trim()),
                None => {
                    return Err(Error::InvalidOptions(format!(
                        "{source}: the query file is not FASTA (no '>' before the first sequence)"
                    )));
                }
            }
        }
    }
    if v.is_empty() {
        return Err(Error::InvalidOptions(format!(
            "{source}: no query sequences"
        )));
    }
    Ok(v)
}

/// Read the query FASTA (`-` = stdin).
pub fn read_queries(path: &Path) -> Result<Vec<QueryRec>, Error> {
    if path.as_os_str() == "-" {
        return parse_queries(std::io::stdin().lock(), "<stdin>");
    }
    let f = std::fs::File::open(path).map_err(|e| Error::io(path, e))?;
    parse_queries(BufReader::new(f), &path.display().to_string())
}

/// Write the records (all, or the indices in `only`) as `>Q<i> <defline>`.
pub fn write_queries(
    recs: &[QueryRec],
    path: &Path,
    only: Option<&HashSet<usize>>,
) -> Result<(), Error> {
    let f = std::fs::File::create(path).map_err(|e| Error::io(path, e))?;
    let mut w = std::io::BufWriter::new(f);
    let io = |e| Error::io(path, e);
    for (i, q) in recs.iter().enumerate() {
        if only.is_some_and(|o| !o.contains(&i)) {
            continue;
        }
        writeln!(w, ">{} {}", q.id, q.defline).map_err(io)?;
        for chunk in q.seq.as_bytes().chunks(80) {
            w.write_all(chunk).map_err(io)?;
            w.write_all(b"\n").map_err(io)?;
        }
    }
    w.flush().map_err(io)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rename_and_subset() {
        let q = parse_queries(">a one\nACGT\nAC\n>b\nGG\n".as_bytes(), "x").unwrap();
        assert_eq!(q[0].id, "Q1");
        assert_eq!(q[0].name, "a");
        assert_eq!(q[0].seq, "ACGTAC");
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("q.fa");
        write_queries(&q, &p, Some(&[1].into_iter().collect())).unwrap();
        assert_eq!(std::fs::read_to_string(&p).unwrap(), ">Q2 b\nGG\n");
        assert!(parse_queries("ACGT\n".as_bytes(), "x").is_err());
    }
}

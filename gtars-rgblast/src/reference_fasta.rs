//! A FASTA that `makeblastdb -parse_seqids` turns into the same database as the plan
//! (feature `reference-fasta`; used by the BLAST+ tests and the rgblastn golden harness).

use crate::error::Error;
use crate::pipeline::run_ordered;
use crate::plan::SubjectPlan;
use crate::store::decode;
use crate::store::fetch::SeqFetch;
use std::io::Write;
use std::path::Path;

/// The FASTA header of one subject (without '>'): its deflines joined by `\x01`, each
/// `ids.join("|") + " " + title`.
pub fn header(plan: &SubjectPlan, oid: usize) -> String {
    plan.subjects[oid]
        .deflines
        .iter()
        .map(|d| {
            let ids: Vec<String> = d.ids.iter().map(|i| i.fasta()).collect();
            format!("{} {}", ids.join("|"), d.title)
        })
        .collect::<Vec<_>>()
        .join("\x01")
}

/// Write one record per subject, in OID order, sequence from `decode::to_ascii`, 80 columns.
pub fn write(plan: &SubjectPlan, fetch: &dyn SeqFetch, path: &Path) -> Result<(), Error> {
    let f = std::fs::File::create(path).map_err(|e| Error::io(path, e))?;
    let mut w = std::io::BufWriter::with_capacity(8 << 20, f);
    let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4)
        .min(16);
    let idx: Vec<usize> = (0..plan.subjects.len()).collect();
    run_ordered(
        &idx,
        threads,
        1 << 30,
        |&i| plan.subjects[i].length as u64 * 2 + 1024,
        |&i| {
            let s = &plan.subjects[i];
            let stored = fetch
                .read(&s.digest)
                .map_err(|e| e.with_digest(&s.digest))?;
            let ascii = decode::to_ascii(plan.store_mode, s.alphabet, &stored, s.length)
                .map_err(|e| e.with_digest(&s.digest))?;
            let h = header(plan, i);
            let mut o = Vec::with_capacity(ascii.len() + ascii.len() / 80 + h.len() + 4);
            o.push(b'>');
            o.extend_from_slice(h.as_bytes());
            o.push(b'\n');
            for line in ascii.chunks(80) {
                o.extend_from_slice(line);
                o.push(b'\n');
            }
            Ok(o)
        },
        |_, o| w.write_all(&o).map_err(|e| Error::io(path, e)),
    )?;
    w.flush().map_err(|e| Error::io(path, e))
}

/// A makeblastdb `-taxid_map` (`<id> <taxid>` lines) that, with [`write`]'s FASTA, gives the
/// plan's taxids: `SQ.<digest>` and each local id by itself, each `gnl|db|tag` id by its bare
/// tag (makeblastdb looks general ids up by the tag alone, and the first line for a tag wins).
/// Taxid 0 (makeblastdb's default) is not written. Plans this map can't express are an error:
/// a bare tag whose deflines need two taxids (`chr1` in human and mouse), or a mixed `SQ.`
/// defline (taxid omitted), which makeblastdb would give a taxid anyway.
pub fn write_taxid_map(plan: &SubjectPlan, path: &Path) -> Result<(), Error> {
    use crate::plan::SeqId;
    let mut seen: std::collections::HashMap<String, i32> = std::collections::HashMap::new();
    let mut out = String::new();
    for (oid, s) in plan.subjects.iter().enumerate() {
        for d in &s.deflines {
            let t = d.taxid.ok_or_else(|| {
                Error::InvalidOptions(format!(
                    "-taxid_map can't express OID {oid} (SQ.{}): its collections have different taxids, so its SQ. defline has none",
                    s.digest
                ))
            })?;
            for id in &d.ids {
                let key = match id {
                    SeqId::Local(x) => x.clone(),
                    SeqId::General { tag, .. } => tag.clone(),
                };
                match seen.get(&key) {
                    Some(&prev) if prev != t => {
                        return Err(Error::InvalidOptions(format!(
                            "-taxid_map can't express id '{key}': makeblastdb looks it up by its bare tag, which needs taxids {prev} and {t}"
                        )));
                    }
                    Some(_) => {}
                    None => {
                        seen.insert(key.clone(), t);
                        if t != 0 {
                            out += &format!("{key} {t}\n");
                        }
                    }
                }
            }
        }
    }
    std::fs::write(path, out).map_err(|e| Error::io(path, e))
}

//! Version 5 side files makeblastdb also writes: `.nos` (OID -> seqids), `.not` (OID ->
//! taxids), `.nto` (taxid -> OIDs) and `.ntf` (LMDB taxid2offset), from per-OID taxid sets
//! (the sorted unique taxids of each OID's deflines that carry one).
//!
//! - `.not`: u64 LE number of OIDs; one u64 LE cumulative end (counted in entries) per OID;
//!   then each OID's taxids as i32 LE.
//! - `.nto`: for each taxid in ascending order, a u32 LE count, then its OIDs (u32 LE, ascending).
//! - `.ntf`: LMDB `taxid2offset`: i32 LE taxid -> u64 LE byte offset of its count in `.nto`.

use super::lmdb;
use super::output::{PartialOutput, with_ext};
use crate::error::Error;
use crate::plan::Subject;
use std::path::{Path, PathBuf};

/// The per-OID taxid sets of a plan's subjects.
pub fn taxid_sets(subjects: &[Subject]) -> Vec<Vec<i32>> {
    subjects
        .iter()
        .map(|s| crate::plan::taxonomy::taxid_set(&s.deflines))
        .collect()
}

/// Distinct taxids over all OIDs (0 included), ascending.
pub fn distinct_taxids(sets: &[Vec<i32>]) -> Vec<i32> {
    let mut v: Vec<i32> = sets.iter().flatten().copied().collect();
    v.sort_unstable();
    v.dedup();
    v
}

/// `.not` bytes.
pub fn not_bytes(sets: &[Vec<i32>]) -> Vec<u8> {
    let entries: usize = sets.iter().map(|s| s.len()).sum();
    let mut w = Vec::with_capacity(8 + 8 * sets.len() + 4 * entries);
    w.extend_from_slice(&(sets.len() as u64).to_le_bytes());
    let mut end = 0u64;
    for s in sets {
        end += s.len() as u64;
        w.extend_from_slice(&end.to_le_bytes());
    }
    for s in sets {
        for t in s {
            w.extend_from_slice(&t.to_le_bytes());
        }
    }
    w
}

/// `.nto` bytes and the `.ntf` entries (taxid, byte offset of its count).
pub fn nto_bytes(sets: &[Vec<i32>]) -> (Vec<u8>, Vec<(i32, u64)>) {
    // counting sort over (taxid, OID) pairs: OIDs come out ascending within each taxid
    let taxa = distinct_taxids(sets);
    let index: std::collections::HashMap<i32, usize> =
        taxa.iter().enumerate().map(|(i, t)| (*t, i)).collect();
    let mut counts = vec![0u32; taxa.len()];
    for s in sets {
        for t in s {
            counts[index[t]] += 1;
        }
    }
    let mut starts = vec![0usize; taxa.len()];
    let mut acc = 0usize;
    for (i, c) in counts.iter().enumerate() {
        starts[i] = acc;
        acc += *c as usize;
    }
    let mut oids = vec![0u32; acc];
    let mut fill = starts.clone();
    for (oid, s) in sets.iter().enumerate() {
        for t in s {
            let k = index[t];
            oids[fill[k]] = oid as u32;
            fill[k] += 1;
        }
    }
    let mut w = Vec::with_capacity(4 * taxa.len() + 4 * acc);
    let mut offsets = Vec::with_capacity(taxa.len());
    for (i, t) in taxa.iter().enumerate() {
        offsets.push((*t, w.len() as u64));
        w.extend_from_slice(&counts[i].to_le_bytes());
        for o in &oids[starts[i]..starts[i] + counts[i] as usize] {
            w.extend_from_slice(&o.to_le_bytes());
        }
    }
    (w, offsets)
}

/// Size of `.nos` for these ids: u64 count, u64 end offset per OID, then each id with a
/// length byte (or 0xFF + u32 when >= 255).
pub fn nos_size(n_oids: u64, id_lens: impl Iterator<Item = usize>) -> u64 {
    8 + 8 * n_oids + id_lens.map(nos_id_size).sum::<u64>()
}
/// Bytes one id takes in `.nos`: its length field and its bytes.
pub fn nos_id_size(len: usize) -> u64 {
    if len >= 255 {
        5 + len as u64
    } else {
        1 + len as u64
    }
}
/// Size of `.not` for these taxid sets.
pub fn not_size(sets: &[Vec<i32>]) -> u64 {
    8 + 8 * sets.len() as u64 + 4 * sets.iter().map(|s| s.len() as u64).sum::<u64>()
}
/// Size of `.nto` for these taxid sets.
pub fn nto_size(sets: &[Vec<i32>]) -> u64 {
    4 * distinct_taxids(sets).len() as u64 + 4 * sets.iter().map(|s| s.len() as u64).sum::<u64>()
}

fn write(out: &PartialOutput, p: PathBuf, bytes: &[u8]) -> Result<PathBuf, Error> {
    out.record(&p);
    std::fs::write(&p, bytes).map_err(|e| Error::io(&p, e))?;
    Ok(p)
}

/// Write `.nos`, `.not`, `.nto` and `.ntf` (named after `base`, the LMDB's stem). `.nos` takes
/// the original-spelling keys of `keys`, which must still be in OID order (before
/// [`KeyTable::sort_dedup`](lmdb::KeyTable::sort_dedup)); `taxids` holds one sorted set per OID.
pub fn write_oid_files(
    base: &Path,
    n_oids: u32,
    keys: &lmdb::KeyTable,
    taxids: &[Vec<i32>],
    out: &PartialOutput,
) -> Result<Vec<PathBuf>, Error> {
    if taxids.len() != n_oids as usize {
        return Err(Error::Store(anyhow::anyhow!(
            "internal: {} taxid sets for {n_oids} OIDs",
            taxids.len()
        )));
    }
    let mut files = Vec::new();
    let n = n_oids as usize;
    let mut ids = keys.nos_ids().peekable();
    if ids.peek().is_some() {
        // per OID its ids sorted
        let mut data = Vec::new();
        let mut ends = Vec::with_capacity(n);
        let mut v: Vec<&[u8]> = Vec::new();
        for oid in 0..n_oids {
            while let Some((_, id)) = ids.next_if(|(o, _)| *o == oid) {
                v.push(id);
            }
            v.sort();
            for id in v.drain(..) {
                if id.len() >= 255 {
                    data.push(0xFF);
                    data.extend_from_slice(&(id.len() as u32).to_le_bytes());
                } else {
                    data.push(id.len() as u8);
                }
                data.extend_from_slice(id);
            }
            ends.push(data.len() as u64);
        }
        let mut w = Vec::with_capacity(8 + 8 * n + data.len());
        w.extend_from_slice(&(n as u64).to_le_bytes());
        for e in ends {
            w.extend_from_slice(&e.to_le_bytes());
        }
        w.extend_from_slice(&data);
        files.push(write(out, with_ext(base, "nos"), &w)?);
    }
    files.push(write(out, with_ext(base, "not"), &not_bytes(taxids))?);
    let (nto, offsets) = nto_bytes(taxids);
    files.push(write(out, with_ext(base, "nto"), &nto)?);
    let p = with_ext(base, "ntf");
    out.record(&p);
    lmdb::write_ntf(&p, &offsets)?;
    files.push(p);
    Ok(files)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn taxonomy_file_bytes() {
        // demo layout from makeblastdb (M4): three OIDs with taxid 0 only
        let sets = vec![vec![0], vec![0], vec![0]];
        let (nto, off) = nto_bytes(&sets);
        assert_eq!(off, vec![(0, 0)]);
        assert_eq!(nto.len() as u64, nto_size(&sets));
        assert_eq!(not_bytes(&sets).len() as u64, not_size(&sets));
        // mixed: OID 0 {0, 9606, 10090, 32630}, OID 1 {9606}, OID 2 {0}
        let sets = vec![vec![0, 9606, 10090, 32630], vec![9606], vec![0]];
        let not = not_bytes(&sets);
        let u64s = |b: &[u8]| -> Vec<u64> {
            b.chunks_exact(8)
                .map(|c| u64::from_le_bytes(c.try_into().unwrap()))
                .collect()
        };
        assert_eq!(u64s(&not[..32]), vec![3, 4, 5, 6]);
        let i32s: Vec<i32> = not[32..]
            .chunks_exact(4)
            .map(|c| i32::from_le_bytes(c.try_into().unwrap()))
            .collect();
        assert_eq!(i32s, vec![0, 9606, 10090, 32630, 9606, 0]);
        let (nto, off) = nto_bytes(&sets);
        let u32s: Vec<u32> = nto
            .chunks_exact(4)
            .map(|c| u32::from_le_bytes(c.try_into().unwrap()))
            .collect();
        // 0: [0, 2]; 9606: [0, 1]; 10090: [0]; 32630: [0]
        assert_eq!(u32s, vec![2, 0, 2, 2, 0, 1, 1, 0, 1, 0]);
        assert_eq!(off, vec![(0, 0), (9606, 12), (10090, 24), (32630, 32)]);
        assert_eq!(nto.len() as u64, nto_size(&sets));
        assert_eq!(distinct_taxids(&sets), vec![0, 9606, 10090, 32630]);
    }
}

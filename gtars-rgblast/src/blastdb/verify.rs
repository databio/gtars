//! The functional contract, without BLAST+: every listed id resolves to its OID, every
//! sequence read back has the store's md5, and the `.nin` totals equal `.njs`.

use super::output::with_ext;
use super::reader::DbReader;
use super::{njs, sidecars};
use crate::convert::nucleotide_bits;
use crate::error::Error;
use crate::pipeline::run_ordered;
use crate::progress::{ProgressEvent, Throttle};
use crate::store::decode;
use crate::store::fetch::SeqFetch;
use gtars_refget::digest::AlphabetType;
use gtars_refget::store::{ReadonlyRefgetStore, StorageMode};
use md5::{Digest, Md5};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;
use std::time::Instant;

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct VerifyReport {
    pub ids_checked: u64,
    pub ids_failed: u64,
    pub sequences_checked: u64,
    pub md5_mismatches: u64,
    /// `.nin` OID count and letters equal the `.njs` numbers.
    pub totals_ok: bool,
    /// Up to 20 failure descriptions.
    pub examples: Vec<String>,
}

impl VerifyReport {
    pub fn ok(&self) -> bool {
        self.ids_failed == 0 && self.md5_mismatches == 0 && self.totals_ok
    }

    fn example(&mut self, s: String) {
        if self.examples.len() < 20 {
            self.examples.push(s);
        }
    }
}

/// The acc2oid key blastdbcmd looks up for an id string.
pub fn lookup_key(id: &str) -> String {
    if let Some(rest) = id.strip_prefix("gnl|")
        && let Some((db, tag)) = rest.split_once('|')
    {
        return format!("{db}:{tag}");
    }
    id.strip_prefix("lcl|").unwrap_or(id).to_string()
}

fn md5_hex(b: &[u8]) -> String {
    let d = Md5::digest(b);
    d.iter().map(|x| format!("{x:02x}")).collect()
}

/// Check the database at `db` against `store`, reading store sequences through `fetch` only
/// for the slow path (U is stored as T in BLAST, so 3-bit and IUPAC mismatches are rechecked
/// against the store sequence with U replaced by T).
pub fn verify_db(
    db: &Path,
    store: &ReadonlyRefgetStore,
    fetch: &dyn SeqFetch,
    threads: usize,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<VerifyReport, Error> {
    let mut r = VerifyReport::default();
    let reader = DbReader::open(db)?;
    let ndb = reader.ndb()?;
    let ids = sidecars::read_ids(&with_ext(db, "ids.tsv"))?;

    // ids
    let mut oid_digest: HashMap<u32, String> = HashMap::new();
    for (id, oid, digest) in &ids {
        r.ids_checked += 1;
        let got = ndb.lookup(&lookup_key(id))?;
        if got != [*oid] {
            r.ids_failed += 1;
            r.example(format!("id {id}: expected OID {oid}, found {got:?}"));
        }
        if id.starts_with("SQ.") {
            oid_digest.insert(*oid, digest.clone());
        }
    }

    // totals
    let n = reader.num_oids();
    let nj = njs::read(&with_ext(db, "njs"));
    r.totals_ok = match &nj {
        Ok(j) => j.number_of_sequences == n as u64 && j.number_of_letters == reader.total_letters(),
        Err(_) => false,
    };
    if !r.totals_ok {
        r.example(format!(
            ".nin has {n} sequences and {} letters; .njs {}",
            reader.total_letters(),
            match &nj {
                Ok(j) => format!("has {} and {}", j.number_of_sequences, j.number_of_letters),
                Err(e) => format!("unreadable: {e}"),
            }
        ));
    }
    if oid_digest.len() as u32 != n {
        r.example(format!(
            "{} OIDs in the database, {} listed with an SQ. id",
            n,
            oid_digest.len()
        ));
        r.ids_failed += (n as u64).abs_diff(oid_digest.len() as u64);
    }

    // sequences
    let mode = store.storage_mode();
    let oids: Vec<u32> = (0..n).collect();
    let lens: Vec<u64> = oids
        .iter()
        .map(|&o| reader.length(o).map(|l| l as u64))
        .collect::<Result<_, _>>()?;
    let bases_total: u64 = lens.iter().sum();
    let check = |&oid: &u32| -> Result<Option<String>, Error> {
        let Some(d) = oid_digest.get(&oid) else {
            return Ok(Some(format!("OID {oid}: no SQ. id")));
        };
        let Some(meta) = store.get_sequence_metadata(d) else {
            return Ok(Some(format!("OID {oid}: SQ.{d} is not in the store")));
        };
        let seq = reader.sequence(oid)?;
        let got = md5_hex(&seq);
        if got == meta.md5 {
            return Ok(None);
        }
        if matches!(
            meta.alphabet,
            AlphabetType::Dna3bit | AlphabetType::DnaIupac
        ) {
            let len = meta.length as u64;
            let stored_len = match mode {
                StorageMode::Encoded => {
                    nucleotide_bits(meta.alphabet).map(|b| (len * b).div_ceil(8))
                }
                StorageMode::Raw => Some(len),
                StorageMode::Zstd => None,
            };
            let stored = fetch.read(d)?;
            let stored = &stored[..];
            if let Some(n) = stored_len
                && stored.len() as u64 != n
            {
                return Err(Error::sequence(
                    d,
                    format!("got {} stored bytes, expected {n}", stored.len()),
                ));
            }
            let mut ascii = decode::to_ascii(mode, meta.alphabet, stored, meta.length as u32)
                .map_err(|e| e.with_digest(d))?;
            for b in ascii.iter_mut() {
                if *b == b'U' {
                    *b = b'T';
                }
            }
            if md5_hex(&ascii) == got {
                return Ok(None);
            }
        }
        Ok(Some(format!(
            "OID {oid} (SQ.{d}): md5 {got}, store says {}",
            meta.md5
        )))
    };
    let t0 = Instant::now();
    let mut throttle = Throttle::new(2000);
    let (mut done, mut bases_done) = (0u64, 0u64);
    run_ordered(
        &oids,
        threads.max(1),
        1 << 30,
        |&o| lens[o as usize] + 64,
        check,
        |i, res| {
            done += 1;
            bases_done += lens[i];
            r.sequences_checked += 1;
            if let Some(msg) = res {
                r.md5_mismatches += 1;
                r.example(msg);
            }
            if throttle.ready()
                && !progress(&ProgressEvent::Sequences {
                    done,
                    total: n as u64,
                    bases_done,
                    bases_total,
                    elapsed_s: t0.elapsed().as_secs_f64(),
                })
            {
                return Err(Error::Cancelled);
            }
            Ok(())
        },
    )?;
    Ok(r)
}

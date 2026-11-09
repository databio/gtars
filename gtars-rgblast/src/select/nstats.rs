//! N counts: how many N (and other ambiguous) bases each sequence holds, computed once and
//! cached by digest. The N filter (`--max-n-fraction`, `--exclude-masked`) and the
//! `n_fraction` sidecar column use them.
//!
//! `n` counts N symbols (dna3bit code `100`, dnaio code `0000`, `N`/`n` in Raw and Zstd
//! stores). `ambiguous` counts every symbol other than A, C, G, T and U. dna2bit holds only
//! ACGT, so its counts are zero without a read.
//!
//! Cache: `<root>/blast/nstats/<first two digest chars>.tsv`, lines `digest length n ambiguous`,
//! appended with `O_APPEND` (one write per batch). Entries are content-addressed, so they never
//! go stale and are valid for any store. Readers tolerate duplicates and a partial last line;
//! an entry whose length disagrees with the store is ignored.

use crate::error::Error;
use crate::pipeline::run_ordered;
use crate::progress::{ProgressEvent, Throttle};
use crate::store::decode;
use crate::store::fetch::SeqFetch;
use crate::store::{OpenedStore, StoreFiles};
use gtars_refget::digest::{AlphabetType, DNA_3BIT_ALPHABET, DNA_IUPAC_ALPHABET};
use gtars_refget::store::StorageMode;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::io::Write;
use std::path::PathBuf;

/// Per-sequence counts.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct NCount {
    pub length: u64,
    pub n: u64,
    pub ambiguous: u64,
}

/// Where N counts come from besides the store.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default)]
pub struct NstatsOptions {
    /// Cache root (`--cache-dir`); None = the shared default ([`crate::cache_root::cache_root`]).
    /// Counts live in `<root>/blast/nstats/`.
    pub cache_dir: Option<PathBuf>,
    /// Append computed counts (and build by-products) to the cache.
    pub write_cache: bool,
}

impl Default for NstatsOptions {
    fn default() -> Self {
        NstatsOptions {
            cache_dir: None,
            write_cache: true,
        }
    }
}

impl NstatsOptions {
    /// The N-count cache (always known now that the root falls back to gtars-refget's default).
    pub fn cache(&self) -> Option<NstatsCache> {
        Some(NstatsCache {
            dir: crate::cache_root::blast_dir(self.cache_dir.as_deref()).join("nstats"),
        })
    }
}

const fn is_n(b: u8) -> bool {
    b == b'N'
}
const fn is_amb(b: u8) -> bool {
    !matches!(b, b'A' | b'C' | b'G' | b'T' | b'U')
}

/// 12-bit dna3bit half-group (4 codes) -> (N count, ambiguous count).
static T3: [(u8, u8); 4096] = {
    let mut t = [(0u8, 0u8); 4096];
    let mut x = 0;
    while x < 4096 {
        let (mut n, mut a) = (0u8, 0u8);
        let mut k = 0;
        while k < 4 {
            let c = (x >> (9 - 3 * k)) & 7;
            let l = DNA_3BIT_ALPHABET.decoding_array[c];
            n += is_n(l) as u8;
            a += is_amb(l) as u8;
            k += 1;
        }
        t[x] = (n, a);
        x += 1;
    }
    t
};

/// dnaio byte (2 codes) -> (N count, ambiguous count).
static T4: [(u8, u8); 256] = {
    let mut t = [(0u8, 0u8); 256];
    let mut b = 0;
    while b < 256 {
        let (l1, l2) = (
            DNA_IUPAC_ALPHABET.decoding_array[b >> 4],
            DNA_IUPAC_ALPHABET.decoding_array[b & 15],
        );
        t[b] = (
            is_n(l1) as u8 + is_n(l2) as u8,
            is_amb(l1) as u8 + is_amb(l2) as u8,
        );
        b += 1;
    }
    t
};

/// True when counting needs the sequence bytes (every nucleotide alphabet but dna2bit).
pub fn needs_read(alphabet: AlphabetType) -> bool {
    matches!(alphabet, AlphabetType::Dna3bit | AlphabetType::DnaIupac)
}

/// Count one stored record. dna2bit needs no bytes (pass an empty slice).
pub fn count(
    mode: StorageMode,
    alphabet: AlphabetType,
    stored: &[u8],
    len: u64,
) -> Result<NCount, Error> {
    let mut c = NCount {
        length: len,
        ..Default::default()
    };
    if alphabet == AlphabetType::Dna2bit {
        return Ok(c);
    }
    match (mode, alphabet) {
        (StorageMode::Encoded, AlphabetType::Dna3bit) => {
            let need = (len * 3).div_ceil(8);
            check_len(stored.len() as u64, need, len, alphabet)?;
            let groups = (len / 8) as usize;
            for g in stored[..3 * groups].chunks_exact(3) {
                let w = (g[0] as usize) << 16 | (g[1] as usize) << 8 | g[2] as usize;
                let (a, b) = (T3[w >> 12], T3[w & 0xFFF]);
                c.n += (a.0 + b.0) as u64;
                c.ambiguous += (a.1 + b.1) as u64;
            }
            for p in 8 * groups as u64..len {
                let bit = (3 * p) as usize;
                let i = bit / 8;
                let hi = stored[i] as u32;
                let lo = *stored.get(i + 1).unwrap_or(&0) as u32;
                let code = (((hi << 8) | lo) >> (13 - (bit % 8))) & 7;
                let l = DNA_3BIT_ALPHABET.decoding_array[code as usize];
                c.n += is_n(l) as u64;
                c.ambiguous += is_amb(l) as u64;
            }
        }
        (StorageMode::Encoded, AlphabetType::DnaIupac) => {
            let need = len.div_ceil(2);
            check_len(stored.len() as u64, need, len, alphabet)?;
            let full = (len / 2) as usize;
            for &b in &stored[..full] {
                let (n, a) = T4[b as usize];
                c.n += n as u64;
                c.ambiguous += a as u64;
            }
            if len % 2 == 1 {
                let l = DNA_IUPAC_ALPHABET.decoding_array[(stored[full] >> 4) as usize];
                c.n += is_n(l) as u64;
                c.ambiguous += is_amb(l) as u64;
            }
        }
        _ => {
            let ascii = decode::to_ascii(mode, alphabet, stored, len as u32)?;
            for &b in &ascii {
                c.n += is_n(b) as u64;
                c.ambiguous += is_amb(b) as u64;
            }
        }
    }
    Ok(c)
}

fn check_len(got: u64, need: u64, len: u64, alphabet: AlphabetType) -> Result<(), Error> {
    if got != need {
        return Err(Error::sequence(
            "",
            format!("stored record is {got} bytes, expected {need} for {len} bp of {alphabet}"),
        ));
    }
    Ok(())
}

/// The on-disk cache under `<root>/nstats/`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NstatsCache {
    pub dir: PathBuf,
}

fn shard(d: &str) -> String {
    d.chars().take(2).collect()
}

/// Parse lines in the cache format; malformed lines (a partial last line) are skipped.
pub fn parse_lines(text: &str, out: &mut HashMap<String, Vec<NCount>>) {
    for l in text.lines() {
        if l.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() != 4 {
            continue;
        }
        let (Ok(length), Ok(n), Ok(ambiguous)) = (f[1].parse(), f[2].parse(), f[3].parse()) else {
            continue;
        };
        out.entry(f[0].to_string()).or_default().push(NCount {
            length,
            n,
            ambiguous,
        });
    }
}

/// The first entry of `d` whose length matches.
fn pick(map: &HashMap<String, Vec<NCount>>, d: &str, len: u64) -> Option<NCount> {
    map.get(d)?.iter().find(|c| c.length == len).copied()
}

impl NstatsCache {
    pub fn path(&self, digest: &str) -> PathBuf {
        self.dir.join(format!("{}.tsv", shard(digest)))
    }

    /// Every cached entry of these digests (all entries; callers check lengths).
    pub fn read(&self, digests: &[&str]) -> HashMap<String, Vec<NCount>> {
        let mut shards: BTreeMap<String, ()> = BTreeMap::new();
        for d in digests {
            shards.insert(shard(d), ());
        }
        let mut all = HashMap::new();
        for s in shards.keys() {
            if let Ok(t) = std::fs::read_to_string(self.dir.join(format!("{s}.tsv"))) {
                parse_lines(&t, &mut all);
            }
        }
        let want: std::collections::HashSet<&str> = digests.iter().copied().collect();
        all.retain(|k, _| want.contains(k.as_str()));
        all
    }

    /// Append entries, one `O_APPEND` write per shard file.
    pub fn append(&self, entries: &[(String, NCount)]) -> Result<(), Error> {
        if entries.is_empty() {
            return Ok(());
        }
        std::fs::create_dir_all(&self.dir).map_err(|e| Error::io(&self.dir, e))?;
        let mut by: BTreeMap<String, String> = BTreeMap::new();
        for (d, c) in entries {
            let b = by.entry(shard(d)).or_default();
            b.push_str(&format!("{d}\t{}\t{}\t{}\n", c.length, c.n, c.ambiguous));
        }
        for (s, text) in by {
            let p = self.dir.join(format!("{s}.tsv"));
            let mut f = std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(&p)
                .map_err(|e| Error::io(&p, e))?;
            f.write_all(text.as_bytes()).map_err(|e| Error::io(&p, e))?;
        }
        Ok(())
    }
}

/// N counts of these sequence digests (the sequence index must be loaded). Filled from the
/// cache, then a scan of the store (announced and reported through
/// `ProgressEvent::NScan`). Non-nucleotide sequences are left out of the result.
pub fn ensure(
    store: &OpenedStore,
    digests: &[String],
    opts: &NstatsOptions,
    threads: usize,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<HashMap<String, NCount>, Error> {
    let ro = &store.store;
    let mut out: HashMap<String, NCount> = HashMap::new();
    let mut need: Vec<(String, AlphabetType, u64)> = Vec::new();
    for d in digests {
        if out.contains_key(d) {
            continue;
        }
        let m = ro.get_sequence_metadata(d).ok_or_else(|| {
            Error::Store(anyhow::anyhow!("sequence {d} is not in sequences.rgsi"))
        })?;
        let len = m.length as u64;
        match m.alphabet {
            AlphabetType::Dna2bit => {
                out.insert(
                    d.clone(),
                    NCount {
                        length: len,
                        n: 0,
                        ambiguous: 0,
                    },
                );
            }
            a if needs_read(a) => need.push((d.clone(), a, len)),
            _ => {}
        }
    }
    need.sort_by(|a, b| a.0.cmp(&b.0));
    need.dedup_by(|a, b| a.0 == b.0);
    if need.is_empty() {
        return Ok(out);
    }
    let cache = opts.cache();
    if let Some(c) = &cache {
        let ds: Vec<&str> = need.iter().map(|x| x.0.as_str()).collect();
        let got = c.read(&ds);
        need.retain(|(d, _, len)| match pick(&got, d, *len) {
            Some(x) => {
                out.insert(d.clone(), x);
                false
            }
            None => true,
        });
    }
    if need.is_empty() {
        return Ok(out);
    }
    let digests: Vec<String> = need.iter().map(|(d, _, _)| d.clone()).collect();
    store.pull(&digests)?;
    let fetch = StoreFiles::new(store);
    let scanned = scan(store, &fetch, &need, threads, progress)?;
    if opts.write_cache
        && let Some(c) = &cache
    {
        c.append(&scanned)?;
    }
    out.extend(scanned);
    Ok(out)
}

/// Read and count each sequence on `threads` workers; returns (digest, counts) in input order.
fn scan(
    store: &OpenedStore,
    fetch: &dyn SeqFetch,
    need: &[(String, AlphabetType, u64)],
    threads: usize,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<Vec<(String, NCount)>, Error> {
    let mode = store.mode();
    let stored_len = |a: AlphabetType, len: u64| -> Option<u64> {
        match mode {
            StorageMode::Encoded => {
                crate::convert::nucleotide_bits(a).map(|b| (len * b).div_ceil(8))
            }
            StorageMode::Raw => Some(len),
            StorageMode::Zstd => None,
        }
    };
    let bytes_total: u64 = need
        .iter()
        .map(|(_, a, l)| stored_len(*a, *l).unwrap_or(*l))
        .sum();
    let total = need.len() as u64;
    let ev = |done, bytes_read| ProgressEvent::NScan {
        done,
        total,
        bytes_read,
        bytes_total,
    };
    if !progress(&ev(0, 0)) {
        return Err(Error::Cancelled);
    }
    let mut out = Vec::with_capacity(need.len());
    let mut bytes_read = 0u64;
    let mut throttle = Throttle::new(1000);
    let mut cancelled = false;
    run_ordered(
        need,
        threads,
        1 << 30,
        |(_, a, l)| stored_len(*a, *l).unwrap_or(*l) + 64,
        |(d, a, l)| {
            let stored = fetch.read(d).map_err(|e| e.with_digest(d))?;
            let stored = &stored[..];
            if let Some(n) = stored_len(*a, *l)
                && stored.len() as u64 != n
            {
                return Err(Error::sequence(
                    d,
                    format!("got {} stored bytes, expected {n}", stored.len()),
                ));
            }
            let c = count(mode, *a, stored, *l).map_err(|e| e.with_digest(d))?;
            Ok((d.clone(), c, stored.len() as u64))
        },
        |_, (d, c, n)| {
            bytes_read += n;
            out.push((d, c));
            if throttle.ready() && !progress(&ev(out.len() as u64, bytes_read)) {
                cancelled = true;
                return Err(Error::Cancelled);
            }
            Ok(())
        },
    )?;
    if cancelled || !progress(&ev(out.len() as u64, bytes_read)) {
        return Err(Error::Cancelled);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use gtars_refget::digest::{encode_sequence, lookup_alphabet};

    struct Rng(u64);
    impl Rng {
        fn below(&mut self, n: u64) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0 % n
        }
    }

    fn naive(ascii: &[u8]) -> NCount {
        NCount {
            length: ascii.len() as u64,
            n: ascii.iter().filter(|&&b| b == b'N' || b == b'n').count() as u64,
            ambiguous: ascii
                .iter()
                .filter(|&&b| !b"ACGTUacgtu".contains(&b))
                .count() as u64,
        }
    }

    #[test]
    fn counts_match_naive_decode() {
        let mut rng = Rng(0xfeed_beef_1234_5678);
        let alphabets: [(AlphabetType, &[u8]); 3] = [
            (AlphabetType::Dna2bit, b"ACGT"),
            (AlphabetType::Dna3bit, b"ACGTNRYU"),
            (AlphabetType::DnaIupac, b"ACGTNMRSKWYUBDHV"),
        ];
        let mut lens: Vec<usize> = (0..=17).collect();
        lens.extend([100, 1001, 4099, 20_000]);
        for (at, symbols) in alphabets {
            for &len in &lens {
                for _ in 0..4 {
                    let mut s = Vec::with_capacity(len);
                    while s.len() < len {
                        if symbols.contains(&b'N') && rng.below(6) == 0 {
                            let run = rng.below(30) as usize + 1;
                            s.extend(std::iter::repeat_n(b'N', run));
                        } else {
                            s.push(symbols[rng.below(symbols.len() as u64) as usize]);
                        }
                    }
                    s.truncate(len);
                    let want = naive(&s);
                    let packed = encode_sequence(&s, lookup_alphabet(&at)).unwrap();
                    let got = count(StorageMode::Encoded, at, &packed, len as u64).unwrap();
                    assert_eq!(got, want, "{at} encoded len {len}");
                    let got = count(StorageMode::Raw, at, &s, len as u64).unwrap();
                    assert_eq!(got, want, "{at} raw len {len}");
                    let lower: Vec<u8> = s.iter().map(|b| b.to_ascii_lowercase()).collect();
                    assert_eq!(
                        count(StorageMode::Raw, at, &lower, len as u64).unwrap(),
                        want
                    );
                    let z = zstd::encode_all(&s[..], 3).unwrap();
                    assert_eq!(count(StorageMode::Zstd, at, &z, len as u64).unwrap(), want);
                }
            }
        }
        // the converter's by-product counts agree
        let s = b"ACGTNNNNRYUacgtnK".to_vec();
        let (mut q, mut a) = (Vec::new(), Vec::new());
        let bc = crate::convert::ascii_to_blast(&s, &mut q, &mut a).unwrap();
        let c = naive(&s);
        assert_eq!((bc.n, bc.ambiguous), (c.n, c.ambiguous));
    }

    #[test]
    fn wrong_record_size_is_an_error() {
        assert!(count(StorageMode::Encoded, AlphabetType::Dna3bit, &[0; 2], 8).is_err());
        assert!(count(StorageMode::Encoded, AlphabetType::DnaIupac, &[0; 2], 5).is_err());
    }

    #[test]
    fn cache_roundtrip_mismatch_and_partial_lines() {
        let d = tempfile::tempdir().unwrap();
        let c = NstatsCache {
            dir: d.path().join("nstats"),
        };
        let e = |l, n| NCount {
            length: l,
            n,
            ambiguous: n + 1,
        };
        c.append(&[
            ("abX".into(), e(10, 2)),
            ("abY".into(), e(20, 3)),
            ("cd1".into(), e(5, 0)),
        ])
        .unwrap();
        // a duplicate and a partial last line
        c.append(&[("abX".into(), e(10, 2))]).unwrap();
        let mut f = std::fs::OpenOptions::new()
            .append(true)
            .open(c.path("abZ"))
            .unwrap();
        f.write_all(b"abZ\t7\t1").unwrap();
        let got = c.read(&["abX", "abY", "abZ", "cd1", "zz9"]);
        assert_eq!(pick(&got, "abX", 10), Some(e(10, 2)));
        assert_eq!(pick(&got, "abY", 20), Some(e(20, 3)));
        assert_eq!(pick(&got, "cd1", 5), Some(e(5, 0)));
        assert_eq!(pick(&got, "abZ", 7), None, "partial line ignored");
        assert_eq!(pick(&got, "zz9", 1), None);
        // a length that disagrees with the store is ignored
        assert_eq!(pick(&got, "abX", 11), None);
    }

    #[test]
    fn cache_concurrent_appends() {
        let d = tempfile::tempdir().unwrap();
        let c = NstatsCache {
            dir: d.path().join("nstats"),
        };
        std::thread::scope(|s| {
            for t in 0..8u64 {
                let c = &c;
                s.spawn(move || {
                    for i in 0..50u64 {
                        let batch: Vec<(String, NCount)> = (0..20)
                            .map(|k| {
                                (
                                    format!("aa{t}_{i}_{k}"),
                                    NCount {
                                        length: t * 1000 + i * 20 + k,
                                        n: k,
                                        ambiguous: k,
                                    },
                                )
                            })
                            .collect();
                        c.append(&batch).unwrap();
                    }
                });
            }
        });
        let text = std::fs::read_to_string(c.path("aa")).unwrap();
        assert_eq!(text.lines().count(), 8 * 50 * 20);
        let mut m = HashMap::new();
        parse_lines(&text, &mut m);
        assert_eq!(m.len(), 8 * 50 * 20, "every line intact");
    }
}

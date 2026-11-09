//! Preflight: everything known before any write (all that `--dry-run` does). Exact sizes for
//! the parts that follow from the plan, an allowance for the parts that need the sequences.

use super::asn;
use super::lmdb::ndb_estimate;
use super::oidfiles::{nos_id_size, nos_size, not_size, nto_size};
use super::writer::{fits, nin_overhead, nin_size};
use crate::error::Error;
use crate::plan::{SkipReason, SubjectPlan};
use gtars_refget::digest::AlphabetType;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

/// Free space kept in reserve on top of the predicted size.
pub const SPACE_MARGIN: u64 = 64 << 20;

/// Per-subject data computed once and reused by the writer. Headers and acc2oid keys are not
/// held: preflight takes their exact sizes from [`asn::encoded_len`] and
/// [`SeqId::lmdb_key_lens`](crate::plan::SeqId::lmdb_key_lens), and the writer encodes each
/// subject's header and keys when it writes it.
pub struct Prepared {
    /// Sorted unique taxids of each subject's deflines (`.not` / `.nto`).
    pub taxids: Vec<Vec<i32>>,
}

pub fn prepare(plan: &SubjectPlan) -> Prepared {
    Prepared {
        taxids: super::oidfiles::taxid_sets(&plan.subjects),
    }
}

/// Sizes that follow exactly from the plan.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct ExactSizes {
    /// `.nsq` without ambiguity tables: Σ(len/4 + 1) + 1 per volume.
    pub nsq_without_ambiguity: u64,
    pub nhr: u64,
    pub nin: u64,
    pub nos: u64,
    pub not: u64,
    pub nto: u64,
}

impl ExactSizes {
    pub fn total(&self) -> u64 {
        self.nsq_without_ambiguity + self.nhr + self.nin + self.nos + self.not + self.nto
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct PreflightSummary {
    pub store: String,
    pub store_label: String,
    pub mode: String,
    pub remote: bool,
    pub collections: u64,
    pub names: u64,
    pub subjects: u64,
    pub deflines: u64,
    pub ids: u64,
    pub total_bases: u64,
    pub max_len: u64,
    pub skipped_oversize: u64,
    pub skipped_non_nucleotide: u64,
    pub skipped_empty: u64,
    pub sanitized: u64,
    pub dash_leading_collections: u64,
    pub ambiguous_aliases: u64,
    pub orphans: u64,
    pub exact: ExactSizes,
    pub exact_bytes: u64,
    /// Ambiguity tables of 3-bit and IUPAC sequences and the LMDB page layout: not knowable
    /// without reading every sequence, so an allowance, never exact.
    pub allowance_bytes: u64,
    /// Volumes predicted from the exact parts (ambiguity tables can add one).
    pub predicted_volumes: u64,
    /// Free bytes in the output directory; None when unknown or not checked.
    pub free_bytes: Option<u64>,
    pub out: PathBuf,
    /// The resolved selection in one line (`ResolvedSelection::summary`); empty when unknown.
    #[serde(default)]
    pub selection: String,
}

pub struct PreflightInput<'a> {
    pub store: String,
    pub store_label: String,
    pub remote: bool,
    pub out: &'a Path,
    pub title: &'a str,
    pub date: &'a str,
    pub max_file_size: u64,
}

pub fn compute(plan: &SubjectPlan, prep: &Prepared, i: &PreflightInput) -> PreflightSummary {
    let lmdb_name = format!("{}.ndb", super::writer::file_name(i.out));
    let ov = nin_overhead(i.title, &lmdb_name, i.date.len()).div_ceil(8) * 8;
    // volume split from the exact parts
    let mut vols: Vec<u64> = Vec::new();
    let (mut n, mut hdr, mut seq) = (0u64, 0u64, 1u64);
    let mut nsq = 0u64;
    let mut nhr = 0u64;
    let mut amb_allow = 0u64;
    // acc2oid keys (every key, as the writer adds them) and the `.nos` ids
    let (mut key_bytes, mut nos_ids) = (0u64, 0u64);
    for s in &plan.subjects {
        let h = asn::encoded_len(&s.deflines) as u64;
        let data = s.length as u64 / 4 + 1;
        if n > 0 && !fits(i.max_file_size, ov, n, hdr, seq, h, data) {
            vols.push(n);
            (n, hdr, seq) = (0, 0, 1);
        }
        n += 1;
        hdr += h;
        seq += data;
        nsq += data;
        nhr += h;
        if s.alphabet != AlphabetType::Dna2bit {
            amb_allow += data / 20 + 4;
        }
        for id in s.deflines.iter().flat_map(|d| &d.ids) {
            id.lmdb_key_lens(|len, save| {
                key_bytes += len as u64 + 24;
                if save {
                    nos_ids += nos_id_size(len);
                }
            });
        }
    }
    vols.push(n);
    let nvol = vols.len() as u64;
    let nin: u64 = vols
        .iter()
        .map(|&n| nin_size(i.title, &lmdb_name, i.date, n))
        .sum();
    let n_oids = plan.subjects.len() as u64;
    let exact = ExactSizes {
        nsq_without_ambiguity: nsq + nvol,
        nhr,
        nin,
        nos: if n_oids > 0 {
            nos_size(n_oids, std::iter::empty()) + nos_ids
        } else {
            0
        },
        not: not_size(&prep.taxids),
        nto: nto_size(&prep.taxids),
    };
    let lmdb_allow = 2 * ndb_estimate(key_bytes, vols.len()) + (48 << 10);
    let exact_bytes = exact.total();
    PreflightSummary {
        store: i.store.clone(),
        store_label: i.store_label.clone(),
        mode: format!("{:?}", plan.store_mode),
        remote: i.remote,
        collections: plan.collections.len() as u64,
        names: plan.n_names(),
        subjects: n_oids,
        deflines: plan.n_deflines(),
        ids: plan.n_lookup_ids(),
        total_bases: plan.totals.total_length,
        max_len: plan.totals.max_len as u64,
        skipped_oversize: plan.n_skipped(SkipReason::Oversize),
        skipped_non_nucleotide: plan.n_skipped(SkipReason::NonNucleotide),
        skipped_empty: plan.n_skipped(SkipReason::Empty),
        sanitized: plan.sanitized.len() as u64,
        dash_leading_collections: plan.n_dash_leading_collections(),
        ambiguous_aliases: plan.ambiguous_aliases.len() as u64,
        orphans: plan.n_orphans(),
        exact,
        exact_bytes,
        allowance_bytes: amb_allow + lmdb_allow,
        predicted_volumes: nvol,
        free_bytes: None,
        out: i.out.to_path_buf(),
        selection: String::new(),
    }
}

/// Free bytes for unprivileged writes in `dir` (or its nearest existing ancestor).
#[cfg(unix)]
pub fn free_space(dir: &Path) -> Option<u64> {
    use std::os::unix::ffi::OsStrExt;
    let mut d = if dir.as_os_str().is_empty() {
        Path::new(".")
    } else {
        dir
    };
    while !d.exists() {
        d = d
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        if d == Path::new(".") {
            break;
        }
    }
    let c = std::ffi::CString::new(d.as_os_str().as_bytes()).ok()?;
    let mut st: libc::statvfs = unsafe { std::mem::zeroed() };
    if unsafe { libc::statvfs(c.as_ptr(), &mut st) } != 0 {
        return None;
    }
    Some(st.f_bavail as u64 * st.f_frsize as u64)
}

#[cfg(not(unix))]
pub fn free_space(_dir: &Path) -> Option<u64> {
    None
}

/// Fill in free space and fail if exact + allowance + margin does not fit.
pub fn check_space(s: &mut PreflightSummary, check: bool) -> Result<(), Error> {
    if !check {
        return Ok(());
    }
    let dir = s.out.parent().map(Path::to_path_buf).unwrap_or_default();
    s.free_bytes = free_space(&dir);
    if let Some(avail) = s.free_bytes {
        let need = s.exact_bytes + s.allowance_bytes + SPACE_MARGIN;
        if need > avail {
            let dir = if dir.as_os_str().is_empty() {
                PathBuf::from(".")
            } else {
                dir
            };
            return Err(Error::InsufficientSpace {
                dir,
                need,
                available: avail,
            });
        }
    }
    Ok(())
}

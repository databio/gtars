//! Tabular output (`-outfmt 6`), computed from the engine's HSPs.
//!
//! Source: BLAST+ 2.17.0 `objtools/align_format/tabular.cpp` (`x_SetFieldsToShow`,
//! `SetFields`, `SetScores`, `x_PrintField`, `x_SetQueryCovSeqalign` and the inline printers in `tabular.hpp`:
//! `pident = DoubleToString(num_ident / align_length * 100, 3)`, `mismatch = length -
//! num_ident - gaps`), `algo/blast/format/blastfmtutil.cpp::InsertSubjectScores` (`qcovs`,
//! `qcovus`), `align_format_util.cpp` (`GetScoreString`, `GetMasterCoverage`,
//! `GetUniqSeqCoverage`, `x_AcknowledgeBlastSequence`), `showdefline.cpp::GetSeqIdListString`
//! and `corelib/ncbistr.cpp::NStr::DoubleToString` (`%.*f` / `%.*e`).
//!
//! Coordinates come out of the engine already in blastn's printed orientation
//! (`engine::search::Hsp`). Ids print as `tabular.cpp` prints them from the deflines: a local
//! id bare, a general id as `gnl|db|tag` (`sseqid`) or `db:tag` (`saccver`). Each defline's
//! ids are all of one kind, so `FindBestChoice(CSeq_id::Score)` is its first id.
//!
//! `qseq`/`sseq` (and the fields derived from them: `btop`, `positive`, `ppos`, and the
//! recount of `nident`/`pident`/`mismatch` that blastn does when any of those is requested)
//! need the alignment strings filled in by `format::alignment::fill`.

use crate::engine::search::{Hit, Hsp, QueryResult, SearchResults};
use crate::error::Error;
use crate::plan::{Defline, SeqId, SubjectPlan};
use std::io::Write;

/// Every nucleotide field blastn 2.17 prints, in `sc_FormatSpecifiers` order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Field {
    Qseqid,
    Qgi,
    Qacc,
    Qaccver,
    Qlen,
    Sseqid,
    Sallseqid,
    Sgi,
    Sallgi,
    Sacc,
    Saccver,
    Sallacc,
    Slen,
    Qstart,
    Qend,
    Sstart,
    Send,
    Qseq,
    Sseq,
    Evalue,
    Bitscore,
    Score,
    Length,
    Pident,
    Nident,
    Mismatch,
    Positive,
    Gapopen,
    Gaps,
    Ppos,
    Frames,
    Qframe,
    Sframe,
    Btop,
    Staxids,
    Staxid,
    Stitle,
    Salltitles,
    Sstrand,
    Qcovs,
    Qcovhsp,
    Qcovus,
}

/// (short name, field).
const FIELDS: [(&str, Field); 42] = [
    ("qseqid", Field::Qseqid),
    ("qgi", Field::Qgi),
    ("qacc", Field::Qacc),
    ("qaccver", Field::Qaccver),
    ("qlen", Field::Qlen),
    ("sseqid", Field::Sseqid),
    ("sallseqid", Field::Sallseqid),
    ("sgi", Field::Sgi),
    ("sallgi", Field::Sallgi),
    ("sacc", Field::Sacc),
    ("saccver", Field::Saccver),
    ("sallacc", Field::Sallacc),
    ("slen", Field::Slen),
    ("qstart", Field::Qstart),
    ("qend", Field::Qend),
    ("sstart", Field::Sstart),
    ("send", Field::Send),
    ("qseq", Field::Qseq),
    ("sseq", Field::Sseq),
    ("evalue", Field::Evalue),
    ("bitscore", Field::Bitscore),
    ("score", Field::Score),
    ("length", Field::Length),
    ("pident", Field::Pident),
    ("nident", Field::Nident),
    ("mismatch", Field::Mismatch),
    ("positive", Field::Positive),
    ("gapopen", Field::Gapopen),
    ("gaps", Field::Gaps),
    ("ppos", Field::Ppos),
    ("frames", Field::Frames),
    ("qframe", Field::Qframe),
    ("sframe", Field::Sframe),
    ("btop", Field::Btop),
    ("staxids", Field::Staxids),
    ("staxid", Field::Staxid),
    ("stitle", Field::Stitle),
    ("salltitles", Field::Salltitles),
    ("sstrand", Field::Sstrand),
    ("qcovs", Field::Qcovs),
    ("qcovhsp", Field::Qcovhsp),
    ("qcovus", Field::Qcovus),
];

/// Fields that need NCBI's `taxdb` (blastn prints `N/A` and warns without it); refused so a
/// user with taxdb installed never gets a silently different column.
const TAXDB_FIELDS: [&str; 8] = [
    "ssciname",
    "scomname",
    "sblastname",
    "sskingdom",
    "sscinames",
    "scomnames",
    "sblastnames",
    "sskingdoms",
];

impl Field {
    pub fn parse(s: &str) -> Option<Field> {
        FIELDS.iter().find(|(n, _)| *n == s).map(|(_, f)| *f)
    }

    /// The `-outfmt` token.
    pub fn name(self) -> &'static str {
        FIELDS.iter().find(|(_, f)| *f == self).unwrap().0
    }
}

/// `std` (`kDfltArgTabularOutputFmt`): note `qaccver`/`saccver`, not `qseqid`/`sseqid`; they
/// print the same string for a local id but differ for a general one.
pub const STD_FIELDS: [Field; 12] = [
    Field::Qaccver,
    Field::Saccver,
    Field::Pident,
    Field::Length,
    Field::Mismatch,
    Field::Gapopen,
    Field::Qstart,
    Field::Qend,
    Field::Sstart,
    Field::Send,
    Field::Evalue,
    Field::Bitscore,
];

#[derive(Debug, Clone, PartialEq)]
pub struct TabularSpec {
    pub fields: Vec<Field>,
    /// Tokens blastn ignores silently (`x_SetFieldsToShow` skips names it does not know);
    /// the CLI warns about them on stderr.
    pub ignored: Vec<String>,
}

impl TabularSpec {
    /// blastn's `-outfmt` grammar for format 6: `6`, `"6 std stitle"`. Other formats
    /// (including the tabular 7 and 10) are refused.
    pub fn parse_outfmt(s: &str) -> Result<TabularSpec, Error> {
        let mut it = s.split_whitespace();
        let code = it
            .next()
            .ok_or_else(|| Error::InvalidOptions("empty -outfmt".into()))?;
        if code != "6" {
            return Err(Error::unsupported(&format!(
                "-outfmt {code} (rgblastn writes -outfmt 6)"
            )));
        }
        let mut fields = Vec::new();
        let mut ignored = Vec::new();
        for tok in it {
            if tok == "std" {
                fields.extend_from_slice(&STD_FIELDS);
            } else if let Some(f) = Field::parse(tok) {
                fields.push(f);
            } else if TAXDB_FIELDS.contains(&tok) {
                return Err(Error::unsupported(&format!(
                    "-outfmt field {tok} (taxonomy names need NCBI's taxdb, which rgblastn does not read; use staxid or staxids)"
                )));
            } else {
                ignored.push(tok.to_string());
            }
        }
        if fields.is_empty() {
            fields.extend_from_slice(&STD_FIELDS);
        }
        Ok(TabularSpec { fields, ignored })
    }

    pub fn has(&self, f: Field) -> bool {
        self.fields.contains(&f)
    }

    /// `qseq`/`sseq` must be filled (`format::alignment::fill`): blastn builds the alignment
    /// strings for these fields and then recounts identities from them.
    pub fn needs_sequences(&self) -> bool {
        [
            Field::Qseq,
            Field::Sseq,
            Field::Btop,
            Field::Positive,
            Field::Ppos,
        ]
        .iter()
        .any(|f| self.has(*f))
    }
}

/// C `%.{prec}e`: Rust's exponent lacks the sign and zero padding.
fn c_sci(x: f64, prec: usize) -> String {
    let s = format!("{:.*e}", prec, x);
    let (mant, exp) = s.split_once('e').unwrap();
    let (sign, digits) = if let Some(d) = exp.strip_prefix('-') {
        ('-', d)
    } else {
        ('+', exp)
    };
    if digits.len() < 2 {
        format!("{mant}e{sign}0{digits}")
    } else {
        format!("{mant}e{sign}{digits}")
    }
}

/// `GetScoreString`'s e-value plus `SetScores`' override of the `[1e-180, 0.0009)` band.
pub fn evalue_string(e: f64) -> String {
    // SetScores overrides GetScoreString's "%2.0le" (< 1e-99) and "%3.0le" bands with
    // NStr::DoubleToString(evalue, 2, fDoubleScientific) for the whole [1e-180, 0.0009) range.
    if e < 1.0e-180 {
        "0.0".into()
    } else if e < 0.0009 {
        c_sci(e, 2)
    } else if e < 0.1 {
        format!("{:.3}", e)
    } else if e < 1.0 {
        format!("{:.2}", e)
    } else if e < 10.0 {
        format!("{:.1}", e)
    } else {
        format!("{:.0}", e)
    }
}

/// `GetScoreString`'s bit score.
pub fn bitscore_string(b: f64) -> String {
    if b > 99999.0 {
        c_sci(b, 3)
    } else if b > 99.9 {
        format!("{:>3}", b as i64)
    } else {
        format!("{:>4.1}", b)
    }
}

/// `GetSeqIdListString(ids, eFullId)`: a local id bare, a general id as FASTA.
fn seqid_string(id: &SeqId) -> String {
    match id {
        SeqId::Local(s) => s.clone(),
        SeqId::General { db, tag } => format!("gnl|{db}|{tag}"),
    }
}

/// `GetLabel(eContent)` of the best id: a local id bare, a general id as `db:tag`.
fn acc_string(id: &SeqId) -> String {
    match id {
        SeqId::Local(s) => s.clone(),
        SeqId::General { db, tag } => format!("{db}:{tag}"),
    }
}

fn first_id(d: &Defline, f: fn(&SeqId) -> String) -> String {
    d.ids.first().map(f).unwrap_or_else(|| "Unknown".into())
}

fn per_defline(deflines: &[Defline], f: fn(&SeqId) -> String) -> String {
    deflines
        .iter()
        .map(|d| first_id(d, f))
        .collect::<Vec<_>>()
        .join(";")
}

fn title_or_na(t: &str) -> &str {
    if t.is_empty() { "N/A" } else { t }
}

/// `x_SetQueryCovSeqalign` / `InsertSubjectScores`: percent with blastn's rounding
/// (`+0.5` below 99, then truncation).
fn coverage_percent(covered: i64, query_len: i64) -> i64 {
    let mut tmp = 100.0 * covered as f64 / query_len as f64;
    if tmp < 99.0 {
        tmp += 0.5;
    }
    tmp as i64
}

/// `GetMasterCoverage`: the merged query ranges of a subject's HSPs (overlap merges,
/// adjacency does not).
fn master_coverage(hsps: &[Hsp]) -> i64 {
    let mut ranges: Vec<(i64, i64)> = hsps
        .iter()
        .map(|h| {
            let (a, b) = (h.q_start as i64, h.q_end as i64);
            if a > b { (b, a) } else { (a, b) }
        })
        .collect();
    ranges.sort_by_key(|r| r.0);
    let mut merged: Vec<(i64, i64)> = Vec::new();
    for r in ranges {
        match merged.last_mut() {
            Some(prev) if prev.0 <= r.1 && r.0 <= prev.1 => {
                prev.0 = prev.0.min(r.0);
                prev.1 = prev.1.max(r.1);
            }
            _ => merged.push(r),
        }
    }
    merged.iter().map(|(a, b)| b - a + 1).sum()
}

/// `CRangeCollection<TSeqPos>`: disjoint sorted closed ranges; `add` merges overlapping and
/// abutting ranges, `subtract` removes a set of ranges.
#[derive(Default, Clone)]
struct RangeColl(Vec<(u32, u32)>);

impl RangeColl {
    fn add(&mut self, r: (u32, u32)) {
        let (mut from, mut to) = r;
        let mut out = Vec::with_capacity(self.0.len() + 1);
        for &(a, b) in &self.0 {
            // x_CombineWith: overlapping or abutting (GetToOpen() >= from) ranges merge
            let touches = a as u64 <= to as u64 + 1 && b as u64 + 1 >= from as u64;
            if touches {
                from = from.min(a);
                to = to.max(b);
            } else {
                out.push((a, b));
            }
        }
        out.push((from, to));
        out.sort_unstable();
        self.0 = out;
    }

    fn subtract(&mut self, other: &RangeColl) {
        for &(a, b) in &other.0 {
            let mut out = Vec::with_capacity(self.0.len() + 1);
            for &(x, y) in &self.0 {
                if y < a || x > b {
                    out.push((x, y));
                } else {
                    if x < a {
                        out.push((x, a - 1));
                    }
                    if y > b {
                        out.push((b + 1, y));
                    }
                }
            }
            self.0 = out;
        }
    }

    fn covered_length(&self) -> u32 {
        self.0.iter().fold(0u32, |acc, &(a, b)| {
            acc.wrapping_add(b.wrapping_sub(a).wrapping_add(1))
        })
    }
}

/// `s_FixMinusStrandRange`.
fn fix_range(r: (u32, u32)) -> (u32, u32) {
    if r.0 > r.1 { (r.1, r.0) } else { r }
}

/// `CAlnMap::GetSeqPosFromAlnPos(0, GetAlnPosFromSeqPos(1, subject_pos))`: the plus-strand
/// query position (0-based) aligned with a subject position (0-based), or -1 at a query gap.
fn query_pos_for_subject_pos(h: &Hsp, query_len: i64, subject_pos: i64) -> i64 {
    let Some(es) = &h.edit_script else {
        // ungapped: a straight diagonal in context coordinates
        let ctx_q = h.query_offset as i64 + (subject_pos - h.subject_offset as i64);
        return if h.minus {
            query_len - 1 - ctx_q
        } else {
            ctx_q
        };
    };
    let (mut q, mut s) = (h.query_offset as i64, h.subject_offset as i64);
    for &(op, n) in &es.ops {
        let n = n as i64;
        match op {
            crate::engine::hsp::EditOp::Sub => {
                if subject_pos < s + n {
                    let ctx_q = q + (subject_pos - s);
                    return if h.minus {
                        query_len - 1 - ctx_q
                    } else {
                        ctx_q
                    };
                }
                q += n;
                s += n;
            }
            crate::engine::hsp::EditOp::Del => {
                if subject_pos < s + n {
                    return -1;
                }
                s += n;
            }
            crate::engine::hsp::EditOp::Ins => q += n,
        }
    }
    -1
}

/// `GetUniqSeqCoverage`: query bases covered by the parts of each HSP whose subject range no
/// earlier HSP of the same subject already covered.
fn uniq_coverage(hsps: &[Hsp], query_len: i64) -> i32 {
    let Some(first) = hsps.first() else { return 0 };
    let subj_range = |h: &Hsp| -> (u32, u32) {
        let (a, b) = (h.s_start as u32 - 1, h.s_end as u32 - 1);
        if a > b { (b, a) } else { (a, b) }
    };
    let query_range =
        |h: &Hsp| -> (u32, u32) { fix_range((h.q_start as u32 - 1, h.q_end as u32 - 1)) };
    let mut subj_coll = RangeColl::default();
    subj_coll.add(subj_range(first));
    let mut query_coll = RangeColl::default();
    query_coll.add(query_range(first));
    for h in &hsps[1..] {
        let sr = subj_range(h);
        let mut coll = RangeColl::default();
        coll.add(sr);
        coll.subtract(&subj_coll);
        if coll.0.is_empty() {
            continue;
        }
        if coll.0[0] == sr {
            query_coll.add(query_range(h));
            subj_coll.add(sr);
        } else {
            for &(a, b) in &coll.0.clone() {
                let qa = query_pos_for_subject_pos(h, query_len, a as i64) as u32;
                let qb = query_pos_for_subject_pos(h, query_len, b as i64) as u32;
                query_coll.add(fix_range((qa, qb)));
                subj_coll.add((a, b));
            }
        }
    }
    query_coll.covered_length() as i32
}

/// Per-hit values shared by its HSPs (`InsertSubjectScores` sets them on the first alignment
/// of a subject; the following ones keep them).
struct HitCoverage {
    qcovs: Option<i64>,
    qcovus: Option<i64>,
}

fn hit_coverage(spec: &TabularSpec, q: &QueryResult, hit: &Hit) -> HitCoverage {
    let qlen = q.length as i64;
    let qcovs = spec
        .has(Field::Qcovs)
        .then(|| master_coverage(&hit.hsps))
        .filter(|c| *c != 0)
        .map(|c| coverage_percent(c, qlen));
    let qcovus = spec
        .has(Field::Qcovus)
        .then(|| uniq_coverage(&hit.hsps, qlen))
        .filter(|c| *c != 0)
        .map(|c| coverage_percent(c as i64, qlen));
    HitCoverage { qcovs, qcovus }
}

/// `SetFields`' identity recount: when the alignment strings are built, `nident`,
/// `positives` and the BTOP come from comparing them character by character.
fn recount(h: &Hsp) -> (i32, i32, String) {
    let (Some(qs), Some(ss)) = (&h.qseq, &h.sseq) else {
        return (h.num_ident, h.num_ident, String::new());
    };
    let (qb, sb) = (qs.as_bytes(), ss.as_bytes());
    let mut btop = String::new();
    let mut ident = 0i32;
    let mut run = 0i64;
    for i in 0..qb.len().min(sb.len()) {
        if qb[i] == sb[i] {
            ident += 1;
            run += 1;
        } else {
            if run > 0 {
                btop.push_str(&run.to_string());
                run = 0;
            }
            btop.push(qb[i] as char);
            btop.push(sb[i] as char);
        }
    }
    if run > 0 {
        btop.push_str(&run.to_string());
    }
    (ident, ident, btop)
}

struct Row<'a> {
    q: &'a QueryResult,
    subj: &'a crate::plan::Subject,
    h: &'a Hsp,
    cov: &'a HitCoverage,
    num_ident: i32,
    positives: i32,
    btop: String,
}

fn field_value(f: Field, r: &Row) -> String {
    let (q, subj, h) = (r.q, r.subj, r.h);
    let d = &subj.deflines;
    match f {
        // the query is a local id whose first title token prints everywhere
        // (s_ReplaceLocalId) and has no gi
        Field::Qseqid | Field::Qacc | Field::Qaccver => q.id.clone(),
        Field::Qgi | Field::Sgi => "0".into(),
        Field::Qlen => q.length.to_string(),
        Field::Sseqid => d
            .first()
            .map(|x| first_id(x, seqid_string))
            .unwrap_or_else(|| "Unknown".into()),
        Field::Sallseqid => per_defline(d, seqid_string),
        Field::Sallgi => d.iter().map(|_| "0").collect::<Vec<_>>().join(";"),
        Field::Sacc | Field::Saccver => d
            .first()
            .map(|x| first_id(x, acc_string))
            .unwrap_or_else(|| "Unknown".into()),
        Field::Sallacc => per_defline(d, acc_string),
        Field::Slen => subj.length.to_string(),
        Field::Qstart => h.q_start.to_string(),
        Field::Qend => h.q_end.to_string(),
        Field::Sstart => h.s_start.to_string(),
        Field::Send => h.s_end.to_string(),
        Field::Qseq => h.qseq.clone().unwrap_or_default(),
        Field::Sseq => h.sseq.clone().unwrap_or_default(),
        Field::Evalue => evalue_string(h.evalue),
        Field::Bitscore => bitscore_string(h.bit_score),
        Field::Score => h.score.to_string(),
        Field::Length => h.align_length.to_string(),
        Field::Pident => {
            let p = if h.align_length > 0 {
                r.num_ident as f64 / h.align_length as f64 * 100.0
            } else {
                0.0
            };
            format!("{:.3}", p)
        }
        Field::Nident => r.num_ident.to_string(),
        Field::Mismatch => (h.align_length - r.num_ident - h.gaps).to_string(),
        Field::Positive => r.positives.to_string(),
        Field::Gapopen => h.gap_opens.to_string(),
        Field::Gaps => h.gaps.to_string(),
        Field::Ppos => {
            let p = if h.align_length > 0 {
                r.positives as f64 / h.align_length as f64 * 100.0
            } else {
                0.0
            };
            format!("{:.2}", p)
        }
        Field::Frames => {
            if h.minus {
                "1/-1".into()
            } else {
                "1/1".into()
            }
        }
        Field::Qframe => "1".into(),
        Field::Sframe => {
            if h.minus {
                "-1".into()
            } else {
                "1".into()
            }
        }
        Field::Btop => r.btop.clone(),
        Field::Staxid => d
            .iter()
            .filter_map(|x| x.taxid)
            .find(|t| *t != 0)
            .map(|t| t.to_string())
            .unwrap_or_else(|| "N/A".into()),
        Field::Staxids => {
            let mut v: Vec<i32> = d.iter().filter_map(|x| x.taxid).collect();
            v.sort_unstable();
            v.dedup();
            if v.is_empty() {
                "N/A".into()
            } else {
                v.iter()
                    .map(|t| t.to_string())
                    .collect::<Vec<_>>()
                    .join(";")
            }
        }
        Field::Stitle => d
            .first()
            .map(|x| title_or_na(&x.title).to_string())
            .unwrap_or_else(|| "N/A".into()),
        Field::Salltitles => {
            if d.is_empty() {
                "N/A".into()
            } else {
                d.iter()
                    .map(|x| title_or_na(&x.title))
                    .collect::<Vec<_>>()
                    .join("<>")
            }
        }
        Field::Sstrand => {
            if h.minus {
                "minus".into()
            } else {
                "plus".into()
            }
        }
        Field::Qcovs => r
            .cov
            .qcovs
            .map(|c| c.to_string())
            .unwrap_or_else(|| "N/A".into()),
        Field::Qcovus => r
            .cov
            .qcovus
            .map(|c| c.to_string())
            .unwrap_or_else(|| "N/A".into()),
        Field::Qcovhsp => {
            coverage_percent((h.q_end - h.q_start).abs() as i64 + 1, q.length as i64).to_string()
        }
    }
}

/// Write every HSP of `results` as tab-separated lines (blastn's hit and HSP order).
pub fn write(
    results: &SearchResults,
    plan: &SubjectPlan,
    spec: &TabularSpec,
    out: &mut dyn Write,
) -> Result<(), Error> {
    let need_recount = spec.needs_sequences();
    let mut text = String::new();
    for q in &results.queries {
        for hit in &q.hits {
            let cov = hit_coverage(spec, q, hit);
            let subj = &plan.subjects[hit.oid as usize];
            for h in &hit.hsps {
                let (num_ident, positives, btop) = if need_recount {
                    recount(h)
                } else {
                    (h.num_ident, h.num_ident, String::new())
                };
                let row = Row {
                    q,
                    subj,
                    h,
                    cov: &cov,
                    num_ident,
                    positives,
                    btop,
                };
                for (i, f) in spec.fields.iter().enumerate() {
                    if i > 0 {
                        text.push('\t');
                    }
                    text.push_str(&field_value(*f, &row));
                }
                text.push('\n');
            }
        }
        if text.len() > (1 << 20) {
            out.write_all(text.as_bytes())
                .map_err(|e| Error::io("<output>", e))?;
            text.clear();
        }
    }
    out.write_all(text.as_bytes())
        .map_err(|e| Error::io("<output>", e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn score_strings_match_c_printf() {
        assert_eq!(evalue_string(0.0), "0.0");
        assert_eq!(evalue_string(3.2e-120), "3.20e-120");
        assert_eq!(evalue_string(1.234e-5), "1.23e-05");
        assert_eq!(evalue_string(0.00089), "8.90e-04");
        assert_eq!(evalue_string(0.0123456), "0.012");
        assert_eq!(evalue_string(0.5), "0.50");
        assert_eq!(evalue_string(2.345), "2.3");
        assert_eq!(evalue_string(10.0), "10");
        assert_eq!(bitscore_string(1234.7), "1234");
        assert_eq!(bitscore_string(99.95), " 99"); // > 99.9 -> %3.0ld of the truncated long
        assert_eq!(bitscore_string(99.9), "99.9");
        assert_eq!(bitscore_string(45.0), "45.0");
        assert_eq!(bitscore_string(9.0), " 9.0");
        assert_eq!(bitscore_string(123456.0), "1.235e+05");
        // round-half-even like glibc
        assert_eq!(format!("{:.0}", 2.5), "2");
        assert_eq!(format!("{:.2}", 0.125), "0.12");
        assert_eq!(format!("{:.3}", 7.0 / 8.0 * 100.0), "87.500");
    }

    #[test]
    fn outfmt_grammar() {
        let s = TabularSpec::parse_outfmt("6").unwrap();
        assert_eq!(s.fields, STD_FIELDS.to_vec());
        assert_eq!(s.fields[0], Field::Qaccver);
        let s = TabularSpec::parse_outfmt("6 std stitle").unwrap();
        assert_eq!(s.fields.len(), 13);
        // blastn ignores unknown tokens; the CLI warns
        let s = TabularSpec::parse_outfmt("6 std bogus").unwrap();
        assert_eq!((s.fields.len(), s.ignored), (12, vec!["bogus".to_string()]));
        for f in ["0", "7 qseqid", "10"] {
            assert!(matches!(
                TabularSpec::parse_outfmt(f),
                Err(Error::InvalidOptions(_))
            ));
        }
        assert!(matches!(
            TabularSpec::parse_outfmt("6 ssciname"),
            Err(Error::InvalidOptions(_))
        ));
        assert!(
            TabularSpec::parse_outfmt("6 std qseq")
                .unwrap()
                .needs_sequences()
        );
        assert!(
            TabularSpec::parse_outfmt("6 std positive")
                .unwrap()
                .needs_sequences()
        );
        assert!(
            !TabularSpec::parse_outfmt("6 std nident")
                .unwrap()
                .needs_sequences()
        );
        for (n, f) in FIELDS {
            assert_eq!(Field::parse(n), Some(f));
            assert_eq!(f.name(), n);
        }
    }

    #[test]
    fn coverage_rounding() {
        // x_SetQueryCovSeqalign: +0.5 below 99, then truncation
        assert_eq!(coverage_percent(996, 1000), 99);
        assert_eq!(coverage_percent(985, 1000), 99);
        assert_eq!(coverage_percent(984, 1000), 98);
        assert_eq!(coverage_percent(1000, 1000), 100);
        assert_eq!(coverage_percent(1, 1000), 0);
    }

    #[test]
    fn range_collection() {
        let mut c = RangeColl::default();
        c.add((10, 20));
        c.add((21, 30)); // abutting merges
        assert_eq!(c.0, vec![(10, 30)]);
        c.add((40, 50));
        c.add((5, 12));
        assert_eq!(c.0, vec![(5, 30), (40, 50)]);
        let mut d = RangeColl::default();
        d.add((0, 100));
        d.subtract(&c);
        assert_eq!(d.0, vec![(0, 4), (31, 39), (51, 100)]);
        assert_eq!(d.covered_length(), 5 + 9 + 50);
    }
}

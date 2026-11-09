//! `qseq`/`sseq`: the aligned sequence strings of each HSP, built from its edit script.
//!
//! Source: BLAST+ 2.17.0 `objtools/align_format/tabular.cpp::SetFields` (the Dense-seg of a
//! minus-strand hit is reversed so the query shows its plus strand and the subject its minus
//! strand; `CAlnVec::GetWholeAlnSeqString` with gap char `-`). The subject letters are the
//! real IUPAC codes (`blastna`, as the traceback saw them), the query letters are upper case
//! with `U` read as `T` (what blastn's FASTA reader stores).

use crate::engine::encoding::{BLASTNA_COMPLEMENT, BLASTNA_TO_IUPAC, reverse_complement};
use crate::engine::hsp::EditOp;
use crate::engine::query::QuerySet;
use crate::engine::search::{Hsp, SearchResults};
use crate::engine::subject::StoreSubjects;
use crate::error::Error;

/// Aligned (query, subject) letters for one HSP, in blastn's printed orientation. `subject`
/// holds the subject's letters from `subject_base` on.
pub fn strings(
    h: &Hsp,
    context_query: &[u8],
    subject: &[u8],
    subject_base: i32,
) -> (String, String) {
    let mut q = Vec::new();
    let mut s = Vec::new();
    let (mut qi, mut si) = (
        h.query_offset as usize,
        (h.subject_offset - subject_base) as usize,
    );
    let mut push = |q: &mut Vec<u8>, s: &mut Vec<u8>, op: EditOp, n: usize| match op {
        EditOp::Sub => {
            q.extend(
                context_query[qi..qi + n]
                    .iter()
                    .map(|&c| BLASTNA_TO_IUPAC[c as usize]),
            );
            s.extend(
                subject[si..si + n]
                    .iter()
                    .map(|&c| BLASTNA_TO_IUPAC[c as usize]),
            );
            qi += n;
            si += n;
        }
        EditOp::Del => {
            q.extend(std::iter::repeat_n(b'-', n));
            s.extend(
                subject[si..si + n]
                    .iter()
                    .map(|&c| BLASTNA_TO_IUPAC[c as usize]),
            );
            si += n;
        }
        EditOp::Ins => {
            q.extend(
                context_query[qi..qi + n]
                    .iter()
                    .map(|&c| BLASTNA_TO_IUPAC[c as usize]),
            );
            s.extend(std::iter::repeat_n(b'-', n));
            qi += n;
        }
    };
    match &h.edit_script {
        Some(es) => {
            for &(op, n) in &es.ops {
                push(&mut q, &mut s, op, n as usize);
            }
        }
        None => push(
            &mut q,
            &mut s,
            EditOp::Sub,
            (h.query_end - h.query_offset) as usize,
        ),
    }
    if h.minus {
        // the context is the reverse complement of the query: flip both rows so the query
        // reads on its plus strand and the subject on its minus strand
        q = complement_letters(&q);
        s = complement_letters(&s);
    }
    (String::from_utf8(q).unwrap(), String::from_utf8(s).unwrap())
}

fn complement_letters(v: &[u8]) -> Vec<u8> {
    v.iter()
        .rev()
        .map(|&c| {
            if c == b'-' {
                b'-'
            } else {
                let code = crate::engine::encoding::ASCII_TO_BLASTNA[c as usize];
                BLASTNA_TO_IUPAC[BLASTNA_COMPLEMENT[code as usize] as usize]
            }
        })
        .collect()
}

/// Fill `qseq`/`sseq` of every HSP, decoding only the subject letters the HSP covers.
pub fn fill(
    results: &mut SearchResults,
    queries: &QuerySet,
    subjects: &StoreSubjects,
) -> Result<(), Error> {
    let mut window = Vec::new();
    for (qi, q) in results.queries.iter_mut().enumerate() {
        let plus = &queries.queries[qi].seq;
        let minus = reverse_complement(plus);
        for hit in &mut q.hits {
            for h in &mut hit.hsps {
                subjects.blastna_window(
                    hit.oid as i32,
                    h.subject_offset..h.subject_end,
                    &mut window,
                )?;
                let ctx: &[u8] = if h.minus { &minus } else { plus };
                let (qs, ss) = strings(h, ctx, &window, h.subject_offset);
                h.qseq = Some(qs);
                h.sseq = Some(ss);
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::encoding::ASCII_TO_BLASTNA;
    use crate::engine::hsp::EditScript;

    fn hsp(minus: bool, qo: i32, so: i32, ops: &[(EditOp, i32)]) -> Hsp {
        Hsp {
            score: 0,
            bit_score: 0.0,
            evalue: 0.0,
            num_ident: 0,
            align_length: 0,
            gaps: 0,
            gap_opens: 0,
            q_start: 0,
            q_end: 0,
            s_start: 0,
            s_end: 0,
            minus,
            query_offset: qo,
            query_end: qo
                + ops
                    .iter()
                    .filter(|(o, _)| *o != EditOp::Del)
                    .map(|(_, n)| n)
                    .sum::<i32>(),
            subject_offset: so,
            subject_end: 0,
            context: 0,
            edit_script: Some(EditScript { ops: ops.to_vec() }),
            qseq: None,
            sseq: None,
        }
    }

    fn code(s: &str) -> Vec<u8> {
        s.bytes().map(|c| ASCII_TO_BLASTNA[c as usize]).collect()
    }

    #[test]
    fn gapped_strings_and_minus_orientation() {
        let q = code("ACGTTGCA");
        let s = code("ACGGTTGCA");
        // ACG-TTGCA vs ACGGTTGCA: a deletion (gap in the query) after 3 matches
        let h = hsp(
            false,
            0,
            0,
            &[(EditOp::Sub, 3), (EditOp::Del, 1), (EditOp::Sub, 5)],
        );
        assert_eq!(
            strings(&h, &q, &s, 0),
            ("ACG-TTGCA".into(), "ACGGTTGCA".into())
        );
        // the same hit seen from the reverse-complemented query context prints the query on
        // its plus strand and the subject reverse-complemented
        let rc = reverse_complement(&q);
        let h = hsp(
            true,
            0,
            0,
            &[(EditOp::Sub, 5), (EditOp::Ins, 1), (EditOp::Sub, 2)],
        );
        let (qs, ss) = strings(&h, &rc, &code("TGCAAGT"), 0);
        // the same letters as a window starting at subject letter 10
        let mut h10 = h.clone();
        h10.subject_offset = 10;
        assert_eq!(
            strings(&h10, &rc, &code("TGCAAGT"), 10),
            (qs.clone(), ss.clone())
        );
        assert_eq!(qs.len(), ss.len());
        assert_eq!(qs, "ACGTTGCA");
        assert_eq!(ss, "AC-TTGCA");
    }
}

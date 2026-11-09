//! Report writers: the TSV summary, JSON and the per-genome detail table. Every sort is total.

use super::rank::Report;
use std::fmt::Write;

pub const TSV_HEADER: &str = "#query\tstatus\ttaxid\ttaxon\trank\ttop_genome\ttop_seq\tpident\tqcov\tbits\tevalue\ttied_genomes\ttied_taxa\tsearched\tflags";
pub const DETAIL_HEADER: &str = "#query\ttaxid\ttaxon\tgenome\tcollection\tseq_name\tseq_digest\tsstart\tsend\tstrand\tpident\tlength\tqcov\tevalue\tbits\ttied\tstage";

/// `--format tsv`: one line per query.
pub fn tsv(r: &Report) -> String {
    let mut s = format!("{TSV_HEADER}\n");
    for x in r.summary_rows() {
        let _ = writeln!(
            s,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            x.query,
            x.status,
            x.taxid,
            x.taxon,
            x.rank,
            x.top_genome,
            x.top_seq,
            x.pident,
            x.qcov,
            x.bits,
            x.evalue,
            x.tied_genomes,
            x.tied_taxa,
            x.searched,
            x.flags
        );
    }
    s
}

/// `--format json`.
pub fn json(r: &Report) -> String {
    serde_json::to_string_pretty(r).unwrap() + "\n"
}

/// `--detail`: one row per (query, genome), grouped by taxon; `top` keeps the best N genomes
/// per taxon.
pub fn detail(r: &Report, top: Option<usize>) -> String {
    let mut s = format!("{DETAIL_HEADER}\n");
    for q in &r.queries {
        for t in &q.taxa {
            for g in t.genomes.iter().take(top.unwrap_or(usize::MAX)) {
                let _ = writeln!(
                    s,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{}\t{}\t{}\t{}",
                    q.query,
                    t.taxid,
                    t.name,
                    g.genome,
                    g.collection,
                    g.seq_name,
                    g.seq_digest,
                    g.sstart,
                    g.send,
                    g.strand,
                    g.pident,
                    g.length,
                    g.qcov,
                    g.evalue,
                    g.bits,
                    if g.tied { "yes" } else { "no" },
                    g.stage
                );
            }
        }
    }
    s
}

#[cfg(test)]
mod column_tests {
    #[test]
    fn summary_columns_match_the_tsv_header() {
        assert_eq!(
            format!("#{}", super::super::rank::SUMMARY_COLUMNS.join("\t")),
            super::TSV_HEADER
        );
    }
}

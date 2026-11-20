//! The representative rule on the jungle-shaped membership table (79 collections, 27,721 rows;
//! `scripts/make_origin_fixtures.sh` writes it from a jungle mirror).

use flate2::read::GzDecoder;
use gtars_rgblast::origin::reps::{GroupKey, RepColl, choose};
use std::io::{BufRead, BufReader};
use std::path::Path;

const DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/fixtures/origin");

fn jungle(max_n: Option<f64>) -> Vec<RepColl> {
    let mut colls: Vec<RepColl> =
        std::fs::read_to_string(Path::new(DIR).join("jungle_collections.tsv"))
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                let taxid: u32 = f[3].parse().unwrap();
                RepColl {
                    digest: f[1].into(),
                    display: f[2].into(),
                    group: GroupKey::Taxid(taxid),
                    group_name: if taxid == 9606 {
                        "Homo sapiens"
                    } else {
                        "Mus musculus"
                    }
                    .into(),
                    n_fraction: f[4].parse().ok(),
                    in_namespace: false,
                    seqs: Vec::new(),
                }
            })
            .collect();
    let gz = GzDecoder::new(
        std::fs::File::open(Path::new(DIR).join("jungle_membership.tsv.gz")).unwrap(),
    );
    for l in BufReader::new(gz).lines().map_while(Result::ok) {
        if l.starts_with('#') {
            continue;
        }
        let f: Vec<u64> = l.split('\t').map(|x| x.parse().unwrap()).collect();
        colls[f[0] as usize].seqs.push((f[1] as u32, f[2]));
    }
    if let Some(m) = max_n {
        colls.retain(|c| c.n_fraction.is_none_or(|f| f < m));
    }
    colls
}

#[test]
fn jungle_representatives_and_costs() {
    let all = jungle(None);
    assert_eq!(all.len(), 79);
    assert_eq!(all.iter().map(|c| c.seqs.len()).sum::<usize>(), 27_721);
    let t = choose(&all, &[], 0).unwrap();
    let reps: Vec<(&str, &str, &str)> = t
        .groups
        .iter()
        .map(|g| (g.rep.as_str(), g.rep_collection.as_str(), g.rule.as_str()))
        .collect();
    assert_eq!(
        reps,
        [
            (
                "GRCh38.p14-fasta-genomic",
                "u1HyLgIlq8M_XvEwy0oGqAvKGHJMGtxH",
                "most-shared"
            ),
            (
                "GRCm38.p6-fasta-genomic-GCA",
                "JL56x8L1q1Fs_-jHvZxBG01Vitac-CmO",
                "most-shared"
            ),
        ]
    );
    let gbp = |b: u64| format!("{:.2}", b as f64 / 1e9);
    assert_eq!(
        (gbp(t.summary.reps_bases), gbp(t.summary.full_bases)),
        ("6.12".into(), "41.81".into())
    );

    // --exclude-masked: the same picks, 6.12 of 21.08 Gbp, worst staged query 15.57 Gbp
    let unmasked = jungle(Some(0.25));
    assert_eq!(unmasked.len(), 67);
    let t2 = choose(&unmasked, &[], 0).unwrap();
    assert_eq!(
        t2.groups
            .iter()
            .map(|g| &g.rep_collection)
            .collect::<Vec<_>>(),
        t.groups
            .iter()
            .map(|g| &g.rep_collection)
            .collect::<Vec<_>>()
    );
    let table = t2.render();
    let lines: Vec<&str> = table.lines().collect();
    assert_eq!(
        lines[1],
        "9606\tHomo sapiens\t45\tGRCh38.p14-fasta-genomic\tu1HyLgIlq8M_XvEwy0oGqAvKGHJMGtxH\tmost-shared\t3.30\t12.75\t3225\t9.45"
    );
    assert_eq!(
        lines[2],
        "10090\tMus musculus\t22\tGRCm38.p6-fasta-genomic-GCA\tJL56x8L1q1Fs_-jHvZxBG01Vitac-CmO\tmost-shared\t2.82\t8.33\t219\t5.51"
    );
    assert_eq!(
        lines[3],
        "# representatives: 2 of 67 collections, 6.12 of 21.08 Gbp (3.4x smaller)"
    );
    assert_eq!(
        lines[4],
        "# worst-case staged query (genome level): 15.57 Gbp = 0.74 of a full search -> staged mode on"
    );
    // all 79: the human worst case is 25.31 of 41.81 Gbp
    assert_eq!(gbp(t.summary.worst_genome_bases), "25.31");
}

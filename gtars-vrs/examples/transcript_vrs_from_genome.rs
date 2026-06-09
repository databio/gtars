//! Demo: transcript-level VRS lookups from genome-level sequences only.
//!
//! The point of this example is that the **transcript sequence is never
//! stored**. The RefgetStore holds a single chromosome. We use the reftx
//! transcript structure (exons + strand) to *derive* the spliced mRNA at
//! runtime (`mature_mrna_for_transcript`, the region-concatenation lookup),
//! compute that derived transcript's refget identity, and anchor a VRS
//! Allele on it — producing a transcript-level `ga4gh:VA.` identifier built
//! entirely out of genome-level sequence.
//!
//! Run with:
//!   cargo run -p gtars-vrs --example transcript_vrs_from_genome --features filesystem

#[cfg(not(feature = "transcripts"))]
fn main() {
    eprintln!(
        "This demo needs the `transcripts` feature. Re-run with:\n  \
         cargo run -p gtars-vrs --example transcript_vrs_from_genome --features filesystem"
    );
}

#[cfg(feature = "transcripts")]
fn main() {
    use gtars_refget::digest::{digest_sequence, sha512t24u};
    use gtars_refget::store::RefgetStore;
    use gtars_refget::transcripts::{build_reftx_bytes_in_memory, ReadonlyTxStore};
    use gtars_refget::{mature_mrna_for_transcript, Exon, ManeStatus, Strand, Transcript};
    use gtars_vrs::{
        allele_identifier, hgvs_str_to_transcript_vrs_id_readonly, Allele, AlleleState,
        SequenceLocation, SequenceReference,
    };

    // ----------------------------------------------------------------------
    // 1. A tiny genome: ONE chromosome. This is the only sequence we store.
    //    Layout (0-based, half-open):
    //      [0,4)   TTTT       upstream
    //      [4,12)  ACGTACGT   << exon 1
    //      [12,20) GGGGGGGG   intron (spliced out — never in the mRNA)
    //      [20,28) TTAACCGG   << exon 2
    //      [28,32) AAAA       downstream
    // ----------------------------------------------------------------------
    let chrom_seq = "TTTTACGTACGTGGGGGGGGTTAACCGGAAAA";
    let chrom_record = digest_sequence("chr_demo", chrom_seq.as_bytes());
    let chrom_key = chrom_record.metadata().sha512t24u.clone(); // base64url, no prefix

    let mut store = RefgetStore::in_memory();
    store.disable_encoding(); // Raw mode: substrings come back as literal ASCII bases
    store.add_sequence_record(chrom_record, true).unwrap();
    let store = store.into_readonly();

    println!("== Stored genome ==");
    println!("  chr_demo  ({} bp)  SQ.{}", chrom_seq.len(), chrom_key);
    println!("  (no transcript sequence is stored — only this chromosome)\n");

    // ----------------------------------------------------------------------
    // 2. A transcript defined purely as STRUCTURE over that chromosome:
    //    two forward-strand exons. chrom_digest is the raw 24 bytes that
    //    base64url-encode to the chromosome's store key.
    // ----------------------------------------------------------------------
    let chrom_digest: [u8; 24] = base64_url::decode(&chrom_key)
        .unwrap()
        .try_into()
        .expect("sha512t24u is 24 bytes");

    let tx = Transcript {
        accession: "NM_DEMO.1".to_string(),
        gene: "DEMO".to_string(),
        chrom_digest,
        strand: Strand::Forward,
        cds_start: None,
        cds_end: None,
        exons: vec![Exon { start: 4, end: 12 }, Exon { start: 20, end: 28 }],
        mane: ManeStatus::default(),
    };

    // ----------------------------------------------------------------------
    // 3. DERIVE the mature mRNA from genome + exon structure (the #4 lookup).
    // ----------------------------------------------------------------------
    let mrna = mature_mrna_for_transcript(&store, &tx).unwrap();
    let tx_key = sha512t24u(&mrna); // the derived transcript's refget identity

    println!("== Derived transcript (NOT stored) ==");
    println!("  exon1 [4,12)  + exon2 [20,28)  spliced (intron GGGGGGGG dropped)");
    println!("  mature mRNA   = {}  ({} bp)", mrna, mrna.len());
    println!("  refget digest = SQ.{}  <- computed from the DERIVED sequence\n", tx_key);

    // ----------------------------------------------------------------------
    // 4. A variant, expressed on the TRANSCRIPT as HGVS: NM_DEMO.1:n.6C>A.
    //    The transcript-anchored bridge does ALL of the work: it derives the
    //    mature mRNA from the genome + exon structure, digests it, projects
    //    the n. coordinate onto the mature-mRNA offset, validates the ref base
    //    against the derived sequence, and emits a transcript-anchored VRS id.
    //    No hand-rolled coordinate math here.
    // ----------------------------------------------------------------------
    let tx_pos: u64 = 5; // interbase start of n.6 (1-based 6 -> 0-based 5)
    let alt = "A";
    let ref_base = &mrna[tx_pos as usize..tx_pos as usize + 1];

    // Build a one-transcript readonly store for the bridge to look up.
    let tx_bytes = build_reftx_bytes_in_memory(&[tx.clone()]).unwrap();
    let tx_store = ReadonlyTxStore::from_bytes(tx_bytes).unwrap();

    let tx_vrs_id =
        hgvs_str_to_transcript_vrs_id_readonly("NM_DEMO.1:n.6C>A", &store, &tx_store)
            .expect("transcript-anchored bridge")
            .value;

    println!("== Transcript-anchored VRS allele ==");
    println!("  NM_DEMO.1 : n.6{}>{}  (interbase [{},{}) on the mature mRNA)", ref_base, alt, tx_pos, tx_pos + 1);
    println!("  reference = SQ.{} (the derived transcript)", tx_key);
    println!("  VRS id    = {}  <- via hgvs_str_to_transcript_vrs_id_readonly\n", tx_vrs_id);

    // ----------------------------------------------------------------------
    // 5. Contrast: the SAME physical variant anchored on the GENOME.
    //    Project the transcript offset onto genomic coordinates (forward
    //    strand: walk exons until the offset lands inside one). reftx's
    //    CoordinateMapper does this generally; the loop here keeps the demo
    //    self-contained.
    // ----------------------------------------------------------------------
    let mut remaining = tx_pos;
    let mut genomic_pos = None;
    for ex in &tx.exons {
        let len = (ex.end - ex.start) as u64;
        if remaining < len {
            genomic_pos = Some(ex.start as u64 + remaining);
            break;
        }
        remaining -= len;
    }
    let genomic_pos = genomic_pos.unwrap();
    let genome_ref_base = &chrom_seq[genomic_pos as usize..genomic_pos as usize + 1];

    let genome_allele = Allele {
        location: SequenceLocation {
            sequence_reference: SequenceReference {
                refget_accession: format!("SQ.{}", chrom_key),
            },
            start: genomic_pos,
            end: genomic_pos + 1,
        },
        state: AlleleState::LiteralSequenceExpression {
            sequence: alt.to_string(),
        },
    };
    let genome_vrs_id = allele_identifier(&genome_allele);

    println!("== Genome-anchored VRS allele (same physical variant) ==");
    println!("  chr_demo  : interbase [{},{}) {}>{}", genomic_pos, genomic_pos + 1, genome_ref_base, alt);
    println!("  reference = SQ.{} (the chromosome)", chrom_key);
    println!("  VRS id    = {}\n", genome_vrs_id);

    // Sanity: the reference base agrees on both coordinate systems.
    assert_eq!(
        ref_base, genome_ref_base,
        "transcript and genome reference base must agree"
    );

    println!("== Takeaway ==");
    println!("  One stored chromosome yielded BOTH a genome-anchored and a");
    println!("  transcript-anchored VRS id for the same variant. The transcript");
    println!("  sequence behind SQ.{} was never stored — it was composed on the", tx_key);
    println!("  fly from the genome via the transcript's exon structure.");
}

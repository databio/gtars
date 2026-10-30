# Real cdot fixtures

Copied unchanged from the cdot repository's own test data
(https://github.com/SACGF/cdot, `tests/test_data/`, MIT license):

| File | cdot version | Build | Transcripts |
|------|--------------|-------|-------------|
| `cdot.refseq.grch37.json` | 0.2.10 | GRCh37 | NM_001637.3 (minus strand, has an alignment gap), NR_023343.1 (non-coding) |
| `cdot.ensembl.grch38.json` | 0.2.26 | GRCh38 | ENST00000617537.5 (minus strand, MANE Select) |

They exist so `TxStoreBuilder::ingest_cdot` is tested against the real file
layout (`genome_builds.<build>`, `"+"`/`"-"` strand, 6-value exons) rather
than a shape we made up ourselves.

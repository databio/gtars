# gtars-reftx tests

## Vendored HGVS corpora

The `tests/data/hgvs/biocommons/gcp/*.tsv` truth tables come from
biocommons/hgvs (Apache-2.0). They reference real GRCh37/GRCh38
transcripts (NM_xxxxxx). Running `project_biocommons_gcp_full` therefore
requires a built reftx binary that contains those transcripts.

## Building a reftx fixture for the gcp tests

```bash
# Pick a cdot transcript JSON release (https://github.com/SACGF/cdot/releases)
# that includes MANE Select for the relevant genes:
#   ADRA2B, BAHCC1, DNAH11, FOLR3, JRK, NEFL, ORAI1, ZCCHC3
# plus everything else the gcp `real.tsv`, `regression.tsv`, and `noncoding.tsv`
# tables reference (see grep below).
#
# Required transcript accessions (extracted from gcp/*.tsv HGVSc column):
grep -h -oE 'NM_[0-9]+\.[0-9]+|NR_[0-9]+\.[0-9]+|XM_[0-9]+\.[0-9]+' \
  tests/data/hgvs/biocommons/gcp/*.tsv | sort -u

# Build reftx with the gtars-cli reftx subcommand once available:
#   gtars reftx build --cdot path/to/cdot.json --output reftx.bin
GTARS_REFTX_FIXTURE=/abs/path/reftx.bin \
  cargo test -p gtars-reftx --test test_hgvs_corpus_mapper -- --ignored
```

Until the reftx builder ships, the gcp test stays `#[ignore]` and CI
runs only the cheap "list & parse" sanity test.

# Refreshing the vendored HGVS coordinate-mapping corpora

These fixtures are vendored verbatim from biocommons/hgvs. Do not edit in
place — re-run the script below to refresh from upstream.

## biocommons/hgvs (Apache-2.0)

```bash
SHA=fa02ca5ec0c9e28a1a6377b551e41ae77b450ce9   # update to current upstream main
BASE=https://raw.githubusercontent.com/biocommons/hgvs/$SHA
DEST=gtars-reftx/tests/data/hgvs/biocommons

# Synthetic NM_999999.x transcript metadata used by the parser harness and
# (eventually) by an in-process SampleDataProvider. Note: upstream calls
# this file `sanity_cp.tsv` and stores it at the top of tests/data/, but it
# is in fact a synthetic transcript table — we put it under sample_data/
# to reflect that.
curl -fsSL "$BASE/tests/data/sanity_cp.tsv" -o "$DEST/sample_data/sanity_cp.tsv"

# Real-transcript HGVSg <-> HGVSc projection truth tables. These reference
# real GRCh37/GRCh38 transcripts (NM_xxxxxx) and require a built reftx
# binary that contains those transcripts in order to run.
for f in ADRA2B-dbSNP.tsv BAHCC1-dbSNP.tsv DNAH11-HGMD.tsv \
         DNAH11-dbSNP-NM_001277115.tsv DNAH11-dbSNP-NM_003777.tsv \
         DNAH11-dbSNP.tsv FOLR3-dbSNP.tsv JRK-dbSNP.tsv NEFL-dbSNP.tsv \
         ORAI1-dbSNP.tsv ZCCHC3-dbSNP.tsv noncoding.tsv real.tsv \
         regression.tsv ; do
  curl -fsSL "$BASE/tests/data/gcp/$f" -o "$DEST/gcp/$f"
done
```

After refreshing, also update the commit SHA in the sibling `NOTICE` file.

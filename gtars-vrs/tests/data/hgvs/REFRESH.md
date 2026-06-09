# Refreshing the vendored HGVS parser corpora

These fixtures are vendored verbatim. Do not edit them in place — re-run
the steps below to refresh from upstream.

## varfish-org/hgvs-rs (Apache-2.0)

```bash
SHA=05e75a9c63849453e538b624ae713c2331831aa3   # update to current upstream main
BASE=https://raw.githubusercontent.com/varfish-org/hgvs-rs/$SHA
DEST=gtars-vrs/tests/data/hgvs/varfish/parser
curl -fsSL "$BASE/tests/data/parser/gauntlet" -o "$DEST/gauntlet"
curl -fsSL "$BASE/tests/data/parser/reject"   -o "$DEST/reject"
```

## biocommons/hgvs (Apache-2.0)

```bash
SHA=fa02ca5ec0c9e28a1a6377b551e41ae77b450ce9   # update to current upstream main
BASE=https://raw.githubusercontent.com/biocommons/hgvs/$SHA
DEST=gtars-vrs/tests/data/hgvs/biocommons
curl -fsSL "$BASE/tests/data/grammar_test.tsv" -o "$DEST/grammar_test.tsv"
```

After refreshing, also update the commit SHA in the sibling `NOTICE` file
and run `cargo test -p gtars-vrs --test test_hgvs_corpus_parser` to spot
any new corpus rows the local parser does not yet handle. Add new
unsupported lines to `known_skips.toml` rather than editing the upstream
files.

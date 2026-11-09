#!/usr/bin/env bash
# Write the jungle query set (for comparing origin runs) to stdout: the four
# small_jungle queries (human chrM, mouse chrM, EBV, random 90 bp) plus 1 kb slices of
# hg38 chr21, an hg38 alt contig, an hs38d1 decoy, hg19 chr21, hg19 chr6 (HLA region),
# mm10 chr19 and GRCm39 chr19.
#   jungle_queries.sh STORE > queries.fa   (WORK, IMAGE as in make_origin_fixtures.sh)
set -euo pipefail
STORE=$1
WORK=${WORK:-$HOME/work/gtars-rgblast}
IMAGE=${IMAGE:-ncbi/blast:2.17.0}
HERE=$(cd "$(dirname "$0")/.." && pwd)
TMP=$(mktemp -d "$WORK/queries-XXXXXX")
trap 'rm -rf "$TMP"' EXIT
# name, sequence digest, range
SLICES="q_hg38_chr21_20M GXCqkr6wwzuwvHIzIV9oSL8FQWHf01gu 20000001-20001000
q_hg38_alt_GL000250 pTtOib7jGZL9lETO295VEWn8jZ9RYCvj 2000001-2001000
q_hs38d1_decoy --gC85hq68PrM8BDnoC-Sonh0NdyXVIv 301-1300
q_hg19_chr21_20M LpTaNW-hwuY_yARP0rtarCnpCQLkgVCg 20000001-20001000
q_hg19_chr6_hla KqaUhJMW3CDjhoVtBetdEKT1n6hM-7Ek 29910001-29911000
q_mm10_chr19_30M mvQ1Z25vmFfx78hVjsgTJjCLadvZm6Ea 30000001-30001000
q_grcm39_chr19_30M MOUSE39 30000001-30001000"
M39=$(awk -F'\t' '$1 == "NC_000085.7" {print $4}' "$STORE/collections/fXjBOJjw-DYsSnnfDBl5vtZu1N7lbnUl.rgsi")
SLICES=${SLICES/MOUSE39/$M39}
echo "$SLICES" | awk '{print $2}' > "$TMP/digests.txt"
"$WORK/target/release/rgmakeblastdb" build --store "$STORE" --all --sequences-file "$TMP/digests.txt" \
  --out "$TMP/db" --quiet > /dev/null 2>&1
cat "$HERE/tests/data/fixtures/origin/small_jungle_queries.fa"
echo "$SLICES" | while read -r name digest range; do
  docker run --rm -u "$(id -u):$(id -g)" -v "$TMP:$TMP" -w "$TMP" "$IMAGE" \
    blastdbcmd -db db -entry "SQ.$digest" -range "$range" | sed "1s/.*/>$name/"
done

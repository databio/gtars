#!/usr/bin/env bash
# Regenerate tests/data/fixtures/origin/small_jungle* from a local jungle mirror:
# a database of the jungle sequences up to 200 kb (all 79 collections in the deflines), four
# queries (human chrM 2001-3000, mouse chrM 5001-6000, EBV 100001-101000, a random 90 bp),
# blastn -outfmt 7 hits through 'rgblast search --engine blastn' on that database (blastn in
# docker through a shim script), the same search under an mm10 gnl| seqid
# list and under an SQ. seqid list, and the sidecar rows those hits need.
#
#   WORK (default ~/work/gtars-rgblast), STORES (default ~/work/refget-blast/stores),
#   IMAGE (default ncbi/blast:2.17.0); binaries from $WORK/target/release.
set -euo pipefail
WORK=${WORK:-$HOME/work/gtars-rgblast}
STORES=${STORES:-$HOME/work/refget-blast/stores}
IMAGE=${IMAGE:-ncbi/blast:2.17.0}
STORE=$STORES/jungle
HERE=$(cd "$(dirname "$0")/.." && pwd)
OUT=$HERE/tests/data/fixtures/origin
RG=$WORK/target/release/rgblast
MK=$WORK/target/release/rgmakeblastdb
TMP=$(mktemp -d "$WORK/fixtures-XXXXXX")
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$OUT"
cd "$TMP"

docker_blast() {
  docker run --rm -u "$(id -u):$(id -g)" -v "$TMP:$TMP" -w "$TMP" "$IMAGE" "$@"
}

# sequences up to 200 kb
awk -F'\t' 'NR > 1 && $2 <= 200000 {print $4}' "$STORE/sequences.rgsi" > small.txt
"$MK" build --store "$STORE" --all --sequences-file small.txt --out "$TMP/db" --quiet --cache-dir "$TMP/cache"

# queries
q() { docker_blast blastdbcmd -db db -entry "$1" -range "$2" | sed "1s/.*/>$3/"; }
{
  q 'gnl|hg38-p14-ucsc|chrM' 2001-3000 q_human_chrM_2001
  q 'gnl|mm10_ucsc|chrM' 5001-6000 q_mouse_chrM_5001
  q 'gnl|hg19_full_analysis|chrEBV' 100001-101000 q_ebv_100001
  echo '>q_random'
  awk 'BEGIN { srand(42); s = ""; for (i = 0; i < 90; i++) s = s substr("ACGT", int(rand() * 4) + 1, 1); print s }'
} > queries.fa

printf '#!/bin/sh\nexec docker run --rm -u %s -v "%s:%s" -w "$PWD" %s blastn "$@"\n' \
  "$(id -u):$(id -g)" "$TMP" "$TMP" "$IMAGE" > blastn
chmod +x blastn
search() {
  "$RG" search --engine blastn --db "$TMP/db" --blastn "$TMP/blastn" --quiet "$@"
}
search -- -query queries.fa -outfmt 7 -out full.tsv
DB=$TMP/db

# the mm10 gnl| view and the SQ. view of every subject hit
awk -F'\t' '$3 == "99TjKCwZJJLjpBqkLpTgC2E_Y3OgKtMz" {print "gnl|mm10-ucsc-initial-soft-masked|" $5}' "$DB.members.tsv" > mm10.txt
grep -v '^#' full.tsv | cut -f2 | sort -u > sq.txt
search -- -query queries.fa -outfmt 7 -seqidlist mm10.txt -out mm10view.tsv
search -- -query queries.fa -outfmt 7 -seqidlist sq.txt -out sqview.tsv

# BLAST version and database paths vary; drop the '# BLASTN' and '# Database' lines
strip() { grep -v -e '^# BLASTN' -e '^# Database' "$1"; }
strip full.tsv > "$OUT/small_jungle_hits.tsv"
strip mm10view.tsv > "$OUT/small_jungle_hits.mm10view.tsv"
strip sqview.tsv > "$OUT/small_jungle_hits.sqview.tsv"
cp queries.fa "$OUT/small_jungle_queries.fa"

# sidecars: every collection; members of the digests in the hits
cp "$DB.collections.tsv" "$OUT/small_jungle.collections.tsv"
sed 's/^SQ\.//' sq.txt > digests.txt
awk -F'\t' 'NR == FNR {d[$1] = 1; next} /^#/ || ($2 in d)' digests.txt "$DB.members.tsv" > "$OUT/small_jungle.members.tsv"
wc -l "$OUT"/small_jungle*

# The jungle-shaped membership table for the representative-rule test: every collection
# (digest, display, taxid, N fraction) and its rows as (collection index, sequence index,
# length), sequence digests replaced by integers.
"$MK" nstats --store "$STORE" --all > nstats.tsv 2> /dev/null
python3 - "$STORE" nstats.tsv "$OUT" <<'EOF'
import gzip, json, os, sys
store, nstats, out = sys.argv[1:]
nf = {}
for l in open(nstats):
    if l.startswith('#'):
        continue
    f = l.rstrip('\n').split('\t')
    nf[f[0]] = (f[1], f[6])
colls = sorted(nf)
seq_ix = {}
with open(f"{out}/jungle_collections.tsv", "w") as c, gzip.open(f"{out}/jungle_membership.tsv.gz", "wt") as m:
    c.write("#coll\tdigest\tdisplay\ttaxid\tn_fraction\n")
    m.write("#coll\tseq\tlength\n")
    for i, d in enumerate(colls):
        fhr = json.load(open(f"{store}/fhr/{d}.fhr.json"))
        taxid = fhr["taxon"]["uri"].rsplit(":", 1)[1]
        c.write(f"{i}\t{d}\t{nf[d][0]}\t{taxid}\t{nf[d][1]}\n")
        for l in open(f"{store}/collections/{d}.rgsi"):
            if l.startswith('#'):
                continue
            f = l.rstrip('\n').split('\t')
            s = seq_ix.setdefault(f[3], len(seq_ix))
            m.write(f"{i}\t{s}\t{f[1]}\n")
EOF

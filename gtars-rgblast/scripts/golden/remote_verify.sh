#!/usr/bin/env bash
# Remote build of hg38 from S3, compared byte for byte with a local build of the same collection,
# then `rgmakeblastdb verify` against the remote store and md5 of every sequence through blastdbcmd.
# usage: remote_verify.sh   (env: WORK, STORES, B, BLAST_IMAGE, CACHE)
# CACHE is the gtars cache root for the remote store's metadata (default $WORK/refget-cache);
# downloaded sequences are not kept (no --resume).
set -u
WORK=${WORK:-$HOME/work/gtars-rgblast}
STORES=${STORES:-$HOME/work/refget-blast/stores}
B=${B:-$WORK/target/release/rgmakeblastdb}
IMG=${BLAST_IMAGE:-ncbi/blast:2.17.0}
export REFGET_CACHE_DIR=${CACHE:-$WORK/refget-cache}
HERE=$(cd "$(dirname "$0")" && pwd)
R=https://refgenie.s3.us-east-1.amazonaws.com/refget-store/jungle
C=EiFob05aCWgVU_B_Ae0cypnQut3cxUP1
mkdir -p $WORK/hg38 && cd $WORK/hg38
rm -rf local remote && mkdir local remote
SOURCE_DATE_EPOCH=1794139200 $B build --store $STORES/jungle --collection $C --title hg38 --out local/hg38 -q
SOURCE_DATE_EPOCH=1794139200 $B build --store $R --collection $C --title hg38 --out remote/hg38 --fetch-jobs 32 --verbose 2>&1 | grep -v '^sequences '
echo "== 1. byte comparison, remote vs local"
for f in $(ls local); do
  if cmp -s local/$f remote/$f; then echo "$f: identical"; else echo "$f: DIFFERS in $(cmp -l local/$f remote/$f 2>/dev/null | wc -l) bytes"; fi
done
echo "== 2. rgmakeblastdb verify against the remote store"
$B verify --db remote/hg38 --store $R
echo "== 3. md5 of every sequence"
python3 - <<'PY'
import os
oids = {}
for l in open('remote/hg38.ids.tsv'):
    if l.startswith('#'): continue
    i, oid, d = l.rstrip('\n').split('\t')
    if i.startswith('SQ.'): oids[d] = oid
with open('expected_md5.tsv', 'w') as f:
    for l in open(os.environ.get('STORES', os.environ['HOME'] + '/work/refget-blast/stores') + '/jungle/sequences.rgsi'):
        if l.startswith('#'): continue
        x = l.rstrip('\n').split('\t')
        if x[3] in oids: f.write(f"{oids[x[3]]}\t{x[4]}\n")
PY
docker run --rm -u $(id -u):$(id -g) -v $PWD:/w -w /w $IMG blastdbcmd -db remote/hg38 -entry all -outfmt '%o %s' | uv run -q --script $HERE/md5check.py expected_md5.tsv
rm -f expected_md5.tsv

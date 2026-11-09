#!/usr/bin/env bash
# hg38 (jungle collection EiFob05aCWgVU_B_Ae0cypnQut3cxUP1): makeblastdb on an exported FASTA vs
# rgmakeblastdb at 1 and 16 threads, 3 interleaved repetitions.
# usage: bench_hg38.sh   (env: WORK, STORES, B, REPS, BLAST_IMAGE; hg38.fa must exist in $WORK/hg38 for the
#        makeblastdb runs: refget store export --collection EiFob05aCWgVU_B_Ae0cypnQut3cxUP1 -o hg38.fa)
set -u
WORK=${WORK:-$HOME/work/gtars-rgblast}
STORES=${STORES:-$HOME/work/refget-blast/stores}
B=${B:-$WORK/target/release/rgmakeblastdb}
IMG=${BLAST_IMAGE:-ncbi/blast:2.17.0}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p $WORK/hg38 && cd $WORK/hg38
U="$(id -u):$(id -g)"
DK="docker run --rm -u $U -v $PWD:/w -w /w -v /usr/bin/time:/usr/bin/time:ro --entrypoint /usr/bin/time $IMG -f 'TIME %e %M'"
S="--store $STORES/jungle --collection EiFob05aCWgVU_B_Ae0cypnQut3cxUP1 --title hg38"
T="/usr/bin/time -f 'TIME %e %M'"
JOBS=(
 "rgmakeblastdb_t1::rm -rf mine::$T $B build $S --out mine/hg38 --threads 1 -q"
 "rgmakeblastdb_t16::rm -rf mine::$T $B build $S --out mine/hg38 --threads 16 -q"
)
if [ -f hg38.fa ]; then
  JOBS+=("mkdb_v5_parse::rm -rf ref; mkdir ref::$DK makeblastdb -in hg38.fa -dbtype nucl -out ref/hg38 -blastdb_version 5 -parse_seqids -title hg38")
fi
python3 $HERE/bench.py bench_hg38.jsonl ${REPS:-3} "${JOBS[@]}"
python3 $HERE/summarize.py bench_hg38.jsonl
echo BENCH_DONE

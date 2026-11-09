#!/usr/bin/env bash
# Whole jungle store: FASTA export + makeblastdb vs rgmakeblastdb --all.
# usage: bench_jungle.sh   (env: WORK, STORES, B, REPS, BLAST_IMAGE)
# The makeblastdb input comes from the refget CLI's whole-store export (export plan):
#   refget store export --all -l '{ids}' (the label template that reproduces rgmakeblastdb's deflines)
set -u
WORK=${WORK:-$HOME/work/gtars-rgblast}
STORES=${STORES:-$HOME/work/refget-blast/stores}
B=${B:-$WORK/target/release/rgmakeblastdb}
IMG=${BLAST_IMAGE:-ncbi/blast:2.17.0}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p $WORK/jungle && cd $WORK/jungle
U="$(id -u):$(id -g)"
DK="docker run --rm -u $U -v $PWD:/w -w /w -v /usr/bin/time:/usr/bin/time:ro --entrypoint /usr/bin/time $IMG -f 'TIME %e %M'"
S="--store $STORES/jungle --all --title jungle"
T="/usr/bin/time -f 'TIME %e %M'"
python3 $HERE/bench.py bench_jungle.jsonl ${REPS:-3} \
 "export_fasta::rm -f jungle.fa::$T refget store export --all --store $STORES/jungle -o jungle.fa" \
 "mkdb_v5_parse::rm -rf ref; mkdir ref::$DK makeblastdb -in jungle.fa -dbtype nucl -parse_seqids -blastdb_version 5 -out ref/jungle -title jungle" \
 "rgmakeblastdb_t16::rm -rf mine::$T $B build $S --out mine/jungle --threads 16 -q"
python3 $HERE/summarize.py bench_jungle.jsonl
echo BENCH_DONE

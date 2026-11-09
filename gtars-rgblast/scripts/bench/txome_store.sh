#!/usr/bin/env bash
# Transcriptome-scale store for manual benchmarks (plan 8 step 34): GENCODE v47 transcripts
# (385,659 records, 384,354 unique sequences, 562 Mbp) as a standalone RefgetStore, and the
# inputs the build and search benchmarks use.
#
# usage: txome_store.sh   (env: WORK, GTARS)
#   WORK   output folder (default ~/work/gtars-rgblast/bench/txome)
#   GTARS  a gtars binary with `refget build` (default: gtars on PATH; build one with
#          `cargo build --release -p gtars-cli --no-default-features --features refget`)
#
# Writes into $WORK:
#   dl/gencode.v47.transcripts.fa.gz  the download (102 MB)
#   txome.fa        the transcripts with deflines cut to the ENST id
#   standalone/     the store (one collection, alias name:gencode_v47)
#   tx_concat.fa    the same bases as 113 records (per-subject overhead comparisons)
#   q_hits.fa       100 exact 500-bp slices of transcripts (seeded)
#   q_nohit.fa      20 random 500-bp sequences (seeded)
#   q_isoform.fa    20 300-bp slices shared by the most isoforms of the 20 genes with most isoforms
# The packed store (packed/) joins when the packed layout lands (plan 8 section C).
#
# The makeblastdb reference database (reference_fasta::write output, then makeblastdb) and the
# byte-identity check come from an ignored golden test:
#   GTARS_RGBLAST_STORE_TXOME=$WORK/standalone GTARS_RGBLAST_WORK=$WORK/golden \
#     cargo test -p gtars-rgblast --release --features golden --test golden -- --ignored golden_txome_gencode
set -euo pipefail
WORK=${WORK:-$HOME/work/gtars-rgblast/bench/txome}
GTARS=${GTARS:-gtars}
URL=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.transcripts.fa.gz
mkdir -p "$WORK/dl"
cd "$WORK"
[ -s dl/gencode.v47.transcripts.fa.gz ] || curl -sS -fL -o dl/gencode.v47.transcripts.fa.gz "$URL"

python3 -I - dl/gencode.v47.transcripts.fa.gz "$WORK" <<'EOF'
import gzip, random, sys
src, out = sys.argv[1], sys.argv[2]
names, genes, seqs = [], [], []
cur = []
with gzip.open(src, "rt") as f:
    for line in f:
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
            cur = []
            fields = line[1:].strip().split("|")
            names.append(fields[0])
            genes.append(fields[1] if len(fields) > 1 else fields[0])
        else:
            cur.append(line.strip().upper())
if cur:
    seqs.append("".join(cur))
assert len(names) == len(seqs) == len(set(names)), "ENST ids must be unique"

def write(path, recs):
    with open(path, "w") as f:
        for name, s in recs:
            f.write(f">{name}\n")
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + "\n")

write(f"{out}/txome.fa", zip(names, seqs))
total = sum(map(len, seqs))
target, recs, buf, acc, k = total / 113, [], [], 0, 1
for s in seqs:
    buf.append(s)
    acc += len(s)
    if acc >= k * target and k < 113:
        recs.append((f"concat{k}", "".join(buf)))
        buf, k = [], k + 1
recs.append((f"concat{k}", "".join(buf)))
write(f"{out}/tx_concat.fa", recs)
rnd = random.Random(20261108)
long = [i for i, s in enumerate(seqs) if len(s) >= 500]
hits = []
for q in range(100):
    i = rnd.choice(long)
    p = rnd.randrange(len(seqs[i]) - 499)
    hits.append((f"hit{q}_{names[i]}_{p}", seqs[i][p:p + 500]))
write(f"{out}/q_hits.fa", hits)
write(f"{out}/q_nohit.fa", [(f"nohit{q}", "".join(rnd.choice("ACGT") for _ in range(500))) for q in range(20)])
by_gene = {}
for i, g in enumerate(genes):
    by_gene.setdefault(g, []).append(i)
iso = []
for g in sorted(by_gene, key=lambda g: (-len(by_gene[g]), g))[:20]:
    members = [seqs[i] for i in by_gene[g]]
    first = max(members, key=len)
    best = max(range(0, max(1, len(first) - 299), 100),
               key=lambda p: (sum(first[p:p + 300] in m for m in members), -p))
    iso.append((f"iso_{g}_{best}", first[best:best + 300]))
write(f"{out}/q_isoform.fa", iso)
print(f"transcripts {len(seqs)}  bases {total}  concat records {len(recs)}")
EOF

if [ ! -d standalone ]; then
  "$GTARS" refget build --collection-alias name:gencode_v47 -o standalone txome.fa
fi
echo TXOME_STORE_DONE

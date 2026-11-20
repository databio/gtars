#!/usr/bin/env bash
# Smoke test of remote stores against the public S3 stores (loose layout) and a local
# Range-capable server (packed layout): export, cache reuse, XDG fallback, ranged downloads and
# resume after an interrupt, rgblastn == blastn on a remote selection, rgblast search/origin, S3's
# 403-for-missing, offline reopen from the cache, and a packed store served over HTTP.
# Downloads about 200 MB (hg38 chr1, chr21, chr22) into a fresh cache under $WORK.
# usage: remote_smoke.sh   (env: WORK, BIN, BLAST_BIN, BLAST_IMAGE, INT_AFTER, EXPORT_INT_AFTER, PY, PORT)
#   BIN        directory with gtars, rgmakeblastdb, rgblastn, rgblast
#   BLAST_BIN  directory with BLAST+ 2.17.0 blastn; without it, docker image BLAST_IMAGE
#   INT_AFTER  seconds after the build's first fetch line before the interrupt in step 4;
#              EXPORT_INT_AFTER the same for the export, counted from its start
#   PY         a python with gtars installed, for step 8 (skipped if it cannot import gtars)
set -euo pipefail
WORK=${WORK:-$HOME/work/rgblast-release/smoke}
BIN=${BIN:-$HOME/work/rgblast-release/target/release}
IMG=${BLAST_IMAGE:-ncbi/blast:2.17.0}
INT_AFTER=${INT_AFTER:-0}
PY=${PY:-python3}
PORT=${PORT:-8089}
S3=https://refgenie.s3.us-east-1.amazonaws.com/refget-store
DEMO=$S3/demo
JUNGLE=$S3/jungle
HG38=EiFob05aCWgVU_B_Ae0cypnQut3cxUP1
G=$BIN/gtars
MK=$BIN/rgmakeblastdb

fail() { echo "FAIL: $*" >&2; exit 1; }
blastn() {
  if [ -n "${BLAST_BIN:-}" ]; then "$BLAST_BIN/blastn" "$@"
  else docker run --rm -u "$(id -u):$(id -g)" -v "$WORK:$WORK" -w "$PWD" "$IMG" blastn "$@"; fi
}
netline() { grep '^network:' "$1" | tail -1; }
cache_hits() { netline "$1" | sed -E 's/.*cache hits ([0-9]+).*/\1/'; }
net_requests() { netline "$1" | sed -E 's/^network: [0-9]+ bytes in ([0-9]+) requests.*/\1/'; }
seqfile() { find "$REFGET_CACHE_DIR/sequences" -name "$1.seq"; }
# interrupt LOG PATTERN CMD...: run CMD, send SIGINT INT_AFTER seconds after LOG shows PATTERN
# (empty PATTERN: after the start), set rc to its exit status. Job control (set -m) keeps
# SIGINT from being ignored in the background job.
interrupt() {
  local log=$1 pat=$2; shift 2
  set -m; "$@" > "$log" 2>&1 & local pid=$!; set +m
  if [ -n "$pat" ]; then
    while kill -0 "$pid" 2>/dev/null && ! grep -qs "$pat" "$log"; do sleep 0.1; done
  fi
  sleep "$INT_AFTER"; kill -INT "$pid" 2>/dev/null || true
  set +e; wait "$pid"; rc=$?; set -e
}

rm -rf "$WORK" && mkdir -p "$WORK/db" && cd "$WORK"
export REFGET_CACHE_DIR=$WORK/cache

echo "== 1. open_remote and export of a small collection"
"$G" refget export --remote "$DEMO" --all --resume --output demo.fa
nrec=$(grep -c '^>' demo.fa)
nseq=$(find "$REFGET_CACHE_DIR/sequences" -name '*.seq' | wc -l)
[ "$nrec" -gt 0 ] && [ "$nrec" -eq "$nseq" ] || fail "export wrote $nrec records, cache has $nseq .seq files"
test -f "$REFGET_CACHE_DIR/CACHEDIR.TAG" || fail "no CACHEDIR.TAG"
ls "$REFGET_CACHE_DIR"/remotes/*/rgstore.json >/dev/null || fail "no remotes/*/rgstore.json"
echo "ok: $nrec records, $nseq cached sequences"

echo "== 2. cache reuse"
"$MK" build --store "$DEMO" --all --out db/demo --verbose > build_demo.log 2>&1
netline build_demo.log
[ "$(cache_hits build_demo.log)" -eq "$nseq" ] || fail "expected $nseq cache hits"
[ "$(net_requests build_demo.log)" -eq 0 ] || fail "sequences were downloaded again"
echo "ok"

echo "== 3. XDG fallback and URL normalization"
env -u REFGET_CACHE_DIR XDG_CACHE_HOME="$WORK/xdg" "$G" refget export --remote "$DEMO/" --all --resume --output demo2.fa
cmp demo.fa demo2.fa || fail "trailing-slash export differs"
test -f "$WORK/xdg/refget/CACHEDIR.TAG" || fail "cache not under XDG_CACHE_HOME/refget"
env -u REFGET_CACHE_DIR XDG_CACHE_HOME="$WORK/xdg" "$G" refget export --remote "$DEMO" --all --resume --output /dev/null
[ "$(ls "$WORK/xdg/refget/remotes" | wc -l)" -eq 1 ] || fail "URL with and without '/' got two remotes/ keys"
echo "ok"

echo "== 4. ranged downloads and resume after an interrupt"
curl -fsS "$JUNGLE/collections/$HG38.rgsi" > hg38.rgsi
dig() { awk -F'\t' -v n="$1" '$1==n {print $4}' hg38.rgsi; }
CHR1=$(dig chr1); CHR21=$(dig chr21); CHR22=$(dig chr22)
printf '%s\n' "$CHR21" > chr21.txt
printf '%s\n%s\n' "$CHR21" "$CHR22" > chr21_22.txt
printf '%s\n%s\n%s\n' "$CHR1" "$CHR21" "$CHR22" > three.txt
"$MK" build --store "$JUNGLE" --collection "$HG38" --sequences-file chr21_22.txt --out db/two --verbose > build_two.log 2>&1
netline build_two.log
[ -n "$(seqfile "$CHR21")" ] && [ -n "$(seqfile "$CHR22")" ] || fail "chr21/chr22 not cached"
interrupt build_three_int.log '^fetch: ' "$MK" build --store "$JUNGLE" --collection "$HG38" \
  --sequences-file three.txt --fetch-jobs 4 --out db/three --verbose
[ "$rc" -ne 0 ] || fail "build finished before the interrupt; lower INT_AFTER"
! ls db/three.* >/dev/null 2>&1 || fail "interrupted build left db/three.* behind"
[ -z "$(seqfile "$CHR1")" ] || fail "chr1 finished before the interrupt; lower INT_AFTER"
echo "interrupted (exit $rc), nothing left behind"
"$MK" build --store "$JUNGLE" --collection "$HG38" --sequences-file three.txt --fetch-jobs 4 \
  --out db/three --verbose > build_three.log 2>&1
netline build_three.log
[ "$(cache_hits build_three.log)" -eq 2 ] || fail "resume: expected 2 cache hits (chr21, chr22)"
[ -n "$(seqfile "$CHR1")" ] || fail "chr1 not cached after resume"
[ -z "$(ls -A "$REFGET_CACHE_DIR/tmp")" ] || fail "stale files in cache/tmp: $(ls -A "$REFGET_CACHE_DIR/tmp")"
echo "ok: resume fetched only chr1"
# The same through `gtars refget export`: drop chr1 from the cache, interrupt, resume.
rm -f "$(seqfile "$CHR1")"
INT_AFTER=${EXPORT_INT_AFTER:-1.5} interrupt export_int.log '' "$G" refget export --remote "$JUNGLE" --collection "$HG38" \
  --names chr1 chr21 chr22 --label '{name}' --resume --fetch-jobs 4 --output three.fa
[ "$rc" -ne 0 ] || fail "export finished before the interrupt; lower EXPORT_INT_AFTER"
[ -z "$(seqfile "$CHR1")" ] || fail "export finished chr1 before the interrupt; lower EXPORT_INT_AFTER"
"$G" refget export --remote "$JUNGLE" --collection "$HG38" --names chr1 chr21 chr22 --label '{name}' \
  --resume --fetch-jobs 4 --output three.fa > export.log 2>&1
netline export.log
[ "$(cache_hits export.log)" -eq 2 ] || fail "export resume: expected 2 cache hits"
[ -z "$(ls -A "$REFGET_CACHE_DIR/tmp")" ] || fail "stale files in cache/tmp after export"
# Every record against the md5 column of the collection's rgsi.
"$PY" -I - three.fa hg38.rgsi <<'PY' || fail "md5 check failed"
import hashlib, sys
want = {}
for line in open(sys.argv[2]):
    if not line.startswith("#"):
        x = line.rstrip("\n").split("\t")
        want[x[0]] = x[4]
def check(name, h):
    ok = h.hexdigest() == want[name]
    print(f"{name} md5 {'ok' if ok else 'MISMATCH'}")
    return ok
good, name, h = True, None, None
with open(sys.argv[1], "rb") as f:
    for line in f:
        if line.startswith(b">"):
            if name:
                good &= check(name, h)
            name, h = line[1:].split()[0].decode(), hashlib.md5()
        else:
            h.update(line.rstrip(b"\n").upper())
if name:
    good &= check(name, h)
sys.exit(0 if good else 1)
PY
[ "$(grep -c '^>' three.fa)" -eq 3 ] || fail "three.fa should hold 3 records"

echo "== 5. rgblastn and rgblast with --store <URL>"
"$PY" -I - three.fa > q.fa <<'PY'
import sys
seqs, name = {}, None
for line in open(sys.argv[1]):
    line = line.strip()
    if line.startswith(">"):
        name = line[1:].split()[0]; seqs[name] = []
    else:
        seqs[name].append(line)
seqs = {k: "".join(v) for k, v in seqs.items()}
q = list(seqs["chr21"][30_000_000:30_000_600])
for i in (100, 300, 500):
    q[i] = {"A": "C", "C": "G", "G": "T", "T": "A"}.get(q[i], "A")
rc = seqs["chr22"][25_000_000:25_000_400][::-1].translate(str.maketrans("ACGTN", "TGCAN"))
print(">q21\n" + "".join(q) + "\n>q22rc\n" + rc)
PY
"$MK" build --store "$JUNGLE" --collection "$HG38" --sequences-file chr21.txt --out db/chr21 --verbose > build_chr21.log 2>&1
[ "$(cache_hits build_chr21.log)" -eq 1 ] && [ "$(net_requests build_chr21.log)" -eq 0 ] || fail "chr21 build downloaded again"
"$BIN/rgblastn" --store "$JUNGLE" --collection "$HG38" --sequences-file chr21.txt -query q.fa -outfmt 6 -out rg.tsv
blastn -db db/chr21 -query q.fa -outfmt 6 -out bl.tsv
cmp rg.tsv bl.tsv || fail "rgblastn and blastn differ"
grep -q '^q21' rg.tsv || fail "no q21 hit"
echo "ok: rgblastn == blastn ($(wc -l < rg.tsv) hits)"
"$BIN/rgblast" search --store "$JUNGLE" --collection "$HG38" --sequences-file chr21.txt --quiet \
  -- -query q.fa -outfmt 6 -out hits.tsv
grep -q '^q21' hits.tsv || fail "rgblast search: no q21 hit"
"$BIN/rgblast" origin --store "$JUNGLE" --collection "$HG38" --sequences-file chr21.txt --quiet \
  --query q.fa --output origin.tsv
grep -q '^q21' origin.tsv || fail "rgblast origin: no q21 line"
echo "ok"

echo "== 6. S3 error semantics (403 for a missing key)"
set +e
"$G" refget export --remote "$JUNGLE" --collection AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA --output x.fa > missing.log 2>&1
rc=$?
"$G" refget export --remote "$S3/no-such-store" --all --output x.fa > missing_store.log 2>&1
rc2=$?
set -e
[ "$rc" -ne 0 ] && grep -qi 'not found' missing.log || fail "missing collection: $(head -c 300 missing.log)"
[ "$rc2" -ne 0 ] || fail "missing store succeeded"
! grep -qi 'retr' missing_store.log || fail "missing store was retried: $(head -c 300 missing_store.log)"
echo "ok: $(head -c 200 missing_store.log)"

echo "== 7. offline reopen from the cache"
docker run --rm --network none -u "$(id -u):$(id -g)" -v "$BIN:/b:ro" -v "$WORK:$WORK" -w "$WORK" \
  -e REFGET_CACHE_DIR ubuntu:24.04 /b/gtars refget export --remote "$DEMO" --all --resume --output demo3.fa
cmp demo.fa demo3.fa || fail "offline export differs"
echo "ok"

echo "== 8. shared cache across bindings"
if "$PY" -I -c 'import gtars.refget' 2>/dev/null; then
  "$PY" -I - "$REFGET_CACHE_DIR" "$JUNGLE" "$CHR21" <<'PY'
import sys
from gtars.refget import RefgetStore
cache, url, d = sys.argv[1:]
s = RefgetStore.open_remote(cache, url, fetch_jobs=8)
r = s.cache_sequences([d])
hits = r["cache_hits"] if isinstance(r, dict) else r.cache_hits
assert hits == 1, r
print("ok: python saw the CLI's cached chr21")
PY
else
  echo "skipped: $PY cannot import gtars"
fi

echo "== 9. packed store over a Range-capable server (local nginx)"
"$G" refget build --packed --output "$WORK/packed" three.fa demo.fa > /dev/null
docker rm -f rg-packed > /dev/null 2>&1 || true
docker run -d --name rg-packed -p "127.0.0.1:$PORT:80" -v "$WORK/packed:/usr/share/nginx/html:ro" nginx:alpine > /dev/null
trap 'docker rm -f rg-packed > /dev/null 2>&1 || true' EXIT
for _ in $(seq 50); do curl -fsS "http://127.0.0.1:$PORT/rgstore.json" > /dev/null 2>&1 && break; sleep 0.2; done
"$G" refget export --store "$WORK/packed" --all --label '{name}' --output local.fa
REFGET_CACHE_DIR=$WORK/cache-packed "$G" refget export --remote "http://127.0.0.1:$PORT" --all --label '{name}' --output remote.fa
cmp local.fa remote.fa || fail "packed remote export differs from local export"
docker logs rg-packed 2>/dev/null > nginx.log
grep -q 'packs/.* 206 ' nginx.log || fail "no 206 responses for packs/"
grep -q 'GET /sequences.pack.idx' nginx.log || fail "no GET of sequences.pack.idx"
printf '%s\n' "$CHR21" > packed_chr21.txt
"$BIN/rgblastn" --store "$WORK/packed" --all --sequences-file packed_chr21.txt \
  -query q.fa -outfmt 6 -out rg_packed_local.tsv
REFGET_CACHE_DIR=$WORK/cache-packed "$BIN/rgblastn" --store "http://127.0.0.1:$PORT" --all \
  --sequences-file packed_chr21.txt -query q.fa -outfmt 6 -out rg_packed.tsv
grep -q '^q21' rg_packed.tsv || fail "rgblastn on the packed remote: no q21 hit"
cmp rg_packed_local.tsv rg_packed.tsv || fail "rgblastn: packed remote differs from packed local"
echo "ok: packed remote export == local, $(grep -c ' 206 ' nginx.log) ranged responses, rgblastn remote == local"

echo "ALL OK"

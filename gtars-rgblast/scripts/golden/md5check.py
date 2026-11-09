# /// script
# dependencies = []
# ///
"""Check sequence md5s read back with blastdbcmd against expected values.

usage: blastdbcmd -db DB -entry all -outfmt '%o %s' | uv run --script md5check.py expected.tsv
expected.tsv: key \\t md5 (key = OID, or whatever the first output column is).
Reports matches, mismatches and expected-but-missing entries.
"""
import hashlib
import sys

exp = {}
for l in open(sys.argv[1]):
    i, m = l.rstrip("\n").split("\t")[:2]
    exp[i] = m
seen = set()
bad = []
n = 0
bases = 0
for line in sys.stdin.buffer:
    i, _, s = line.rstrip(b"\n").partition(b" ")
    i = i.decode()
    if i in seen:
        continue
    seen.add(i)
    n += 1
    bases += len(s)
    m = hashlib.md5(s.upper()).hexdigest()
    if exp.get(i) != m:
        bad.append((i, m, exp.get(i)))
missing = [i for i in exp if i not in seen]
print(f"checked {n} sequences, {bases} bases; md5 mismatches {len(bad)}; expected-but-missing {len(missing)}")
for b in bad[:10]:
    print("MISMATCH", b)
for m in missing[:10]:
    print("MISSING", m)

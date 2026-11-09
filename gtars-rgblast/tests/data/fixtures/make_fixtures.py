#!/usr/bin/env python3
"""Write the gtars-rgblast test fixtures (deterministic; run from this directory).

usage: python3 make_fixtures.py
Every output is a pure function of the seeds below; re-running gives identical files.
asmD's single name is found by a search loop so its collection digest starts with '-'.
"""
import base64, hashlib, json, random

rnd = random.Random(20261108)


def acgt(n):
    return "".join(rnd.choice("ACGT") for _ in range(n))


def sprinkle(seq, letters, every):
    s = list(seq)
    for i in range(every // 2, len(s), every):
        s[i] = letters[(i // every) % len(letters)]
    return "".join(s)


def fasta(path, records):
    with open(path, "w") as f:
        for name, desc, seq in records:
            f.write(f">{name}{' ' + desc if desc else ''}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")


def sha512t24u(b):
    return base64.urlsafe_b64encode(hashlib.sha512(b).digest()[:24]).decode()


def canon(x):
    return json.dumps(x, separators=(",", ":"), sort_keys=True, ensure_ascii=False).encode()


def seqcol_digest(records):
    """Top-level seqcol digest: the inherent attributes are names and sequences."""
    names = [r[0] for r in records]
    seqs = ["SQ." + sha512t24u(r[2].upper().encode()) for r in records]
    level1 = {k: sha512t24u(canon(v)) for k, v in (("names", names), ("sequences", seqs))}
    return sha512t24u(canon(level1))


chr1 = acgt(1000)
chr2 = ""
for run in [1, 15, 16, 4095, 4096, 4097]:
    chr2 += acgt(50) + "N" * run
chr2 += acgt(50)
chrM = sprinkle(acgt(400), "RYU", 37)
iupac1 = sprinkle(acgt(300), "MKSWBDHVN", 23)

fasta("asmA.fa", [
    ("chr1", "", chr1),
    ("chr2", "N runs 1 15 16 4095 4096 4097", chr2),
    ("chrM", "", chrM),
    ("iupac1", "", iupac1),
    ("ENA|KI270740|KI270740.1", "Homo sapiens test contig", acgt(200)),
    ("chrX;alt", "", acgt(150)),
    ("HLA-A*01:01:01:01", "", acgt(120)),
])
fasta("asmB.fa", [
    ("chr1", "", chr1),
    ("chrB1", "", acgt(1001)),
    ("chrB2", "", acgt(1002)),
    ("chrB3", "", acgt(1003)),
    ("chrUn_x", "", chr2),
])
fasta("asmC.fa", [
    ("c1", "", acgt(500)),
    ("prot1", "a protein", "MKVLAAGIVGLLLAQEEFWPRTYHHSDKLMNQPRSTVWYACDEFGHIKLMNPQRSTVWY"),
    ("c2", "", acgt(64)),
])
fasta("numeric.fa", [(n, "", acgt(40)) for n in
                     ["0", "1", "22", "00", "007", "0012", "2147483647", "2147483648", "4294967296", "+1", "-1"]])
dseq = acgt(80)
i = 0
while not seqcol_digest([(f"dash{i}", "", dseq)]).startswith("-"):
    i += 1
asm_d = [(f"dash{i}", "", dseq)]
fasta("asmD.fa", asm_d)
fasta("orphan.fa", [("orphan1", "", acgt(90))])
fasta("collide.fa", [("a|b", "", acgt(50)), ("a_b", "", acgt(50))])
with open("aliases.tsv", "w") as f:
    f.write("#kind\tnamespace\talias\tcollection\tname\n")
    for row in [
        ("sequence", "insdc", "CM000001.1", "asmA", "chr1"),
        ("sequence", "refseq", "NC_TEST.1", "asmA", "chr1"),
        ("sequence", "insdc", "123456", "asmA", "chrM"),
        ("sequence", "insdc", "0099", "asmB", "chrB1"),
        ("collection", "name", "asmC-alt", "asmC", ""),
        ("collection", "refseq", "GCF_TEST.1", "asmC", ""),
    ]:
        f.write("\t".join(row) + "\n")
with open("asmA.fhr.json", "w") as f:
    json.dump({"schema": "https://raw.githubusercontent.com/FAIRtracks/fhr-spec/main/fhr.json",
               "genome": "Testus fixturus", "version": "asmA.1"}, f, indent=2)
    f.write("\n")
print(f"asmD name dash{i}, digest {seqcol_digest(asm_d)}")

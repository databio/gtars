#!/usr/bin/env python3
"""Write the `taxa` test fixture (deterministic; run from this directory).

usage: python3 make_taxa.py
Every output is a pure function of the seed below; re-running gives identical files.
"""
import random

rnd = random.Random(20261109)


def acgt(n):
    return "".join(rnd.choice("ACGT") for _ in range(n))


def fasta(path, records):
    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")


h_chr1 = acgt(300)
h_chr2 = acgt(280)
h_chr3 = acgt(260)
h_chrm = acgt(200)
decoy = acgt(250)
b_chr2b = acgt(240)
# human_masked chr1: chr1 with two thirds hard-masked (N runs), so 3-bit
masked = list(h_chr1)
for start in (10, 120, 220):
    for i in range(start, start + 66):
        masked[i] = "N"
h_chr1_masked = "".join(masked)
m_chr1 = acgt(300)
m_7 = acgt(220)
# one 4-bit (IUPAC) sequence, no N
m_iupac = "".join(rnd.choice("ACGT") if i % 9 else rnd.choice("MKSWBDHV") for i in range(150))
u1 = acgt(180)
orphan = acgt(120)

fasta("human_a.fa", [("chr1", h_chr1), ("chr2", h_chr2), ("chr3", h_chr3), ("chrM", h_chrm)])
fasta("human_b.fa", [("chr1", h_chr1), ("chrM", h_chrm), ("decoy", decoy), ("chr2b", b_chr2b)])
fasta("human_masked.fa", [("chr1", h_chr1_masked), ("chr2", h_chr2)])
fasta("mouse_a.fa", [("chr1", m_chr1), ("decoy", decoy), ("7", m_7), ("chrUn_iupac", m_iupac)])
fasta("nofhr.fa", [("u1", u1)])
fasta("orphan.fa", [("o1", orphan)])

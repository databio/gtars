#!/usr/bin/env python3
"""Demo of FAI indexing with gtars-refget"""

import gtars.refget as refget

print("=== Demo 1: Fast FAI-only computation ===")
test_file = "tests/data/fasta/base.fa"
print(f"Computing FAI for: {test_file}")

fai_records = refget.compute_fai(test_file)
print(f"Found {len(fai_records)} sequences\n")

for record in fai_records:
    print(f"{record.name}: {record.length} bp")
    if record.fai:
        print(f"  offset: {record.fai.offset}")
        print(f"  line_bases: {record.fai.line_bases}")
        print(f"  line_bytes: {record.fai.line_bytes}")
    else:
        print(f"  FAI: None (gzipped file)")
    print()

print("\n=== Demo 2: digest_fasta() also computes FAI ===")
print(f"Processing: {test_file}")
seqcol = refget.digest_fasta(test_file)
print(f"Found {len(seqcol.sequences)} sequences\n")

for seq in seqcol.sequences[:2]:  # Show first 2
    print(f"{seq.metadata.name}:")
    print(f"  Length: {seq.metadata.length}")
    print(f"  SHA512: {seq.metadata.sha512t24u}")
    if seq.fai:
        print(f"  FAI offset: {seq.fai.offset}")
        print(f"  FAI line_bases: {seq.fai.line_bases}")
        print(f"  FAI line_bytes: {seq.fai.line_bytes}")
    print()

print("\n=== Demo 3: Gzipped files return None for FAI ===")
gzipped_file = "tests/data/fasta/base.fa.gz"
print(f"Computing FAI for: {gzipped_file}")
fai_records = refget.compute_fai(gzipped_file)
for record in fai_records:
    print(f"{record.name}: fai={record.fai}")

print("\n=== Done ===")

#!/usr/bin/env python3
"""Compute VRS Allele IDs from a gnomAD VCF using a real reference genome.

Usage:
    python gnomad_vrs.py <fasta_path> <vcf_path>
"""

import sys
import time

from gtars.refget import RefgetStore


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <fasta_path> <vcf_path>", file=sys.stderr)
        sys.exit(1)

    fasta_path = sys.argv[1]
    vcf_path = sys.argv[2]

    # Step 1: Load reference genome
    print(f"Loading FASTA: {fasta_path}", file=sys.stderr)
    t0 = time.time()
    store = RefgetStore.in_memory()
    store.add_sequence_collection_from_fasta(fasta_path)
    load_time = time.time() - t0
    print(f"Loaded in {load_time:.1f}s", file=sys.stderr)

    # Get collection digest
    collections = store.list_collections()
    assert len(collections) > 0, "No collections found"
    collection_digest = collections[0].digest
    print(
        f"Collection: {collection_digest} ({collections[0].n_sequences} sequences)",
        file=sys.stderr,
    )

    # Step 2: Compute VRS IDs
    print(f"Processing VCF: {vcf_path}", file=sys.stderr)
    t1 = time.time()
    results = store.compute_vrs_ids(collection_digest, vcf_path)
    vrs_time = time.time() - t1

    # Step 3: Print results
    print("chrom\tpos\tref\talt\tvrs_id")
    for r in results:
        print(f"{r['chrom']}\t{r['pos']}\t{r['ref']}\t{r['alt']}\t{r['vrs_id']}")

    total_time = time.time() - t0
    print(
        f"\nComputed {len(results)} VRS IDs in {vrs_time:.3f}s "
        f"({len(results) / vrs_time:.0f} variants/sec)",
        file=sys.stderr,
    )
    print(f"Total time: {total_time:.1f}s", file=sys.stderr)


if __name__ == "__main__":
    main()

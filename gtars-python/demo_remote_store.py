#!/usr/bin/env python3
"""
Demo: Remote RefgetStore

This script demonstrates loading a RefgetStore from a remote URL (S3)
and fetching sequences on-demand with local caching.

This complements the local RefgetStore tutorial which shows how to create
stores from FASTA files. Here we focus on accessing pre-built remote stores.

See also: refgenie-docs/docs/refget/notebooks/refgetstore.ipynb
"""

import os
import tempfile
from pathlib import Path

# Import from refget wrapper (consistent with other tutorials)
# Can also use: from gtars.refget import GlobalRefgetStore
from refget.refget_store import GlobalRefgetStore

# Remote store URL (2023 Human Pangenome Reference - 47 haplotype-resolved assemblies)
REMOTE_URL = "https://refgenie.s3.us-east-1.amazonaws.com/pangenome_refget_store"

# Persistent cache directory
CACHE_DIR = Path.home() / ".refget" / "pangenome_cache"

# Example collection from the pangenome (HG03540.pri.mat.f1_v2)
EXAMPLE_COLLECTION = "0aHV7I-94paL9Z1H4LNlqsW3WxJhlou5"
EXAMPLE_SEQ_NAME = "JAGYVX010000006.1 unmasked:primary_assembly HG03540.pri.mat.f1_v2:JAGYVX010000006.1:1:96320881:1"


def main():
    print("=" * 60)
    print("Remote RefgetStore Demo")
    print("=" * 60)

    # 1. Load remote store (fetches metadata only, ~1.5 MB)
    print(f"\n1. Loading remote store from:\n   {REMOTE_URL}")
    print(f"   Cache directory: {CACHE_DIR}\n")

    store = GlobalRefgetStore.load_remote(
        cache_path=str(CACHE_DIR),
        remote_url=REMOTE_URL
    )

    print(f"   Loaded! {len(store)} sequences available (metadata only)")

    # 2. Show store statistics
    print(f"\n2. Store statistics:")
    stats = store.stats()
    for key, value in stats.items():
        print(f"   {key}: {value}")

    # 3. List sequences (first 5)
    print(f"\n3. Listing sequences (first 5 of {len(store)}):")
    records = store.sequence_records()
    for i, rec in enumerate(records[:5]):
        m = rec.metadata
        print(f"   {i+1}. {m.name[:50]}...")
        print(f"      sha512t24u: {m.sha512t24u}")
        print(f"      length: {m.length:,} bp")

    # 4. Fetch a sequence by ID (downloads sequence data on first access)
    seq_digest = "du4GiRD_OcmdmCn_RmImyb71YZ4XoCdk"
    print(f"\n4. Get sequence record by ID (fetches from remote):")
    record = store.get_sequence_by_id(seq_digest)
    if record:
        print(f"   Name: {record.metadata.name}")
        print(f"   Length: {record.metadata.length:,} bp")
        print(f"   MD5: {record.metadata.md5}")

    # 5. Get a substring (sequence already cached from step 4)
    print(f"\n5. Get substring [0:100] (uses cached sequence):")
    sub_seq = store.get_substring(seq_digest, 0, 100)
    print(f"   Result: {sub_seq}")

    # 6. Get another substring from same sequence
    print(f"\n6. Get substring [1000:1100]:")
    sub_seq2 = store.get_substring(seq_digest, 1000, 1100)
    print(f"   Result: {sub_seq2}")

    # 7. Export sequences to FASTA file by digest
    print(f"\n7. Exporting sequences to FASTA by digest:")
    output_fasta = os.path.join(tempfile.gettempdir(), "demo_export.fa")
    digests_to_export = [
        "du4GiRD_OcmdmCn_RmImyb71YZ4XoCdk",
        "cPD3x19YSSfB_TzCKAnp1tzjOKlQVu7l",
        "d8BGSj_irEbXexaV7pStsWf9mFEZFL-8",
    ]
    print(f"   Exporting {len(digests_to_export)} sequences to: {output_fasta}")
    store.export_fasta_by_digests(digests_to_export, output_fasta, 80)

    with open(output_fasta) as f:
        for i, line in enumerate(f):
            if i >= 6:
                print("   ...")
                break
            print(f"   {line.rstrip()[:70]}")

    # ============================================================
    # Collection-based features (using remote pangenome collections)
    # ============================================================
    print("\n" + "=" * 60)
    print("Collection-Based Features")
    print("=" * 60)

    # 8. Get sequence by collection + name (triggers lazy-load of collection)
    print(f"\n8. Get sequence by collection + name:")
    print(f"   Collection: {EXAMPLE_COLLECTION}")
    print(f"   Sequence: {EXAMPLE_SEQ_NAME[:50]}...")

    record = store.get_sequence_by_collection_and_name(EXAMPLE_COLLECTION, EXAMPLE_SEQ_NAME)
    if record:
        print(f"   Found! Length: {record.metadata.length:,} bp")
        print(f"   Digest: {record.metadata.sha512t24u}")

    # 9. Batch retrieval with substrings_from_regions()
    print(f"\n9. Batch retrieval with substrings_from_regions():")
    temp_dir = tempfile.mkdtemp(prefix="refget_demo_")
    bed_path = os.path.join(temp_dir, "regions.bed")

    # Create BED file with regions from the example sequence
    bed_content = f"""{EXAMPLE_SEQ_NAME}\t0\t100
{EXAMPLE_SEQ_NAME}\t1000\t1100
{EXAMPLE_SEQ_NAME}\t50000\t50050
"""
    with open(bed_path, "w") as f:
        f.write(bed_content)

    sequences = store.substrings_from_regions(EXAMPLE_COLLECTION, bed_path)
    for seq in sequences:
        print(f"   {seq.start:,}-{seq.end:,}: {seq.sequence[:40]}...")

    # 10. Export BED regions to FASTA with export_fasta_from_regions()
    print(f"\n10. Export BED regions to FASTA with export_fasta_from_regions():")
    output_regions = os.path.join(temp_dir, "regions.fa")
    store.export_fasta_from_regions(EXAMPLE_COLLECTION, bed_path, output_regions)

    print(f"    Output: {output_regions}")
    with open(output_regions) as f:
        for i, line in enumerate(f):
            if i >= 6:
                print("    ...")
                break
            print(f"    {line.rstrip()[:65]}")

    # Summary
    print("\n" + "=" * 60)
    print("Demo complete!")
    print("=" * 60)
    print(f"\nCache directory: {CACHE_DIR}")
    print(f"Temp files: {temp_dir}")
    print("\nKey features demonstrated:")
    print("  - load_remote(): Load store from URL, fetch sequences on-demand")
    print("  - get_sequence_by_id(): Lookup by SHA-512/24u or MD5 digest")
    print("  - get_sequence_by_collection_and_name(): Lookup by sequence name")
    print("  - substrings_from_regions(): Batch retrieval from BED file")
    print("  - export_fasta_by_digests(): Export sequences by digest")
    print("  - export_fasta_from_regions(): Export BED regions to FASTA")


if __name__ == "__main__":
    main()

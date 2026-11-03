# Gtars digests using python bindings

## Introduction

You can use this to compute md5 or GA4GH sha512t24u digests of strings, or FASTA files.

## Tutorial

Computing digests for all sequences in a FASTA file:


```python
import os
import tempfile
from gtars.refget import digest_fasta

with tempfile.TemporaryDirectory() as temp_dir:
    fasta_content = (
        ">chr1 description\n"
        "ATGCATGCATGC\n"
        ">chr2\n"
        "GGGGAAAA\n"
    )
    fasta_path = os.path.join(temp_dir, "example.fa")
    with open(fasta_path, "w") as f:
        f.write(fasta_content)

    collection = digest_fasta(fasta_path)
    print(f"Collection-level digest: {collection.digest}")
    print(f"Number of sequences in collection: {len(collection)}")
    print(f"Metadata for first sequence: {collection[0].metadata.name}, Length: {collection[0].metadata.length}")
```

Compute a digest for a given sequence:

```python
from gtars.refget import sha512t24u_digest

sequence_data = "AGCT"
digest = sha512t24u_digest(sequence_data)
print(f"SHA512t24u digest for '{sequence_data}': {digest}")

```
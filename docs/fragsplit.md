# gtars-fragsplit

Split fragment files according to pseudobulk assignments for single-cell ATAC-seq data processing.

## Features

- Efficient fragment file splitting
- Pseudobulk generation from single-cell data
- Support for cell type annotations
- Gzip compression support
- Memory-efficient streaming

## Usage

### Command Line
```bash
gtars fragsplit \
  --fragments fragments.tsv.gz \
  --barcodes barcodes.txt \
  --output-dir pseudobulk/
```

### From Python
```python
from gtars import fragsplit

fragsplit.split_by_clusters(
    fragments="fragments.tsv.gz",
    clusters="clusters.csv",
    output_dir="pseudobulk/"
)
```

## Input Formats

### Fragment File
Standard 10x Genomics format:

```
chr1    1000    1500    AAACCCAAGAAACACT-1    1
```

### Barcode Assignment
```
AAACCCAAGAAACACT-1    cluster1
AAACCCAAGAAACAGC-1    cluster2
```

## Output

Creates one fragment file per cluster/pseudobulk group with appropriate cell barcodes.
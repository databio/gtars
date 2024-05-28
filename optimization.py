import scanpy as sc
import numpy as np

from rich.progress import track
from genimtools.tokenizers import TreeTokenizer
# from genimtools.utils import write_tokens_to_gtok
from geniml.io import Region

adata = sc.read_h5ad("luecken2021_parsed.h5ad")
# adata = sc.pp.subsample(adata, n_obs=10_000, copy=True)
t = TreeTokenizer("screen.bed")

adata_features = [
    Region(chr, int(start), int(end))
    for chr, start, end in track(
        zip(adata.var["chr"], adata.var["start"], adata.var["end"]),
        total=adata.var.shape[0],
        description="Extracting regions from AnnData",
    )
]

features = np.ndarray(len(adata_features), dtype=object)

for i, region in enumerate(adata_features):
    features[i] = region

del adata_features

tokenized = []
for row in track(
    range(adata.shape[0]),
    total=adata.shape[0],
    description="Tokenizing",
):
    _, non_zeros = adata.X[row].nonzero()
    regions = features[non_zeros]
    tokenized.append(t(regions))




# the issue is that
t(regions)
# is, for some reason much slower and less efficient than
t.encode(regions)
# why?
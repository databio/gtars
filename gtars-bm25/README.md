# gtars-bm25
This crate implements a BM25 sparse embedding for genomic intervals, motivated by Qdrants own [bm25 implementation](https://github.com/qdrant/fastembed/blob/main/fastembed/sparse/bm25.py) within fastembed. The implementation is actually BM25-_like_, in that it assumes a constant, prior-known average document length. This enables us to compute the BM25 scores for a query interval without needing to know teh distribution of document lengths in the corpus. Moreover, sparse BM25 embeddings for documents need not be recomputed as the corpus grows.

This method is designed to be used in conjunction with one of our dense embedding models, such as Atacformer or Region2Vec to enable hybrid search. The sparse BM25 embedding can be used perform "key-word" search (look for specific regions), while the dense embedding can be used to perform "semantic search" (look for similar biology). By combining the two with a fusion strategy, we can achieve better recall and precision than either method alone.

## Example usage
Here is an example usage of the BM25 embedding:

```python
from gtars.bm25 import Bm25
from gtars.models import RegionSet

model = Bm25(
    tokenizer="/path/to/vocab.bed",
    k=1.5,
    b=0.75,
    avg_doc_length=1_000
)

query = RegionSet("path/to/query.bed")
embedding = model.embed(query)

print(embedding.indices) # [1, 5, 10]
print(embedding.values) # [0.5, 1.0, 0.75]
```

## Use with Atacformer and Qdrant
BM25 can be used with dense embedding models like Atacformer to perform hybrid search in Qdrant.


First, we need to create a Qdrant collection with both dense and sparse vector configurations:
```python
from geniml.atacformer import AtacformerForCellClustering
from gtars.bm25 import Bm25
from gtars.models import RegionSet
from gtars.tokenizers import Tokenizer

from qdrant_client import models as qdrant_models
from qdrant_client import QdrantClient

# instantiate the qdrant collection
client = QdrantClient("http://localhost:6333")
client.recreate_collection(
    collection_name="bedbase",
    # atacformer embeddings
    vectors_config={
        "dense": qdrant_models.VectorParams(
            size=384,
            distance=qdrant_models.Distance.COSINE
        ),
    },
    # bm25 sparse embeddings
    sparse_vectors_config={
        "sparse": qdrant_models.SparseVectorsConfig(
            modifier=qdrant_models.Modifier.IDF
        )
    }
)
```

Then we can instantiate our Atacformer and BM25 models, and insert some data into the collection:
```python
# instantiate the models
tokenizer = Tokenizer.from_pretrained("databio/atacformer-ctft-hg38")
atacformer = AtacformerForCellClustering.from_pretrained("databio/atacformer-ctft-hg38")
bm25 = Bm25(
    tokenizer=tokenizer,
    k=1.5,
    b=0.75,
    avg_doc_length=1_000 # bed files are usually very large
)

documents = [
    RegionSet("path/to/document1.bed"),
    RegionSet("path/to/document2.bed"),
    RegionSet("path/to/document3.bed"),
    RegionSet("path/to/document4.bed"),
    RegionSet("path/to/document5.bed"),
]

for i, document in enumerate(documents):
    dense_embedding = atacformer.embed(document)
    sparse_embedding = bm25.embed(document)

    client.upsert(
        collection_name="bedbase",
        points=[
            qdrant_models.PointStruct(
                id=i,
                vector=dense_embedding,
                sparse_vector=sparse_embedding
            )
        ]
    )
```

Finally, we can perform a hybrid search using both the dense and sparse embeddings:
```python
query = RegionSet("path/to/query.bed")
dense_query_embedding = atacformer.embed(query)
sparse_query_embedding = bm25.embed(query)

response = client.query_points(
    collection_name="bedbase",
    prefetch=[
        qdrant_models.Prefetch(
            query=sparse_query_embedding,
            using="sparse",
            limit=3,
        ),
        qdrant_models.Prefetch(
            query=dense_query_embedding,
            using="dense",
            limit=3,
        )
    ],
    query=qdrant_models.FusionQuery(fusion=qdrant_models.Fusion.RRF),
    limit=3,
)
```
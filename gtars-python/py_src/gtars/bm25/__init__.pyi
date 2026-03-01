from typing import Any, List, Union

from gtars.tokenizers import Tokenizer

class SparseVector:
    """
    A sparse vector with indices and values, compatible with Qdrant's sparse vector format.
    """

    @property
    def indices(self) -> List[int]:
        """
        The non-zero dimension indices.
        """

    @property
    def values(self) -> List[float]:
        """
        The values at each non-zero dimension.
        """

    def __len__(self) -> int: ...
    def __repr__(self) -> str: ...

class Bm25:
    """
    BM25 sparse embedding model for genomic intervals.

    Computes BM25-like term frequency scores over a vocabulary of genomic regions.
    Designed for hybrid search with dense models like Atacformer via Qdrant.
    """

    def __init__(
        self,
        tokenizer: Union[str, Tokenizer],
        k: float = 1.0,
        b: float = 0.75,
        avg_doc_length: float = 1000.0,
    ) -> None:
        """
        Create a new BM25 model.

        Args:
            tokenizer: A path to a BED/BED.GZ/TOML vocabulary file, or an existing Tokenizer object.
            k: Term frequency saturation parameter. Higher values increase the impact of term frequency.
            b: Document length normalization parameter (0 to 1). Higher values penalize longer documents more.
            avg_doc_length: Assumed average document length (constant across the corpus).
        """

    def embed(self, regions: Any) -> SparseVector:
        """
        Embed a set of regions into a BM25 sparse vector.

        Args:
            regions: A file path to a BED file, or an iterable of objects with chr, start, end attributes.

        Returns:
            SparseVector: Sparse vector with token ID indices and BM25 term-frequency scores as values.
        """

    def tokenize(self, regions: Any) -> List[int]:
        """
        Tokenize regions into token IDs without computing BM25 scores.

        Args:
            regions: A file path to a BED file, or an iterable of objects with chr, start, end attributes.

        Returns:
            List[int]: Token IDs of overlapping vocabulary regions.
        """

    @property
    def vocab_size(self) -> int:
        """
        The number of regions in the vocabulary.
        """

    @property
    def k(self) -> float:
        """
        The term frequency saturation parameter.
        """

    @property
    def b(self) -> float:
        """
        The document length normalization parameter.
        """

    @property
    def avg_doc_length(self) -> float:
        """
        The assumed average document length.
        """

    def __repr__(self) -> str: ...
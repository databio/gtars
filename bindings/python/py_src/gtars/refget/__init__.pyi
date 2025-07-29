from typing import Union, Optional, List, Type
from enum import Enum
from os import PathLike

class AlphabetType(Enum):
    """
    Represents the type of alphabet for a sequence.
    """

    Dna2bit: int
    Dna3bit: int
    DnaIupac: int
    Protein: int
    Ascii: int
    Unknown: int

    def __str__(self) -> str: ...

class SequenceMetadata:
    """
    Metadata for a biological sequence.
    """

    name: str
    length: int
    sha512t24u: str
    md5: str
    alphabet: AlphabetType

    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

class SequenceRecord:
    """
    A record representing a biological sequence, including its metadata and optional data.
    """

    metadata: SequenceMetadata
    data: Optional[bytes]  # Vec<u8> maps to bytes in Python

    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

class SeqColDigestLvl1:
    """
    Level 1 digests for a sequence collection.
    """

    sequences_digest: str
    names_digest: str
    lengths_digest: str

    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

class SequenceCollection:
    """
    A collection of biological sequences.
    """

    sequences: List[SequenceRecord]
    digest: str
    lvl1: SeqColDigestLvl1
    file_path: Optional[str]
    has_data: bool

    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

class RetrievedSequence:
    """
    Represents a retrieved sequence segment with its metadata.
    Exposed from the Rust `PyRetrievedSequence` struct.
    """

    sequence: str
    chrom_name: str
    start: int
    end: int

    def __init__(
        self, sequence: str, chrom_name: str, start: int, end: int
    ) -> None: ...
    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

class StorageMode(Enum):
    """
    Defines how sequence data is stored in the Refget store.
    """

    Raw: int
    Encoded: int

class GlobalRefgetStore:
    """
    A global store for refget sequences, allowing import, retrieval, and storage operations.
    """

    def __init__(self, mode: StorageMode) -> None: ...
    def import_fasta(self, file_path: Union[str, PathLike]) -> None:
        """
        Import a fasta into the GlobalRefgetStore
        """
        ...
    def get_sequence_by_id(self, digest: str) -> Optional[SequenceRecord]:
        """
        Retrieves a sequence record by its SHA512t24u or MD5 digest.
        """
        ...
    def get_sequence_by_collection_and_name(
        self, collection_digest: str, sequence_name: str
    ) -> Optional[SequenceRecord]:
        """
        Retrieve a SequenceRecord from the store by its collection digest and name
        """

        ...
    def get_substring(self, seq_digest: str, start: int, end: int) -> Optional[str]:
        """
        Retrieves a substring from an encoded sequence by its SHA512t24u digest.
        Args:
                seq_digest - str - the path to import from
                start - int - The start index of the substring (inclusive)
                end - int - The end index of the substring (exclusive)
        Returns:
                substring - str - returns substring if found, None if not found

        """
        ...
    def write_store_to_directory(
        self, root_path: Union[str, PathLike], seqdata_path_template: str
    ) -> None:
        """
        Write a GlobalRefgetStore object to a directory
        """
        ...
    @classmethod
    def load_from_directory(
        cls: Type["GlobalRefgetStore"], root_path: Union[str, PathLike]
    ) -> "GlobalRefgetStore":
        """
        Load a GlobalRefgetStore from a directory path
        """

        ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

def sha512t24u_digest(readable: Union[str, bytes]) -> str:
    """
    Computes the GA4GH SHA512t24u digest for a given string or bytes.
    """
    ...

def md5_digest(readable: Union[str, bytes]) -> str:
    """
    Computes the MD5 digest for a given string or bytes.
    """
    ...

def digest_fasta(fasta: Union[str, PathLike]) -> SequenceCollection:
    """
    Digests a FASTA file and returns a SequenceCollection object.
    """
    ...

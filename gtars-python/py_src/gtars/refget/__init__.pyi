"""Type stubs and documentation for the gtars.refget module.

This file serves two purposes:

1. **Type Hints**: Provides type annotations for IDE autocomplete and static
   type checking tools like mypy.

2. **Documentation**: Contains Google-style docstrings that mkdocstrings uses
   to generate the API reference documentation website.

Note: The actual implementation is in Rust (gtars-python/src/refget/mod.rs)
and compiled via PyO3. This stub file provides the Python interface definition
and structured documentation that tools can parse properly.
"""

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
    """A global store for GA4GH refget sequences with lazy-loading support.

    GlobalRefgetStore provides content-addressable storage for reference genome
    sequences following the GA4GH refget specification. Supports both local and
    remote stores with on-demand sequence loading.

    Attributes:
        cache_path: Local directory path where the store is located or cached.
            None for in-memory stores.
        remote_url: Remote URL of the store if loaded remotely, None otherwise.

    Examples:
        Create a new store and import sequences::

            from gtars.refget import GlobalRefgetStore, StorageMode
            store = GlobalRefgetStore(StorageMode.Encoded)
            store.import_fasta("genome.fa")

        Load an existing local store::

            store = GlobalRefgetStore.load_local("/data/hg38")
            seq = store.get_substring("chr1_digest", 0, 1000)

        Load a remote store with caching::

            store = GlobalRefgetStore.load_remote(
                "/local/cache",
                "https://example.com/hg38"
            )
    """

    cache_path: Optional[str]
    remote_url: Optional[str]

    def __init__(self, mode: StorageMode) -> None:
        """Create a new empty GlobalRefgetStore.

        Args:
            mode: Storage mode - StorageMode.Raw (uncompressed) or
                StorageMode.Encoded (bit-packed, space-efficient).

        Example::

            store = GlobalRefgetStore(StorageMode.Encoded)
        """
        ...

    @classmethod
    def load_local(cls, cache_path: Union[str, PathLike]) -> "GlobalRefgetStore":
        """Load a local RefgetStore from a directory.

        Loads metadata from the local store immediately; sequence data is loaded
        on-demand when first accessed. This is efficient for large genomes where
        you may only need specific sequences.

        Args:
            cache_path: Local directory containing the refget store (must have
                index.json and sequences.farg files).

        Returns:
            GlobalRefgetStore with metadata loaded, sequences lazy-loaded.

        Raises:
            IOError: If the store directory or index files cannot be read.

        Example::

            store = GlobalRefgetStore.load_local("/data/hg38_store")
            seq = store.get_substring("chr1_digest", 0, 1000)
        """
        ...

    @classmethod
    def load_remote(
        cls, cache_path: Union[str, PathLike], remote_url: str
    ) -> "GlobalRefgetStore":
        """Load a remote RefgetStore with local caching.

        Fetches metadata (index.json, sequences.farg) from a remote URL immediately.
        Sequence data (.seq files) are downloaded on-demand when first accessed and
        cached locally. This is ideal for working with large remote genomes where
        you only need specific sequences.

        By default, persistence is enabled (sequences are cached to disk).
        Call `disable_persistence()` after loading to keep only in memory.

        Args:
            cache_path: Local directory to cache downloaded metadata and sequences.
                Created if it doesn't exist.
            remote_url: Base URL of the remote refget store (e.g.,
                "https://example.com/hg38" or "s3://bucket/hg38").

        Returns:
            GlobalRefgetStore with metadata loaded, sequences fetched on-demand.

        Raises:
            IOError: If remote metadata cannot be fetched or cache cannot be written.

        Example::

            store = GlobalRefgetStore.load_remote(
                "/data/cache/hg38",
                "https://refget-server.com/hg38"
            )
            # First access fetches from remote and caches
            seq = store.get_substring("chr1_digest", 0, 1000)
            # Second access uses cache
            seq2 = store.get_substring("chr1_digest", 1000, 2000)
        """
        ...

    def set_encoding_mode(self, mode: StorageMode) -> None:
        """Change the storage mode, re-encoding/decoding existing sequences as needed.

        When switching from Raw to Encoded, all Full sequences in memory are
        encoded (2-bit packed). When switching from Encoded to Raw, all Full
        sequences in memory are decoded back to raw bytes.

        Args:
            mode: The storage mode to switch to (StorageMode.Raw or StorageMode.Encoded).

        Example::

            store = GlobalRefgetStore.in_memory()
            store.set_encoding_mode(StorageMode.Raw)
        """
        ...

    def enable_persistence(self, path: Union[str, PathLike]) -> None:
        """Enable disk persistence for this store.

        Sets up the store to write sequences to disk. Any in-memory Full sequences
        are flushed to disk and converted to Stubs.

        Args:
            path: Directory for storing sequences and metadata.

        Raises:
            IOError: If the directory cannot be created or written to.

        Example::

            store = GlobalRefgetStore.in_memory()
            store.add_sequence_collection_from_fasta("genome.fa")
            store.enable_persistence("/data/store")  # Flush to disk
        """
        ...

    def disable_persistence(self) -> None:
        """Disable disk persistence for this store.

        New sequences will be kept in memory only. Existing Stub sequences
        can still be loaded from disk if local_path is set.

        Example::

            store = GlobalRefgetStore.load_remote("/cache", "https://example.com")
            store.disable_persistence()  # Stop caching new sequences
        """
        ...

    def import_fasta(self, file_path: Union[str, PathLike]) -> None:
        """Import sequences from a FASTA file into the store.

        Reads all sequences from a FASTA file and adds them to the store.
        Computes GA4GH digests and creates a sequence collection.

        Args:
            file_path: Path to the FASTA file.

        Raises:
            IOError: If the file cannot be read or parsed.

        Example::

            store = GlobalRefgetStore(StorageMode.Encoded)
            store.import_fasta("genome.fa")
        """
        ...

    def get_sequence_by_id(self, digest: str) -> Optional[SequenceRecord]:
        """Retrieve a sequence record by its digest (SHA-512/24u or MD5).

        Searches for a sequence by its GA4GH SHA-512/24u digest. If not found
        and the input looks like an MD5 digest (32 hex characters), tries MD5 lookup.

        Args:
            digest: Sequence digest (SHA-512/24u base64url or MD5 hex string).

        Returns:
            The sequence record if found, None otherwise.

        Example::

            record = store.get_sequence_by_id("aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2")
            if record:
                print(f"Found: {record.metadata.name}")
        """
        ...

    def get_sequence_by_collection_and_name(
        self, collection_digest: str, sequence_name: str
    ) -> Optional[SequenceRecord]:
        """Retrieve a sequence by collection digest and sequence name.

        Looks up a sequence within a specific collection using its name
        (e.g., "chr1", "chrM"). This is useful when you know the genome assembly
        (collection) and chromosome name.

        Args:
            collection_digest: The collection's SHA-512/24u digest.
            sequence_name: Name of the sequence within that collection.

        Returns:
            The sequence record if found, None otherwise.

        Example::

            record = store.get_sequence_by_collection_and_name(
                "uC_UorBNf3YUu1YIDainBhI94CedlNeH",
                "chr1"
            )
        """
        ...

    def get_substring(self, seq_digest: str, start: int, end: int) -> Optional[str]:
        """Extract a substring from a sequence.

        Retrieves a specific region from a sequence using 0-based, half-open
        coordinates [start, end). Automatically loads sequence data if not
        already cached (for lazy-loaded stores).

        Args:
            seq_digest: Sequence digest (SHA-512/24u).
            start: Start position (0-based, inclusive).
            end: End position (0-based, exclusive).

        Returns:
            The substring sequence if found, None otherwise.

        Example::

            # Get first 1000 bases of chr1
            seq = store.get_substring("chr1_digest", 0, 1000)
            print(f"First 50bp: {seq[:50]}")
        """
        ...

    def list_sequences(self) -> List[SequenceMetadata]:
        """List all sequence metadata in the store.

        Returns:
            List of metadata for all sequences in the store.
        """
        ...

    def list_collections(self) -> List[SequenceCollection]:
        """List all sequence collections in the store.

        Returns:
            List of all sequence collections.
        """
        ...

    def write_store_to_directory(
        self, root_path: Union[str, PathLike], seqdata_path_template: str
    ) -> None:
        """Write the store to a directory on disk.

        Persists the store with all sequences and metadata to disk using the
        RefgetStore directory format.

        Args:
            root_path: Directory path to write the store to.
            seqdata_path_template: Path template for sequence files (e.g.,
                "sequences/%s2/%s.seq" where %s2 = first 2 chars of digest,
                %s = full digest).

        Example::

            store.write_store_to_directory(
                "/data/my_store",
                "sequences/%s2/%s.seq"
            )
        """
        ...

    def get_seqs_bed_file(
        self,
        collection_digest: str,
        bed_file_path: Union[str, PathLike],
        output_fasta_path: Union[str, PathLike],
    ) -> None:
        """Extract sequences for BED regions and write to FASTA.

        Args:
            collection_digest: Collection digest to look up sequence names.
            bed_file_path: Path to BED file with regions.
            output_fasta_path: Path to write output FASTA file.
        """
        ...

    def get_seqs_bed_file_to_vec(
        self, collection_digest: str, bed_file_path: Union[str, PathLike]
    ) -> List[RetrievedSequence]:
        """Extract sequences for BED regions and return as list.

        Args:
            collection_digest: Collection digest to look up sequence names.
            bed_file_path: Path to BED file with regions.

        Returns:
            List of retrieved sequence segments.
        """
        ...

    def export_fasta(
        self,
        collection_digest: str,
        output_path: Union[str, PathLike],
        sequence_names: Optional[List[str]] = None,
        line_width: Optional[int] = None,
    ) -> None:
        """Export sequences from a collection to a FASTA file.

        Args:
            collection_digest: Collection to export from.
            output_path: Path to write FASTA file.
            sequence_names: Optional list of sequence names to export. If None,
                exports all sequences in the collection.
            line_width: Optional line width for wrapping sequences. If None,
                uses default of 80.
        """
        ...

    def export_fasta_by_digests(
        self,
        digests: List[str],
        output_path: Union[str, PathLike],
        line_width: Optional[int] = None,
    ) -> None:
        """Export sequences by their digests to a FASTA file.

        Args:
            digests: List of sequence digests to export.
            output_path: Path to write FASTA file.
            line_width: Optional line width for wrapping sequences. If None,
                uses default of 80.
        """
        ...

    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

def sha512t24u_digest(readable: Union[str, bytes]) -> str:
    """Compute the GA4GH SHA-512/24u digest for a sequence.

    This function computes the GA4GH refget standard digest (truncated SHA-512,
    base64url encoded) for a given sequence string or bytes.

    Args:
        readable: Input sequence as str or bytes.

    Returns:
        The SHA-512/24u digest (32 character base64url string).

    Raises:
        TypeError: If input is not str or bytes.

    Example::
        from gtars.refget import sha512t24u_digest
        digest = sha512t24u_digest("ACGT")
        print(digest)
        # Output: 'aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2'
    """
    ...

def md5_digest(readable: Union[str, bytes]) -> str:
    """Compute the MD5 digest for a sequence.

    This function computes the MD5 hash for a given sequence string or bytes.
    MD5 is supported for backward compatibility with legacy systems.

    Args:
        readable: Input sequence as str or bytes.

    Returns:
        The MD5 digest (32 character hexadecimal string).

    Raises:
        TypeError: If input is not str or bytes.

    Example::
        from gtars.refget import md5_digest
        digest = md5_digest("ACGT")
        print(digest)
        # Output: 'f1f8f4bf413b16ad135722aa4591043e'
    """
    ...

def digest_fasta(fasta: Union[str, PathLike]) -> SequenceCollection:
    """Digest all sequences in a FASTA file and compute collection-level digests.

    This function reads a FASTA file and computes GA4GH-compliant digests for
    each sequence, as well as collection-level digests (Level 1 and Level 2)
    following the GA4GH refget specification.

    Args:
        fasta: Path to FASTA file (str or PathLike).

    Returns:
        Collection containing all sequences with their metadata and computed digests.

    Raises:
        IOError: If the FASTA file cannot be read or parsed.

    Example::
        from gtars.refget import digest_fasta
        collection = digest_fasta("genome.fa")
        print(f"Collection digest: {collection.digest}")
        print(f"Number of sequences: {len(collection)}")
    """
    ...

def compute_fai(fasta: Union[str, PathLike]) -> List["FaiRecord"]:
    """Compute FASTA index (FAI) metadata for all sequences in a FASTA file.

    This function computes the FAI index metadata (offset, line_bases, line_bytes)
    for each sequence in a FASTA file, compatible with samtools faidx format.
    Only works with uncompressed FASTA files.

    Args:
        fasta: Path to FASTA file (str or PathLike). Must be uncompressed.

    Returns:
        List of FAI records, one per sequence, containing name, length,
        and FAI metadata (offset, line_bases, line_bytes).

    Raises:
        IOError: If the FASTA file cannot be read or is compressed.

    Example::
        from gtars.refget import compute_fai
        fai_records = compute_fai("genome.fa")
        for record in fai_records:
        ...     print(f"{record.name}: {record.length} bp")
    """
    ...

def load_fasta(fasta: Union[str, PathLike]) -> SequenceCollection:
    """Load a FASTA file with sequence data into a SequenceCollection.

    This function reads a FASTA file and loads all sequences with their data
    into memory. Unlike digest_fasta(), this includes the actual sequence data,
    not just metadata.

    Args:
        fasta: Path to FASTA file (str or PathLike).

    Returns:
        Collection containing all sequences with their metadata and sequence data loaded.

    Raises:
        IOError: If the FASTA file cannot be read or parsed.

    Example::
        from gtars.refget import load_fasta
        collection = load_fasta("genome.fa")
        first_seq = collection[0]
        print(f"Sequence: {first_seq.data[:50]}...")
    """
    ...

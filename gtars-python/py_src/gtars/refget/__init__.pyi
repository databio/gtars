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
    sequence: Optional[bytes]  # Raw sequence bytes when loaded (Full), None when stub-only

    def decode(self) -> Optional[str]:
        """Decode and return the sequence data as a string.

        For Full records with sequence data, returns the decoded sequence.
        For Stub records without sequence data, returns None.

        Returns:
            Decoded sequence string if data is available, None otherwise.
        """
        ...

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

class SequenceCollectionMetadata:
    """
    Metadata for a sequence collection.

    Contains the collection digest and level 1 digests for names, sequences, and lengths.
    This is a lightweight representation of a collection without the actual sequence list.

    Attributes:
        digest: The collection's SHA-512/24u digest.
        n_sequences: Number of sequences in the collection.
        names_digest: Level 1 digest of the names array.
        sequences_digest: Level 1 digest of the sequences array.
        lengths_digest: Level 1 digest of the lengths array.
    """

    digest: str
    n_sequences: int
    names_digest: str
    sequences_digest: str
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

class SequenceCollectionRecord:
    """
    A record representing a sequence collection, which may be a Stub or Full.

    Stub records contain only metadata (digest, n_sequences, level 1 digests).
    Full records contain metadata plus the list of SequenceRecord objects.
    """

    metadata: SequenceCollectionMetadata

    @property
    def sequences(self) -> Optional[List[SequenceRecord]]:
        """Get the sequences if loaded (Full), None if stub-only."""
        ...

    def has_sequences(self) -> bool:
        """Check if this record has sequences loaded (is Full, not Stub)."""
        ...

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

class RefgetStore:
    """A global store for GA4GH refget sequences with lazy-loading support.

    RefgetStore provides content-addressable storage for reference genome
    sequences following the GA4GH refget specification. Supports both local and
    remote stores with on-demand sequence loading.

    Attributes:
        cache_path: Local directory path where the store is located or cached.
            None for in-memory stores.
        remote_url: Remote URL of the store if loaded remotely, None otherwise.

    Note:
        **Boolean evaluation**: RefgetStore follows Python container semantics,
        meaning ``bool(store)`` is ``False`` for empty stores (like ``list``,
        ``dict``, etc.). To check if a store variable is initialized (not None),
        use ``if store is not None:`` rather than ``if store:``.

        Example::

            store = RefgetStore.in_memory()  # Empty store
            bool(store)  # False (empty container)
            len(store)   # 0

            # Wrong: checks emptiness, not initialization
            if store:
                process(store)

            # Right: checks if variable is set
            if store is not None:
                process(store)

    Examples:
        Create a new store and import sequences::

            from gtars.refget import RefgetStore, StorageMode
            store = RefgetStore(StorageMode.Encoded)
            store.import_fasta("genome.fa")

        Open an existing local store::

            store = RefgetStore.open_local("/data/hg38")
            seq = store.get_substring("chr1_digest", 0, 1000)

        Open a remote store with caching::

            store = RefgetStore.open_remote(
                "/local/cache",
                "https://example.com/hg38"
            )
    """

    cache_path: Optional[str]
    remote_url: Optional[str]

    def __init__(self, mode: StorageMode) -> None:
        """Create a new empty RefgetStore.

        Args:
            mode: Storage mode - StorageMode.Raw (uncompressed) or
                StorageMode.Encoded (bit-packed, space-efficient).

        Example::

            store = RefgetStore(StorageMode.Encoded)
        """
        ...

    @classmethod
    def in_memory(cls) -> "RefgetStore":
        """Create a new in-memory RefgetStore.

        Creates a store that keeps all sequences in memory. Use this for
        temporary processing or when you don't need disk persistence.

        Returns:
            New empty RefgetStore with Encoded storage mode.

        Example::

            store = RefgetStore.in_memory()
            store.import_fasta("genome.fa")
        """
        ...

    @classmethod
    def on_disk(cls, cache_path: Union[str, PathLike]) -> "RefgetStore":
        """Create or load a disk-backed RefgetStore.

        If the directory contains an existing store (rgstore.json),
        loads it. Otherwise creates a new store with Encoded mode.

        Args:
            cache_path: Directory path for the store. Created if it doesn't exist.

        Returns:
            RefgetStore (new or loaded from disk).

        Example::

            store = RefgetStore.on_disk("/data/my_store")
            store.import_fasta("genome.fa")
            # Store is automatically persisted to disk
        """
        ...

    @classmethod
    def open_local(cls, path: Union[str, PathLike]) -> "RefgetStore":
        """Open a local RefgetStore from a directory.

        Loads only lightweight metadata and stubs. Collections and sequences
        remain as stubs until explicitly accessed with get_collection()/get_sequence().

        Expects: rgstore.json, sequences.rgsi, collections.rgci, collections/*.rgsi

        Args:
            path: Local directory containing the refget store.

        Returns:
            RefgetStore with metadata loaded, sequences lazy-loaded.

        Raises:
            IOError: If the store directory or index files cannot be read.

        Example::

            store = RefgetStore.open_local("/data/hg38_store")
            seq = store.get_substring("chr1_digest", 0, 1000)
        """
        ...

    @classmethod
    def open_remote(
        cls, cache_path: Union[str, PathLike], remote_url: str
    ) -> "RefgetStore":
        """Open a remote RefgetStore with local caching.

        Loads only lightweight metadata and stubs from the remote URL.
        Data is fetched on-demand when get_collection()/get_sequence() is called.

        By default, persistence is enabled (sequences are cached to disk).
        Call `disable_persistence()` after loading to keep only in memory.

        Args:
            cache_path: Local directory to cache downloaded metadata and sequences.
                Created if it doesn't exist.
            remote_url: Base URL of the remote refget store (e.g.,
                "https://example.com/hg38" or "s3://bucket/hg38").

        Returns:
            RefgetStore with metadata loaded, sequences fetched on-demand.

        Raises:
            IOError: If remote metadata cannot be fetched or cache cannot be written.

        Example::

            store = RefgetStore.open_remote(
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

            store = RefgetStore.in_memory()
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

            store = RefgetStore.in_memory()
            store.add_sequence_collection_from_fasta("genome.fa")
            store.enable_persistence("/data/store")  # Flush to disk
        """
        ...

    def disable_persistence(self) -> None:
        """Disable disk persistence for this store.

        New sequences will be kept in memory only. Existing Stub sequences
        can still be loaded from disk if local_path is set.

        Example::

            store = RefgetStore.open_remote("/cache", "https://example.com")
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

            store = RefgetStore(StorageMode.Encoded)
            store.import_fasta("genome.fa")
        """
        ...

    # =========================================================================
    # Collection API
    # =========================================================================

    def list_collections(self) -> List[SequenceCollectionMetadata]:
        """List all collection metadata in the store.

        Returns metadata for all collections without loading full collection data.
        Use this for browsing/inventory operations.

        Returns:
            List of metadata for all collections.

        Example::

            for meta in store.list_collections():
                print(f"Collection {meta.digest}: {meta.n_sequences} sequences")
        """
        ...

    def get_collection_metadata(self, collection_digest: str) -> Optional[SequenceCollectionMetadata]:
        """Get metadata for a collection by digest.

        Returns lightweight metadata without loading the full collection.
        Use this for quick lookups of collection information.

        Args:
            collection_digest: The collection's SHA-512/24u digest.

        Returns:
            Collection metadata if found, None otherwise.

        Example::

            meta = store.get_collection_metadata("uC_UorBNf3YUu1YIDainBhI94CedlNeH")
            if meta:
                print(f"Collection has {meta.n_sequences} sequences")
        """
        ...

    def get_collection(self, collection_digest: str) -> SequenceCollection:
        """Get a collection by digest with all sequences loaded.

        Loads the collection and all its sequence data into memory.
        Use this when you need full access to sequence content.

        Args:
            collection_digest: The collection's SHA-512/24u digest.

        Returns:
            The collection with all sequence data loaded.

        Raises:
            IOError: If the collection cannot be loaded.

        Example::

            collection = store.get_collection("uC_UorBNf3YUu1YIDainBhI94CedlNeH")
            for seq in collection.sequences:
                print(f"{seq.metadata.name}: {seq.decode()[:20]}...")
        """
        ...

    def iter_collections(self) -> List[SequenceCollection]:
        """Iterate over all collections with their sequences loaded.

        This loads all collection data upfront and returns a list of
        SequenceCollection objects with full sequence data.

        For browsing without loading data, use list_collections() instead.

        Returns:
            List of all collections with loaded sequences.

        Example::

            for coll in store.iter_collections():
                print(f"{coll.digest}: {len(coll.sequences)} sequences")
        """
        ...

    def is_collection_loaded(self, collection_digest: str) -> bool:
        """Check if a collection is fully loaded.

        Returns True if the collection's sequence list is loaded in memory,
        False if it's only metadata (stub).

        Args:
            collection_digest: The collection's SHA-512/24u digest.

        Returns:
            True if loaded, False otherwise.
        """
        ...

    # =========================================================================
    # Sequence API
    # =========================================================================

    def list_sequences(self) -> List[SequenceMetadata]:
        """List all sequence metadata in the store.

        Returns metadata for all sequences without loading sequence data.
        Use this for browsing/inventory operations.

        Returns:
            List of metadata for all sequences in the store.

        Example::

            for meta in store.list_sequences():
                print(f"{meta.name}: {meta.length} bp")
        """
        ...

    def get_sequence_metadata(self, seq_digest: str) -> Optional[SequenceMetadata]:
        """Get metadata for a sequence by digest (no data loaded).

        Use this for lightweight lookups when you don't need the actual sequence.

        Args:
            seq_digest: The sequence's SHA-512/24u digest.

        Returns:
            Sequence metadata if found, None otherwise.
        """
        ...

    def get_sequence(self, digest: str) -> Optional[SequenceRecord]:
        """Retrieve a sequence record by its digest (SHA-512/24u or MD5).

        Loads the sequence data if not already in memory. Supports lookup
        by either SHA-512/24u (preferred) or MD5 digest.

        Args:
            digest: Sequence digest (SHA-512/24u base64url or MD5 hex string).

        Returns:
            The sequence record with data if found, None otherwise.

        Example::

            record = store.get_sequence("aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2")
            if record:
                print(f"Found: {record.metadata.name}")
                print(f"Sequence: {record.decode()[:50]}...")
        """
        ...

    def get_sequence_by_name(
        self, collection_digest: str, sequence_name: str
    ) -> Optional[SequenceRecord]:
        """Retrieve a sequence by collection digest and sequence name.

        Looks up a sequence within a specific collection using its name
        (e.g., "chr1", "chrM"). Loads the sequence data if needed.

        Args:
            collection_digest: The collection's SHA-512/24u digest.
            sequence_name: Name of the sequence within that collection.

        Returns:
            The sequence record with data if found, None otherwise.

        Example::

            record = store.get_sequence_by_name(
                "uC_UorBNf3YUu1YIDainBhI94CedlNeH",
                "chr1"
            )
            if record:
                print(f"Sequence: {record.decode()[:50]}...")
        """
        ...

    def iter_sequences(self) -> List[SequenceRecord]:
        """Iterate over all sequences with their data loaded.

        This ensures all sequence data is loaded and returns a list of
        SequenceRecord objects with full sequence data.

        For browsing without loading data, use list_sequences() instead.

        Returns:
            List of all sequences with loaded data.

        Example::

            for seq in store.iter_sequences():
                print(f"{seq.metadata.name}: {seq.decode()[:20]}...")
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

    # =========================================================================
    # Store Management
    # =========================================================================

    def stats(self) -> dict:
        """Returns statistics about the store.

        Returns:
            dict with keys:
                - 'n_sequences': Total number of sequences (Stub + Full)
                - 'n_sequences_loaded': Number of sequences with data loaded (Full)
                - 'n_collections': Total number of collections (Stub + Full)
                - 'n_collections_loaded': Number of collections with sequences loaded (Full)
                - 'storage_mode': Storage mode ('Raw' or 'Encoded')
                - 'total_disk_size': Total size of all files on disk in bytes

        Note:
            n_collections_loaded only reflects collections fully loaded in memory.
            For remote stores, collections are loaded on-demand when accessed.

        Example::

            stats = store.stats()
            print(f"Store has {stats['n_sequences']} sequences")
            print(f"Collections: {stats['n_collections']} total, {stats['n_collections_loaded']} loaded")
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

    # =========================================================================
    # BED/FASTA Export
    # =========================================================================

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

def digest_sequence(
    name: str,
    data: bytes,
    description: Optional[str] = None,
) -> SequenceRecord:
    """Create a SequenceRecord from raw data, computing all metadata.

    This is the sequence-level parallel to digest_fasta() for collections.
    It computes the GA4GH sha512t24u digest, MD5 digest, detects the alphabet,
    and returns a SequenceRecord with computed metadata and the original data.

    The input data is automatically uppercased to ensure consistent digest
    computation (matching FASTA processing behavior).

    Args:
        name: The sequence name (e.g., "chr1").
        data: The raw sequence bytes (e.g., b"ACGTACGT").
        description: Optional description text for the sequence.

    Returns:
        A SequenceRecord with computed metadata and the original data (uppercased).

    Example::
        from gtars.refget import digest_sequence
        seq = digest_sequence("chr1", b"ACGTACGT")
        print(seq.metadata.name, seq.metadata.length)
        # Output: chr1 8

        # With description
        seq2 = digest_sequence("chr1", b"ACGT", description="Chromosome 1")
        print(seq2.metadata.description)
        # Output: Chromosome 1
    """
    ...

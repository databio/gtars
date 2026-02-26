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

from typing import Any, Union, Optional, List, Iterator
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

class FaiMetadata:
    """FASTA index (FAI) metadata for a sequence.

    Contains the information needed to quickly seek to a sequence
    in a FASTA file, compatible with samtools faidx format.

    Attributes:
        offset: Byte offset of the first base in the FASTA file.
        line_bases: Number of bases per line.
        line_bytes: Number of bytes per line (including newline).
    """

    offset: int
    line_bases: int
    line_bytes: int

    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

class FaiRecord:
    """A FASTA index record for a single sequence.

    Represents one line of a .fai index file with sequence name,
    length, and FAI metadata for random access.

    Attributes:
        name: Sequence name.
        length: Sequence length in bases.
        fai: FAI metadata (None for gzipped files).
    """

    name: str
    length: int
    fai: Optional[FaiMetadata]

    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

class SequenceMetadata:
    """Metadata for a biological sequence.

    Contains identifying information and computed digests for a sequence,
    without the actual sequence data.

    Attributes:
        name: Sequence name (first word of FASTA header).
        description: Description from FASTA header (text after first whitespace).
        length: Length of the sequence in bases.
        sha512t24u: GA4GH SHA-512/24u digest (32-char base64url).
        md5: MD5 digest (32-char hex string).
        alphabet: Detected alphabet type (DNA, protein, etc.).
        fai: FASTA index metadata if available.
    """

    name: str
    description: Optional[str]
    length: int
    sha512t24u: str
    md5: str
    alphabet: AlphabetType
    fai: Optional[FaiMetadata]

    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

class SequenceRecord:
    """A record representing a biological sequence, including its metadata and optional data.

    SequenceRecord can be either a "stub" (metadata only) or "full" (metadata + data).
    Stubs are used for lazy-loading where sequence data is fetched on demand.

    Attributes:
        metadata: Sequence metadata (name, length, digests).
        sequence: Raw sequence data if loaded, None for stubs.
        is_loaded: Whether sequence data is loaded (True) or just metadata (False).
    """

    metadata: SequenceMetadata
    sequence: Optional[bytes]

    @property
    def is_loaded(self) -> bool:
        """Whether sequence data is loaded (true) or just metadata (false)."""
        ...

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
    """Metadata for a sequence collection.

    Contains the collection digest and level 1 digests for names, sequences, and lengths.
    This is a lightweight representation of a collection without the actual sequence list.

    Attributes:
        digest: The collection's SHA-512/24u digest.
        n_sequences: Number of sequences in the collection.
        names_digest: Level 1 digest of the names array.
        sequences_digest: Level 1 digest of the sequences array.
        lengths_digest: Level 1 digest of the lengths array.
        name_length_pairs_digest: Ancillary digest (if computed).
        sorted_name_length_pairs_digest: Ancillary digest (if computed).
        sorted_sequences_digest: Ancillary digest (if computed).
    """

    digest: str
    n_sequences: int
    names_digest: str
    sequences_digest: str
    lengths_digest: str
    name_length_pairs_digest: Optional[str]
    sorted_name_length_pairs_digest: Optional[str]
    sorted_sequences_digest: Optional[str]

    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

class SequenceCollection:
    """A collection of biological sequences (e.g., a genome assembly).

    SequenceCollection represents a set of sequences with collection-level
    digests following the GA4GH seqcol specification. Supports iteration,
    indexing, and len().

    Attributes:
        sequences: List of sequence records.
        digest: Collection-level SHA-512/24u digest (Level 2).
        lvl1: Level 1 digests for names, lengths, sequences.
        file_path: Source file path if loaded from FASTA.

    Examples:
        Iterate over sequences::

            for seq in collection:
                print(f"{seq.metadata.name}: {seq.metadata.length} bp")

        Access by index::

            first_seq = collection[0]
            last_seq = collection[-1]

        Get length::

            n = len(collection)
    """

    sequences: List[SequenceRecord]
    digest: str
    lvl1: SeqColDigestLvl1
    file_path: Optional[str]

    def write_fasta(self, file_path: str, line_width: Optional[int] = None) -> None:
        """Write the collection to a FASTA file.

        Args:
            file_path: Path to the output FASTA file.
            line_width: Number of bases per line (default: 70).

        Raises:
            IOError: If any sequence doesn't have data loaded.

        Example::

            collection = load_fasta("genome.fa")
            collection.write_fasta("output.fa")
            collection.write_fasta("output.fa", line_width=60)
        """
        ...

    def __len__(self) -> int: ...
    def __getitem__(self, idx: int) -> SequenceRecord: ...
    def __iter__(self) -> Iterator[SequenceRecord]: ...
    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

class RetrievedSequence:
    """Represents a retrieved sequence segment with its metadata.

    Returned by methods that extract subsequences from specific regions,
    such as substrings_from_regions().

    Attributes:
        sequence: The extracted sequence string.
        chrom_name: Chromosome/sequence name (e.g., "chr1").
        start: Start position (0-based, inclusive).
        end: End position (0-based, exclusive).
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
    """Defines how sequence data is stored in the Refget store.

    Variants:
        Raw: Store sequences as raw bytes (1 byte per base).
        Encoded: Store sequences with 2-bit encoding (4 bases per byte).
    """

    Raw: int
    Encoded: int

class FhrMetadata:
    """FAIR Headers Reference genome metadata for a sequence collection.

    Fields match the FHR 1.0 specification. All fields are optional.
    Note: ``schema_version`` is a number (int or float) per spec, passed as
    a Python numeric type and stored as ``serde_json::Number`` internally.
    """

    genome: Optional[str]
    version: Optional[str]
    masking: Optional[str]
    genome_synonym: Optional[list[str]]
    voucher_specimen: Optional[str]
    documentation: Optional[str]
    identifier: Optional[list[str]]
    scholarly_article: Optional[str]
    funding: Optional[str]

    def __init__(self, **kwargs: Any) -> None: ...

    @staticmethod
    def from_json(path: str) -> "FhrMetadata": ...

    def to_dict(self) -> dict[str, Any]: ...
    def to_json(self, path: str) -> None: ...

    def __repr__(self) -> str: ...

class RefgetStore:
    """A global store for GA4GH refget sequences with lazy-loading support.

    RefgetStore provides content-addressable storage for reference genome
    sequences following the GA4GH refget specification. Supports both local and
    remote stores with on-demand sequence loading.

    Attributes:
        cache_path: Local directory path where the store is located or cached.
            None for in-memory stores.
        remote_url: Remote URL of the store if loaded remotely, None otherwise.
        quiet: Whether the store suppresses progress output.
        storage_mode: Current storage mode (Raw or Encoded).

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

            from gtars.refget import RefgetStore
            store = RefgetStore.in_memory()
            store.add_sequence_collection_from_fasta("genome.fa")

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

    @property
    def quiet(self) -> bool:
        """Whether the store is in quiet mode."""
        ...

    @property
    def storage_mode(self) -> StorageMode:
        """Current storage mode (Raw or Encoded)."""
        ...

    @property
    def is_persisting(self) -> bool:
        """Whether the store is currently persisting to disk.

        Example::

            store = RefgetStore.in_memory()
            print(store.is_persisting)  # False
            store.enable_persistence("/data/store")
            print(store.is_persisting)  # True
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
            store.add_sequence_collection_from_fasta("genome.fa")
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
            store.add_sequence_collection_from_fasta("genome.fa")
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

    def enable_encoding(self) -> None:
        """Enable 2-bit encoding for space efficiency.

        Re-encodes any existing Raw sequences in memory.

        Example::

            store = RefgetStore.in_memory()
            store.disable_encoding()  # Switch to Raw
            store.enable_encoding()   # Back to Encoded
        """
        ...

    def disable_encoding(self) -> None:
        """Disable encoding, use raw byte storage.

        Decodes any existing Encoded sequences in memory.

        Example::

            store = RefgetStore.in_memory()
            store.disable_encoding()  # Switch to Raw mode
        """
        ...

    def set_quiet(self, quiet: bool) -> None:
        """Set whether to suppress progress output.

        When quiet is True, operations like add_sequence_collection_from_fasta
        will not print progress messages.

        Args:
            quiet: Whether to suppress progress output.

        Example::

            store = RefgetStore.in_memory()
            store.set_quiet(True)
            store.add_sequence_collection_from_fasta("genome.fa")  # No output
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

    def add_sequence_collection_from_fasta(
        self,
        file_path: Union[str, PathLike],
        force: bool = False,
        namespaces: Optional[List[str]] = None,
    ) -> tuple[SequenceCollectionMetadata, bool]:
        """Add a sequence collection from a FASTA file.

        Reads a FASTA file, digests the sequences, creates a SequenceCollection,
        and adds it to the store along with all its sequences.

        Args:
            file_path: Path to the FASTA file to import.
            force: If True, overwrite existing collections/sequences.
                If False (default), skip duplicates.
            namespaces: Optional list of namespace prefixes to extract aliases from
                FASTA headers. For example, ["ncbi", "refseq"] will scan headers
                for tokens like ``ncbi:NC_000001.11`` and register them as aliases.

        Returns:
            A tuple containing:
                - SequenceCollectionMetadata: Metadata for the collection.
                - bool: True if the collection was newly added, False if it already existed.

        Raises:
            IOError: If the file cannot be read or processed.

        Example::

            store = RefgetStore.in_memory()
            metadata, was_new = store.add_sequence_collection_from_fasta("genome.fa")
            print(f"{'Added' if was_new else 'Skipped'}: {metadata.digest}")

            # Extract aliases from FASTA headers
            metadata, was_new = store.add_sequence_collection_from_fasta(
                "genome.fa", namespaces=["ncbi", "refseq"]
            )
        """
        ...

    def add_sequence_collection(
        self, collection: SequenceCollection, force: bool = False
    ) -> None:
        """Add a pre-built SequenceCollection to the store.

        Adds a SequenceCollection (created via ``digest_fasta()`` or programmatically)
        directly to the store without reading from a FASTA file.

        Args:
            collection: A SequenceCollection to add.
            force: If True, overwrite existing collections/sequences.
                If False (default), skip duplicates.

        Raises:
            IOError: If the collection cannot be stored.

        Example::

            from gtars.refget import RefgetStore, digest_fasta
            store = RefgetStore.in_memory()
            collection = digest_fasta("genome.fa")
            store.add_sequence_collection(collection)
        """
        ...

    def add_sequence(self, sequence: SequenceRecord, force: bool = False) -> None:
        """Add a sequence to the store without collection association.

        The sequence can be created using ``digest_sequence()`` and later
        retrieved by its digest via ``get_sequence()``.

        Args:
            sequence: A SequenceRecord created by ``digest_sequence()``.
            force: If True, overwrite existing. If False (default), skip duplicates.

        Raises:
            IOError: If the sequence cannot be stored.

        Example::

            from gtars.refget import RefgetStore, digest_sequence
            store = RefgetStore.in_memory()
            seq = digest_sequence(b"ACGTACGT")
            store.add_sequence(seq)
            retrieved = store.get_sequence(seq.metadata.sha512t24u)
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
        Automatically strips "SQ." prefix from digest if present.

        Args:
            seq_digest: The sequence's SHA-512/24u digest, optionally with "SQ." prefix.

        Returns:
            Sequence metadata if found, None otherwise.
        """
        ...

    def get_sequence(self, digest: str) -> SequenceRecord:
        """Retrieve a sequence record by its digest (SHA-512/24u or MD5).

        Loads the sequence data if not already in memory. Supports lookup
        by either SHA-512/24u (preferred) or MD5 digest. Automatically
        strips "SQ." prefix if present (case-insensitive).

        Args:
            digest: Sequence digest (SHA-512/24u base64url or MD5 hex string),
                optionally with "SQ." prefix.

        Returns:
            The sequence record with data.

        Raises:
            KeyError: If the sequence is not found.

        Example::

            record = store.get_sequence("aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2")
            print(f"Found: {record.metadata.name}")
            # Also works with SQ. prefix
            record = store.get_sequence("SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2")
        """
        ...

    def get_sequence_by_name(
        self, collection_digest: str, sequence_name: str
    ) -> SequenceRecord:
        """Retrieve a sequence by collection digest and sequence name.

        Looks up a sequence within a specific collection using its name
        (e.g., "chr1", "chrM"). Loads the sequence data if needed.
        Automatically strips "SQ." prefix from collection digest if present.

        Args:
            collection_digest: The collection's SHA-512/24u digest, optionally
                with "SQ." prefix.
            sequence_name: Name of the sequence within that collection.

        Returns:
            The sequence record with data.

        Raises:
            KeyError: If the sequence is not found.

        Example::

            record = store.get_sequence_by_name(
                "uC_UorBNf3YUu1YIDainBhI94CedlNeH",
                "chr1"
            )
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

    def get_substring(self, seq_digest: str, start: int, end: int) -> str:
        """Extract a substring from a sequence.

        Retrieves a specific region from a sequence using 0-based, half-open
        coordinates [start, end). Automatically loads sequence data if not
        already cached (for lazy-loaded stores). Automatically strips "SQ."
        prefix from digest if present.

        Args:
            seq_digest: Sequence digest (SHA-512/24u), optionally with "SQ." prefix.
            start: Start position (0-based, inclusive).
            end: End position (0-based, exclusive).

        Returns:
            The substring sequence.

        Raises:
            KeyError: If the sequence is not found.

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

    def write(self) -> None:
        """Write the store using its configured paths.

        Convenience method for disk-backed stores. Uses the store's own
        local_path and seqdata_path_template.

        Raises:
            IOError: If the store cannot be written.
        """
        ...

    def write_store_to_dir(
        self,
        root_path: Union[str, PathLike],
        seqdata_path_template: Optional[str] = None,
    ) -> None:
        """Write the store to a directory on disk.

        Persists the store with all sequences and metadata to disk using the
        RefgetStore directory format.

        Args:
            root_path: Directory path to write the store to.
            seqdata_path_template: Optional path template for sequence files (e.g.,
                "sequences/%s2/%s.seq" where %s2 = first 2 chars of digest,
                %s = full digest). Uses default if not specified.

        Example::

            store.write_store_to_dir("/data/my_store")
            store.write_store_to_dir("/data/my_store", "sequences/%s2/%s.seq")
        """
        ...

    # =========================================================================
    # Seqcol Features
    # =========================================================================

    def get_collection_level1(self, digest: str) -> dict:
        """Get level 1 representation (attribute digests) for a collection.

        Args:
            digest: Collection digest.

        Returns:
            dict with spec-compliant field names (names, lengths, sequences,
            plus optional name_length_pairs, sorted_name_length_pairs, sorted_sequences).
        """
        ...

    def get_collection_level2(self, digest: str) -> dict:
        """Get level 2 representation (full arrays, spec format) for a collection.

        Args:
            digest: Collection digest.

        Returns:
            dict with names (list[str]), lengths (list[int]), sequences (list[str]).
        """
        ...

    def compare(self, digest_a: str, digest_b: str) -> dict:
        """Compare two collections by digest.

        Args:
            digest_a: First collection digest.
            digest_b: Second collection digest.

        Returns:
            dict with keys: digests, attributes, array_elements.
        """
        ...

    def find_collections_by_attribute(self, attr_name: str, attr_digest: str) -> List[str]:
        """Find collections by attribute digest.

        Args:
            attr_name: Attribute name (names, lengths, sequences,
                name_length_pairs, sorted_name_length_pairs, sorted_sequences).
            attr_digest: The digest to search for.

        Returns:
            List of collection digests that have the matching attribute.
        """
        ...

    def get_attribute(self, attr_name: str, attr_digest: str) -> Optional[list]:
        """Get attribute array by digest.

        Args:
            attr_name: Attribute name (names, lengths, or sequences).
            attr_digest: The digest to search for.

        Returns:
            The attribute array, or None if not found.
        """
        ...

    def enable_ancillary_digests(self) -> None:
        """Enable computation of ancillary digests."""
        ...

    def disable_ancillary_digests(self) -> None:
        """Disable computation of ancillary digests."""
        ...

    def has_ancillary_digests(self) -> bool:
        """Returns whether ancillary digests are enabled."""
        ...

    def has_attribute_index(self) -> bool:
        """Returns whether the on-disk attribute index is enabled."""
        ...

    def enable_attribute_index(self) -> None:
        """Enable indexed attribute lookup (not yet implemented)."""
        ...

    def disable_attribute_index(self) -> None:
        """Disable indexed attribute lookup, using brute-force scan instead."""
        ...

    # =========================================================================
    # BED/FASTA Export
    # =========================================================================

    def export_fasta_from_regions(
        self,
        collection_digest: str,
        bed_file_path: Union[str, PathLike],
        output_file_path: Union[str, PathLike],
    ) -> None:
        """Export sequences from BED file regions to a FASTA file.

        Reads a BED file defining genomic regions and exports the sequences
        for those regions to a FASTA file.

        Args:
            collection_digest: The collection's SHA-512/24u digest.
            bed_file_path: Path to BED file defining regions.
            output_file_path: Path to write the output FASTA file.

        Raises:
            IOError: If files cannot be read/written or sequences not found.

        Example::

            store.export_fasta_from_regions(
                "uC_UorBNf3YUu1YIDainBhI94CedlNeH",
                "regions.bed",
                "output.fa"
            )
        """
        ...

    def substrings_from_regions(
        self, collection_digest: str, bed_file_path: Union[str, PathLike]
    ) -> List[RetrievedSequence]:
        """Get substrings for BED file regions as a list.

        Reads a BED file and returns a list of sequences for each region.

        Args:
            collection_digest: The collection's SHA-512/24u digest.
            bed_file_path: Path to BED file defining regions.

        Returns:
            List of retrieved sequence segments.

        Raises:
            IOError: If files cannot be read or sequences not found.

        Example::

            sequences = store.substrings_from_regions(
                "uC_UorBNf3YUu1YIDainBhI94CedlNeH",
                "regions.bed"
            )
            for seq in sequences:
                print(f"{seq.chrom_name}:{seq.start}-{seq.end}")
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
        seq_digests: List[str],
        output_path: Union[str, PathLike],
        line_width: Optional[int] = None,
    ) -> None:
        """Export sequences by their digests to a FASTA file.

        Args:
            seq_digests: List of sequence digests to export.
            output_path: Path to write FASTA file.
            line_width: Optional line width for wrapping sequences. If None,
                uses default of 80.
        """
        ...

    # =========================================================================
    # Alias API
    # =========================================================================

    # Sequence aliases
    def add_sequence_alias(self, namespace: str, alias: str, digest: str) -> None:
        """Add a sequence alias: namespace/alias maps to sequence digest."""
        ...
    def get_sequence_by_alias(self, namespace: str, alias: str) -> Optional[SequenceRecord]:
        """Resolve a sequence alias to the sequence record."""
        ...
    def get_aliases_for_sequence(self, digest: str) -> list[tuple[str, str]]:
        """Reverse lookup: find all (namespace, alias) pairs pointing to this sequence digest."""
        ...
    def list_sequence_alias_namespaces(self) -> list[str]:
        """List all sequence alias namespaces."""
        ...
    def list_sequence_aliases(self, namespace: str) -> Optional[list[str]]:
        """List all aliases in a sequence alias namespace."""
        ...
    def remove_sequence_alias(self, namespace: str, alias: str) -> bool:
        """Remove a single sequence alias. Returns True if it existed."""
        ...
    def load_sequence_aliases(self, namespace: str, path: str) -> int:
        """Load sequence aliases from a TSV file (alias\\tdigest per line)."""
        ...

    # Collection aliases
    def add_collection_alias(self, namespace: str, alias: str, digest: str) -> None:
        """Add a collection alias: namespace/alias maps to collection digest."""
        ...
    def get_collection_by_alias(self, namespace: str, alias: str) -> Optional[SequenceCollectionMetadata]:
        """Resolve a collection alias to the collection metadata."""
        ...
    def get_aliases_for_collection(self, digest: str) -> list[tuple[str, str]]:
        """Reverse lookup: find all (namespace, alias) pairs pointing to this collection digest."""
        ...
    def list_collection_alias_namespaces(self) -> list[str]:
        """List all collection alias namespaces."""
        ...
    def list_collection_aliases(self, namespace: str) -> Optional[list[str]]:
        """List all aliases in a collection alias namespace."""
        ...
    def remove_collection_alias(self, namespace: str, alias: str) -> bool:
        """Remove a single collection alias. Returns True if it existed."""
        ...
    def load_collection_aliases(self, namespace: str, path: str) -> int:
        """Load collection aliases from a TSV file (alias\\tdigest per line)."""
        ...

    # FHR metadata
    def set_fhr_metadata(self, collection_digest: str, metadata: FhrMetadata) -> None:
        """Set FHR metadata for a collection."""
        ...
    def get_fhr_metadata(self, collection_digest: str) -> Optional[FhrMetadata]:
        """Get FHR metadata for a collection. Returns None if missing."""
        ...
    def remove_fhr_metadata(self, collection_digest: str) -> bool:
        """Remove FHR metadata for a collection."""
        ...
    def list_fhr_metadata(self) -> list[str]:
        """List all collection digests that have FHR metadata."""
        ...
    def load_fhr_metadata(self, collection_digest: str, path: str) -> None:
        """Load FHR metadata from a JSON file and attach it to a collection."""
        ...

    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[SequenceMetadata]: ...
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

def compute_fai(fasta: Union[str, PathLike]) -> List[FaiRecord]:
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
            print(f"{record.name}: {record.length} bp")
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
        print(f"Sequence: {first_seq.decode()[:50]}...")
    """
    ...

def digest_sequence(
    data: bytes,
    name: Optional[str] = None,
    description: Optional[str] = None,
) -> SequenceRecord:
    """Create a SequenceRecord from raw data, computing all metadata.

    This is the sequence-level parallel to digest_fasta() for collections.
    It computes the GA4GH sha512t24u digest, MD5 digest, detects the alphabet,
    and returns a SequenceRecord with computed metadata and the original data.

    The input data is automatically uppercased to ensure consistent digest
    computation (matching FASTA processing behavior).

    Args:
        data: The raw sequence bytes (e.g., b"ACGTACGT").
        name: Optional sequence name (e.g., "chr1"). Defaults to "" if not provided.
        description: Optional description text for the sequence.

    Returns:
        A SequenceRecord with computed metadata and the original data (uppercased).

    Example::
        from gtars.refget import digest_sequence
        seq = digest_sequence(b"ACGTACGT")
        print(seq.metadata.length)
        # Output: 8

        seq = digest_sequence(b"ACGT", name="chr1")
        print(seq.metadata.name, seq.metadata.length)
        # Output: chr1 4

        # With description
        seq2 = digest_sequence(b"ACGT", name="chr1", description="Chromosome 1")
        print(seq2.metadata.description)
        # Output: Chromosome 1
    """
    ...

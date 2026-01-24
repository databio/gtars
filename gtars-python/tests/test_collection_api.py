"""
Tests for RefgetStore Collection/Sequence Retrieval API

Tests the 9-method API:
- list_collections(), list_sequences()
- get_collection_metadata(), get_collection()
- get_sequence_metadata(), get_sequence(), get_sequence_by_name()
- iter_collections(), iter_sequences()
"""

import pytest
from gtars.refget import RefgetStore, digest_fasta


class TestCollectionAPI:
    """Test suite for the RefgetStore collection/sequence retrieval API."""

    @pytest.fixture
    def store_with_data(self):
        """Create a store with test data loaded."""
        store = RefgetStore.in_memory()
        store.add_sequence_collection_from_fasta("../tests/data/fasta/base.fa")
        return store

    @pytest.fixture
    def expected_data(self):
        """Get expected data from digesting the FASTA directly."""
        return digest_fasta("../tests/data/fasta/base.fa")

    # =========================================================================
    # list_collections() tests
    # =========================================================================

    def test_list_collections_returns_metadata_list(self, store_with_data):
        """list_collections() returns a list of SequenceCollectionMetadata."""
        collections = store_with_data.list_collections()

        assert isinstance(collections, list)
        assert len(collections) == 1  # One FASTA = one collection

        meta = collections[0]
        # Check it's metadata, not full collection
        assert hasattr(meta, 'digest')
        assert hasattr(meta, 'n_sequences')
        assert hasattr(meta, 'names_digest')
        assert hasattr(meta, 'sequences_digest')
        assert hasattr(meta, 'lengths_digest')
        # Should NOT have sequences attribute (it's metadata only)
        assert not hasattr(meta, 'sequences')

    def test_list_collections_n_sequences_correct(self, store_with_data):
        """list_collections() reports correct n_sequences."""
        collections = store_with_data.list_collections()
        meta = collections[0]

        assert meta.n_sequences == 3  # base.fa has 3 sequences

    # =========================================================================
    # list_sequences() tests
    # =========================================================================

    def test_list_sequences_returns_metadata_list(self, store_with_data):
        """list_sequences() returns a list of SequenceMetadata."""
        sequences = store_with_data.list_sequences()

        assert isinstance(sequences, list)
        assert len(sequences) == 3  # base.fa has 3 sequences

        for meta in sequences:
            assert hasattr(meta, 'name')
            assert hasattr(meta, 'length')
            assert hasattr(meta, 'sha512t24u')
            assert hasattr(meta, 'md5')
            assert hasattr(meta, 'alphabet')
            # Should NOT have sequence data
            assert not hasattr(meta, 'sequence')

    # =========================================================================
    # get_collection_metadata() tests
    # =========================================================================

    def test_get_collection_metadata_returns_metadata(self, store_with_data, expected_data):
        """get_collection_metadata() returns SequenceCollectionMetadata."""
        digest = expected_data.digest
        meta = store_with_data.get_collection_metadata(digest)

        assert meta is not None
        assert meta.digest == digest
        assert meta.n_sequences == 3
        assert hasattr(meta, 'names_digest')
        assert hasattr(meta, 'sequences_digest')
        assert hasattr(meta, 'lengths_digest')
        # Should NOT have sequences
        assert not hasattr(meta, 'sequences')

    def test_get_collection_metadata_not_found(self, store_with_data):
        """get_collection_metadata() returns None for non-existent digest."""
        result = store_with_data.get_collection_metadata("nonexistent_digest")
        assert result is None

    # =========================================================================
    # get_collection() tests
    # =========================================================================

    def test_get_collection_returns_full_collection(self, store_with_data, expected_data):
        """get_collection() returns SequenceCollection with loaded sequences."""
        digest = expected_data.digest
        collection = store_with_data.get_collection(digest)

        assert collection is not None
        assert collection.digest == digest
        assert hasattr(collection, 'sequences')
        assert hasattr(collection, 'lvl1')

    def test_get_collection_has_correct_sequence_count(self, store_with_data, expected_data):
        """get_collection() returns collection with correct number of sequences."""
        digest = expected_data.digest
        collection = store_with_data.get_collection(digest)

        # The key test: sequences list should match n_sequences
        assert len(collection.sequences) == 3

    def test_get_collection_sequences_have_data(self, store_with_data, expected_data):
        """get_collection() returns sequences that can be decoded."""
        digest = expected_data.digest
        collection = store_with_data.get_collection(digest)

        for seq in collection.sequences:
            decoded = seq.decode()
            assert decoded is not None, f"Sequence {seq.metadata.name} has no data"
            assert len(decoded) > 0, f"Sequence {seq.metadata.name} is empty"
            assert len(decoded) == seq.metadata.length

    # =========================================================================
    # get_sequence_metadata() tests
    # =========================================================================

    def test_get_sequence_metadata_returns_metadata(self, store_with_data, expected_data):
        """get_sequence_metadata() returns SequenceMetadata."""
        seq_digest = expected_data.sequences[0].metadata.sha512t24u
        meta = store_with_data.get_sequence_metadata(seq_digest)

        assert meta is not None
        assert meta.sha512t24u == seq_digest
        assert hasattr(meta, 'name')
        assert hasattr(meta, 'length')
        assert hasattr(meta, 'md5')

    def test_get_sequence_metadata_not_found(self, store_with_data):
        """get_sequence_metadata() returns None for non-existent digest."""
        result = store_with_data.get_sequence_metadata("nonexistent_digest")
        assert result is None

    # =========================================================================
    # get_sequence() tests
    # =========================================================================

    def test_get_sequence_returns_record_with_data(self, store_with_data, expected_data):
        """get_sequence() returns SequenceRecord with loadable data."""
        seq_digest = expected_data.sequences[0].metadata.sha512t24u
        record = store_with_data.get_sequence(seq_digest)

        assert record is not None
        assert record.metadata.sha512t24u == seq_digest

        decoded = record.decode()
        assert decoded is not None
        assert len(decoded) == record.metadata.length

    def test_get_sequence_by_md5(self, store_with_data, expected_data):
        """get_sequence() works with MD5 digest."""
        md5_digest = expected_data.sequences[0].metadata.md5
        record = store_with_data.get_sequence(md5_digest)

        assert record is not None
        assert record.metadata.md5 == md5_digest
        assert record.decode() is not None

    def test_get_sequence_not_found(self, store_with_data):
        """get_sequence() returns None for non-existent digest."""
        result = store_with_data.get_sequence("nonexistent_digest")
        assert result is None

    # =========================================================================
    # get_sequence_by_name() tests
    # =========================================================================

    def test_get_sequence_by_name_returns_record_with_data(self, store_with_data, expected_data):
        """get_sequence_by_name() returns SequenceRecord with loadable data."""
        collection_digest = expected_data.digest
        seq_name = expected_data.sequences[0].metadata.name

        record = store_with_data.get_sequence_by_name(collection_digest, seq_name)

        assert record is not None
        assert record.metadata.name == seq_name

        decoded = record.decode()
        assert decoded is not None
        assert len(decoded) == record.metadata.length

    def test_get_sequence_by_name_not_found(self, store_with_data, expected_data):
        """get_sequence_by_name() returns None for non-existent name."""
        collection_digest = expected_data.digest
        result = store_with_data.get_sequence_by_name(collection_digest, "nonexistent_seq")
        assert result is None

    # =========================================================================
    # iter_collections() tests
    # =========================================================================

    def test_iter_collections_returns_full_collections(self, store_with_data):
        """iter_collections() returns list of SequenceCollection with data."""
        collections = store_with_data.iter_collections()

        assert isinstance(collections, list)
        assert len(collections) == 1

        coll = collections[0]
        assert hasattr(coll, 'digest')
        assert hasattr(coll, 'sequences')
        assert len(coll.sequences) == 3

    def test_iter_collections_sequences_decodable(self, store_with_data):
        """iter_collections() returns collections with decodable sequences."""
        collections = store_with_data.iter_collections()

        for coll in collections:
            for seq in coll.sequences:
                decoded = seq.decode()
                assert decoded is not None, f"Sequence {seq.metadata.name} not decodable"

    # =========================================================================
    # iter_sequences() tests
    # =========================================================================

    def test_iter_sequences_returns_records_with_data(self, store_with_data):
        """iter_sequences() returns list of SequenceRecord with data."""
        sequences = store_with_data.iter_sequences()

        assert isinstance(sequences, list)
        assert len(sequences) == 3

    def test_iter_sequences_all_decodable(self, store_with_data):
        """iter_sequences() returns sequences that can all be decoded."""
        sequences = store_with_data.iter_sequences()

        for seq in sequences:
            decoded = seq.decode()
            assert decoded is not None, f"Sequence {seq.metadata.name} not decodable"
            assert len(decoded) == seq.metadata.length

    # =========================================================================
    # Cross-method consistency tests
    # =========================================================================

    def test_list_vs_get_collection_consistency(self, store_with_data):
        """list_collections() and get_collection() report same n_sequences."""
        meta_list = store_with_data.list_collections()
        meta = meta_list[0]

        full_collection = store_with_data.get_collection(meta.digest)

        # n_sequences from metadata should equal actual sequence count
        assert meta.n_sequences == len(full_collection.sequences)

    def test_list_vs_iter_sequences_consistency(self, store_with_data):
        """list_sequences() and iter_sequences() return same count."""
        metadata_list = store_with_data.list_sequences()
        full_list = store_with_data.iter_sequences()

        assert len(metadata_list) == len(full_list)

    def test_metadata_vs_full_sequence_consistency(self, store_with_data, expected_data):
        """Metadata from get_sequence_metadata matches get_sequence."""
        seq_digest = expected_data.sequences[0].metadata.sha512t24u

        meta = store_with_data.get_sequence_metadata(seq_digest)
        full = store_with_data.get_sequence(seq_digest)

        assert meta.name == full.metadata.name
        assert meta.length == full.metadata.length
        assert meta.sha512t24u == full.metadata.sha512t24u
        assert meta.md5 == full.metadata.md5

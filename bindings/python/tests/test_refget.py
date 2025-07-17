import pytest
from gtars.refget import (
    GlobalRefgetStore,
    StorageMode,
    digest_fasta,
    sha512t24u_digest,
    md5_digest
)

import os
import tempfile

class TestRefget:
    
    def test_digest_fasta_basic(self):
        """Test basic functionality of digest_fasta"""
        fasta_path = "../../gtars/tests/data/fasta/base.fa"
        
        # Test that digest_fasta returns a PySequenceCollection
        result = digest_fasta(fasta_path)
        
        # Test that it has the expected structure
        assert hasattr(result, 'sequences')
        assert hasattr(result, 'digest')
        assert hasattr(result, 'lvl1')
        assert hasattr(result, 'file_path')
        assert hasattr(result, 'has_data')
        
        # Test that sequences is a list
        assert isinstance(result.sequences, list)
        
        # Test that we have the expected number of sequences (from the Rust test)
        assert len(result.sequences) == 3
        
        # Test the first sequence
        seq0 = result.sequences[0]
        assert hasattr(seq0, 'metadata')
        assert hasattr(seq0, 'data')
        
        # Test metadata
        metadata = seq0.metadata
        assert hasattr(metadata, 'name')
        assert hasattr(metadata, 'length')
        assert hasattr(metadata, 'sha512t24u')
        assert hasattr(metadata, 'md5')
        assert hasattr(metadata, 'alphabet')
        
        # Test specific values from the Rust test
        assert metadata.length == 8
        assert metadata.sha512t24u == "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        assert metadata.md5 == "5f63cfaa3ef61f88c9635fb9d18ec945"
        
        # Test lvl1 digests
        lvl1 = result.lvl1
        assert hasattr(lvl1, 'sequences_digest')
        assert hasattr(lvl1, 'names_digest')
        assert hasattr(lvl1, 'lengths_digest')
        
        # Test specific digest values from the Rust test
        assert lvl1.sequences_digest == "0uDQVLuHaOZi1u76LjV__yrVUIz9Bwhr"
        assert lvl1.names_digest == "Fw1r9eRxfOZD98KKrhlYQNEdSRHoVxAG"
        assert lvl1.lengths_digest == "cGRMZIb3AVgkcAfNv39RN7hnT5Chk7RX"
        
        # Test string representations
        assert str(result).startswith("SequenceCollection with 3 sequences")
        assert str(metadata).startswith("SequenceMetadata for sequence")
        assert str(lvl1).startswith("SeqColDigestLvl1:")
        
    def test_digest_fasta_nonexistent_file(self):
        """Test that digest_fasta raises an error for non-existent files"""
        with pytest.raises(Exception):
            digest_fasta("nonexistent.fa")

    def test_global_refget_store_basic(self):
            """Test basic creation and operations of GlobalRefgetStore"""
            # Test store creation with both modes
            store_raw = GlobalRefgetStore(StorageMode.Raw)
            store_encoded = GlobalRefgetStore(StorageMode.Encoded)

            # Test string representations
            #print(store_raw.__repr__)
            assert repr(store_raw).startswith("<GlobalRefgetStore")
            assert repr(store_encoded).startswith("<GlobalRefgetStore")

    def test_store_import_and_retrieve(self):
        """Test importing FASTA and retrieving sequences"""
        store = GlobalRefgetStore(StorageMode.Raw)
        fasta_path = "../../gtars/tests/data/fasta/base.fa"
        # fasta_path = os.path.abspath("../../gtars/tests/data/fasta/base.fa")
        # Import FASTA
        store.import_fasta(fasta_path)

        # Get known sequence by SHA512t24u digest
        sha512 = "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        seq = store.get_sequence_by_id(sha512)
        assert seq is not None
        assert seq.metadata.length == 8

        # TODO This doesnt work
        # Get known sequence by MD5
        # md5 = "5f63cfaa3ef61f88c9635fb9d18ec945"
        # seq = store.get_sequence_by_id(md5)
        # assert seq is not None
        # assert seq.metadata.length == 8

    def test_store_substring(self):
        """Test substring retrieval"""
        store = GlobalRefgetStore(StorageMode.Raw)
        fasta_path = "../../gtars/tests/data/fasta/base.fa"
        store.import_fasta(fasta_path)

        # Get a substring from a known sequence
        sha512 = "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        substr = store.get_substring(sha512, 0, 4)
        assert substr is not None
        assert len(substr) == 4
#
    def test_store_persistence(self):
        """Test saving and loading store"""
        store = GlobalRefgetStore(StorageMode.Raw)
        fasta_path = "../../gtars/tests/data/fasta/base.fa"
        store.import_fasta(fasta_path)

        # Save store to temporary directory
        with tempfile.TemporaryDirectory() as tmpdir:
            store.write_store_to_directory(tmpdir, "sequences/%s2/%s.seq")

            # Load store back
            loaded_store = GlobalRefgetStore.load_from_directory(tmpdir)

            # Verify same sequences exist
            sha512 = "iYtREV555dUFKg2_agSJW6suquUyPpMw"
            seq1 = store.get_sequence_by_id(sha512)
            seq2 = loaded_store.get_sequence_by_id(sha512)

            assert seq1 is not None
            assert seq2 is not None
            assert seq1.metadata.length == seq2.metadata.length
            assert seq1.metadata.sha512t24u == seq2.metadata.sha512t24u
#
    def test_store_errors(self):
        """Test error conditions"""
        store = GlobalRefgetStore(StorageMode.Raw)

        # Test importing non-existent file
        with pytest.raises(Exception):
            store.import_fasta("nonexistent.fa")

        # Test getting non-existent sequence
        bogus_digest = "not_a_sequence"
        assert store.get_sequence_by_id(bogus_digest) is None

        # TODO fix md5
        #assert store.get_sequence_by_md5(bogus_digest) is None

        # Test invalid substring parameters
        sha512 = "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        store.import_fasta("../../gtars/tests/data/fasta/base.fa")
        assert store.get_substring(sha512, 10, 5) is None  # end < start
        assert store.get_substring(sha512, 0, 100) is None  # end > length
#
    def test_store_collection_operations(self):
        """Test collection-related operations"""
        store = GlobalRefgetStore(StorageMode.Raw)
        fasta_path = "../../gtars/tests/data/fasta/base.fa"

        # Import sequences and get sequence by collection and name
        store.import_fasta(fasta_path)

        # Get sequence from default collection
        result = digest_fasta(fasta_path)
        seq = store.get_sequence_by_collection_and_name(
            result.digest,
            result.sequences[0].metadata.name
        )

        assert seq is not None
        assert seq.metadata.length == 8
        assert seq.metadata.sha512t24u == "iYtREV555dUFKg2_agSJW6suquUyPpMw"
import pytest
from gtars.refget import digest_fasta

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

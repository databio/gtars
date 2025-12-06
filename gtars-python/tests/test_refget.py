import pytest
from gtars.refget import (
    GlobalRefgetStore,
    StorageMode,
    digest_fasta,
    load_fasta,
    sha512t24u_digest,
    md5_digest,
    RetrievedSequence,
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
        assert hasattr(result, "sequences")
        assert hasattr(result, "digest")
        assert hasattr(result, "lvl1")
        assert hasattr(result, "file_path")

        # Test that sequences is a list
        assert isinstance(result.sequences, list)

        # Test that we have the expected number of sequences (from the Rust test)
        assert len(result.sequences) == 3

        # Test the first sequence
        seq0 = result.sequences[0]
        assert hasattr(seq0, "metadata")
        assert hasattr(seq0, "data")

        # Test metadata
        metadata = seq0.metadata
        assert hasattr(metadata, "name")
        assert hasattr(metadata, "length")
        assert hasattr(metadata, "sha512t24u")
        assert hasattr(metadata, "md5")
        assert hasattr(metadata, "alphabet")

        # Test specific values from the Rust test
        assert metadata.length == 8
        assert metadata.sha512t24u == "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        assert metadata.md5 == "5f63cfaa3ef61f88c9635fb9d18ec945"

        # Test lvl1 digests
        lvl1 = result.lvl1
        assert hasattr(lvl1, "sequences_digest")
        assert hasattr(lvl1, "names_digest")
        assert hasattr(lvl1, "lengths_digest")

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
        # New repr format shows n_sequences and location
        assert repr(store_raw).startswith("GlobalRefgetStore")
        assert repr(store_encoded).startswith("GlobalRefgetStore")
        assert "n_sequences=0" in repr(store_raw)  # Empty store
        assert "memory-only" in repr(store_raw)  # No local path

    def test_store_import_and_retrieve(self):
        """Test importing FASTA and retrieving sequences"""
        store = GlobalRefgetStore(StorageMode.Raw)
        fasta_path = "../../gtars/tests/data/fasta/base.fa"
        # fasta_path = os.path.abspath("../../gtars/tests/data/fasta/base.fa")
        # Import FASTA
        store.add_sequence_collection_from_fasta(fasta_path)

        # Get known sequence by SHA512t24u digest
        sha512 = "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        seq = store.get_sequence_by_id(sha512)
        assert seq is not None
        assert seq.metadata.length == 8

        md5 = "5f63cfaa3ef61f88c9635fb9d18ec945"
        seq = store.get_sequence_by_id(md5)
        assert seq is not None
        assert seq.metadata.length == 8

    def test_store_substring(self):
        """Test substring retrieval"""
        store = GlobalRefgetStore(StorageMode.Raw)
        fasta_path = "../../gtars/tests/data/fasta/base.fa"
        store.add_sequence_collection_from_fasta(fasta_path)

        # Get a substring from a known sequence
        sha512 = "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        substr = store.get_substring(sha512, 0, 4)
        assert substr is not None
        assert len(substr) == 4

    def test_store_persistence(self):
        """Test saving and loading store"""
        store = GlobalRefgetStore(StorageMode.Raw)
        fasta_path = "../../gtars/tests/data/fasta/base.fa"
        store.add_sequence_collection_from_fasta(fasta_path)

        # Save store to temporary directory
        with tempfile.TemporaryDirectory() as tmpdir:
            store.write_store_to_dir(tmpdir, "sequences/%s2/%s.seq")

            # Load store back
            loaded_store = GlobalRefgetStore.load_local(tmpdir)

            # Verify same sequences exist
            sha512 = "iYtREV555dUFKg2_agSJW6suquUyPpMw"
            seq1 = store.get_sequence_by_id(sha512)
            seq2 = loaded_store.get_sequence_by_id(sha512)

            assert seq1 is not None
            assert seq2 is not None
            assert seq1.metadata.length == seq2.metadata.length
            assert seq1.metadata.sha512t24u == seq2.metadata.sha512t24u

    def test_store_errors(self):
        """Test error conditions"""
        store = GlobalRefgetStore(StorageMode.Raw)

        # Test importing non-existent file
        with pytest.raises(Exception):
            store.add_sequence_collection_from_fasta("nonexistent.fa")

        # Test getting non-existent sequence
        bogus_digest = "not_a_sequence"
        assert store.get_sequence_by_id(bogus_digest) is None

        # Test invalid substring parameters
        sha512 = "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        store.add_sequence_collection_from_fasta("../../gtars/tests/data/fasta/base.fa")
        assert store.get_substring(sha512, 10, 5) is None  # end < start
        assert store.get_substring(sha512, 0, 100) is None  # end > length

    def test_store_collection_operations(self):
        """Test collection-related operations"""
        store = GlobalRefgetStore(StorageMode.Raw)
        fasta_path = "../../gtars/tests/data/fasta/base.fa"

        # Import sequences and get sequence by collection and name
        store.add_sequence_collection_from_fasta(fasta_path)

        # Get sequence from default collection
        result = digest_fasta(fasta_path)
        seq = store.get_sequence_by_collection_and_name(
            result.digest, result.sequences[0].metadata.name
        )

        assert seq is not None
        assert seq.metadata.length == 8
        assert seq.metadata.sha512t24u == "iYtREV555dUFKg2_agSJW6suquUyPpMw"

    def test_get_seqs_from_bed_file_bindings(self):
        """
        Test both export_fasta_from_regions (writes to file)
        and substrings_from_regions (returns list) bindings.
        """

        with tempfile.TemporaryDirectory() as temp_dir_str:
            temp_dir = temp_dir_str

            fasta_content = ">chr1\n" "ATGCATGCATGC\n" ">chr2\n" "GGGGAAAA\n"
            temp_fasta_path = os.path.join(temp_dir, "test.fa")
            with open(temp_fasta_path, "w") as f:
                f.write(fasta_content)


            store = GlobalRefgetStore(StorageMode.Encoded)
            imported_collection = store.add_sequence_collection_from_fasta(temp_fasta_path)
            result = digest_fasta(temp_fasta_path)

            collection_digest = result.digest


            chr1_sha = sha512t24u_digest(b"ATGCATGCATGC")
            chr1_md5 = md5_digest(b"ATGCATGCATGC")
            chr2_sha = sha512t24u_digest(b"GGGGAAAA")
            chr2_md5 = md5_digest(b"GGGGAAAA")

            bed_content = (
                "chr1\t0\t5\n"
                "chr1\t8\t12\n"
                "chr2\t0\t4\n"
                "chr_nonexistent\t10\t20\n"
                "chr1\t-5\t100"
            )
            temp_bed_path = os.path.join(temp_dir, "test.bed")
            with open(temp_bed_path, "w") as f:
                f.write(bed_content)

            temp_output_fa_path = os.path.join(temp_dir, "output.fa")
            store.export_fasta_from_regions(
                collection_digest, temp_bed_path, temp_output_fa_path
            )

            with open(temp_output_fa_path, "r") as f:
                output_fa_content = f.read()

            expected_fa_content = f""">chr1 12 dna2bit {chr1_sha} {chr1_md5}
ATGCAATGC
>chr2 8 dna2bit {chr2_sha} {chr2_md5}
GGGG
"""

            assert (
                output_fa_content.strip() == expected_fa_content.strip()
            ), "Output FASTA file content mismatch"
            print("✓ export_fasta_from_regions binding test passed.")

            vec_result = store.substrings_from_regions(
                collection_digest, temp_bed_path
            )

            expected_vec = [
                RetrievedSequence(sequence="ATGCA", chrom_name="chr1", start=0, end=5),
                RetrievedSequence(sequence="ATGC", chrom_name="chr1", start=8, end=12),
                RetrievedSequence(sequence="GGGG", chrom_name="chr2", start=0, end=4),
            ]

            assert len(vec_result) == len(
                expected_vec
            ), "Length of retrieved sequence list mismatch"
            for i in range(len(vec_result)):
                assert vec_result[i].sequence == expected_vec[i].sequence
                assert vec_result[i].chrom_name == expected_vec[i].chrom_name
                assert vec_result[i].start == expected_vec[i].start
                assert vec_result[i].end == expected_vec[i].end

            print("✓ substrings_from_regions binding test passed.")

    def test_decode_with_no_data(self):
        """Test that decode() returns None when sequence data is not loaded"""
        fasta_path = "../../gtars/tests/data/fasta/base.fa"

        # digest_fasta should not load sequence data
        result = digest_fasta(fasta_path)

        for seq_record in result.sequences:
            assert seq_record.data is None, "digest_fasta should not load sequence data"
            assert seq_record.decode() is None, "decode() should return None when data is None"

    def test_decode_with_loaded_data(self):
        """Test that decode() returns correct sequences when data is loaded"""
        fasta_path = "../../gtars/tests/data/fasta/base.fa"

        # load_fasta should load sequence data
        result = load_fasta(fasta_path)

        # Expected sequences from base.fa
        expected_sequences = [
            ("chrX", "TTGGGGAA"),
            ("chr1", "GGAA"),
            ("chr2", "GCGC"),
        ]

        assert len(result.sequences) == len(expected_sequences)

        for seq_record, (expected_name, expected_seq) in zip(result.sequences, expected_sequences):
            assert seq_record.metadata.name == expected_name
            assert seq_record.data is not None, "load_fasta should load sequence data"

            decoded = seq_record.decode()
            assert decoded is not None, "decode() should return Some when data is present"
            assert decoded == expected_seq, f"Decoded sequence for {expected_name} should match expected"

    def test_decode_with_store_sequences(self):
        """Test decode() with sequences retrieved from a store"""
        store = GlobalRefgetStore(StorageMode.Encoded)
        fasta_path = "../../gtars/tests/data/fasta/base.fa"
        store.add_sequence_collection_from_fasta(fasta_path)

        # Get sequence by ID
        sha512 = "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        seq = store.get_sequence_by_id(sha512)

        assert seq is not None
        assert seq.data is not None, "Store should have loaded sequence data"

        decoded = seq.decode()
        assert decoded is not None
        assert decoded == "TTGGGGAA", "Should correctly decode sequence from encoded store"

    def test_decode_raw_vs_encoded_storage(self):
        """Test that decode() works with both Raw and Encoded storage modes"""
        fasta_path = "../../gtars/tests/data/fasta/base.fa"

        # Test with Raw storage mode
        store_raw = GlobalRefgetStore(StorageMode.Raw)
        store_raw.add_sequence_collection_from_fasta(fasta_path)
        sha512 = "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        seq_raw = store_raw.get_sequence_by_id(sha512)
        decoded_raw = seq_raw.decode()

        # Test with Encoded storage mode
        store_encoded = GlobalRefgetStore(StorageMode.Encoded)
        store_encoded.add_sequence_collection_from_fasta(fasta_path)
        seq_encoded = store_encoded.get_sequence_by_id(sha512)
        decoded_encoded = seq_encoded.decode()

        # Both should produce the same decoded sequence
        assert decoded_raw == decoded_encoded
        assert decoded_raw == "TTGGGGAA"

    def test_load_fasta_function(self):
        """Test the new load_fasta() function"""
        fasta_path = "../../gtars/tests/data/fasta/base.fa"

        # Test that load_fasta returns a SequenceCollection with data
        result = load_fasta(fasta_path)

        # Should have same structure as digest_fasta
        assert hasattr(result, "sequences")
        assert hasattr(result, "digest")
        assert hasattr(result, "lvl1")
        assert hasattr(result, "file_path")

        # But sequences should have data loaded
        assert len(result.sequences) == 3
        for seq_record in result.sequences:
            assert seq_record.data is not None, "load_fasta should load sequence data"
            assert seq_record.decode() is not None, "Should be able to decode loaded data"

        # Verify digests match digest_fasta
        digest_result = digest_fasta(fasta_path)
        assert result.digest == digest_result.digest
        assert result.lvl1.sequences_digest == digest_result.lvl1.sequences_digest
        assert result.lvl1.names_digest == digest_result.lvl1.names_digest
        assert result.lvl1.lengths_digest == digest_result.lvl1.lengths_digest

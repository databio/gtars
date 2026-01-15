import pytest
from gtars.refget import (
    GlobalRefgetStore,
    StorageMode,
    digest_fasta,
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
        assert hasattr(result, "has_data")

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
        # print(store_raw.__repr__)
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

        md5 = "5f63cfaa3ef61f88c9635fb9d18ec945"
        seq = store.get_sequence_by_id(md5)
        assert seq is not None
        assert seq.metadata.length == 8

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

    def test_store_errors(self):
        """Test error conditions"""
        store = GlobalRefgetStore(StorageMode.Raw)

        # Test importing non-existent file
        with pytest.raises(Exception):
            store.import_fasta("nonexistent.fa")

        # Test getting non-existent sequence
        bogus_digest = "not_a_sequence"
        assert store.get_sequence_by_id(bogus_digest) is None

        # Test invalid substring parameters
        sha512 = "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        store.import_fasta("../../gtars/tests/data/fasta/base.fa")
        assert store.get_substring(sha512, 10, 5) is None  # end < start
        assert store.get_substring(sha512, 0, 100) is None  # end > length

    def test_store_collection_operations(self):
        """Test collection-related operations"""
        store = GlobalRefgetStore(StorageMode.Raw)
        fasta_path = "../../gtars/tests/data/fasta/base.fa"

        # Import sequences and get sequence by collection and name
        store.import_fasta(fasta_path)

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
        Test both get_seqs_bed_file (writes to file)
        and get_seqs_bed_file_to_vec (returns list) bindings.
        """

        with tempfile.TemporaryDirectory() as temp_dir_str:
            temp_dir = temp_dir_str

            fasta_content = ">chr1\n" "ATGCATGCATGC\n" ">chr2\n" "GGGGAAAA\n"
            temp_fasta_path = os.path.join(temp_dir, "test.fa")
            with open(temp_fasta_path, "w") as f:
                f.write(fasta_content)

            store = GlobalRefgetStore(StorageMode.Encoded)
            imported_collection = store.import_fasta(temp_fasta_path)
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
            store.get_seqs_bed_file(
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
            print("✓ get_seqs_bed_file binding test passed.")

            vec_result = store.get_seqs_bed_file_to_vec(
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

            print("✓ get_seqs_bed_file_to_vec binding test passed.")

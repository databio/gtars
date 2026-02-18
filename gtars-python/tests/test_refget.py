import pytest
from gtars.refget import (
    RefgetStore,
    StorageMode,
    digest_fasta,
    digest_sequence,
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
        fasta_path = "../tests/data/fasta/base.fa"

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
        assert hasattr(seq0, "sequence")

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
        """Test basic creation and operations of RefgetStore"""
        # Test store creation with both modes
        store_raw = RefgetStore.in_memory()
        store_encoded = RefgetStore.in_memory()

        # Test string representations
        # New repr format shows n_sequences and location
        assert repr(store_raw).startswith("RefgetStore")
        assert repr(store_encoded).startswith("RefgetStore")
        assert "n_sequences=0" in repr(store_raw)  # Empty store
        assert "memory-only" in repr(store_raw)  # No local path

    def test_store_import_and_retrieve(self):
        """Test importing FASTA and retrieving sequences"""
        store = RefgetStore.in_memory()
        fasta_path = "../tests/data/fasta/base.fa"
        # fasta_path = os.path.abspath("../tests/data/fasta/base.fa")
        # Import FASTA
        store.add_sequence_collection_from_fasta(fasta_path)

        # Get known sequence by SHA512t24u digest
        sha512 = "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        seq = store.get_sequence(sha512)
        assert seq is not None
        assert seq.metadata.length == 8

        md5 = "5f63cfaa3ef61f88c9635fb9d18ec945"
        seq = store.get_sequence(md5)
        assert seq is not None
        assert seq.metadata.length == 8

    def test_store_substring(self):
        """Test substring retrieval"""
        store = RefgetStore.in_memory()
        fasta_path = "../tests/data/fasta/base.fa"
        store.add_sequence_collection_from_fasta(fasta_path)

        # Get a substring from a known sequence
        sha512 = "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        substr = store.get_substring(sha512, 0, 4)
        assert substr is not None
        assert len(substr) == 4

    def test_store_persistence(self):
        """Test saving and loading store"""
        store = RefgetStore.in_memory()
        fasta_path = "../tests/data/fasta/base.fa"
        store.add_sequence_collection_from_fasta(fasta_path)

        # Save store to temporary directory - test with explicit template
        with tempfile.TemporaryDirectory() as tmpdir:
            store.write_store_to_dir(tmpdir, "sequences/%s2/%s.seq")

            # Load store back
            loaded_store = RefgetStore.open_local(tmpdir)

            # Verify same sequences exist
            sha512 = "iYtREV555dUFKg2_agSJW6suquUyPpMw"
            seq1 = store.get_sequence(sha512)
            seq2 = loaded_store.get_sequence(sha512)

            assert seq1 is not None
            assert seq2 is not None
            assert seq1.metadata.length == seq2.metadata.length
            assert seq1.metadata.sha512t24u == seq2.metadata.sha512t24u

        # Test with default template (no second argument)
        with tempfile.TemporaryDirectory() as tmpdir:
            store.write_store_to_dir(tmpdir)

            # Load store back
            loaded_store = RefgetStore.open_local(tmpdir)

            # Verify same sequences exist
            seq2 = loaded_store.get_sequence(sha512)
            assert seq2 is not None
            assert seq1.metadata.length == seq2.metadata.length
            assert seq1.metadata.sha512t24u == seq2.metadata.sha512t24u

    def test_store_errors(self):
        """Test error conditions"""
        store = RefgetStore.in_memory()

        # Test importing non-existent file
        with pytest.raises(Exception):
            store.add_sequence_collection_from_fasta("nonexistent.fa")

        # Test getting non-existent sequence
        bogus_digest = "not_a_sequence"
        assert store.get_sequence(bogus_digest) is None

        # Test invalid substring parameters
        sha512 = "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        store.add_sequence_collection_from_fasta("../tests/data/fasta/base.fa")
        assert store.get_substring(sha512, 10, 5) is None  # end < start
        assert store.get_substring(sha512, 0, 100) is None  # end > length

    def test_store_collection_operations(self):
        """Test collection-related operations"""
        store = RefgetStore.in_memory()
        fasta_path = "../tests/data/fasta/base.fa"

        # Import sequences and get sequence by collection and name
        store.add_sequence_collection_from_fasta(fasta_path)

        # Get sequence from default collection
        result = digest_fasta(fasta_path)
        seq = store.get_sequence_by_name(
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

            store = RefgetStore.in_memory()
            store.add_sequence_collection_from_fasta(temp_fasta_path)
            result = digest_fasta(temp_fasta_path)

            collection_digest = result.digest

            chr1_sha = sha512t24u_digest(b"ATGCATGCATGC")
            chr1_md5 = md5_digest(b"ATGCATGCATGC")
            chr2_sha = sha512t24u_digest(b"GGGGAAAA")
            chr2_md5 = md5_digest(b"GGGGAAAA")

            # Use only valid entries - errors now propagate instead of being silently skipped
            bed_content = (
                "chr1\t0\t5\n"
                "chr1\t8\t12\n"
                "chr2\t0\t4"
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
        fasta_path = "../tests/data/fasta/base.fa"

        # digest_fasta should not load sequence data
        result = digest_fasta(fasta_path)

        for seq_record in result.sequences:
            assert seq_record.sequence is None, "digest_fasta should not load sequence data"
            assert seq_record.decode() is None, "decode() should return None when data is None"

    def test_decode_with_loaded_data(self):
        """Test that decode() returns correct sequences when data is loaded"""
        fasta_path = "../tests/data/fasta/base.fa"

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
            assert seq_record.sequence is not None, "load_fasta should load sequence data"

            decoded = seq_record.decode()
            assert decoded is not None, "decode() should return Some when data is present"
            assert decoded == expected_seq, f"Decoded sequence for {expected_name} should match expected"

    def test_decode_with_store_sequences(self):
        """Test decode() with sequences retrieved from a store"""
        store = RefgetStore.in_memory()
        fasta_path = "../tests/data/fasta/base.fa"
        store.add_sequence_collection_from_fasta(fasta_path)

        # Get sequence by ID
        sha512 = "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        seq = store.get_sequence(sha512)

        assert seq is not None
        assert seq.sequence is not None, "Store should have loaded sequence data"

        decoded = seq.decode()
        assert decoded is not None
        assert decoded == "TTGGGGAA", "Should correctly decode sequence from encoded store"

    def test_decode_raw_vs_encoded_storage(self):
        """Test that decode() works with both Raw and Encoded storage modes"""
        fasta_path = "../tests/data/fasta/base.fa"

        # Test with Raw storage mode
        store_raw = RefgetStore.in_memory()
        store_raw.set_encoding_mode(StorageMode.Raw)
        store_raw.add_sequence_collection_from_fasta(fasta_path)
        sha512 = "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        seq_raw = store_raw.get_sequence(sha512)
        decoded_raw = seq_raw.decode()

        # Test with Encoded storage mode
        store_encoded = RefgetStore.in_memory()
        store_encoded.add_sequence_collection_from_fasta(fasta_path)
        seq_encoded = store_encoded.get_sequence(sha512)
        decoded_encoded = seq_encoded.decode()

        # Both should produce the same decoded sequence
        assert decoded_raw == decoded_encoded
        assert decoded_raw == "TTGGGGAA"

    def test_load_fasta_function(self):
        """Test the new load_fasta() function"""
        fasta_path = "../tests/data/fasta/base.fa"

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
            assert seq_record.sequence is not None, "load_fasta should load sequence data"
            assert seq_record.decode() is not None, "Should be able to decode loaded data"

        # Verify digests match digest_fasta
        digest_result = digest_fasta(fasta_path)
        assert result.digest == digest_result.digest
        assert result.lvl1.sequences_digest == digest_result.lvl1.sequences_digest
        assert result.lvl1.names_digest == digest_result.lvl1.names_digest
        assert result.lvl1.lengths_digest == digest_result.lvl1.lengths_digest

    def test_enable_persistence(self):
        """Test enable_persistence() flushes in-memory store to disk"""
        fasta_path = "../tests/data/fasta/base.fa"

        # Create in-memory store and add sequences
        store = RefgetStore.in_memory()
        store.add_sequence_collection_from_fasta(fasta_path)

        # Enable persistence to a temp directory
        with tempfile.TemporaryDirectory() as tmpdir:
            store.enable_persistence(tmpdir)

            # Check that rgstore.json was created (new format)
            assert os.path.exists(os.path.join(tmpdir, "rgstore.json"))

            # Load the store back and verify sequences are accessible
            loaded_store = RefgetStore.open_local(tmpdir)

            sha512 = "iYtREV555dUFKg2_agSJW6suquUyPpMw"
            seq1 = store.get_sequence(sha512)
            seq2 = loaded_store.get_sequence(sha512)

            assert seq1 is not None
            assert seq2 is not None
            assert seq1.metadata.sha512t24u == seq2.metadata.sha512t24u

    def test_disable_persistence(self):
        """Test disable_persistence() stops writing to disk"""
        fasta_path = "../tests/data/fasta/base.fa"

        with tempfile.TemporaryDirectory() as tmpdir:
            # Create disk-backed store
            store = RefgetStore.on_disk(tmpdir)

            # Disable persistence
            store.disable_persistence()

            # Store should still work in memory
            store.add_sequence_collection_from_fasta(fasta_path)

            sha512 = "iYtREV555dUFKg2_agSJW6suquUyPpMw"
            seq = store.get_sequence(sha512)
            assert seq is not None
            assert seq.metadata.length == 8

    def test_compute_fai(self):
        """Test compute_fai() returns correct FAI records"""
        from gtars.refget import compute_fai

        fasta_path = "../tests/data/fasta/base.fa"
        fai_records = compute_fai(fasta_path)

        # Should have 3 records for base.fa
        assert len(fai_records) == 3

        # Check first record structure
        rec = fai_records[0]
        assert rec.name == "chrX"
        assert rec.length == 8
        assert rec.fai is not None
        assert rec.fai.offset > 0  # Byte offset after header
        assert rec.fai.line_bases == 8  # All bases on one line
        assert rec.fai.line_bytes == 9  # 8 bases + newline

        # Check other records
        assert fai_records[1].name == "chr1"
        assert fai_records[1].length == 4
        assert fai_records[2].name == "chr2"
        assert fai_records[2].length == 4

    def test_compute_fai_gzipped(self):
        """Test compute_fai() returns fai=None for gzipped files"""
        from gtars.refget import compute_fai

        fasta_path = "../tests/data/fasta/base.fa.gz"
        fai_records = compute_fai(fasta_path)

        # Should still have 3 records
        assert len(fai_records) == 3

        # But FAI should be None for gzipped files
        for rec in fai_records:
            assert rec.fai is None

    def test_sequence_collection_pythonic_interface(self):
        """Test SequenceCollection supports len(), indexing, and iteration"""
        result = digest_fasta("../tests/data/fasta/base.fa")

        # Test __len__
        assert len(result) == 3

        # Test __getitem__ with positive index
        assert result[0].metadata.name == "chrX"
        assert result[1].metadata.name == "chr1"
        assert result[2].metadata.name == "chr2"

        # Test __getitem__ with negative index
        assert result[-1].metadata.name == "chr2"
        assert result[-3].metadata.name == "chrX"

        # Test index out of range
        with pytest.raises(IndexError):
            _ = result[10]

        # Test iteration
        names = [seq.metadata.name for seq in result]
        assert names == ["chrX", "chr1", "chr2"]

    def test_refget_store_pythonic_interface(self):
        """Test RefgetStore supports len() and iteration"""
        store = RefgetStore.in_memory()
        store.add_sequence_collection_from_fasta("../tests/data/fasta/base.fa")

        # Test __len__
        assert len(store) == 3

        # Test iteration yields SequenceMetadata
        count = 0
        for seq_meta in store:
            assert hasattr(seq_meta, 'name')
            assert hasattr(seq_meta, 'length')
            assert hasattr(seq_meta, 'sha512t24u')
            count += 1
        assert count == 3

    def test_collection_inspection_methods(self):
        """Test collection listing and metadata retrieval"""
        store = RefgetStore.in_memory()
        fasta_path = "../tests/data/fasta/base.fa"
        store.add_sequence_collection_from_fasta(fasta_path)

        # Get expected digest
        result = digest_fasta(fasta_path)
        expected_digest = result.digest

        # Test list_collections - now returns metadata objects
        collections = store.list_collections()
        assert len(collections) == 1
        assert collections[0].digest == expected_digest

        # Test get_collection_metadata
        meta = store.get_collection_metadata(expected_digest)
        assert meta is not None
        assert meta.digest == expected_digest
        assert meta.n_sequences == 3
        assert meta.names_digest == result.lvl1.names_digest
        assert meta.sequences_digest == result.lvl1.sequences_digest
        assert meta.lengths_digest == result.lvl1.lengths_digest

        # Test str/repr
        assert expected_digest in str(meta)
        assert "n_sequences=3" in repr(meta)

        # Test is_collection_loaded (in-memory store should be loaded)
        assert store.is_collection_loaded(expected_digest)

        # Test non-existent collection
        assert store.get_collection_metadata("nonexistent") is None

    def test_sequence_enumeration_methods(self):
        """Test list_sequences() method"""
        store = RefgetStore.in_memory()
        store.add_sequence_collection_from_fasta("../tests/data/fasta/base.fa")

        # Test list_sequences - returns metadata
        metadata_list = store.list_sequences()
        assert len(metadata_list) == 3
        for meta in metadata_list:
            assert hasattr(meta, 'name')
            assert hasattr(meta, 'length')
            assert hasattr(meta, 'sha512t24u')
            assert hasattr(meta, 'md5')

    def test_export_fasta(self):
        """Test export_fasta() exports full collection or subset"""
        store = RefgetStore.in_memory()
        fasta_path = "../tests/data/fasta/base.fa"
        store.add_sequence_collection_from_fasta(fasta_path)

        result = digest_fasta(fasta_path)
        collection_digest = result.digest

        with tempfile.TemporaryDirectory() as tmpdir:
            # Export all sequences
            output_path = os.path.join(tmpdir, "all.fa")
            store.export_fasta(collection_digest, output_path, None, None)

            with open(output_path) as f:
                content = f.read()
            assert ">chrX" in content
            assert ">chr1" in content
            assert ">chr2" in content

            # Export subset
            subset_path = os.path.join(tmpdir, "subset.fa")
            store.export_fasta(collection_digest, subset_path, ["chr1", "chr2"], 60)

            with open(subset_path) as f:
                subset_content = f.read()
            assert ">chrX" not in subset_content
            assert ">chr1" in subset_content
            assert ">chr2" in subset_content

    def test_export_fasta_by_digests(self):
        """Test export_fasta_by_digests() exports specific sequences by digest"""
        store = RefgetStore.in_memory()
        store.add_sequence_collection_from_fasta("../tests/data/fasta/base.fa")

        # Get digests for chr1 and chr2
        sha_chr1 = sha512t24u_digest(b"GGAA")  # chr1 sequence
        sha_chr2 = sha512t24u_digest(b"GCGC")  # chr2 sequence

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "by_digest.fa")
            store.export_fasta_by_digests([sha_chr1, sha_chr2], output_path, None)

            with open(output_path) as f:
                content = f.read()

            # Should contain chr1 and chr2 sequences
            assert "GGAA" in content
            assert "GCGC" in content
            # chrX sequence should not be present
            assert "TTGGGGAA" not in content

    def test_sequence_collection_write_fasta(self):
        """Test SequenceCollection.write_fasta() method"""
        # load_fasta returns collection with data
        collection = load_fasta("../tests/data/fasta/base.fa")

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "written.fa")
            collection.write_fasta(output_path, 80)  # default line width

            with open(output_path) as f:
                content = f.read()

            assert ">chrX" in content
            assert "TTGGGGAA" in content
            assert ">chr1" in content
            assert "GGAA" in content

            # Test with custom line width
            output_path2 = os.path.join(tmpdir, "written2.fa")
            collection.write_fasta(output_path2, 4)

            with open(output_path2) as f:
                content2 = f.read()
            # With line_width=4, TTGGGGAA should be split
            assert "TTGG\n" in content2 or "GGAA\n" in content2

    def test_quiet_mode(self):
        """Test store quiet mode suppresses output"""
        store = RefgetStore.in_memory()

        # Test getter
        assert store.quiet == False

        # Test setter
        store.set_quiet(True)
        assert store.quiet == True

        store.set_quiet(False)
        assert store.quiet == False

    def test_get_collection_method(self):
        """Test store.get_collection() returns SequenceCollection with loaded sequences"""
        store = RefgetStore.in_memory()
        fasta_path = "../tests/data/fasta/base.fa"
        store.add_sequence_collection_from_fasta(fasta_path)

        # Get the collection digest
        collections = store.list_collections()
        assert len(collections) == 1
        digest = collections[0].digest

        # Get full collection with sequences
        coll = store.get_collection(digest)
        assert hasattr(coll, 'digest')
        assert hasattr(coll, 'sequences')
        assert hasattr(coll, 'lvl1')
        assert len(coll.sequences) == 3
        # Sequences should have data loaded
        for seq in coll.sequences:
            assert seq.decode() is not None

    def test_iter_collections(self):
        """Test iter_collections() returns loaded collections"""
        store = RefgetStore.in_memory()
        fasta_path = "../tests/data/fasta/base.fa"
        store.add_sequence_collection_from_fasta(fasta_path)

        collections = store.iter_collections()
        assert len(collections) == 1

        coll = collections[0]
        assert hasattr(coll, 'digest')
        assert hasattr(coll, 'sequences')
        assert len(coll.sequences) == 3
        # All sequences should have data loaded
        for seq in coll.sequences:
            assert seq.decode() is not None

    def test_iter_sequences(self):
        """Test iter_sequences() returns loaded sequences"""
        store = RefgetStore.in_memory()
        fasta_path = "../tests/data/fasta/base.fa"
        store.add_sequence_collection_from_fasta(fasta_path)

        sequences = store.iter_sequences()
        assert len(sequences) == 3

        # All sequences should have data loaded
        for seq in sequences:
            assert hasattr(seq, 'metadata')
            assert seq.decode() is not None

    def test_string_representations(self):
        """Test __str__ and __repr__ for all types"""
        from gtars.refget import compute_fai

        fasta_path = "../tests/data/fasta/base.fa"

        # SequenceCollection
        coll = digest_fasta(fasta_path)
        assert "3 sequences" in str(coll)
        assert "SequenceCollection" in repr(coll)

        # SequenceRecord
        rec = coll.sequences[0]
        assert "chrX" in str(rec)
        assert "SequenceRecord" in repr(rec)

        # SequenceMetadata
        meta = rec.metadata
        assert "chrX" in str(meta)
        assert "SequenceMetadata" in repr(meta)

        # SeqColDigestLvl1
        lvl1 = coll.lvl1
        assert "SeqColDigestLvl1" in str(lvl1)
        assert "SeqColDigestLvl1" in repr(lvl1)

        # FaiRecord
        fai_records = compute_fai(fasta_path)
        fai = fai_records[0]
        assert "chrX" in str(fai)
        assert "FaiRecord" in repr(fai)

        # FaiMetadata
        if fai.fai:
            assert "FaiMetadata" in str(fai.fai)
            assert "FaiMetadata" in repr(fai.fai)

        # AlphabetType - just check it has a string representation
        assert str(meta.alphabet) is not None

        # RetrievedSequence
        rs = RetrievedSequence(sequence="ATGC", chrom_name="chr1", start=0, end=4)
        assert "chr1" in str(rs)
        assert "RetrievedSequence" in repr(rs)

        # RefgetStore
        store = RefgetStore.in_memory()
        assert "RefgetStore" in repr(store)
        assert "memory-only" in repr(store)

    def test_add_standalone_sequence(self):
        """Test adding a sequence without a collection."""
        store = RefgetStore.in_memory()

        # Nameless sequence
        seq = digest_sequence(b"ACGTACGTACGT")
        store.add_sequence(seq)
        retrieved = store.get_sequence(seq.metadata.sha512t24u)
        assert retrieved is not None
        assert retrieved.metadata.length == 12

        # Named sequence
        seq2 = digest_sequence(b"TTTTAAAA", name="my_seq")
        store.add_sequence(seq2)
        retrieved2 = store.get_sequence(seq2.metadata.sha512t24u)
        assert retrieved2 is not None

        # Adding again should not error (skip duplicate)
        store.add_sequence(seq)

        # Force add should also work
        store.add_sequence(seq, force=True)

    def test_digest_sequence_name_optional(self):
        """Test that digest_sequence works with and without name."""
        # No name
        seq1 = digest_sequence(b"ACGT")
        assert seq1.metadata.length == 4
        assert seq1.metadata.name == ""

        # With name
        seq2 = digest_sequence(b"ACGT", name="chr1")
        assert seq2.metadata.name == "chr1"

        # Same data should produce same digest regardless of name
        assert seq1.metadata.sha512t24u == seq2.metadata.sha512t24u

    def test_sequence_aliases(self):
        """Test sequence alias add, forward lookup, reverse lookup, and remove."""
        store = RefgetStore.in_memory()
        seq = digest_sequence(b"ACGTACGT", name="chr1")
        store.add_sequence(seq)

        digest = seq.metadata.sha512t24u
        store.add_sequence_alias("ncbi", "NC_000001.11", digest)
        store.add_sequence_alias("ucsc", "chr1", digest)

        # Forward lookup
        found = store.get_sequence_by_alias("ncbi", "NC_000001.11")
        assert found is not None
        assert found.metadata.name == "chr1"

        # Reverse lookup
        aliases = store.get_aliases_for_sequence(digest)
        assert len(aliases) == 2
        assert ("ncbi", "NC_000001.11") in aliases
        assert ("ucsc", "chr1") in aliases

        # List namespaces
        ns = store.list_sequence_alias_namespaces()
        assert "ncbi" in ns
        assert "ucsc" in ns

        # List aliases in namespace
        alias_list = store.list_sequence_aliases("ncbi")
        assert "NC_000001.11" in alias_list

        # Remove alias
        assert store.remove_sequence_alias("ncbi", "NC_000001.11")
        assert store.get_sequence_by_alias("ncbi", "NC_000001.11") is None

    def test_collection_aliases(self):
        """Test collection alias add, forward lookup, and reverse lookup."""
        store = RefgetStore.in_memory()
        meta, _ = store.add_sequence_collection_from_fasta("../tests/data/fasta/base.fa")

        store.add_collection_alias("ucsc", "hg38", meta.digest)

        coll = store.get_collection_by_alias("ucsc", "hg38")
        assert coll is not None
        assert coll.digest == meta.digest

        aliases = store.get_aliases_for_collection(meta.digest)
        assert ("ucsc", "hg38") in aliases

    def test_alias_persistence(self, tmp_path):
        """Test that aliases persist across store save/reload."""
        store_path = tmp_path / "store"
        store = RefgetStore.on_disk(str(store_path))
        meta, _ = store.add_sequence_collection_from_fasta("../tests/data/fasta/base.fa")

        store.add_collection_alias("ucsc", "hg38", meta.digest)

        # Reload
        store2 = RefgetStore.open_local(str(store_path))
        assert store2.get_collection_by_alias("ucsc", "hg38") is not None


# =========================================================================
# FHR Metadata Tests
# =========================================================================


def test_fhr_metadata_empty_by_default():
    from gtars.refget import RefgetStore, FhrMetadata

    store = RefgetStore.in_memory()
    meta, _ = store.add_sequence_collection_from_fasta("../tests/data/fasta/base.fa")

    # FHR metadata always works -- no enable step needed
    assert store.get_fhr_metadata(meta.digest) is None
    assert store.list_fhr_metadata() == []


def test_fhr_metadata_set_get():
    from gtars.refget import RefgetStore, FhrMetadata

    store = RefgetStore.in_memory()
    meta, _ = store.add_sequence_collection_from_fasta("../tests/data/fasta/base.fa")

    fhr = FhrMetadata(genome="Homo sapiens", version="GRCh38", masking="soft-masked")
    store.set_fhr_metadata(meta.digest, fhr)

    retrieved = store.get_fhr_metadata(meta.digest)
    assert retrieved is not None
    assert retrieved.genome == "Homo sapiens"
    assert retrieved.version == "GRCh38"
    assert retrieved.masking == "soft-masked"


def test_fhr_metadata_none_for_missing():
    from gtars.refget import RefgetStore

    store = RefgetStore.in_memory()
    meta, _ = store.add_sequence_collection_from_fasta("../tests/data/fasta/base.fa")
    assert store.get_fhr_metadata(meta.digest) is None


def test_fhr_metadata_to_dict():
    from gtars.refget import FhrMetadata

    fhr = FhrMetadata(genome="Test", version="1.0")
    d = fhr.to_dict()
    assert d["genome"] == "Test"
    assert d["version"] == "1.0"


def test_fhr_metadata_from_json(tmp_path):
    import json
    from gtars.refget import FhrMetadata

    fhr_file = tmp_path / "test.fhr.json"
    fhr_file.write_text(
        json.dumps(
            {
                "genome": "Homo sapiens",
                "version": "GRCh38.p14",
                "taxon": {
                    "name": "Homo sapiens",
                    "uri": "https://identifiers.org/taxonomy:9606",
                },
                "masking": "soft-masked",
            }
        )
    )

    fhr = FhrMetadata.from_json(str(fhr_file))
    assert fhr.genome == "Homo sapiens"
    assert fhr.version == "GRCh38.p14"


def test_fhr_metadata_persistence(tmp_path):
    import json
    from gtars.refget import RefgetStore, FhrMetadata

    store_path = tmp_path / "store"
    store = RefgetStore.on_disk(str(store_path))
    meta, _ = store.add_sequence_collection_from_fasta("../tests/data/fasta/base.fa")

    fhr = FhrMetadata(genome="Test", version="1.0")
    store.set_fhr_metadata(meta.digest, fhr)

    # Verify sidecar file was written
    fhr_path = store_path / "collections" / f"{meta.digest}.fhr.json"
    assert fhr_path.exists()
    data = json.loads(fhr_path.read_text())
    assert data["genome"] == "Test"

    # Reload -- FHR metadata always loaded
    store2 = RefgetStore.open_local(str(store_path))
    retrieved = store2.get_fhr_metadata(meta.digest)
    assert retrieved is not None
    assert retrieved.genome == "Test"


def test_fhr_metadata_spec_fields():
    """Test creating FhrMetadata with all FHR 1.0 spec fields including new ones."""
    import json
    from gtars.refget import FhrMetadata

    fhr = FhrMetadata(
        genome="Bombas huntii",
        version="0.0.1",
        schemaVersion=1.0,
        masking="soft-masked",
        voucherSpecimen="Located in Freezer 33",
        documentation="Built assembly from...",
        identifier=["beetlebase:TC010103"],
        scholarlyArticle="10.1371/journal.pntd.0008755",
        funding="NIH R01",
        vitalStats={
            "L50": 42,
            "N50": 1000000,
            "totalBasePairs": 3000000000,
            "readTechnology": "hifi",
        },
        accessionID={"name": "PBARC", "url": "https://example.com"},
    )

    assert fhr.genome == "Bombas huntii"
    assert fhr.voucher_specimen == "Located in Freezer 33"
    assert fhr.documentation == "Built assembly from..."
    assert fhr.identifier == ["beetlebase:TC010103"]
    assert fhr.scholarly_article == "10.1371/journal.pntd.0008755"
    assert fhr.funding == "NIH R01"

    # Round-trip through dict
    d = fhr.to_dict()
    assert d["voucherSpecimen"] == "Located in Freezer 33"
    assert d["documentation"] == "Built assembly from..."
    assert d["scholarlyArticle"] == "10.1371/journal.pntd.0008755"
    assert d["funding"] == "NIH R01"
    assert d["vitalStats"]["L50"] == 42
    assert d["vitalStats"]["N50"] == 1000000
    assert d["accessionID"]["name"] == "PBARC"
    # seqcol_digest should NOT appear
    assert "seqcolDigest" not in d
    assert "seqcol_digest" not in d


def test_add_sequence_collection_basic():
    """Test adding a pre-built SequenceCollection to a store."""
    fasta_path = "../tests/data/fasta/base.fa"

    # Build a collection using digest_fasta (no FASTA file write needed)
    collection = digest_fasta(fasta_path)

    store = RefgetStore.in_memory()
    store.add_sequence_collection(collection)

    # Collection should appear in list_collections
    collections = store.list_collections()
    assert len(collections) == 1
    assert collections[0].digest == collection.digest


def test_add_sequence_collection_sequences_retrievable():
    """Test that sequences from add_sequence_collection are retrievable by metadata.

    digest_fasta() produces a collection with DIGEST_ONLY sequences (no raw data),
    so the sequences are stored as Stubs and get_sequence returns metadata but
    get_sequence(digest) returns None for unloaded Stubs. Verify metadata is present
    via list_collections and collection digest lookup.
    """
    fasta_path = "../tests/data/fasta/base.fa"

    collection = digest_fasta(fasta_path)
    store = RefgetStore.in_memory()
    store.add_sequence_collection(collection)

    # The collection metadata should be present and match
    meta = store.get_collection_metadata(collection.digest)
    assert meta is not None, "Collection metadata should be present in store"
    assert meta.digest == collection.digest
    assert meta.n_sequences == len(collection.sequences)

    # Use load_fasta to get a collection with Full sequences, verify those ARE retrievable
    from gtars.refget import load_fasta

    full_collection = load_fasta(fasta_path)
    store2 = RefgetStore.in_memory()
    store2.add_sequence_collection(full_collection)

    for seq_record in full_collection.sequences:
        digest = seq_record.metadata.sha512t24u
        retrieved = store2.get_sequence(digest)
        assert retrieved is not None, f"Sequence {digest} should be retrievable from full collection"
        assert retrieved.metadata.length == seq_record.metadata.length


def test_add_sequence_collection_skip_duplicate():
    """Test that add_sequence_collection skips duplicates by default."""
    fasta_path = "../tests/data/fasta/base.fa"

    collection = digest_fasta(fasta_path)
    store = RefgetStore.in_memory()

    # Add the same collection twice - should not raise
    store.add_sequence_collection(collection)
    store.add_sequence_collection(collection)

    # Should still have exactly one collection
    assert len(store.list_collections()) == 1


def test_add_sequence_collection_force():
    """Test that add_sequence_collection with force=True overwrites existing."""
    fasta_path = "../tests/data/fasta/base.fa"

    collection = digest_fasta(fasta_path)
    store = RefgetStore.in_memory()

    # Add once, then force-overwrite - should not raise
    store.add_sequence_collection(collection)
    store.add_sequence_collection(collection, force=True)

    # Should still have exactly one collection
    assert len(store.list_collections()) == 1


def test_is_persisting_false_for_memory_store():
    """Test that in-memory stores report is_persisting as False."""
    store = RefgetStore.in_memory()
    assert store.is_persisting is False


def test_is_persisting_true_after_enable(tmp_path):
    """Test that is_persisting becomes True after enable_persistence."""
    store = RefgetStore.in_memory()
    assert store.is_persisting is False

    store.enable_persistence(str(tmp_path))
    assert store.is_persisting is True


def test_is_persisting_false_after_disable(tmp_path):
    """Test that is_persisting becomes False after disable_persistence."""
    store = RefgetStore.in_memory()
    store.enable_persistence(str(tmp_path))
    assert store.is_persisting is True

    store.disable_persistence()
    assert store.is_persisting is False

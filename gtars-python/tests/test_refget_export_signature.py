"""Binding tests for export_fasta's optional-argument signature.

Regression test for a pyo3 binding bug: export_fasta (and friends) took
Option<...> parameters in Rust but had no #[pyo3(signature = ...)] attribute,
so Python required all arguments positionally even though the .pyi stub
declared them optional. See:
    TypeError: RefgetStore.export_fasta() missing 1 required positional
    argument: 'line_width'
"""

from gtars.refget import RefgetStore


def _collection(store, path, body):
    path.write_text(body)
    metadata, _ = store.add_sequence_collection_from_fasta(str(path))
    return metadata.digest


def test_export_fasta_with_no_optional_args(tmp_path):
    store = RefgetStore.in_memory()
    digest = _collection(store, tmp_path / "a.fa", ">chr1\nACGT\n>chr2\nGGCC\n")

    out_path = tmp_path / "out_default.fa"
    store.export_fasta(digest, str(out_path))

    assert out_path.exists()
    assert out_path.stat().st_size > 0


def test_export_fasta_with_explicit_none(tmp_path):
    store = RefgetStore.in_memory()
    digest = _collection(store, tmp_path / "a.fa", ">chr1\nACGT\n>chr2\nGGCC\n")

    out_path = tmp_path / "out_none.fa"
    store.export_fasta(digest, str(out_path), None)

    assert out_path.exists()
    assert out_path.stat().st_size > 0


def test_export_fasta_with_line_width_keyword(tmp_path):
    store = RefgetStore.in_memory()
    digest = _collection(store, tmp_path / "a.fa", ">chr1\nACGT\n>chr2\nGGCC\n")

    out_path = tmp_path / "out_kw.fa"
    store.export_fasta(digest, str(out_path), line_width=60)

    assert out_path.exists()
    assert out_path.stat().st_size > 0


def test_export_fasta_with_sequence_names_positional(tmp_path):
    store = RefgetStore.in_memory()
    digest = _collection(store, tmp_path / "a.fa", ">chr1\nACGT\n>chr2\nGGCC\n")

    out_path = tmp_path / "out_names.fa"
    store.export_fasta(digest, str(out_path), ["chr1"])

    assert out_path.exists()
    content = out_path.read_text()
    assert "chr1" in content
    assert "chr2" not in content


def test_export_fasta_by_digests_with_no_optional_args(tmp_path):
    store = RefgetStore.in_memory()
    digest = _collection(store, tmp_path / "a.fa", ">chr1\nACGT\n")
    collection = store.get_collection(digest)
    seq_digest = collection.sequences[0].metadata.sha512t24u

    out_path = tmp_path / "out_by_digest.fa"
    store.export_fasta_by_digests([seq_digest], str(out_path))

    assert out_path.exists()
    assert out_path.stat().st_size > 0

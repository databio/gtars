"""Round-trip tests for RefgetStore: decoded output must equal the original
(uppercased) FASTA input, through every read path, in every storage mode the
Python binding exposes.

Mirrors the in-memory (step 9) and on-disk (step 11) cases of
gtars-refget/tests/roundtrip_store.rs -- the fixture table below is kept in
sync by eye with FIXTURE_RECORDS there, so Python and Rust behavior can be
compared directly.

`PyStorageMode` (gtars-python/src/refget/mod.rs) has no `Zstd` variant, so
Zstd is covered in Rust only; that gap is not introduced by this file.

NOTE: as of this writing, `gtars-py` does not compile on this branch. The
`From<StorageMode> for PyStorageMode` match in
gtars-python/src/refget/mod.rs (~line 978) has no arm for the (separately,
recently added) `StorageMode::Zstd` variant:

    error[E0004]: non-exhaustive patterns: `StorageMode::Zstd` not covered

That is a pre-existing compile break unrelated to this plan (a different,
concurrently-developed plan added `StorageMode::Zstd`). This file is written
test-first anyway, per the plan; it cannot actually be run (`maturin develop`
fails) until that match is made exhaustive elsewhere.
"""

import pytest

from gtars.refget import RefgetStore, StorageMode, sha512t24u_digest

# Keep in sync by eye with `fixture_records()` in
# gtars-refget/tests/roundtrip_store.rs.
FIXTURE_RECORDS = [
    ("dna2bit", "ACGTACGTACacgt"),
    ("dna3bit", "ACGTNNNRYRYACGT"),
    ("iupac_all", "ACGTURYSWKMBDHVNACGTURYSWKMBDHVNACGTURYSWKMBDHVNA"),
    ("iupac_dh", "DDDHHHVVVDHVDHV"),
    ("rna", "ACGUACGUUUUAGCU"),
    ("rna_lower", "acguacguuuu"),
    ("protein_std", "MEFILPQACDEGHKNRSTVWY*"),
    ("selenoprotein", "MPRLLUGSEEAUVLLK*"),
    ("protein_bzoj", "MEBZOJXBZOJK*"),
    ("protein_b_only", "MEFBBBILPQ"),
    ("protein_iupac_letters", "MCVTYHNGTGYC"),
    ("protein_gaps", "MEF-ILP.QX*"),
    ("ascii_misc", "MEF123ILP"),
]


def expected_sequences():
    """Uppercased sequence for each fixture: what the ingest pipeline
    uppercases every sequence byte to before it digests, guesses or encodes.
    """
    return {name: seq.upper() for name, seq in FIXTURE_RECORDS}


def write_fixture_fasta(path):
    """Write the fixture records wrapped at 10 chars/line, matching the Rust
    fixture file so the multi-line join path runs in both languages."""
    with open(path, "w") as f:
        for name, seq in FIXTURE_RECORDS:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 10):
                f.write(seq[i : i + 10] + "\n")


def assert_lossless(label, name, got, meta, expected, failures):
    if got != expected:
        failures.append(
            f"{label}/{name}: got {got!r} (len {len(got)}), "
            f"expected {expected!r} (len {len(expected)})"
        )
        return
    got_sha = sha512t24u_digest(got)
    if got_sha != meta.sha512t24u:
        failures.append(
            f"{label}/{name}: sha512t24u(got)={got_sha} != "
            f"metadata.sha512t24u={meta.sha512t24u}"
        )
    if len(got) != meta.length:
        failures.append(
            f"{label}/{name}: len(got)={len(got)} != metadata.length={meta.length}"
        )


def read_paths(store, meta, expected, label, failures):
    digest = meta.sha512t24u
    length = len(expected)

    # Full stream_sequence().read_all()
    got = store.stream_sequence(digest).read_all()
    assert_lossless(
        f"{label} stream_sequence(full)", meta.name, got, meta, expected, failures
    )

    # Chunked iteration, to cover the O(1)-memory chunked path.
    got_chunked = "".join(store.stream_sequence(digest, chunk_size=3))
    if got_chunked != expected:
        failures.append(
            f"{label}/{meta.name}: chunked stream_sequence = {got_chunked!r}, "
            f"want {expected!r}"
        )

    # get_substring, full range.
    got_sub = store.get_substring(digest, 0, length)
    assert_lossless(
        f"{label} get_substring(full)", meta.name, got_sub, meta, expected, failures
    )

    # get_substring, a sub-range, when the sequence is long enough.
    if length >= 5:
        s, e = 3, length - 2
        got_range = store.get_substring(digest, s, e)
        want = expected[s:e]
        if got_range != want:
            failures.append(
                f"{label}/{meta.name}: get_substring({s},{e}) = {got_range!r}, "
                f"want {want!r}"
            )


def build_in_memory(mode, fasta_path):
    store = RefgetStore.in_memory()
    store.set_encoding_mode(mode)
    store.add_sequence_collection_from_fasta(str(fasta_path))
    return store


def run_all(store, label, failures):
    expected_map = expected_sequences()
    metas = store.list_sequences()
    assert len(metas) == len(FIXTURE_RECORDS), (
        f"{label}: expected {len(FIXTURE_RECORDS)} sequences, found {len(metas)}"
    )
    for meta in metas:
        expected = expected_map[meta.name]
        read_paths(store, meta, expected, label, failures)


@pytest.mark.parametrize("mode", [StorageMode.Encoded, StorageMode.Raw])
def test_in_memory_is_lossless(tmp_path, mode):
    fasta_path = tmp_path / "fixtures.fa"
    write_fixture_fasta(fasta_path)

    store = build_in_memory(mode, fasta_path)

    failures = []
    run_all(store, f"in_memory/{mode}", failures)
    assert not failures, "\n".join(failures)


@pytest.mark.parametrize("mode", [StorageMode.Encoded, StorageMode.Raw])
def test_on_disk_is_lossless(tmp_path, mode):
    fasta_path = tmp_path / "fixtures.fa"
    write_fixture_fasta(fasta_path)

    store_dir = tmp_path / "store"
    store = RefgetStore.on_disk(str(store_dir))
    store.set_encoding_mode(mode)
    store.add_sequence_collection_from_fasta(str(fasta_path))

    failures = []

    # Reopen fresh: the mode must survive the round trip through
    # rgstore.json, and reads go through the disk-backed stub path.
    reopened = RefgetStore.open_local(str(store_dir))
    assert str(reopened.storage_mode) == str(mode), (
        "storage_mode did not survive reopen"
    )
    run_all(reopened, f"on_disk/{mode}/stub", failures)

    # Load everything resident and re-run, to cover the resident-after-load
    # branch.
    reopened.load_all_sequences()
    run_all(reopened, f"on_disk/{mode}/resident", failures)

    assert not failures, "\n".join(failures)

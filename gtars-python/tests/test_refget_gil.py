"""gtars.refget releases the GIL during slow Rust work.

A background thread runs a pure-Python loop and records the longest pause
between two of its steps. If a call on the main thread held the GIL, that
thread could not run at all during the call, so the longest pause would be
about as long as the call itself.
"""

import functools
import http.server
import random
import sys
import threading
import time

import pytest

from gtars.refget import RefgetStore, digest_fasta


class Ticker(threading.Thread):
    """Loop in Python and record the longest pause between steps."""

    def __init__(self):
        super().__init__(daemon=True)
        self.longest = 0.0
        self.stop = threading.Event()

    def run(self):
        last = time.perf_counter()
        while not self.stop.is_set():
            now = time.perf_counter()
            self.longest = max(self.longest, now - last)
            last = now


def assert_releases_gil(fn):
    """Run fn on this thread and check the ticker thread kept running."""
    old_interval = sys.getswitchinterval()
    sys.setswitchinterval(0.001)
    ticker = Ticker()
    ticker.start()
    try:
        time.sleep(0.02)
        ticker.longest = 0.0
        start = time.perf_counter()
        result = fn()
        elapsed = time.perf_counter() - start
    finally:
        ticker.stop.set()
        ticker.join()
        sys.setswitchinterval(old_interval)
    assert elapsed > 0.02, f"call too fast to measure ({elapsed:.3f}s)"
    assert ticker.longest < elapsed / 2, (
        f"other threads froze for {ticker.longest:.3f}s of a {elapsed:.3f}s call"
    )
    return result


@pytest.fixture(scope="module")
def big_fasta(tmp_path_factory):
    """A 64 Mbp FASTA, big enough that digesting it takes a few tenths of a second."""
    path = tmp_path_factory.mktemp("gil") / "big.fa"
    rng = random.Random(42)
    to_acgt = bytes.maketrans(bytes(range(256)), b"ACGT" * 64)
    with open(path, "wb") as f:
        for i in range(4):
            f.write(b">chr%d\n" % i + rng.randbytes(16_000_000).translate(to_acgt) + b"\n")
    return str(path)


def test_digest_import_and_export_release_gil(big_fasta, tmp_path):
    coll = assert_releases_gil(lambda: digest_fasta(big_fasta))
    assert len(coll) == 4

    store = RefgetStore.on_disk(str(tmp_path / "store"))
    meta, was_new = assert_releases_gil(
        lambda: store.add_sequence_collection_from_fasta(big_fasta)
    )
    assert was_new and meta.digest == coll.digest

    out = tmp_path / "out.fa"
    assert_releases_gil(
        lambda: store.export_fasta(str(out), collections=meta.digest, label="name")
    )
    assert digest_fasta(str(out)).digest == coll.digest


class SlowHandler(http.server.SimpleHTTPRequestHandler):
    """Serve static files, waiting before each response like a slow server."""

    def do_GET(self):
        time.sleep(0.2)
        super().do_GET()

    def log_message(self, *args):
        pass


def test_remote_fetch_releases_gil(tmp_path):
    """The server runs in a thread of this process, so it can only answer while
    the waiting call has released the GIL."""
    served = tmp_path / "served"
    fasta = tmp_path / "a.fa"
    fasta.write_text(">chr1\nACGTACGTAC\n>chr2\nGGCCGGCC\n")
    writer = RefgetStore.on_disk(str(served))
    meta, _ = writer.add_sequence_collection_from_fasta(str(fasta))
    writer.write()

    handler = functools.partial(SlowHandler, directory=str(served))
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    try:
        url = f"http://127.0.0.1:{server.server_address[1]}"
        store = assert_releases_gil(
            lambda: RefgetStore.open_remote(str(tmp_path / "cache"), url)
        )
        seq = assert_releases_gil(lambda: store.get_sequence_by_name(meta.digest, "chr1"))
        assert seq.decode() == "ACGTACGTAC"
    finally:
        server.shutdown()
        server.server_close()


def test_threads_share_one_store(big_fasta, tmp_path):
    """Threads using one store take turns: no deadlock, no 'Already borrowed'
    error, and every read sees correct data while a slow import runs."""
    store = RefgetStore.on_disk(str(tmp_path / "store"))
    small = tmp_path / "small.fa"
    small.write_text(">s1\nACGTACGT\n")
    meta, _ = store.add_sequence_collection_from_fasta(str(small))

    errors = []
    reads = []

    def run(fn):
        try:
            fn()
        except Exception as e:
            errors.append(e)

    def import_big():
        store.add_sequence_collection_from_fasta(big_fasta)

    def read_many():
        for _ in range(20):
            reads.append(store.get_sequence_by_name(meta.digest, "s1").decode())
            store.stats()  # a quick call that keeps the GIL

    threads = [threading.Thread(target=run, args=(import_big,))]
    threads += [threading.Thread(target=run, args=(read_many,)) for _ in range(3)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=120)

    assert not any(t.is_alive() for t in threads), "threads deadlocked"
    assert not errors, errors
    assert reads == ["ACGTACGT"] * 60
    assert store.stats()["n_collections"] == "2"

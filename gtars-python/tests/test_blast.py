"""gtars.blast: build, verify, view and nstats over an on-disk RefgetStore."""

import os
import shutil
import subprocess
import tempfile
import threading
from types import SimpleNamespace
from pathlib import Path

import pytest

from gtars.refget import RefgetStore

blast = pytest.importorskip("gtars.blast")

FASTA = Path(__file__).resolve().parents[2] / "tests" / "data" / "fasta"
GTARS_BIN = Path(os.environ.get("GTARS_BIN", ""))


@pytest.fixture(autouse=True)
def isolated_env(tmp_path, monkeypatch):
    # N counts are cached; keep tests off the real cache, and pin the database date
    monkeypatch.setenv("GTARS_RGBLAST_CACHE", str(tmp_path / "nstats-cache"))
    monkeypatch.setenv("SOURCE_DATE_EPOCH", "0")


@pytest.fixture
def store(tmp_path):
    """A store with `base.fa`, `different_names.fa` (same sequences) and a collection of its own."""
    path = tmp_path / "store"
    s = RefgetStore.on_disk(str(path))
    digests = {}
    uniq = tmp_path / "uniq.fa"
    uniq.write_text(">u1\nACGTACGTACGTTTGACCA\n>u2\nCCCGGGATATATTT\n")
    for name, fa in [("base", "base.fa"), ("other", "different_names.fa"), ("uniq", None)]:
        if fa is not None:  # copies: importing writes an .rgsi next to the FASTA
            shutil.copy(FASTA / fa, tmp_path / fa)
            fa = tmp_path / fa
        else:
            fa = uniq
        meta, _ = s.add_sequence_collection_from_fasta(str(fa))
        digests[name] = meta.digest
    del s  # the on-disk store is complete once the collections are added
    return SimpleNamespace(path=str(path), digests=digests)


def files(prefix: Path):
    return sorted(p for p in prefix.parent.iterdir() if p.name.startswith(prefix.name))


def test_build_returns_summary_and_writes_files(store, tmp_path):
    out = tmp_path / "db"
    r = blast.build(store.path, out, all=True, threads=2)
    assert isinstance(r, dict)
    assert r["preflight"]["subjects"] > 0 and r["files"]
    assert any(p.suffix == ".nin" for p in files(out))
    assert blast.verify(out, store.path, threads=2)["ok"] is True


def test_dry_run_writes_nothing(store, tmp_path):
    out = tmp_path / "db"
    r = blast.build(store.path, out, all=True, dry_run=True)
    assert isinstance(r, dict) and r
    assert files(out) == []


def test_progress_events(store, tmp_path):
    events = []

    def cb(ev):
        events.append((ev["event"], threading.current_thread() is threading.main_thread(), "message" in ev))

    blast.build(store.path, tmp_path / "db", all=True, threads=2, progress=cb)
    names = [e[0] for e in events]
    assert names[0] == "opened" and names[-1] == "finished"
    assert "planned" in names
    assert all(e[1] and e[2] for e in events)


def test_progress_true_prints_to_stderr(store, tmp_path, capsys):
    blast.build(store.path, tmp_path / "db", all=True, threads=2, progress=True)
    err = capsys.readouterr().err
    assert err.strip()


def test_callback_exception_cancels_and_propagates(store, tmp_path):
    out = tmp_path / "db"

    def cb(ev):
        if ev["event"] == "planned":
            raise RuntimeError("stop")

    with pytest.raises(RuntimeError, match="stop"):
        blast.build(store.path, out, all=True, progress=cb)
    assert files(out) == []


def test_keyboard_interrupt_from_callback(store, tmp_path):
    def cb(ev):
        if ev["event"] == "planned":
            raise KeyboardInterrupt

    with pytest.raises(KeyboardInterrupt):
        blast.build(store.path, tmp_path / "db", all=True, progress=cb)


def test_errors_carry_kind(store, tmp_path):
    out = tmp_path / "db"
    with pytest.raises(ValueError, match="cannot be used with") as e:
        blast.build(store.path, out, all=True, collection="x")
    assert e.value.kind == "invalid_input"
    with pytest.raises(FileNotFoundError) as e:
        blast.build(store.path, out, collection="no-such-collection")
    assert e.value.kind == "not_found"
    blast.build(store.path, out, all=True, threads=2)
    with pytest.raises(FileExistsError) as e:
        blast.build(store.path, out, all=True, threads=2)
    assert e.value.kind == "conflict"
    blast.build(store.path, out, all=True, threads=2, force=True)
    with pytest.raises(ValueError, match="only applies to remote stores"):
        blast.build(store.path, tmp_path / "db2", all=True, fetch_jobs=4)


def test_selection_accepts_one_value_or_a_list_and_ints(store, tmp_path):
    d = store.digests
    a = blast.build(store.path, tmp_path / "a", collection=d["base"], threads=2)
    b = blast.build(store.path, tmp_path / "b", collection=[d["base"], d["uniq"]], threads=2)
    assert a["preflight"]["subjects"] < b["preflight"]["subjects"]
    with pytest.raises(ValueError):
        blast.build(store.path, tmp_path / "c", taxon=[9606, "mouse"], collection="x")


def test_pathlike_store_and_out(store, tmp_path):
    blast.build(Path(store.path), tmp_path / "p", all=True, threads=2)
    assert (tmp_path / "p.nin").exists() or files(tmp_path / "p")


def test_views_and_view_verb(store, tmp_path):
    out = tmp_path / "db"
    blast.build(store.path, out, all=True, threads=2)
    v = blast.view(out, per_collection=True)
    assert len(v["views"]) == 3
    with pytest.raises(ValueError, match="per-collection"):
        blast.view(out)


def test_nstats_one_row_per_collection(store):
    r = blast.nstats(store.path, all=True, threads=2)
    assert len(r["collections"]) == 3
    assert all(0 <= c["n_fraction"] <= 1 for c in r["collections"])


def test_exception_classes():
    assert issubclass(blast.EngineNotFoundError, FileNotFoundError)
    assert issubclass(blast.EngineError, RuntimeError)


@pytest.mark.skipif(not GTARS_BIN.is_file(), reason="set GTARS_BIN to the gtars binary built from this tree")
def test_build_matches_the_cli_byte_for_byte(store, tmp_path):
    (tmp_path / "cli").mkdir()
    (tmp_path / "py").mkdir()
    subprocess.run(
        [str(GTARS_BIN), "blast", "build", "--store", store.path, "--all", "--out",
         str(tmp_path / "cli" / "db"), "--threads", "2", "-q"],
        check=True,
    )
    blast.build(store.path, tmp_path / "py" / "db", all=True, threads=2)
    a = {p.name: p.read_bytes() for p in (tmp_path / "cli").iterdir()}
    b = {p.name: p.read_bytes() for p in (tmp_path / "py").iterdir()}
    assert a.keys() == b.keys() and a == b


@pytest.mark.skipif(os.name != "posix", reason="needs SIGINT")
def test_sigint_cancels_and_removes_output(store, tmp_path):
    import signal

    out = tmp_path / "db"

    def cb(ev):
        if ev["event"] == "planned":
            os.kill(os.getpid(), signal.SIGINT)

    with pytest.raises(KeyboardInterrupt):
        blast.build(store.path, out, all=True, progress=cb)
    assert files(out) == []


# ---------------------------------------------------------------------------------------------
# search, report, origin, reps

ORIGIN = Path(__file__).resolve().parents[2] / "gtars-rgblast" / "tests" / "data" / "fixtures" / "origin"
COLUMNS = ["query", "status", "taxid", "taxon", "rank", "top_genome", "top_seq", "pident", "qcov",
           "bits", "evalue", "tied_genomes", "tied_taxa", "searched", "flags"]


def find_blastn():
    """How to reach BLAST+: PATH, GTARS_RGBLAST_BIN, or docker image ncbi/blast:2.17.0 (through
    a test-only shim script); None when there is none."""
    if shutil.which("blastn"):
        return shutil.which("blastn")
    d = os.environ.get("GTARS_RGBLAST_BIN")
    if d and (Path(d) / "blastn").is_file():
        return str(Path(d) / "blastn")
    image = os.environ.get("GTARS_RGBLAST_DOCKER", "ncbi/blast:2.17.0")
    if shutil.which("docker") and subprocess.run(
        ["docker", "image", "inspect", image], capture_output=True
    ).returncode == 0:
        return "docker:" + image
    return None


needs_blast = pytest.mark.skipif(find_blastn() is None, reason="no BLAST+ (PATH, GTARS_RGBLAST_BIN or docker)")


@pytest.fixture
def blastn_db(dna, tmp_path):
    """A database built from the dna store, and the blastn keyword arguments to search it."""
    db = tmp_path / "db"
    blast.build(dna.path, db, all=True)
    found = find_blastn()
    if found.startswith("docker:"):
        mounts = {os.path.realpath(p) for p in (tempfile.gettempdir(), tmp_path, os.getcwd())}
        shim = tmp_path / "blastn"
        shim.write_text(
            "#!/bin/sh\nexec docker run --rm -i -u %d:%d %s -w \"$PWD\" %s blastn \"$@\"\n"
            % (os.getuid(), os.getgid(), " ".join(f"-v '{m}:{m}'" for m in sorted(mounts)), found[7:])
        )
        shim.chmod(0o755)
        found = str(shim)
    return dict(engine="blastn", db=db, blastn=found)


def rows_as_tsv(rows):
    return "".join("\t".join(str(r[c]) for c in COLUMNS) + "\n" for r in rows)


def test_report_fields_constant():
    assert blast.REPORT_FIELDS == ["qlen", "score", "sallseqid", "salltitles", "staxids"]


def test_report_matches_the_fixture(tmp_path):
    hits = ORIGIN / "small_jungle_hits.tsv"
    r = blast.report(hits, db=ORIGIN / "small_jungle", query=ORIGIN / "small_jungle_queries.fa")
    want = (ORIGIN / "expected_summary.tsv").read_text().splitlines(keepends=True)[1:]
    assert all(list(row) == COLUMNS for row in r["summary"])
    assert rows_as_tsv(r["summary"]) == "".join(want)
    assert isinstance(r["report"], dict) and r["report"]["queries"]
    out = tmp_path / "r.json"
    blast.report([hits], db=ORIGIN / "small_jungle", output=out, format="json")
    assert out.read_text().lstrip().startswith("{")


def test_report_errors(store):
    hits = ORIGIN / "small_jungle_hits.tsv"
    with pytest.raises(ValueError) as e:
        blast.report(hits)
    assert e.value.kind == "invalid_input"
    with pytest.raises(ValueError):
        blast.report(hits, db=ORIGIN / "small_jungle", store=store.path)
    with pytest.raises(ValueError):
        blast.report(hits, db=ORIGIN / "small_jungle", format="pdf")


def test_reps_on_a_store(store):
    r = blast.reps(store.path, all=True)
    assert set(r) == {"groups", "summary"}
    assert r["groups"]
    with pytest.raises(ValueError):
        blast.reps(store.path)


def test_search_without_an_engine_raises(store, tmp_path):
    q = ORIGIN / "small_jungle_queries.fa"
    with pytest.raises(blast.EngineNotFoundError) as e:
        blast.search(query=q, engine="blastn", db=tmp_path / "db", blastn="/nonexistent/blastn")
    assert isinstance(e.value, FileNotFoundError) and e.value.kind == "not_found"
    with pytest.raises(ValueError):
        blast.search(store.path, all=True, query=q, engine="nope")
    with pytest.raises(ValueError, match="needs --db BASE"):
        blast.search(query=q, engine="blastn")
    with pytest.raises(ValueError, match="select collections when building the database"):
        blast.search(store.path, all=True, query=q, engine="blastn", db=tmp_path / "db")


@pytest.fixture
def dna(tmp_path):
    """A store of two collections of random 2 kb sequences, and a query cut from the first."""
    import random

    rng = random.Random(7)
    seqs = ["".join(rng.choice("ACGT") for _ in range(2000)) for _ in range(3)]
    path = tmp_path / "dnastore"
    s = RefgetStore.on_disk(str(path))
    for i, recs in enumerate([seqs[:2], seqs[1:]]):
        fa = tmp_path / f"c{i}.fa"
        fa.write_text("".join(f">s{i}_{j}\n{x}\n" for j, x in enumerate(recs)))
        s.add_sequence_collection_from_fasta(str(fa))
    del s
    q = tmp_path / "q.fa"
    q.write_text(">q1\n" + seqs[0][500:900] + "\n")
    return SimpleNamespace(path=str(path), query=q)


@needs_blast
def test_blastn_search_captures_output(dna, blastn_db, tmp_path):
    kw = dict(query=dna.query, outfmt="6", **blastn_db)
    events = []
    r1 = blast.search(progress=events.append, **kw)
    assert r1["exit_code"] == 0 and isinstance(r1["output"], str) and r1["output"]
    assert any(e["event"] == "engine_running" or e["event"] == "engine_start" for e in events)
    assert r1["db"] == str(blastn_db["db"]) and r1["store"] is None
    out = tmp_path / "hits.tsv"
    r3 = blast.search(out=out, **kw)
    assert r3["output"] is None and out.read_text() == r1["output"]
    assert (tmp_path / "hits.tsv.meta.json").is_file()
    with pytest.raises(blast.EngineError) as e:
        blast.search(query=dna.query, args=["-task", "nonsense"], **blastn_db)
    assert e.value.returncode != 0


@needs_blast
def test_blastn_origin_returns_a_summary(dna, blastn_db):
    r = blast.origin(query=dna.query, **blastn_db)
    assert set(r) == {"summary", "report", "meta"}
    assert r["summary"] and r["summary"][0]["query"] == "q1"
    assert r["meta"]["mode"] == "exhaustive"
    with pytest.raises(ValueError, match="needs the default engine, rgblastn"):
        blast.origin(query=dna.query, staged=True, **blastn_db)


# ---------------------------------------------------------------------------------------------
# the default engine, rgblastn: in process, no BLAST+ needed

def test_rgblastn_search_and_origin(dna, tmp_path):
    r = blast.search(dna.path, all=True, query=dna.query, outfmt="6")
    assert r["exit_code"] == 0 and r["engine"]["name"] == "rgblastn"
    assert r["output"].startswith("q1\t")
    out = tmp_path / "hits.tsv"
    r2 = blast.search(dna.path, all=True, query=dna.query, outfmt="6", out=out)
    assert r2["output"] is None and out.read_text() == r["output"]
    o = blast.origin(dna.path, all=True, query=dna.query)
    assert o["summary"][0]["query"] == "q1" and o["summary"][0]["status"] == "taxon"


def test_rgblastn_flag_errors_are_value_errors(dna, tmp_path):
    for kw in ({"blastn": "/x/blastn"}, {"db": "x"}):
        with pytest.raises(ValueError) as e:
            blast.search(dna.path, all=True, query=dna.query, **kw)
        assert e.value.kind == "invalid_input"
    with pytest.raises(ValueError):
        blast.search(dna.path, all=True, query=dna.query, engine="rgblastn", args=["-bogus", "1"])

"""BLAST over a RefgetStore: build BLAST databases from a store and search it.

Each function runs the ``gtars blast`` verb of the same name. Keyword options are that verb's
long flags with ``-`` written as ``_`` (``gtars blast build --help`` lists them): ``True``
passes a flag, ``None`` and ``False`` leave it out, a list repeats it, and any other value is
passed as text. Results are plain dicts.

``progress`` is ``None``, ``True`` (print one line per event to stderr) or a callable that gets
one dict per event: ``event`` (for example ``"sequences"``), the event's fields, and
``message``, the line the command line prints. It runs on the calling thread; an exception it
raises cancels the work, removes partial output and is re-raised. Ctrl-C does the same.

Errors carry a ``kind`` attribute: ``"not_found"``, ``"invalid_input"``, ``"conflict"``,
``"io"``, ``"network"``, ``"data"`` or ``"cancelled"``.
"""

import os

from .. import gtars as _ext  # noqa: F401  (registers gtars._blast)
from .._blast import (  # noqa: F401
    REPORT_FIELDS,
    SUMMARY_COLUMNS,
    EngineError,
    EngineNotFoundError,
    _main,
    _run,
)

__all__ = [
    "build",
    "verify",
    "view",
    "nstats",
    "search",
    "report",
    "origin",
    "reps",
    "REPORT_FIELDS",
    "SUMMARY_COLUMNS",
    "EngineError",
    "EngineNotFoundError",
]


def _text(v):
    return os.fspath(v) if isinstance(v, os.PathLike) else str(v)


def _flags(options):
    argv = []
    for key, v in options.items():
        if v is None or v is False:
            continue
        flag = "--" + key.replace("_", "-")
        if v is True:
            argv.append(flag)
        elif isinstance(v, (list, tuple)):
            argv += [f"{flag}={_text(x)}" for x in v]
        else:
            argv.append(f"{flag}={_text(v)}")
    return argv


def _engine_args(query, out, outfmt, args):
    argv = []
    for flag, v in (("-query", query), ("-out", out), ("-outfmt", outfmt)):
        if v is not None:
            argv += [flag, _text(v)]
    return argv + [_text(a) for a in (args or [])]


def _ordered(result):
    """Summary rows with their columns in the TSV's order."""
    result["summary"] = [{c: row[c] for c in SUMMARY_COLUMNS} for row in result["summary"]]
    return result


def _call(verb, lead, options, progress, trailing=()):
    argv = [verb, *lead, *_flags(options)]
    if trailing:
        argv += ["--", *trailing]
    return _run(argv, progress)


def build(store, out, *, progress=None, **options):
    """Write a BLAST version 5 nucleotide database from a store (``gtars blast build``)."""
    return _call("build", ["--store", _text(store), "--out", _text(out)], options, progress)


def verify(db, store, *, progress=None, **options):
    """Check a database against its store (``gtars blast verify``)."""
    return _call("verify", ["--db", _text(db), "--store", _text(store)], options, progress)


def view(db, *, progress=None, **options):
    """Write subset views of a database (``gtars blast view``)."""
    return _call("view", ["--db", _text(db)], options, progress)


def nstats(store, *, progress=None, **options):
    """Count N bases per collection (``gtars blast nstats``)."""
    return _call("nstats", ["--store", _text(store)], options, progress)


def _store(store):
    return [] if store is None else ["--store", _text(store)]


def search(store=None, query=None, *, out=None, outfmt=None, args=None, progress=None, **options):
    """Search the selected collections of a store with rgblastn (``gtars blast search``), or,
    with ``engine="blastn", db=...`` and no store, a database you built with stock blastn.

    ``query``, ``out`` and ``outfmt`` become the engine arguments ``-query``, ``-out`` and
    ``-outfmt``, ahead of ``args`` (more engine arguments, as a list). Without ``out`` the
    engine's output is returned as ``"output"``.
    """
    engine = _engine_args(query, out, outfmt, args)
    return _call("search", _store(store), options, progress, engine)


def report(hits, *, progress=None, **options):
    """Best hit per genome, grouped by species (``gtars blast report``).

    ``hits`` is one tabular hits file or a list of them. Returns ``summary`` (one row per
    query) and ``report``.
    """
    paths = [hits] if isinstance(hits, (str, os.PathLike)) else list(hits)
    return _ordered(_call("report", [_text(p) for p in paths], options, progress))


def origin(store=None, query=None, *, args=None, progress=None, **options):
    """Where did these sequences come from (``gtars blast origin``): a store with rgblastn, or,
    with ``engine="blastn", db=...`` and no store, a database you built (one exhaustive search).

    Returns ``summary``, ``report`` and ``meta``. ``args`` are more engine arguments.
    """
    lead = _store(store)
    if query is not None:
        lead += ["--query", _text(query)]
    return _ordered(_call("origin", lead, options, progress, _engine_args(None, None, None, args)))


def reps(store, *, progress=None, **options):
    """The representative of each species (``gtars blast reps``)."""
    return _call("reps", ["--store", _text(store)], options, progress)

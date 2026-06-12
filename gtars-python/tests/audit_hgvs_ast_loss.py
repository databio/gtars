"""Audit test for finding #8 (MEDIUM): silent AST data loss in PosEditPy::from_rs.

Source under test: gtars-python/src/vrs/hgvs.rs, PosEditPy::from_rs (lines ~275-313).

The HGVS Rust AST (gtars-vrs/src/hgvs/ast.rs) models several "uncertain" and
"whole sequence" location ranges that carry MORE information than the flattened
Python `PosEdit` object can express:

  * LocationRange::WholeSequence  (e.g. `g.=`, `g.?`)
      -> from_rs substitutes a HARDCODED placeholder Position{base:1, offset:0,
         datum:SeqStart} for `start` and sets `end = None`. The fact that the
         expression refers to the whole sequence is silently lost, and the
         caller sees a concrete-looking base=1 position that was never in the
         input. Worse, `uncertain` is False here (parser sets it False for
         whole-sequence), so nothing flags the approximation.

  * LocationRange::UncertainEnd { start, end_low, end_high }
      (e.g. `c.4_(?_246)del`)
      -> from_rs computes `end = end_low.map(...)`. When the LOW bound is `?`
         (end_low == None) but the HIGH bound is present (end_high == Some(246)),
         the end position 246 lives only in `end_high`, which from_rs never
         reads. Result: `pos_edit.end is None` even though the HGVS string
         explicitly bounds the end at 246. The end position is silently dropped.

  * LocationRange::UncertainBoth { start_high, end_low, .. }
      (e.g. `c.(?_6)_(?_246)del`)
      -> same `end = end_low.map(...)` flaw drops the end when end_low is `?`;
         and start falls back to a base=1 placeholder when start_high is `?`.

This test asserts the EXPECTED-CORRECT behavior: the public Python API should
faithfully reflect the position information present in the input (or expose an
explicit approximation flag). A failing assertion CONFIRMS the data-loss bug.

If `gtars.vrs.hgvs` cannot be imported (extension not built with the `vrs`
feature), the whole module is skipped with a clear reason -- verdict NOT-RUN.
"""

import pytest

# Skip the entire module cleanly if the vrs.hgvs binding is not available.
hgvs = pytest.importorskip(
    "gtars.vrs.hgvs",
    reason="gtars.vrs.hgvs extension not built/importable in this environment",
)

# Audit reproduction (PR #256): these assert the correct behavior and currently
# fail, documenting the data-loss bug. Marked xfail (non-strict) so they don't
# gate CI; remove this marker once the bug is fixed and they pass.
pytestmark = pytest.mark.xfail(
    reason="audit reproduction (PR #256): HGVS AST drops uncertain/whole-sequence position info",
    strict=False,
)

parse_hgvs = hgvs.parse_hgvs


def test_uncertain_end_low_unknown_keeps_end_bound():
    """`c.4_(?_246)del` bounds the end at 246 (high bound); it must not vanish.

    The parser stores 246 in `end_high` and leaves `end_low = None` (the `?`).
    PosEditPy::from_rs does `end = end_low.map(...)`, so it silently returns
    `end is None`, dropping the 246 end bound. A faithful binding would expose
    the end position (246) and/or flag it as approximate.
    """
    v = parse_hgvs("NM_000546.6:c.4_(?_246)del")
    pe = v.pos_edit

    # Sanity: this expression IS a range with an uncertain component.
    assert pe.start.base == 4, "start should be the certain position 4"

    # EXPECTED-CORRECT: the end bound 246 from the input must be represented.
    assert pe.end is not None, (
        "end position silently dropped: HGVS specifies an end bound (246) via "
        "the high estimate, but PosEditPy::from_rs reads only end_low and "
        "returns end=None (data-loss bug, hgvs.rs ~line 295)"
    )
    assert pe.end.base == 246, (
        f"end bound should reflect the input (246), got base={pe.end.base}"
    )


def test_whole_sequence_not_misrepresented_as_base_1():
    """`g.=` is a whole-sequence identity; it must not masquerade as base=1.

    PosEditPy::from_rs substitutes Position{base:1, offset:0} and end=None for
    WholeSequence, fabricating a concrete-looking position that was never in the
    input, with no flag. A faithful binding would either expose a whole-sequence
    marker or at least set `uncertain`/some indicator -- not silently report a
    real-looking base=1.
    """
    v = parse_hgvs("NC_000017.11:g.=")
    pe = v.pos_edit

    # EXPECTED-CORRECT: a caller must be able to tell this is NOT a concrete
    # base-1 single position. Either an explicit approximation/whole-sequence
    # flag is set, or the fabricated base=1 placeholder is not presented as a
    # genuine position. The current code reports base=1, offset=0, end=None,
    # uncertain=False -- indistinguishable from a real `g.1` identity.
    fabricated_placeholder = (
        pe.start.base == 1 and pe.start.offset == 0 and pe.end is None
    )
    assert not (fabricated_placeholder and not pe.uncertain), (
        "whole-sequence `g.=` is silently flattened to a fabricated concrete "
        "position base=1/offset=0/end=None with uncertain=False, "
        "indistinguishable from a real single-base identity (hgvs.rs ~line 285)"
    )


def test_uncertain_both_low_unknown_keeps_end_bound():
    """`c.(?_6)_(?_246)del`: end bound 246 lives in end_high and is dropped."""
    v = parse_hgvs("NM_000546.6:c.(?_6)_(?_246)del")
    pe = v.pos_edit

    assert pe.end is not None, (
        "UncertainBoth end bound (246) dropped: from_rs uses end_low.map(...) "
        "and ignores end_high (hgvs.rs ~line 302)"
    )
    assert pe.end.base == 246

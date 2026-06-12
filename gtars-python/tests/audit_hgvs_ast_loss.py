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

# Audit reproduction (PR #256): the bug is now fixed. PosEditPy faithfully
# projects every LocationRange variant via the structured PositionBound API, so
# these tests now assert the correct behavior and must PASS (no xfail).

parse_hgvs = hgvs.parse_hgvs


def test_uncertain_end_low_unknown_keeps_end_bound():
    """`c.4_(?_246)del` bounds the end at 246 (high bound); it must not vanish.

    The parser stores 246 in `end_high` and leaves `end_low = None` (the `?`).
    The faithful binding exposes the end as an uncertain PositionBound whose
    `high` carries 246 and whose `low` is None (the `?`).
    """
    v = parse_hgvs("NM_000546.6:c.4_(?_246)del")
    pe = v.pos_edit

    # This expression IS a range with an uncertain end component.
    assert pe.location_kind == "range"
    assert pe.start.kind == "certain"
    assert pe.start.position.base == 4, "start should be the certain position 4"

    # The end bound 246 from the input must be represented (not dropped).
    assert pe.end is not None, "end endpoint must be present for this range"
    assert pe.end.kind == "uncertain"
    assert pe.end.low is None, "the `?` low bound is represented as None"
    assert pe.end.high is not None, "the high bound (246) must be preserved"
    assert pe.end.high.base == 246, (
        f"end high bound should reflect the input (246), got base={pe.end.high.base}"
    )


def test_whole_sequence_not_misrepresented_as_base_1():
    """`g.=` is a whole-sequence identity; it must not masquerade as base=1.

    The faithful binding marks `location_kind == "whole_sequence"` and exposes
    no fabricated positions: both `start` and `end` are None.
    """
    v = parse_hgvs("NC_000017.11:g.=")
    pe = v.pos_edit

    assert pe.location_kind == "whole_sequence"
    assert pe.start is None, "no fabricated base=1 start position should exist"
    assert pe.end is None, "whole-sequence has no end endpoint"


def test_uncertain_both_low_unknown_keeps_end_bound():
    """`c.(?_6)_(?_246)del`: both endpoints uncertain; bounds preserved."""
    v = parse_hgvs("NM_000546.6:c.(?_6)_(?_246)del")
    pe = v.pos_edit

    assert pe.location_kind == "range"

    assert pe.start is not None
    assert pe.start.kind == "uncertain"
    assert pe.start.low is None, "the `?` start low bound is represented as None"
    assert pe.start.high is not None
    assert pe.start.high.base == 6

    assert pe.end is not None, "UncertainBoth end bound (246) must not be dropped"
    assert pe.end.kind == "uncertain"
    assert pe.end.low is None, "the `?` end low bound is represented as None"
    assert pe.end.high is not None
    assert pe.end.high.base == 246

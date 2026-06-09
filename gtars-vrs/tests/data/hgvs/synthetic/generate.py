#!/usr/bin/env python3
"""
Synthetic HGVS test-set generator (gtars plan 10).

============================================================================
THIS GENERATOR MUST NOT IMPORT, CALL, OR RE-IMPLEMENT BY COPYING ANY CODE
FROM gtars_reftx::CoordinateMapper OR gtars_vrs::hgvs::bridge.

It walks exon tables directly to compute expected genomic coordinates.
Sharing coordinate-transformation code with the Rust crates under test
turns layer-2 of the synthetic test set into a tautology and silently
masks real bugs. If you find yourself wanting to import a Rust helper
to "avoid duplication" -- STOP. Duplication is the entire point.

Layer-1 VRS-ID generation MAY use the gtars PyO3 module because layer 1
is documented as a self-consistency check, not an independent oracle.
============================================================================

Outputs (in --out-dir, default = this script's dir):
  synthetic.fa
  synthetic.fa.fai
  synthetic_transcripts.json   (cdot 0.2.x-shape)
  refget_collection.json       (collection digest + per-sequence digests)
  cases.tsv                    (the test corpus)

Run:
    python3 generate.py            # regenerate fixtures
    python3 generate.py --check    # regenerate to tmpdir, diff against checked-in
"""
from __future__ import annotations

import argparse
import base64
import filecmp
import hashlib
import json
import os
import random
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

# ---------------------------------------------------------------------------
# Module constants -- bumping FIXTURE_VERSION marks a deliberate corpus change.
# ---------------------------------------------------------------------------

FIXTURE_VERSION = 1
GENOME_SEED = 0xDA51
GENOME_LENGTH = 30_000
CHROM_NAME = "chr_synth"


# ---------------------------------------------------------------------------
# Transcript layout (1-based inclusive, cdot convention).
# ---------------------------------------------------------------------------

@dataclass
class Tx:
    accession: str
    strand: int  # +1 or -1
    exons: list[tuple[int, int]]  # (g_start, g_end), genomic-sorted ascending
    cds: Optional[tuple[int, int]]  # (g_start, g_end) of CDS (inclusive)

TRANSCRIPTS: list[Tx] = [
    Tx("NM_synth_simple", +1,
       [(1001, 1300)],
       (1051, 1249)),
    Tx("NM_synth_fwd", +1,
       [(2001, 2100), (2301, 2450), (2701, 2850), (3101, 3200)],
       (2051, 3150)),
    Tx("NM_synth_rev", -1,
       [(5001, 5100), (5301, 5450), (5701, 5850), (6101, 6200)],
       (5051, 6150)),
    Tx("NR_synth_nc", +1,
       [(8001, 8050), (8201, 8300), (8501, 8550)],
       None),
    Tx("NM_synth_tiny", +1,
       [(9001, 9010), (9015, 9030), (9035, 9044)],
       (9003, 9042)),
]
TX_BY_NAME = {t.accession: t for t in TRANSCRIPTS}


# ---------------------------------------------------------------------------
# Synthetic genome
# ---------------------------------------------------------------------------

def build_genome() -> bytes:
    rng = random.Random(GENOME_SEED)
    return "".join(rng.choices("ACGT", k=GENOME_LENGTH)).encode("ascii")


def write_fasta(path: Path, name: str, seq: bytes) -> None:
    with path.open("wb") as fh:
        fh.write(b">" + name.encode("ascii") + b"\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60])
            fh.write(b"\n")


def write_fai(fai_path: Path, name: str, length: int) -> None:
    # samtools .fai columns: name, length, offset (byte after header line),
    # linebases (60), linewidth (61, includes \n).
    header_offset = len(b">" + name.encode("ascii") + b"\n")
    fai_path.write_text(f"{name}\t{length}\t{header_offset}\t60\t61\n")


# ---------------------------------------------------------------------------
# refget digest helpers (independent of gtars; verified to match gtars-refget)
# ---------------------------------------------------------------------------

def sha512t24u(data: bytes) -> str:
    """GA4GH sha512t24u: base64url(no padding) of first 24 bytes of SHA-512."""
    digest = hashlib.sha512(data).digest()[:24]
    return base64.urlsafe_b64encode(digest).rstrip(b"=").decode("ascii")


def md5_hex(data: bytes) -> str:
    return hashlib.md5(data).hexdigest()


# ---------------------------------------------------------------------------
# Brute-force coordinate walker (Layer-2 oracle).
# Independent of gtars_reftx::CoordinateMapper.
# ---------------------------------------------------------------------------

def exons_in_tx_order(tx: Tx) -> list[tuple[int, int]]:
    """Return exons in 5'->3' transcript order."""
    if tx.strand == +1:
        return list(tx.exons)
    return list(reversed(tx.exons))


def transcript_pos_to_genomic(tx: Tx, tx_pos: int) -> Optional[int]:
    """1-based transcript position -> 1-based genomic position. None if OOB."""
    if tx_pos < 1:
        return None
    acc = 0
    for (g_start, g_end) in exons_in_tx_order(tx):
        L = g_end - g_start + 1
        if acc + 1 <= tx_pos <= acc + L:
            within = tx_pos - acc  # 1-based offset into exon (transcript order)
            if tx.strand == +1:
                return g_start + within - 1
            else:
                return g_end - within + 1
        acc += L
    return None


def transcript_length(tx: Tx) -> int:
    return sum(e - s + 1 for (s, e) in tx.exons)


def cds_tx_bounds(tx: Tx) -> Optional[tuple[int, int]]:
    """Return (cds_tx_start, cds_tx_end) in 1-based transcript coords, or None."""
    if tx.cds is None:
        return None
    g_lo, g_hi = tx.cds
    if tx.strand == +1:
        # 5' end of CDS in transcript = genomic g_lo
        cds_5p_g = g_lo
        cds_3p_g = g_hi
    else:
        cds_5p_g = g_hi
        cds_3p_g = g_lo
    s = genomic_to_tx_pos_in_exon(tx, cds_5p_g)
    e = genomic_to_tx_pos_in_exon(tx, cds_3p_g)
    assert s is not None and e is not None, f"CDS endpoints must lie in exons of {tx.accession}"
    return (s, e)


def genomic_to_tx_pos_in_exon(tx: Tx, g_pos: int) -> Optional[int]:
    """If g_pos is inside an exon, return its 1-based transcript position."""
    acc = 0
    for (g_start, g_end) in exons_in_tx_order(tx):
        L = g_end - g_start + 1
        if g_start <= g_pos <= g_end:
            if tx.strand == +1:
                within = g_pos - g_start + 1
            else:
                within = g_end - g_pos + 1
            return acc + within
        acc += L
    return None


def intron_length_after_exon(tx: Tx, exon_idx_tx_order: int) -> Optional[int]:
    """Length (bp) of the intron 3' of the given exon, in transcript order.

    Returns None if there is no following exon.
    """
    exons_to = exons_in_tx_order(tx)
    if exon_idx_tx_order + 1 >= len(exons_to):
        return None
    cur = exons_to[exon_idx_tx_order]
    nxt = exons_to[exon_idx_tx_order + 1]
    if tx.strand == +1:
        # intron between cur.end+1 .. nxt.start-1 inclusive
        return nxt[0] - cur[1] - 1
    else:
        # transcript order reversed; in genomic order the next exon is to the LEFT.
        # cur is the rightmost of the pair, nxt is the leftmost.
        # intron = nxt.end+1 .. cur.start-1 (genomic)
        return cur[0] - nxt[1] - 1


@dataclass
class CResult:
    g_pos: Optional[int] = None  # 1-based genomic
    error: Optional[str] = None  # error tag


def c_pos_to_genomic(tx: Tx, c_pos: int, intron_offset: int,
                     is_cds_end: bool = False) -> CResult:
    """Convert HGVS c. coordinate to 1-based genomic position.

    c_pos: signed integer (negative = 5'UTR, positive = CDS)
    intron_offset: 0 = exonic; +N = N bases into intron 3' of c_pos;
                   -N = N bases into intron 5' of c_pos (i.e. 5' of c_pos+1).
    is_cds_end: if True, c_pos is interpreted as 3'UTR (c.*N) where N>=1.
    """
    if tx.cds is None:
        return CResult(error="NonCoding")
    if c_pos == 0:
        return CResult(error="ParseError")  # not valid HGVS
    bounds = cds_tx_bounds(tx)
    assert bounds is not None
    cds_tx_start, cds_tx_end = bounds

    if is_cds_end:
        # c.*N: N bases past cds_tx_end (1-based offset).
        if c_pos < 1:
            return CResult(error="ParseError")
        tx_pos = cds_tx_end + c_pos
    elif c_pos < 0:
        # c.-N: N bases 5' of CDS start.
        tx_pos = cds_tx_start + c_pos  # c_pos is negative
        if tx_pos < 1:
            return CResult(error="OutOfBounds")
    else:
        # c.N (positive): N-th coding base.
        tx_pos = cds_tx_start + c_pos - 1

    return _tx_pos_with_offset_to_genomic(tx, tx_pos, intron_offset)


def n_pos_to_genomic(tx: Tx, n_pos: int, intron_offset: int) -> CResult:
    """HGVS n.N (+offset / -offset) for non-coding (or coding) transcripts."""
    if n_pos < 1:
        return CResult(error="ParseError")
    return _tx_pos_with_offset_to_genomic(tx, n_pos, intron_offset)


def _tx_pos_with_offset_to_genomic(tx: Tx, tx_pos: int, intron_offset: int) -> CResult:
    """Common tail: an integer transcript position plus an intronic offset."""
    tx_len = transcript_length(tx)
    # The base anchor must lie inside an exon for offset semantics.
    if tx_pos < 1 or tx_pos > tx_len:
        # Only allowed if anchor is the conceptual position just before exon 1
        # (HGVS does not really place anchors outside the transcript for n.,
        #  so we treat these as out-of-bounds).
        return CResult(error="OutOfBounds")

    if intron_offset == 0:
        g = transcript_pos_to_genomic(tx, tx_pos)
        if g is None:
            return CResult(error="OutOfBounds")
        return CResult(g_pos=g)

    # Intronic. Find which exon the anchor sits in (transcript order) and
    # determine which intron the offset shifts into.
    exons_to = exons_in_tx_order(tx)
    acc = 0
    anchor_exon_idx = None
    anchor_within = None  # 1-based offset into exon (transcript order)
    for i, (g_start, g_end) in enumerate(exons_to):
        L = g_end - g_start + 1
        if acc + 1 <= tx_pos <= acc + L:
            anchor_exon_idx = i
            anchor_within = tx_pos - acc
            break
        acc += L
    assert anchor_exon_idx is not None

    if intron_offset > 0:
        # +offset: anchor must be at the 3' end of its exon (last base in tx order).
        cur = exons_to[anchor_exon_idx]
        L = cur[1] - cur[0] + 1
        if anchor_within != L:
            return CResult(error="InvalidIntronicOffset")
        intron_len = intron_length_after_exon(tx, anchor_exon_idx)
        if intron_len is None:
            return CResult(error="OutOfBounds")
        if intron_offset > intron_len:
            return CResult(error="IntronOffsetTooLarge")
        # genomic shift: +offset in tx orientation
        if tx.strand == +1:
            g_anchor = cur[1]  # genomic last base of exon (3' in tx)
            return CResult(g_pos=g_anchor + intron_offset)
        else:
            g_anchor = cur[0]  # in tx order rightmost->leftmost; on minus strand 3' in tx is leftmost genomic
            return CResult(g_pos=g_anchor - intron_offset)
    else:
        # -offset: anchor must be at the 5' end of its exon (first base in tx order).
        if anchor_within != 1:
            return CResult(error="InvalidIntronicOffset")
        if anchor_exon_idx == 0:
            return CResult(error="OutOfBounds")  # no preceding intron
        intron_len = intron_length_after_exon(tx, anchor_exon_idx - 1)
        assert intron_len is not None
        if -intron_offset > intron_len:
            return CResult(error="IntronOffsetTooLarge")
        cur = exons_to[anchor_exon_idx]
        if tx.strand == +1:
            g_anchor = cur[0]  # 5' in tx = leftmost genomic
            return CResult(g_pos=g_anchor + intron_offset)  # intron_offset negative
        else:
            g_anchor = cur[1]  # 5' in tx = rightmost genomic on minus
            return CResult(g_pos=g_anchor - intron_offset)  # subtract negative = add


# ---------------------------------------------------------------------------
# Sequence helpers
# ---------------------------------------------------------------------------

COMP = {ord('A'): ord('T'), ord('T'): ord('A'),
        ord('C'): ord('G'), ord('G'): ord('C'),
        ord('N'): ord('N')}


def revcomp(seq: bytes) -> bytes:
    return bytes(COMP[b] for b in reversed(seq))


def get_genomic_ref(genome: bytes, g_pos_1based: int, length: int) -> bytes:
    """Return forward-strand genomic ref bytes [g_pos_1based, g_pos_1based+length)."""
    start_0 = g_pos_1based - 1
    return genome[start_0:start_0 + length]


# ---------------------------------------------------------------------------
# cdot JSON emission
# ---------------------------------------------------------------------------

def emit_cdot(path: Path) -> None:
    """Write a minimal cdot 0.2.x-shaped JSON document.

    Our internal `Tx` table stores exons and CDS as 1-based inclusive
    (s, e) tuples. cdot / gtars-reftx use 0-based half-open intervals,
    so we convert here: 1-based inclusive [s, e] -> 0-based half-open
    [s-1, e].
    """
    transcripts = {}
    for tx in TRANSCRIPTS:
        transcripts[tx.accession] = {
            "id": tx.accession,
            "gene_name": tx.accession.replace("NM_", "GENE_").replace("NR_", "GENE_"),
            "contig": CHROM_NAME,
            "strand": tx.strand,  # cdot uses +1/-1 ints (per gtars-reftx loader)
            "cds_start": (tx.cds[0] - 1) if tx.cds else None,
            "cds_end": tx.cds[1] if tx.cds else None,
            "exons": [[s - 1, e] for (s, e) in tx.exons],
        }
    doc = {"transcripts": transcripts}
    path.write_text(json.dumps(doc, indent=2, sort_keys=True))


# ---------------------------------------------------------------------------
# refget collection JSON
# ---------------------------------------------------------------------------

def emit_refget_collection(path: Path, seq: bytes) -> str:
    """Write refget collection JSON; return the collection digest."""
    seq_digest = sha512t24u(seq)
    md5 = md5_hex(seq)

    # Compute the collection digest (a sha512t24u of the canonical
    # name-and-length-and-sequence triplet, GA4GH-spec). For our purposes
    # what matters is that the test loads the FASTA into gtars-refget and
    # uses *its* collection digest. We still write a useful manifest.
    collection_payload = json.dumps({
        "names": [CHROM_NAME],
        "lengths": [len(seq)],
        "sequences": [seq_digest],
    }, sort_keys=True, separators=(",", ":")).encode("ascii")
    collection_digest = "SQ." + sha512t24u(collection_payload)

    manifest = {
        "fixture_version": FIXTURE_VERSION,
        "collection_digest_advisory": collection_digest,
        "sequences": {
            CHROM_NAME: {
                "ga4gh": "SQ." + seq_digest,
                "sha512t24u": seq_digest,
                "length": len(seq),
                "md5": md5,
            }
        },
    }
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    return collection_digest


# ---------------------------------------------------------------------------
# Optional: import gtars for layer-1 expected VRS IDs
# ---------------------------------------------------------------------------

def _try_import_gtars():
    try:
        from gtars.vrs import vrs_id  # noqa: F401
        from gtars.refget import RefgetStore  # noqa: F401
        return True
    except Exception:
        return False


def compute_vrs_ids_for_cases(cases: list["Case"], fasta_path: Path) -> dict[str, str]:
    """For each non-error case, compute the VRS ID using gtars.

    Returns a dict {case_id: vrs_id_or_empty}. If gtars cannot be imported,
    raises -- per the plan, no silent skipping.
    """
    try:
        from gtars.refget import RefgetStore
        from gtars.vrs import vrs_id, normalize_allele
    except Exception as e:
        raise RuntimeError(
            "gtars PyO3 module not importable; cannot compute layer-1 expected "
            "VRS IDs. Install with `pip install -e gtars-python/` or run via "
            "PYTHONPATH=gtars-python/py_src.\n"
            f"Underlying error: {e}"
        )

    store = RefgetStore.in_memory()
    store.disable_encoding()
    store.set_quiet(True)
    store.add_sequence_collection_from_fasta(str(fasta_path))
    # Find the seq digest for our chromosome.
    # Use digest_sequence on the in-memory FASTA contents to be safe.
    from gtars.refget import digest_sequence
    seq_bytes = fasta_path.read_bytes()
    # Reconstruct the raw sequence (strip headers / newlines).
    raw = bytearray()
    for line in seq_bytes.splitlines():
        if line.startswith(b">"):
            continue
        raw.extend(line)
    rec = digest_sequence(bytes(raw), name=CHROM_NAME)
    sha = rec.metadata.sha512t24u
    sq = "SQ." + sha

    result = {}
    for case in cases:
        if case.expected_error or not case.expected_chrom:
            result[case.case_id] = ""
            continue
        # Normalize then digest. Use the same path the bridge uses.
        ref_str = case.expected_ref or ""
        alt_str = case.expected_alt or ""
        # 0-based interbase start
        start_ib = case.expected_pos - 1  # case.expected_pos is 1-based
        # For a sub/identity, the VCF-equivalent reference span is len(ref_str).
        # For an ins (ref empty), span is 0; for del (alt empty), span is len(ref).
        # We let normalize handle the rest.
        norm = normalize_allele(bytes(raw).decode("ascii"), start_ib, ref_str, alt_str)
        result[case.case_id] = vrs_id(sq, norm["start"], norm["end"], norm["allele"])
    return result


# ---------------------------------------------------------------------------
# Case generation
# ---------------------------------------------------------------------------

@dataclass
class Case:
    case_id: str = ""
    hgvs_string: str = ""
    ref_type: str = ""
    edit_type: str = ""
    location_class: str = ""
    strand: str = ""
    transcript: str = ""
    expected_chrom: str = ""
    expected_pos: int = 0  # 1-based; 0 if unset
    expected_ref: str = ""
    expected_alt: str = ""
    expected_vrs_id: str = ""
    expected_error: str = ""
    notes: str = ""


def _strand_str(s: int) -> str:
    return "+" if s == +1 else "-"


def _hgvs_c_pos(c_pos: int, offset: int, is_cds_end: bool) -> str:
    """Render the position part of an HGVS c. position (no edit suffix)."""
    if is_cds_end:
        base = f"*{c_pos}"
    elif c_pos < 0:
        base = f"{c_pos}"
    else:
        base = f"{c_pos}"
    if offset > 0:
        return f"{base}+{offset}"
    if offset < 0:
        return f"{base}{offset}"  # negative sign present
    return base


def _hgvs_n_pos(n_pos: int, offset: int) -> str:
    base = f"{n_pos}"
    if offset > 0:
        return f"{base}+{offset}"
    if offset < 0:
        return f"{base}{offset}"
    return base


def _gen_substitution_cases(genome: bytes, cases: list[Case]) -> None:
    """Generate substitution (SNV) cases across reference types and locations."""

    # ---- g. SNV on chr_synth at a known position (forward) ----
    for g_pos in (1500, 4000, 7000):
        ref = chr(genome[g_pos - 1])
        alt = "T" if ref != "T" else "A"
        cases.append(Case(
            hgvs_string=f"{CHROM_NAME}:g.{g_pos}{ref}>{alt}",
            ref_type="g", edit_type="sub", location_class="genomic", strand=".",
            transcript="",
            expected_chrom=CHROM_NAME, expected_pos=g_pos,
            expected_ref=ref, expected_alt=alt,
            notes=f"plain g. SNV at g.{g_pos}",
        ))

    # ---- c. SNV in CDS for each coding tx ----
    coding = [t for t in TRANSCRIPTS if t.cds is not None]
    for tx in coding:
        bounds = cds_tx_bounds(tx)
        assert bounds
        cds_tx_start, cds_tx_end = bounds
        # exercise c.1 (start codon), middle, last
        for label, c_pos in [("cds_start", 1),
                             ("cds", 50 if cds_tx_end - cds_tx_start >= 50 else 2),
                             ("cds_stop", cds_tx_end - cds_tx_start + 1)]:
            r = c_pos_to_genomic(tx, c_pos, 0)
            if r.error or r.g_pos is None:
                continue
            g_ref_byte = chr(genome[r.g_pos - 1])
            if tx.strand == +1:
                hgvs_ref = g_ref_byte
            else:
                hgvs_ref = COMP_C(g_ref_byte)
            hgvs_alt = "T" if hgvs_ref != "T" else "A"
            g_alt = hgvs_alt if tx.strand == +1 else COMP_C(hgvs_alt)
            cases.append(Case(
                hgvs_string=f"{tx.accession}:c.{c_pos}{hgvs_ref}>{hgvs_alt}",
                ref_type="c", edit_type="sub", location_class=label,
                strand=_strand_str(tx.strand), transcript=tx.accession,
                expected_chrom=CHROM_NAME, expected_pos=r.g_pos,
                expected_ref=g_ref_byte, expected_alt=g_alt,
                notes=f"c. SNV at {label}",
            ))

    # ---- 5'UTR (c.-N) ----
    for tx in coding:
        bounds = cds_tx_bounds(tx)
        cds_tx_start = bounds[0]
        if cds_tx_start < 2:
            continue  # no 5'UTR available
        c_pos = -1  # one base before CDS start
        r = c_pos_to_genomic(tx, c_pos, 0)
        if r.error or r.g_pos is None:
            continue
        g_ref_byte = chr(genome[r.g_pos - 1])
        hgvs_ref = g_ref_byte if tx.strand == +1 else COMP_C(g_ref_byte)
        hgvs_alt = "T" if hgvs_ref != "T" else "A"
        g_alt = hgvs_alt if tx.strand == +1 else COMP_C(hgvs_alt)
        cases.append(Case(
            hgvs_string=f"{tx.accession}:c.{c_pos}{hgvs_ref}>{hgvs_alt}",
            ref_type="c", edit_type="sub", location_class="5utr",
            strand=_strand_str(tx.strand), transcript=tx.accession,
            expected_chrom=CHROM_NAME, expected_pos=r.g_pos,
            expected_ref=g_ref_byte, expected_alt=g_alt,
            notes="c.-1 5'UTR SNV",
        ))

    # ---- 3'UTR (c.*N) ----
    for tx in coding:
        # Place at *1 (one base after stop codon).
        r = c_pos_to_genomic(tx, 1, 0, is_cds_end=True)
        if r.error or r.g_pos is None:
            continue
        g_ref_byte = chr(genome[r.g_pos - 1])
        hgvs_ref = g_ref_byte if tx.strand == +1 else COMP_C(g_ref_byte)
        hgvs_alt = "T" if hgvs_ref != "T" else "A"
        g_alt = hgvs_alt if tx.strand == +1 else COMP_C(hgvs_alt)
        cases.append(Case(
            hgvs_string=f"{tx.accession}:c.*1{hgvs_ref}>{hgvs_alt}",
            ref_type="c", edit_type="sub", location_class="3utr",
            strand=_strand_str(tx.strand), transcript=tx.accession,
            expected_chrom=CHROM_NAME, expected_pos=r.g_pos,
            expected_ref=g_ref_byte, expected_alt=g_alt,
            notes="c.*1 3'UTR SNV",
        ))

    # ---- intronic positive offset (c.X+5) ----
    for tx in coding:
        if len(tx.exons) < 2:
            continue
        # pick exon 0 (in tx order). last base of that exon in tx coords.
        exons_to = exons_in_tx_order(tx)
        L0 = exons_to[0][1] - exons_to[0][0] + 1
        # convert tx position L0 -> c_pos (relative to CDS start)
        bounds = cds_tx_bounds(tx); cds_tx_start, _ = bounds
        c_pos = L0 - cds_tx_start + 1  # may be 0 or negative; if so, skip
        if c_pos == 0:
            continue
        intron_len = intron_length_after_exon(tx, 0)
        for off in (1, intron_len):
            if off is None or off < 1:
                continue
            r = c_pos_to_genomic(tx, c_pos, off)
            if r.error or r.g_pos is None:
                continue
            g_ref_byte = chr(genome[r.g_pos - 1])
            hgvs_ref = g_ref_byte if tx.strand == +1 else COMP_C(g_ref_byte)
            hgvs_alt = "T" if hgvs_ref != "T" else "A"
            g_alt = hgvs_alt if tx.strand == +1 else COMP_C(hgvs_alt)
            cases.append(Case(
                hgvs_string=f"{tx.accession}:c.{_hgvs_c_pos(c_pos, off, False)}{hgvs_ref}>{hgvs_alt}",
                ref_type="c", edit_type="sub", location_class="intron_pos",
                strand=_strand_str(tx.strand), transcript=tx.accession,
                expected_chrom=CHROM_NAME, expected_pos=r.g_pos,
                expected_ref=g_ref_byte, expected_alt=g_alt,
                notes=f"intron+{off} (intron len {intron_len})",
            ))

    # ---- intronic negative offset (c.Y-3) ----
    for tx in coding:
        if len(tx.exons) < 2:
            continue
        exons_to = exons_in_tx_order(tx)
        # First base of exon 1 (the exon AFTER the first intron) in tx coords:
        L0 = exons_to[0][1] - exons_to[0][0] + 1
        tx_pos_anchor = L0 + 1  # first base of 2nd exon in tx order
        bounds = cds_tx_bounds(tx); cds_tx_start, _ = bounds
        c_pos = tx_pos_anchor - cds_tx_start + 1
        if c_pos == 0:
            continue
        intron_len = intron_length_after_exon(tx, 0)
        for off in (1, intron_len):
            if off is None or off < 1:
                continue
            r = c_pos_to_genomic(tx, c_pos, -off)
            if r.error or r.g_pos is None:
                continue
            g_ref_byte = chr(genome[r.g_pos - 1])
            hgvs_ref = g_ref_byte if tx.strand == +1 else COMP_C(g_ref_byte)
            hgvs_alt = "T" if hgvs_ref != "T" else "A"
            g_alt = hgvs_alt if tx.strand == +1 else COMP_C(hgvs_alt)
            cases.append(Case(
                hgvs_string=f"{tx.accession}:c.{_hgvs_c_pos(c_pos, -off, False)}{hgvs_ref}>{hgvs_alt}",
                ref_type="c", edit_type="sub", location_class="intron_neg",
                strand=_strand_str(tx.strand), transcript=tx.accession,
                expected_chrom=CHROM_NAME, expected_pos=r.g_pos,
                expected_ref=g_ref_byte, expected_alt=g_alt,
                notes=f"intron-{off} (intron len {intron_len})",
            ))

    # ---- exon boundary first/last base ----
    for tx in coding:
        # First base of CDS-containing transcript in tx coords (the 5'-most exonic base).
        # That maps to the genomic edge of the first exon (in tx order).
        bounds = cds_tx_bounds(tx); cds_tx_start, cds_tx_end = bounds
        # tx_pos = 1 (5' UTR boundary if applicable)
        for label, tx_pos in [("exon_boundary_first", 1),
                              ("exon_boundary_last", transcript_length(tx))]:
            c_pos = tx_pos - cds_tx_start + 1
            if c_pos == 0:
                # land at c.1 (covered elsewhere); shift slightly
                c_pos = 1
            # Check c_pos is non-CDS-end-style (c.-N or c.N or c.*N)?
            # If tx_pos > cds_tx_end, this should be c.*N
            is_end = tx_pos > cds_tx_end
            if is_end:
                c_pos_render = tx_pos - cds_tx_end
                r = c_pos_to_genomic(tx, c_pos_render, 0, is_cds_end=True)
                pos_str = f"*{c_pos_render}"
            elif tx_pos < cds_tx_start:
                c_pos_render = tx_pos - cds_tx_start
                r = c_pos_to_genomic(tx, c_pos_render, 0)
                pos_str = f"{c_pos_render}"
            else:
                c_pos_render = tx_pos - cds_tx_start + 1
                r = c_pos_to_genomic(tx, c_pos_render, 0)
                pos_str = f"{c_pos_render}"
            if r.error or r.g_pos is None:
                continue
            g_ref_byte = chr(genome[r.g_pos - 1])
            hgvs_ref = g_ref_byte if tx.strand == +1 else COMP_C(g_ref_byte)
            hgvs_alt = "T" if hgvs_ref != "T" else "A"
            g_alt = hgvs_alt if tx.strand == +1 else COMP_C(hgvs_alt)
            cases.append(Case(
                hgvs_string=f"{tx.accession}:c.{pos_str}{hgvs_ref}>{hgvs_alt}",
                ref_type="c", edit_type="sub", location_class=label,
                strand=_strand_str(tx.strand), transcript=tx.accession,
                expected_chrom=CHROM_NAME, expected_pos=r.g_pos,
                expected_ref=g_ref_byte, expected_alt=g_alt,
                notes=f"{label} via tx_pos={tx_pos}",
            ))

    # ---- n. SNV on the noncoding transcript ----
    nc = TX_BY_NAME["NR_synth_nc"]
    for n_pos in (1, 50, transcript_length(nc)):
        r = n_pos_to_genomic(nc, n_pos, 0)
        if r.error or r.g_pos is None:
            continue
        g_ref_byte = chr(genome[r.g_pos - 1])
        hgvs_ref = g_ref_byte
        hgvs_alt = "T" if hgvs_ref != "T" else "A"
        cases.append(Case(
            hgvs_string=f"{nc.accession}:n.{n_pos}{hgvs_ref}>{hgvs_alt}",
            ref_type="n", edit_type="sub", location_class="cds",
            strand="+", transcript=nc.accession,
            expected_chrom=CHROM_NAME, expected_pos=r.g_pos,
            expected_ref=g_ref_byte, expected_alt=hgvs_alt,
            notes=f"n.{n_pos} SNV on noncoding tx",
        ))


def COMP_C(b: str) -> str:
    return {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}[b]


def _gen_del_cases(genome: bytes, cases: list[Case]) -> None:
    """Generate del cases for various lengths."""
    coding = [t for t in TRANSCRIPTS if t.cds is not None]
    # g. del at fixed positions, lengths 1, 2, 3, 5
    for length in (1, 2, 3, 5):
        g_pos = 1600 + length * 10
        ref = genome[g_pos - 1: g_pos - 1 + length].decode("ascii")
        if length == 1:
            hgvs = f"{CHROM_NAME}:g.{g_pos}del"
        else:
            hgvs = f"{CHROM_NAME}:g.{g_pos}_{g_pos + length - 1}del"
        cases.append(Case(
            hgvs_string=hgvs,
            ref_type="g", edit_type="del", location_class="genomic", strand=".",
            transcript="",
            expected_chrom=CHROM_NAME, expected_pos=g_pos,
            expected_ref=ref, expected_alt="",
            notes=f"g. del len {length}",
        ))
    # c. del on NM_synth_fwd within CDS
    tx = TX_BY_NAME["NM_synth_fwd"]
    for length in (1, 2, 3):
        c_start = 10
        c_end = c_start + length - 1
        r0 = c_pos_to_genomic(tx, c_start, 0)
        r1 = c_pos_to_genomic(tx, c_end, 0)
        if r0.error or r1.error:
            continue
        g_lo = min(r0.g_pos, r1.g_pos)
        g_hi = max(r0.g_pos, r1.g_pos)
        ref = genome[g_lo - 1: g_hi].decode("ascii")
        if length == 1:
            hgvs = f"{tx.accession}:c.{c_start}del"
        else:
            hgvs = f"{tx.accession}:c.{c_start}_{c_end}del"
        cases.append(Case(
            hgvs_string=hgvs,
            ref_type="c", edit_type="del", location_class="cds",
            strand="+", transcript=tx.accession,
            expected_chrom=CHROM_NAME, expected_pos=g_lo,
            expected_ref=ref, expected_alt="",
            notes=f"c. del len {length} (fwd)",
        ))
    # c. del on NM_synth_rev (minus strand)
    tx = TX_BY_NAME["NM_synth_rev"]
    for length in (1, 3):
        c_start = 10
        c_end = c_start + length - 1
        r0 = c_pos_to_genomic(tx, c_start, 0)
        r1 = c_pos_to_genomic(tx, c_end, 0)
        if r0.error or r1.error:
            continue
        g_lo = min(r0.g_pos, r1.g_pos)
        g_hi = max(r0.g_pos, r1.g_pos)
        ref = genome[g_lo - 1: g_hi].decode("ascii")
        if length == 1:
            hgvs = f"{tx.accession}:c.{c_start}del"
        else:
            hgvs = f"{tx.accession}:c.{c_start}_{c_end}del"
        cases.append(Case(
            hgvs_string=hgvs,
            ref_type="c", edit_type="del", location_class="cds",
            strand="-", transcript=tx.accession,
            expected_chrom=CHROM_NAME, expected_pos=g_lo,
            expected_ref=ref, expected_alt="",
            notes=f"c. del len {length} (rev)",
        ))


def _gen_ins_cases(genome: bytes, cases: list[Case]) -> None:
    """Generate ins cases."""
    # g. ins between two adjacent positions
    for length, alt in [(1, "T"), (2, "AT"), (5, "ACGTA")]:
        g_pos = 1700 + length * 10
        # HGVS notation: g.X_X+1insSEQ
        hgvs = f"{CHROM_NAME}:g.{g_pos}_{g_pos + 1}ins{alt}"
        cases.append(Case(
            hgvs_string=hgvs,
            ref_type="g", edit_type="ins", location_class="genomic", strand=".",
            transcript="",
            expected_chrom=CHROM_NAME, expected_pos=g_pos + 1,  # interbase ins site
            expected_ref="", expected_alt=alt,
            notes=f"g. ins len {length}",
        ))
    # c. ins on fwd transcript
    tx = TX_BY_NAME["NM_synth_fwd"]
    for length, alt in [(1, "G"), (2, "GA"), (5, "GATCA")]:
        c_pos = 20
        r0 = c_pos_to_genomic(tx, c_pos, 0)
        r1 = c_pos_to_genomic(tx, c_pos + 1, 0)
        if r0.error or r1.error:
            continue
        # On forward strand, insertion site after r0 (between r0 and r0+1).
        g_pos = max(r0.g_pos, r1.g_pos)
        hgvs = f"{tx.accession}:c.{c_pos}_{c_pos + 1}ins{alt}"
        cases.append(Case(
            hgvs_string=hgvs,
            ref_type="c", edit_type="ins", location_class="cds",
            strand="+", transcript=tx.accession,
            expected_chrom=CHROM_NAME, expected_pos=g_pos,
            expected_ref="", expected_alt=alt,
            notes=f"c. ins len {length} fwd",
        ))


def _gen_dup_cases(genome: bytes, cases: list[Case]) -> None:
    """Generate dup cases."""
    for length in (1, 3):
        g_pos = 1800 + length * 10
        ref = genome[g_pos - 1: g_pos - 1 + length].decode("ascii")
        if length == 1:
            hgvs = f"{CHROM_NAME}:g.{g_pos}dup"
        else:
            hgvs = f"{CHROM_NAME}:g.{g_pos}_{g_pos + length - 1}dup"
        cases.append(Case(
            hgvs_string=hgvs,
            ref_type="g", edit_type="dup", location_class="genomic", strand=".",
            transcript="",
            expected_chrom=CHROM_NAME, expected_pos=g_pos,
            expected_ref=ref, expected_alt=ref + ref,
            notes=f"g. dup len {length}",
        ))


def _gen_delins_cases(genome: bytes, cases: list[Case]) -> None:
    for ref_len, alt in [(1, "T"), (1, "ACG"), (3, "T"), (5, "AGT")]:
        g_pos = 1900 + ref_len * 10 + len(alt)
        ref = genome[g_pos - 1: g_pos - 1 + ref_len].decode("ascii")
        # Avoid edge case where ref==alt for 1->1
        if ref == alt:
            alt = "G" if alt != "G" else "C"
        if ref_len == 1:
            hgvs = f"{CHROM_NAME}:g.{g_pos}delins{alt}"
        else:
            hgvs = f"{CHROM_NAME}:g.{g_pos}_{g_pos + ref_len - 1}delins{alt}"
        cases.append(Case(
            hgvs_string=hgvs,
            ref_type="g", edit_type="delins", location_class="genomic", strand=".",
            transcript="",
            expected_chrom=CHROM_NAME, expected_pos=g_pos,
            expected_ref=ref, expected_alt=alt,
            notes=f"g. delins {ref_len}->{len(alt)}",
        ))


def _gen_identity_cases(genome: bytes, cases: list[Case]) -> None:
    # g.X= identity
    for g_pos in (2500, 4500):
        ref = chr(genome[g_pos - 1])
        cases.append(Case(
            hgvs_string=f"{CHROM_NAME}:g.{g_pos}=",
            ref_type="g", edit_type="identity", location_class="genomic", strand=".",
            transcript="",
            expected_chrom=CHROM_NAME, expected_pos=g_pos,
            expected_ref=ref, expected_alt=ref,
            notes="g. identity",
        ))


def _gen_negative_cases(genome: bytes, cases: list[Case]) -> None:
    # c. position past CDS bounds
    cases.append(Case(
        hgvs_string="NM_synth_simple:c.500A>T",
        ref_type="c", edit_type="sub", location_class="cds", strand="+",
        transcript="NM_synth_simple",
        expected_error="OutOfBounds",
        notes="CDS only 199bp; c.500 is past end",
    ))
    # tiny intron offset too large
    tiny = TX_BY_NAME["NM_synth_tiny"]
    # CDS layout for tiny: g 9003..9042. tx_len=36, exons 10/16/10.
    # cds_tx_start = ? exon0 has 10 bp (g 9001..9010). g 9003 -> tx_pos 3.
    # So c.1 is g 9003, c.8 is g 9010 (last base of exon 0). intron is 4bp (9011-9014).
    cases.append(Case(
        hgvs_string="NM_synth_tiny:c.8+50A>T",
        ref_type="c", edit_type="sub", location_class="intron_pos", strand="+",
        transcript="NM_synth_tiny",
        expected_error="IntronOffsetTooLarge",
        notes="intron is 4bp, offset is 50",
    ))
    cases.append(Case(
        hgvs_string="NM_synth_tiny:c.9-50A>T",
        ref_type="c", edit_type="sub", location_class="intron_neg", strand="+",
        transcript="NM_synth_tiny",
        expected_error="IntronOffsetTooLarge",
        notes="intron is 4bp, offset -50",
    ))
    # c.0 / c.-0 / c.*0
    for s in ("NM_synth_fwd:c.0A>T",):
        cases.append(Case(
            hgvs_string=s,
            ref_type="c", edit_type="sub", location_class="cds", strand="+",
            transcript="NM_synth_fwd",
            expected_error="ParseError",
            notes="HGVS forbids c.0",
        ))
    # c.*999 past end
    cases.append(Case(
        hgvs_string="NM_synth_fwd:c.*999A>T",
        ref_type="c", edit_type="sub", location_class="3utr", strand="+",
        transcript="NM_synth_fwd",
        expected_error="OutOfBounds",
        notes="past 3'UTR end",
    ))
    # g. past chromosome end
    cases.append(Case(
        hgvs_string=f"{CHROM_NAME}:g.{GENOME_LENGTH + 100}A>T",
        ref_type="g", edit_type="sub", location_class="genomic", strand=".",
        transcript="",
        expected_error="OutOfBounds",
        notes="past chromosome end",
    ))


def generate_cases(genome: bytes) -> list[Case]:
    cases: list[Case] = []
    _gen_substitution_cases(genome, cases)
    _gen_del_cases(genome, cases)
    _gen_ins_cases(genome, cases)
    _gen_dup_cases(genome, cases)
    _gen_delins_cases(genome, cases)
    _gen_identity_cases(genome, cases)
    _gen_negative_cases(genome, cases)

    # Assign deterministic case_ids in stable order.
    cases.sort(key=lambda c: (c.ref_type, c.transcript, c.hgvs_string))
    for i, c in enumerate(cases, 1):
        c.case_id = f"synth_{i:04d}"

    # Generator self-check on minus-strand cases:
    # for every minus-strand row, expected_ref reverse-complemented should
    # match the HGVS-displayed allele. (Skip rows where parser-rejected.)
    for c in cases:
        if c.strand == "-" and c.expected_ref and c.edit_type == "sub":
            # Extract the HGVS ref (the letter just before '>')
            try:
                ref_part, _, _ = c.hgvs_string.partition(">")
                hgvs_ref = ref_part[-1]
                rc = COMP_C(c.expected_ref)
                assert rc == hgvs_ref, (
                    f"GENERATOR BUG: case {c.case_id} {c.hgvs_string}: "
                    f"revcomp(expected_ref)={rc} != hgvs_ref={hgvs_ref}"
                )
            except KeyError:
                pass
    return cases


# ---------------------------------------------------------------------------
# TSV emission
# ---------------------------------------------------------------------------

CASES_HEADER = (
    "case_id\thgvs_string\tref_type\tedit_type\tlocation_class\tstrand\t"
    "transcript\texpected_chrom\texpected_pos\texpected_ref\texpected_alt\t"
    "expected_vrs_id\texpected_error\tnotes\n"
)

EQUIV_HEADER = (
    "group_id\tmember_kind\texpression\ttranscript\texpected_vrs_id\tnotes\n"
)


def emit_cases_tsv(path: Path, cases: list[Case]) -> None:
    with path.open("w") as fh:
        fh.write(f"# FIXTURE_VERSION={FIXTURE_VERSION}\n")
        fh.write(f"# GENOME_SEED=0x{GENOME_SEED:X}\n")
        fh.write(f"# Generated by gtars-vrs/tests/data/hgvs/synthetic/generate.py — do not hand-edit.\n")
        fh.write(CASES_HEADER)
        for c in cases:
            pos_str = "" if c.expected_pos == 0 else str(c.expected_pos)
            fh.write("\t".join([
                c.case_id, c.hgvs_string, c.ref_type, c.edit_type,
                c.location_class, c.strand, c.transcript,
                c.expected_chrom, pos_str, c.expected_ref, c.expected_alt,
                c.expected_vrs_id, c.expected_error, c.notes,
            ]) + "\n")


# ---------------------------------------------------------------------------
# Equivalence-group generation (plan 12).
#
# Coverage targets (only what the synthetic transcript layout permits):
#   - g + c.fwd + vcf 3-way groups: sub, del, ins, dup, delins, identity
#   - g + c.rev + vcf 3-way groups: sub, del, ins, dup, delins, identity
#   - g + n. + vcf 3-way groups (NR_synth_nc): a few subs
#   - g + vcf 2-way smoke groups
#
# 4-way (g + c + n + vcf) and 5-way (g + c.fwd + c.rev + n + vcf) groups are
# IMPOSSIBLE with the current synthetic transcript layout because no genomic
# region is overlapped by both a coding tx and the noncoding tx, and no region
# is covered by both NM_synth_fwd and NM_synth_rev. This is documented in the
# README. If the layout changes (e.g. an overlapping noncoding tx is added),
# extend this function.
# ---------------------------------------------------------------------------

@dataclass
class EquivMember:
    member_kind: str  # one of: hgvs_g, hgvs_c, hgvs_n, vcf
    expression: str
    transcript: str = ""
    notes: str = ""


@dataclass
class EquivGroup:
    group_id: str
    members: list[EquivMember] = field(default_factory=list)
    expected_vrs_id: str = ""
    # The canonical (chrom, pos, ref, alt) used to compute the VRS ID.
    chrom: str = ""
    pos: int = 0
    ref: str = ""
    alt: str = ""


def _tx_position_in_cds(tx: Tx, g_pos: int) -> Optional[int]:
    """If g_pos lies in CDS of `tx`, return c. coordinate (1-based, positive
    integer; CDS region only — not 5'/3'UTR / intronic). Else None."""
    if tx.cds is None:
        return None
    g_lo, g_hi = tx.cds
    if not (g_lo <= g_pos <= g_hi):
        return None
    tx_pos = genomic_to_tx_pos_in_exon(tx, g_pos)
    if tx_pos is None:
        return None
    bounds = cds_tx_bounds(tx)
    assert bounds is not None
    cds_tx_start, cds_tx_end = bounds
    c_pos = tx_pos - cds_tx_start + 1
    if 1 <= c_pos <= (cds_tx_end - cds_tx_start + 1):
        return c_pos
    return None


def _tx_n_position_in_exon(tx: Tx, g_pos: int) -> Optional[int]:
    """If g_pos lies in any exon of `tx`, return its n. coordinate (= 1-based
    transcript position). Else None."""
    return genomic_to_tx_pos_in_exon(tx, g_pos)


def _hgvs_g_for(edit_type: str, g_pos: int, ref: str, alt: str) -> str:
    """Render the genomic (g.) HGVS string for one of the supported edit types."""
    if edit_type == "sub":
        assert len(ref) == 1 and len(alt) == 1
        return f"{CHROM_NAME}:g.{g_pos}{ref}>{alt}"
    if edit_type == "identity":
        assert len(ref) == 1 and ref == alt
        return f"{CHROM_NAME}:g.{g_pos}="
    if edit_type == "del":
        if len(ref) == 1:
            return f"{CHROM_NAME}:g.{g_pos}del"
        return f"{CHROM_NAME}:g.{g_pos}_{g_pos + len(ref) - 1}del"
    if edit_type == "ins":
        # ref is empty; the cases-tsv convention uses pos = right anchor of the gap.
        # HGVS ins is X_X+1insSEQ where X = pos-1 in our convention.
        return f"{CHROM_NAME}:g.{g_pos - 1}_{g_pos}ins{alt}"
    if edit_type == "dup":
        # alt = ref + ref; the duplicated span starts at g_pos with len(ref).
        L = len(ref)
        if L == 1:
            return f"{CHROM_NAME}:g.{g_pos}dup"
        return f"{CHROM_NAME}:g.{g_pos}_{g_pos + L - 1}dup"
    if edit_type == "delins":
        L = len(ref)
        if L == 1:
            return f"{CHROM_NAME}:g.{g_pos}delins{alt}"
        return f"{CHROM_NAME}:g.{g_pos}_{g_pos + L - 1}delins{alt}"
    raise ValueError(f"unsupported edit_type for g.: {edit_type}")


def _revcomp_str(s: str) -> str:
    return "".join(COMP_C(c) for c in reversed(s))


def _hgvs_c_for(tx: Tx, edit_type: str, g_pos: int, ref: str, alt: str) -> Optional[str]:
    """Render a c. HGVS string for `tx` describing the edit at genomic g_pos.

    Returns None if the edit cannot be expressed (e.g. spans tx boundary or
    is not in CDS for the relevant kind). The transcript orientation is
    handled here: on minus-strand transcripts ref/alt are revcomp'd and the
    c. position(s) reflect transcript order.
    """
    if edit_type in ("sub", "identity"):
        # Single base.
        c_pos = _tx_position_in_cds(tx, g_pos)
        if c_pos is None:
            return None
        if tx.strand == +1:
            c_ref, c_alt = ref, alt
        else:
            c_ref, c_alt = COMP_C(ref), COMP_C(alt)
        if edit_type == "identity":
            return f"{tx.accession}:c.{c_pos}="
        return f"{tx.accession}:c.{c_pos}{c_ref}>{c_alt}"
    if edit_type == "del":
        # Span [g_pos .. g_pos+len(ref)-1] must be entirely in CDS.
        L = len(ref)
        c_lo = _tx_position_in_cds(tx, g_pos)
        c_hi = _tx_position_in_cds(tx, g_pos + L - 1)
        if c_lo is None or c_hi is None:
            return None
        # Must be exonic / contiguous (not crossing intron). For the picks
        # below we'll only ever choose tuples that lie strictly within one
        # exon; this assert documents that assumption.
        if abs(c_lo - c_hi) + 1 != L:
            return None
        c_start = min(c_lo, c_hi)
        c_end = max(c_lo, c_hi)
        if L == 1:
            return f"{tx.accession}:c.{c_start}del"
        return f"{tx.accession}:c.{c_start}_{c_end}del"
    if edit_type == "ins":
        # In our cases convention: pos points to the genomic base just AFTER
        # the insertion site (1-based). Anchors are pos-1 and pos.
        c_left = _tx_position_in_cds(tx, g_pos - 1)
        c_right = _tx_position_in_cds(tx, g_pos)
        if c_left is None or c_right is None:
            return None
        if abs(c_left - c_right) != 1:
            return None
        if tx.strand == +1:
            c_lo, c_hi = min(c_left, c_right), max(c_left, c_right)
            c_alt = alt
        else:
            # On minus strand the "left" tx anchor is the higher genomic pos,
            # and the alt sequence is revcomp'd because tx orientation is
            # opposite of genomic.
            c_lo, c_hi = min(c_left, c_right), max(c_left, c_right)
            c_alt = _revcomp_str(alt)
        return f"{tx.accession}:c.{c_lo}_{c_hi}ins{c_alt}"
    if edit_type == "dup":
        # ref is the duplicated block at [g_pos .. g_pos+len(ref)-1].
        L = len(ref)
        c_lo = _tx_position_in_cds(tx, g_pos)
        c_hi = _tx_position_in_cds(tx, g_pos + L - 1)
        if c_lo is None or c_hi is None:
            return None
        if abs(c_lo - c_hi) + 1 != L:
            return None
        c_start = min(c_lo, c_hi)
        c_end = max(c_lo, c_hi)
        if L == 1:
            return f"{tx.accession}:c.{c_start}dup"
        return f"{tx.accession}:c.{c_start}_{c_end}dup"
    if edit_type == "delins":
        L = len(ref)
        c_lo = _tx_position_in_cds(tx, g_pos)
        c_hi = _tx_position_in_cds(tx, g_pos + L - 1)
        if c_lo is None or c_hi is None:
            return None
        if abs(c_lo - c_hi) + 1 != L:
            return None
        c_start = min(c_lo, c_hi)
        c_end = max(c_lo, c_hi)
        c_alt = alt if tx.strand == +1 else _revcomp_str(alt)
        if L == 1:
            return f"{tx.accession}:c.{c_start}delins{c_alt}"
        return f"{tx.accession}:c.{c_start}_{c_end}delins{c_alt}"
    return None


def _hgvs_n_for(tx: Tx, edit_type: str, g_pos: int, ref: str, alt: str) -> Optional[str]:
    """Render an n. HGVS string for `tx` describing the edit at g_pos.

    Used for the noncoding transcript (NR_synth_nc), which is forward strand,
    so no revcomp.
    """
    assert tx.strand == +1, "minus-strand n. not exercised in this fixture"
    if edit_type in ("sub", "identity"):
        n_pos = _tx_n_position_in_exon(tx, g_pos)
        if n_pos is None:
            return None
        if edit_type == "identity":
            return f"{tx.accession}:n.{n_pos}="
        return f"{tx.accession}:n.{n_pos}{ref}>{alt}"
    if edit_type == "del":
        L = len(ref)
        n_lo = _tx_n_position_in_exon(tx, g_pos)
        n_hi = _tx_n_position_in_exon(tx, g_pos + L - 1)
        if n_lo is None or n_hi is None or abs(n_lo - n_hi) + 1 != L:
            return None
        if L == 1:
            return f"{tx.accession}:n.{n_lo}del"
        return f"{tx.accession}:n.{n_lo}_{n_hi}del"
    if edit_type == "ins":
        n_left = _tx_n_position_in_exon(tx, g_pos - 1)
        n_right = _tx_n_position_in_exon(tx, g_pos)
        if n_left is None or n_right is None or abs(n_left - n_right) != 1:
            return None
        n_lo, n_hi = min(n_left, n_right), max(n_left, n_right)
        return f"{tx.accession}:n.{n_lo}_{n_hi}ins{alt}"
    if edit_type == "dup":
        L = len(ref)
        n_lo = _tx_n_position_in_exon(tx, g_pos)
        n_hi = _tx_n_position_in_exon(tx, g_pos + L - 1)
        if n_lo is None or n_hi is None or abs(n_lo - n_hi) + 1 != L:
            return None
        if L == 1:
            return f"{tx.accession}:n.{n_lo}dup"
        return f"{tx.accession}:n.{n_lo}_{n_hi}dup"
    if edit_type == "delins":
        L = len(ref)
        n_lo = _tx_n_position_in_exon(tx, g_pos)
        n_hi = _tx_n_position_in_exon(tx, g_pos + L - 1)
        if n_lo is None or n_hi is None or abs(n_lo - n_hi) + 1 != L:
            return None
        if L == 1:
            return f"{tx.accession}:n.{n_lo}delins{alt}"
        return f"{tx.accession}:n.{n_lo}_{n_hi}delins{alt}"
    return None


def _vcf_member(g_pos: int, ref: str, alt: str) -> EquivMember:
    return EquivMember(
        member_kind="vcf",
        expression=f"{CHROM_NAME}:{g_pos}:{ref}:{alt}",
        transcript="",
    )


# Each pick: (edit_type, g_pos, ref_len_or_alt_seq_for_ins, label_suffix).
# The actual ref bytes and computed alt come from the synthetic genome — we
# never hand-bake biological bytes here, only positions.
#
# For sub:      we read 1 byte at g_pos, pick alt = "T" (or "A" if ref=="T").
# For identity: ref == alt = the genome byte at g_pos.
# For del:      ref = genome[g_pos .. g_pos+ref_len), alt = "".
# For ins:      ref = "", alt = the spec'd insertion sequence; pos = g_pos
#               (interbase site between g_pos-1 and g_pos), per cases.tsv
#               convention.
# For dup:      ref = genome[g_pos .. g_pos+ref_len), alt = ref + ref.
# For delins:   ref = genome[g_pos .. g_pos+ref_len), alt = the spec'd
#               replacement sequence (must differ from ref).
#
# All chosen g_pos values were verified to lie strictly within one exon of
# their target transcript so c./n. forms are renderable.
# Picks are the SAME (chrom, pos, ref, alt) tuples already present in
# cases.tsv (per the plan: "pick a curated subset of the same tuples
# already used in cases.tsv"). This guarantees the cross-link sanity check
# in the harness — every group's `expected_vrs_id` already appears in
# cases.tsv. To verify a pick: find its synth_NNNN row in cases.tsv and
# confirm (expected_pos, expected_ref, expected_alt, edit_type) match.
_EQ_PICKS_RAW: list[tuple[str, int, object, str]] = [
    # ---- Forward-strand CDS picks reusing NM_synth_fwd cases.tsv tuples ----
    ("sub",      2051, 1,         "fwd_sub_at_cds_start"),    # synth_0010 c.1C>T
    ("sub",      2100, 1,         "fwd_sub_in_cds"),          # synth_0017 c.50C>T
    ("sub",      3150, 1,         "fwd_sub_at_cds_stop"),     # synth_0014 c.400A>T
    ("del",      2060, 1,         "fwd_del_len1"),            # synth_0009 c.10del
    ("del",      2060, 2,         "fwd_del_len2"),            # synth_0007 c.10_11del
    ("del",      2060, 3,         "fwd_del_len3"),            # synth_0008 c.10_12del
    ("ins",      2071, "G",       "fwd_ins_len1"),            # synth_0011 c.20_21insG
    # ---- Reverse-strand CDS picks reusing NM_synth_rev cases.tsv tuples ----
    ("sub",      6101, 1,         "rev_sub_in_cds"),          # synth_0030 c.50T>A
    ("sub",      6150, 1,         "rev_sub_at_cds_start"),    # synth_0026 c.1A>T
    ("del",      6141, 1,         "rev_del_len1"),            # synth_0025 c.10del
    ("del",      6139, 3,         "rev_del_len3"),            # synth_0024 c.10_12del
    # ---- Noncoding transcript picks reusing NR_synth_nc cases.tsv tuples ----
    ("sub",      8001, 1,         "nc_sub_at_n1"),            # synth_0073 n.1A>T
    ("sub",      8050, 1,         "nc_sub_at_n50"),           # synth_0075 n.50G>T
    ("sub",      8550, 1,         "nc_sub_at_n200"),          # synth_0074 n.200T>A
    # ---- Genomic-only smoke groups (g + vcf 2-way) reusing intergenic cases.tsv tuples ----
    ("sub",      1500, 1,         "g_only_sub_1500"),         # synth_0054
    ("sub",      4000, 1,         "g_only_sub_4000"),         # synth_0070
    ("sub",      7000, 1,         "g_only_sub_7000"),         # synth_0072
    ("del",      1610, 1,         "g_only_del_len1"),         # synth_0055
    ("del",      1620, 2,         "g_only_del_len2"),         # synth_0056
    ("dup",      1810, 1,         "g_only_dup_len1"),         # synth_0062
    ("dup",      1830, 3,         "g_only_dup_len3"),         # synth_0063
    ("identity", 2500, 1,         "g_only_identity_2500"),    # synth_0068
    ("identity", 4500, 1,         "g_only_identity_4500"),    # synth_0071
]


def _resolve_pick(genome: bytes, edit_type: str, g_pos: int, spec: object
                  ) -> tuple[str, str]:
    """Translate a pick spec into concrete (ref, alt) strings using the genome."""
    def alt_for(ref_byte: str) -> str:
        return "T" if ref_byte != "T" else "A"

    if edit_type == "sub":
        assert isinstance(spec, int) and spec == 1
        ref = chr(genome[g_pos - 1])
        return ref, alt_for(ref)
    if edit_type == "identity":
        assert isinstance(spec, int) and spec == 1
        ref = chr(genome[g_pos - 1])
        return ref, ref
    if edit_type == "del":
        assert isinstance(spec, int) and spec >= 1
        ref = genome[g_pos - 1: g_pos - 1 + spec].decode("ascii")
        return ref, ""
    if edit_type == "ins":
        assert isinstance(spec, str) and len(spec) >= 1
        return "", spec
    if edit_type == "dup":
        assert isinstance(spec, int) and spec >= 1
        ref = genome[g_pos - 1: g_pos - 1 + spec].decode("ascii")
        return ref, ref + ref
    if edit_type == "delins":
        assert isinstance(spec, tuple) and len(spec) == 2
        ref_len, alt = spec
        ref = genome[g_pos - 1: g_pos - 1 + ref_len].decode("ascii")
        # Avoid degenerate ref==alt for 1->1 (would be an identity).
        if ref == alt:
            alt = "G" if alt != "G" else "C"
        return ref, alt
    raise ValueError(f"unknown edit_type {edit_type}")


_EQ_PICKS: list[tuple[str, int, str, str, Optional[str], str]] = []  # filled lazily


def _build_equiv_groups_from_picks(genome: bytes) -> list[EquivGroup]:
    """Walk picks; for each, derive (ref, alt) from the genome and build an
    EquivGroup with all applicable member kinds."""
    groups: list[EquivGroup] = []
    counters: dict[str, int] = {}
    for (edit_type, g_pos, spec, label) in _EQ_PICKS_RAW:
        ref, alt = _resolve_pick(genome, edit_type, g_pos, spec)

        # Stable group id: edit_type + label-suffix + counter (per label).
        key = f"eq_{label}"
        counters[key] = counters.get(key, 0) + 1
        group_id = f"{key}_{counters[key]:03d}"

        members: list[EquivMember] = []

        # Always emit the g. form.
        g_str = _hgvs_g_for(edit_type, g_pos, ref, alt)
        members.append(EquivMember(
            member_kind="hgvs_g",
            expression=g_str,
            transcript="",
            notes=f"genomic {edit_type}",
        ))

        # Walk all coding transcripts; emit c. members where the edit lies
        # entirely in CDS. (target_tx is just documentation; we still emit
        # any other tx that happens to overlap. None overlap each other, so
        # at most one c. emitted per pick in practice.)
        for tx in TRANSCRIPTS:
            if tx.cds is None:
                continue
            c_str = _hgvs_c_for(tx, edit_type, g_pos, ref, alt)
            if c_str is None:
                continue
            members.append(EquivMember(
                member_kind="hgvs_c",
                expression=c_str,
                transcript=tx.accession,
                notes=f"strand={_strand_str(tx.strand)}",
            ))

        # Walk noncoding (forward-strand) transcripts; emit n. members.
        for tx in TRANSCRIPTS:
            if tx.cds is not None:
                continue
            n_str = _hgvs_n_for(tx, edit_type, g_pos, ref, alt)
            if n_str is None:
                continue
            members.append(EquivMember(
                member_kind="hgvs_n",
                expression=n_str,
                transcript=tx.accession,
                notes="noncoding",
            ))

        # Always emit a vcf member.
        members.append(_vcf_member(g_pos, ref, alt))

        # If only g. + vcf made it (no transcript form), the group is still
        # a useful 2-way smoke test.
        if len(members) < 2:
            continue

        groups.append(EquivGroup(
            group_id=group_id,
            members=members,
            chrom=CHROM_NAME,
            pos=g_pos,
            ref=ref,
            alt=alt,
        ))
    return groups


def _compute_equiv_vrs_id(group: EquivGroup, fasta_path: Path) -> str:
    """Compute the canonical VRS ID for a group via the SAME path that
    populates `expected_vrs_id` in cases.tsv (gtars normalize+digest)."""
    from gtars.refget import RefgetStore, digest_sequence
    from gtars.vrs import vrs_id, normalize_allele

    seq_bytes = fasta_path.read_bytes()
    raw = bytearray()
    for line in seq_bytes.splitlines():
        if line.startswith(b">"):
            continue
        raw.extend(line)
    rec = digest_sequence(bytes(raw), name=CHROM_NAME)
    sq = "SQ." + rec.metadata.sha512t24u
    start_ib = group.pos - 1
    norm = normalize_allele(bytes(raw).decode("ascii"), start_ib, group.ref, group.alt)
    return vrs_id(sq, norm["start"], norm["end"], norm["allele"])


def generate_equivalence_groups(genome: bytes, fasta_path: Path,
                                with_vrs: bool = True) -> list[EquivGroup]:
    groups = _build_equiv_groups_from_picks(genome)
    # Sort for deterministic output.
    groups.sort(key=lambda g: g.group_id)
    if with_vrs:
        for g in groups:
            g.expected_vrs_id = _compute_equiv_vrs_id(g, fasta_path)
    return groups


def emit_equivalence_groups_json(path: Path, groups: list[EquivGroup]) -> None:
    rows = []
    for g in groups:
        for m in g.members:
            rows.append({
                "group_id": g.group_id,
                "member_kind": m.member_kind,
                "expression": m.expression,
                "transcript": m.transcript,
                "expected_vrs_id": g.expected_vrs_id,
                "notes": m.notes,
            })
    with path.open("w") as fh:
        json.dump(rows, fh, indent=2)
        fh.write("\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def regenerate(out_dir: Path, *, with_vrs: bool = True) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    seq = build_genome()
    write_fasta(out_dir / "synthetic.fa", CHROM_NAME, seq)
    write_fai(out_dir / "synthetic.fa.fai", CHROM_NAME, len(seq))
    emit_cdot(out_dir / "synthetic_transcripts.json")
    emit_refget_collection(out_dir / "refget_collection.json", seq)
    cases = generate_cases(seq)
    if with_vrs:
        vrs_ids = compute_vrs_ids_for_cases(cases, out_dir / "synthetic.fa")
        for c in cases:
            c.expected_vrs_id = vrs_ids.get(c.case_id, "")
    emit_cases_tsv(out_dir / "cases.tsv", cases)

    # Equivalence groups (plan 12). Built from a curated subset of biological
    # tuples; each group asserts that all member representations collapse to
    # the same VRS ID. expected_vrs_id is computed via the same VCF→VRS path
    # used for cases.tsv, never re-derived in the harness.
    # Output goes to fixtures/ alongside other test corpora.
    groups = generate_equivalence_groups(seq, out_dir / "synthetic.fa",
                                         with_vrs=with_vrs)
    fixtures_dir = out_dir.parent.parent / "fixtures"
    fixtures_dir.mkdir(parents=True, exist_ok=True)
    emit_equivalence_groups_json(fixtures_dir / "equivalence_groups.json", groups)


def diff_dirs(a: Path, b: Path, names: list[str]) -> list[str]:
    diffs = []
    for n in names:
        if not filecmp.cmp(a / n, b / n, shallow=False):
            diffs.append(n)
    return diffs


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", type=Path,
                   default=Path(__file__).resolve().parent,
                   help="output directory (default: this script's dir)")
    p.add_argument("--check", action="store_true",
                   help="regenerate to a tmpdir and diff against the checked-in dir")
    p.add_argument("--no-vrs", action="store_true",
                   help="skip layer-1 VRS-ID computation (development convenience)")
    args = p.parse_args(argv)

    if args.check:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            tmp_synthetic = tmp / "data" / "hgvs" / "synthetic"
            tmp_synthetic.mkdir(parents=True, exist_ok=True)
            regenerate(tmp_synthetic, with_vrs=not args.no_vrs)
            committed = Path(__file__).resolve().parent
            fixtures_committed = committed.parent.parent / "fixtures"
            # Check synthetic files
            synthetic_names = ["synthetic.fa", "synthetic.fa.fai",
                               "synthetic_transcripts.json", "refget_collection.json",
                               "cases.tsv"]
            diffs = diff_dirs(tmp_synthetic, committed, synthetic_names)
            # Check equivalence_groups.json in fixtures/
            tmp_fixtures = tmp / "fixtures"
            if not filecmp.cmp(tmp_fixtures / "equivalence_groups.json",
                               fixtures_committed / "equivalence_groups.json",
                               shallow=False):
                diffs.append("equivalence_groups.json")
            if diffs:
                print(f"FIXTURE DRIFT: {diffs}", file=sys.stderr)
                return 1
            print("OK: synthetic fixtures match generator output")
            return 0
    regenerate(args.out_dir, with_vrs=not args.no_vrs)
    print(f"Wrote synthetic fixtures to {args.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

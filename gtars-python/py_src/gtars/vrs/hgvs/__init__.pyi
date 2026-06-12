"""Type stubs for the gtars.vrs.hgvs module.

The real implementation lives in `gtars-python/src/vrs/hgvs.rs`.
"""

from enum import Enum
from typing import Any, Dict, Optional

from gtars.refget import RefgetStore
from gtars.reftx import ReftxProvider


class ReferenceType(Enum):
    Coding: int
    NonCoding: int
    Genomic: int
    Mitochondrial: int
    Protein: int
    Rna: int


class Datum(Enum):
    SeqStart: int
    Cds: int
    CdsStop: int


class Position:
    base: int
    offset: int
    datum: Datum

    def to_dict(self) -> Dict[str, Any]: ...


class Edit:
    """Edit payload with discriminator `kind` of:
    `"substitution"`, `"deletion"`, `"insertion"`, `"delins"`,
    `"duplication"`, `"inversion"`, `"identity"`, `"unknown"`.
    """

    kind: str
    ref: Optional[str]
    alt: Optional[str]

    def to_dict(self) -> Dict[str, Any]: ...


class PositionBound:
    kind: str          # "certain" | "uncertain"
    position: Optional[Position]   # set when kind == "certain"
    low: Optional[Position]        # uncertain low bound; None == "?"
    high: Optional[Position]       # uncertain high bound; None == "?"

    def to_dict(self) -> Dict[str, Any]: ...


class PosEdit:
    location_kind: str   # "single" | "range" | "whole_sequence"
    start: Optional[PositionBound]
    end: Optional[PositionBound]
    edit: Edit
    uncertain: bool

    def to_dict(self) -> Dict[str, Any]: ...


class HgvsVariant:
    accession: str
    gene: Optional[str]
    reference_type: ReferenceType
    pos_edit: PosEdit

    def to_dict(self) -> Dict[str, Any]: ...


def parse_hgvs(s: str) -> HgvsVariant: ...


def hgvs_to_vrs_id(
    hgvs_str: str,
    provider: ReftxProvider,
    refget: RefgetStore,
    collection_digest: str,
) -> str: ...


class HgvsError(Exception): ...

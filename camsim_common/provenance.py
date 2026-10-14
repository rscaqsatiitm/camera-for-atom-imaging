"""Provenance-tagged values.

Every number that enters a simulation carries a tag stating where it came from:

USER         supplied by the PI in the project prompt or in a USER-supplied document
SOURCE       retrieved and read in a session, with a pointer (arXiv/DOI/URL and location)
SIM          produced by our own code (pointer = config file and output path)
DATA         produced by an analysis of data (incl. figure read-off), with method
DERIVED      computed from other tagged values (pointer = formula or report section)
PLACEHOLDER  recalled or assumed, not yet verified; rendered as [VALUE?: approx. X, verify]
"""
from __future__ import annotations

from enum import Enum
from typing import Any

from pydantic import BaseModel, ConfigDict, field_validator


class Prov(str, Enum):
    USER = "USER"
    SOURCE = "SOURCE"
    SIM = "SIM"
    DATA = "DATA"
    DERIVED = "DERIVED"
    PLACEHOLDER = "PLACEHOLDER"


class ProvValue(BaseModel):
    """A scalar or array value with unit and provenance."""

    model_config = ConfigDict(extra="forbid")

    value: Any
    unit: str = ""
    prov: Prov
    pointer: str = ""
    note: str = ""

    @field_validator("pointer")
    @classmethod
    def _pointer_required(cls, v: str, info):
        prov = info.data.get("prov")
        if prov in (Prov.SOURCE, Prov.USER, Prov.SIM, Prov.DATA, Prov.DERIVED) and not v:
            raise ValueError(f"provenance {prov} requires a pointer")
        return v

    @property
    def is_placeholder(self) -> bool:
        return self.prov == Prov.PLACEHOLDER

    def tag(self) -> str:
        """Human-readable rendering used in reports."""
        if self.is_placeholder:
            return f"[VALUE?: approx. {self.value} {self.unit}, verify]".replace("  ", " ")
        return f"{self.value} {self.unit} ({self.prov.value}: {self.pointer})".replace("  ", " ")

    def __float__(self) -> float:
        return float(self.value)


def pv(value: Any, unit: str, prov: str | Prov, pointer: str = "", note: str = "") -> ProvValue:
    """Short constructor: pv(0.44, '', 'SOURCE', 'arXiv:2403.12021v4 Methods')."""
    return ProvValue(value=value, unit=unit, prov=Prov(prov), pointer=pointer, note=note)


def placeholder(approx: Any, unit: str, note: str = "") -> ProvValue:
    """A value we only recall or assume. Swept or flagged downstream."""
    return ProvValue(value=approx, unit=unit, prov=Prov.PLACEHOLDER, pointer="", note=note)


def val(x: Any) -> Any:
    """Strip provenance: return the raw value of a ProvValue, else x itself."""
    return x.value if isinstance(x, ProvValue) else x

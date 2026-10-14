"""Helpers shared by the validation scripts."""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy import optimize

from camsim_common.provenance import val

ROOT = Path(__file__).resolve().parents[2]


def v(cfg_entry):
    """Value of a config entry written as {value, unit, prov, pointer} or a bare number."""
    return cfg_entry["value"] if isinstance(cfg_entry, dict) and "value" in cfg_entry else cfg_entry


def solve_monotone(f, target, lo, hi, **kw):
    """Root of f(x) = target on [lo, hi]; returns None if not bracketed."""
    flo, fhi = f(lo) - target, f(hi) - target
    if flo * fhi > 0:
        return None
    return optimize.brentq(lambda x: f(x) - target, lo, hi, **kw)


def fmt(x, sig=3):
    if x is None:
        return "n/a"
    if isinstance(x, str):
        return x
    if x == 0:
        return "0"
    ax = abs(x)
    if ax < 1e-3 or ax >= 1e5:
        m, e = f"{x:.{sig - 1}e}".split("e")
        return f"{m}e{int(e)}"
    return f"{x:.{sig}g}"


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    for r in rows:
        out.append("| " + " | ".join(fmt(c) if not isinstance(c, str) else c for c in r) + " |")
    return "\n".join(out)

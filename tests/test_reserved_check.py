"""Checker for the reserved validation points: predicates and their negative controls (R5 b, c).

Criteria are those of configs/production/reserved_validation_2026-11-19.yaml (committed 2026-11-19);
the constants in p1_atomdet/sweep/reserved_check.py transcribe them. These tests were written and
committed before the checker's first run on any reserved point.
"""
import math

import numpy as np
import pytest

from p1_atomdet.counts.lattice import Lattice
from p1_atomdet.sweep import reserved_check as rc


def test_reserved_file_loads_twenty_points():
    pts = rc.load_reserved(rc.ROOT / "configs/production/reserved_validation_2026-11-19.yaml")
    assert len(pts) == 20 and len({p.key for p in pts}) == 20
    assert sum(p.design == "camera" for p in pts) == 10
    assert all(len(p.values) == (6 if p.design == "camera" else 8) for p in pts)


def test_criteria_constants_match_the_registered_file():
    assert (rc.R1_TIMES_S, rc.R1_N_PER_NODE, rc.R1_MIN_EXPECTED, rc.R1_FAMILY_ALPHA) == ((0.020, 0.080), 400_000, 100, 0.05)
    assert (rc.R2_DS, rc.R2_N_TAU, rc.R2_N_OFFSETS, rc.R2_TOL_LOG10) == (0.01, 48, 5, 0.0043)
    assert (rc.R3_MIN_ST, rc.R5_DELTA_FRAC) == (0.05, 0.10)


def _r2(base_log10, refined_log10, reached=(1.0, 1.0)):
    mk = lambda lt, r: {"reached_1e-03": r, "log10t_1e-03": lt}
    return {"targets": [1e-3], "base": mk(base_log10, reached[0]), "refined": mk(refined_log10, reached[1])}


def test_r2_clean_twin_and_negative_control():
    res = _r2(-1.5, -1.5)
    assert rc.r2_verdict(res)["1e-03"]["pass"]                                    # clean twin
    assert not rc.r2_verdict(res, add_to_refined=0.01)["1e-03"]["pass"]           # R5(b)
    assert not rc.r2_verdict(_r2(-1.5, -1.5, reached=(1.0, 0.0)))["1e-03"]["pass"]  # reached flag differs


def _rows(st, conf):
    return [{"param": n, "ST": s, "ST_conf": c} for n, s, c in zip("abc", st, conf)]


def test_r3_clean_twin_and_negative_control():
    a = _rows([0.6, 0.3, 0.02], [0.05, 0.04, 0.01])
    assert rc.r3_verdict(a, a)["pass"]                                            # clean twin
    b = _rows([0.6, 0.3 + 2 * 0.04, 0.02], [0.05, 0.04, 0.01])                    # R5(c): 2 x half-width
    v = rc.r3_verdict(a, b)
    assert not v["pass"] and not v["indices"]["b"]["pass"] and v["n_checked"] == 2


def test_r1_predicate_detects_a_wrong_rate():
    n = 3_600_000
    band = (1.0e-3, 1.0e-3)
    k_ok = int(round(1.0e-3 * n))
    assert rc.r1_agrees(band, k_ok, n, 0.95)                                      # clean twin
    assert not rc.r1_agrees(band, int(round(1.3e-3 * n)), n, 0.95)                # rate off by 30%
    assert rc.testable(1e-3, n) and not rc.testable(1e-6, n)


def test_r1_band_brackets_the_reported_rate():
    p0 = Lattice(np.array([0.5, 0.3, 0.2]), 0.0, 1.0)
    p1 = Lattice(np.array([0.1, 0.3, 0.6]), 0.0, 1.0)
    nodes = [{"weight": 1.0, "threshold": 1.0, "p0": p0, "p1": p1}]
    (flo, fhi), (nlo, nhi) = rc.r1_band(nodes)
    assert flo <= p0.sf(1.0) <= fhi and nlo <= p1.cdf_lt(1.0) <= nhi


def test_same_result_treats_nan_as_equal():
    assert rc.same_result({"t": math.nan, "x": 1.0}, {"t": math.nan, "x": 1.0})
    assert not rc.same_result({"t": math.nan, "x": 1.0}, {"t": math.nan, "x": 1.0000001})

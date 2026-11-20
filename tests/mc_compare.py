"""Comparators shared by the positive checks and their negative controls.

Keeping the pass predicate in one place means a negative control exercises the same code as the
check it controls, not a copy of it.
"""
from p1_atomdet.fidelity import clopper_pearson

# Negative-control perturbation, fixed on 2026-11-19 before the first run (CLAUDE.md rule 10):
# the classifier applies T_opt +/- DELTA_FRAC x (mean score with atom - mean score without atom)
# while the rates are reported at T_opt. Do not change after seeing results.
DELTA_FRAC = 0.10


def mc_agrees(kind, reported, p0, p1, k_fp, k_fn, n):
    """True if MC counts (k_fp, k_fn out of n each) agree with the lattice rates reported at reported.threshold.

    Box: the reported FP and FN must each lie inside the 95% Clopper-Pearson interval of the MC.
    Matched / Bayes: the statistic is rebinned on the lattice while the MC uses exact values, so the
    lattice rates over +/- one lattice step around the threshold must overlap the MC interval.
    (Same predicate as the original test_lattice_vs_bruteforce_mc_at_1em3.)
    """
    lo_fp, hi_fp = clopper_pearson(k_fp, n)
    lo_fn, hi_fn = clopper_pearson(k_fn, n)
    if kind == "box":
        return (lo_fp <= reported.fp <= hi_fp) and (lo_fn <= reported.fn <= hi_fn)
    T, h = reported.threshold, p0.h
    fp_band = (p0.sf(T + h), p0.sf(T - h))
    fn_band = (p1.cdf_lt(T - h), p1.cdf_lt(T + h))
    return (fp_band[0] <= hi_fp and fp_band[1] >= lo_fp) and (fn_band[0] <= hi_fn and fn_band[1] >= lo_fn)


def wrong_thresholds(p0, p1, t_opt):
    """The two deliberately wrong thresholds of the negative control."""
    d = DELTA_FRAC * (p1.mean() - p0.mean())
    return t_opt - d, t_opt + d

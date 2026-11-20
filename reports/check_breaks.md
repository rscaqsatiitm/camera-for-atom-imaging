# Deliberate breaks of new checks (CLAUDE.md rule 11)

Each new check is broken once on purpose, run, seen to fail, and restored. Append new entries; do
not delete. Criteria of every check were committed before its first run (commit 0515d23, then
a1c2742 for the reserved-point criteria).

## 2026-11-19: negative controls and placeholder-list check

Checks added in commit 0515d23: `tests/test_negative_controls.py` (wrong-threshold controls with
clean twins), `tests/test_placeholder_lists.py` (report lists equal scanner output). The positive
MC checks in `tests/test_counts.py` were refactored to call the shared predicate
`tests/mc_compare.py::mc_agrees` (same logic as before).

### First run, unbroken (before any break)

12 passed, 1 failed. All five CMOS and CF controls fired at both wrong thresholds and their clean
twins passed. Wall time 124 s.

**Finding F1 (the registered control does not fire in one case).**
`test_wrong_threshold_fails_mc_check_emccd_box` failed: the low wrong threshold passed the check.

    reported at T_opt=2.000: FP=2.341e-02 FN=1.134e-01; MC:
    T_opt:  T=2.000 FP=2.300e-02 FN=1.137e-01 agrees=True
    T_low:  T=1.710 FP=2.300e-02 FN=1.137e-01 agrees=True
    T_high: T=2.290 FP=1.410e-03 FN=3.330e-01 agrees=False

Root cause: the EMCCD case is in photon-counting mode, so the score is an integer click count.
Every threshold in (1, 2] makes the same decision; the MC counts at T_low are identical to those at
T_opt. The registered shift, 0.10 x separation = 0.29 click, is below the one-click resolution of
the statistic, so the "wrong" threshold is not wrong in effect. The upward shift reaches the next
integer (effective threshold 3) and the check rejects it. Hence the comparator works; the
perturbation rule is degenerate for integer-valued scores.

Action: the test is kept unchanged and marked `xfail(strict=True)` with this diagnosis (a pass
would then be reported as an error). The perturbation was not enlarged after the result (rule 10).
A replacement control for integer-valued scores (for example, a shift of at least one statistic
step) is a decision for the PI and would be a new, separately registered check.

### B1: perturbation set to zero (`DELTA_FRAC = 0.0` in `tests/mc_compare.py`)

Expected: every control reports that a wrong threshold passed; the registration guard fails.
Observed: 6 failed, 1 xfailed (82 s). The five controls failed with "wrong threshold (low) passed
the check" (T_low = T_opt, identical counts), and `test_delta_frac_is_the_registered_value` failed.
Restored with `git checkout tests/mc_compare.py`.

### B2: comparators reject everything (`mc_agrees` returns False; `CF_REL_TOL = 0.0`)

Expected: every clean twin fails; the existing positive MC checks fail.
Observed: 9 failed, 1 passed, 1 xfailed (116 s). All five clean twins failed with "clean twin
failed" (for CF: lattice eps 7.0102e-04 against CF 7.0107e-04, box; 2.1374e-04 equal to 4 digits,
matched, rejected at zero tolerance), and all four positive checks in `tests/test_counts.py`
(box, matched, Bayes, EMCCD) failed. The guard test passed, as it should. Restored by hand; no
break marker remains (`grep "DELIBERATE BREAK"` empty).

### B3: placeholder lists out of step with the scanner

(a) One entry (`dark.Ea_eV` of the ORCA-Quest preset) deleted from the list in
`reports/p1_validation.md`: `test_validation_report_list_matches_scanner` failed.
(b) First entry deleted from the list in `reports/p1_production_note.md`:
`test_production_note_list_matches_full_scan` failed.
(c) A PLACEHOLDER entry appended to `configs/validation/v3_analytic.yaml`: both tests failed.
Each restored with `git checkout`; `git status` showed no change to reports or configs.

## 2026-11-20: step-aware control for integer scores (resolves F1)

Registered in commit dc14c24 before its first run: `tests/mc_compare.py::wrong_thresholds_step_aware`,
shift = max(0.10 x separation, one step of the no-atom score lattice), used by
`test_wrong_threshold_fails_mc_check_emccd_photon_counting`. The strict-xfail test of F1 was removed;
its diagnosis above stands.

First run (unbroken): passed.

    reported T_opt=2.000 FP=2.341e-02 FN=1.134e-01
    T_opt:  T=2.000 FP=2.300e-02 FN=1.137e-01 agrees=True   (clean twin)
    T_low:  T=1.000 FP=2.185e-01 FN=1.795e-02 agrees=False
    T_high: T=3.000 FP=1.410e-03 FN=3.330e-01 agrees=False

### B4: one-step floor removed (`d = DELTA_FRAC x separation` only)

Expected: the F1 failure returns. Observed: 1 failed, "wrong threshold (low) passed the check",
T_low = 1.710 with counts identical to T_opt, as in F1. Restored with `git checkout tests/mc_compare.py`.

The B2 break (comparator rejects everything) already covers the clean-twin assertion, which this
test shares through `_assert_control`.

## Validation-report SNR cut, 2026-11-20

`report.py` now flags Table I cells above `SNR_ROUNDING = 0.05` (was 0.06). This sets the cut to the
criterion the text stated, so it tightens rather than loosens. It is a report classification, not a
pass/fail check; its effect is visible in the report (3 of 18 cells flagged instead of 2).

## Pending

The checker for the reserved-point criteria (`configs/production/reserved_validation_2026-11-19.yaml`)
is not written; its deliberate breaks are criterion R5 of that file and must be recorded here
before its output on the reserved points is read.

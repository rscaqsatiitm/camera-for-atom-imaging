"""Checker for the reserved validation points of the production Sobol sweeps.

    python -m p1_atomdet.sweep.reserved_check --reserved configs/production/reserved_validation_2026-11-19.yaml \
        [--criteria R1,R2] [--points camera:0,full:3] [--workers 4] [--out runs/production/reserved_check]
    python -m p1_atomdet.sweep.reserved_check ... --perturb r1      # negative control R5(a), verdicts only

Implements criteria R1, R2 and R4 of the reserved-points file as written there (committed 2026-11-19,
before this checker existed); R3 is `r3_verdict`, to be applied to the full N = 1024 results.
Nothing here may change a criterion: tolerances are read from the definitions below, which
transcribe the file, and a change to them needs a dated amendment of the file first.

R1 compares the production CF engine (`production.eps_at`, matched filter, as called by
`sweep.evaluate.sobol_point`) with brute-force Monte Carlo (`counts.mc`) at 20 ms and 80 ms. For
each sub-pixel offset node the MC uses the same SiteModel, the CF matched-filter weights and the
node's reported threshold; counts are pooled over the equally weighted nodes. Pooling Bernoulli
trials with unequal success probabilities gives a variance below the binomial one, so the
Clopper-Pearson interval on the pooled count is conservative. Shared layers (not tested by R1):
optics and pixel fractions, SiteModel inputs, TauDistribution, offset nodes and weights, the
classifier weights and the thresholds.
"""
from __future__ import annotations

import argparse
import json
import time
from dataclasses import dataclass
from multiprocessing import get_context
from pathlib import Path

import numpy as np
import yaml

from camsim_common.rng import stream

from ..validation.common import ROOT

# ---- criteria, transcribed from configs/production/reserved_validation_2026-11-19.yaml -------------
R1_TIMES_S = (0.020, 0.080)
R1_N_PER_NODE = 400_000          # samples per class per offset node
R1_MIN_EXPECTED = 100            # testable if samples x production rate >= 100
R1_FAMILY_ALPHA = 0.05           # Bonferroni: per-comparison level 1 - alpha / M
R2_DS, R2_N_TAU, R2_N_OFFSETS = 0.01, 48, 5
R2_TOL_LOG10 = 0.0043            # 1% in t
R3_MIN_ST = 0.05
R5_DELTA_FRAC = 0.10             # rule of tests/mc_compare.py
MC_SEED = 20261119               # base seed of the R1 Monte Carlo streams (set with this checker)


@dataclass
class ReservedPoint:
    design: str                  # "camera" | "full"
    label: str                   # "random_0".. | corner name
    config: str
    values: dict

    @property
    def key(self) -> str:
        return f"{self.design}:{self.label}"


def load_reserved(path) -> list[ReservedPoint]:
    d = yaml.safe_load(Path(path).read_text())
    out = []
    for design, block in (("camera", d["camera_design"]), ("full", d["full_design"])):
        cols = block["columns"]
        rows = [(f"random_{i}", r) for i, r in enumerate(block["random"])] + list(block["corners"].items())
        for label, r in rows:
            out.append(ReservedPoint(design, label, block["config"], dict(zip(cols, map(float, r)))))
    return out


def _point_dict(rp: ReservedPoint) -> dict:
    cfg = yaml.safe_load((ROOT / rp.config).read_text())
    return {**cfg["fixed"], **rp.values}, cfg


def _explicit(p: dict, scenario_overrides: dict | None = None) -> dict:
    """The point with its scenario file merged in explicitly (as sweep.evaluate._with_file does),
    optionally with Scenario attributes overridden."""
    f = yaml.safe_load((ROOT / p["scenario_file"]).read_text())
    q = {**f, **{k: x for k, x in p.items() if k != "scenario_file"}}
    q["scenario"] = {**f.get("scenario", {}), **(scenario_overrides or {})}
    return q


def production_state(p: dict):
    """Scenario, Camera and Geometry exactly as sobol_point builds them."""
    from ..production import Camera, geometry, read_classes
    from .evaluate import _scenario, _with_file
    p = _with_file(p)
    sc = _scenario(p)
    sc.n_offsets = int(p.get("n_offsets", 2))
    sc.offset_rule = p.get("offset_rule", "midpoint")
    sc.ds = float(p.get("ds", 0.02))
    cam = Camera(qe=p["qe"], read_classes=read_classes(p["read_median_e"], p.get("tail_fraction", 0.0),
                                                         p.get("tail_factor", 5.0)),
                 pixels_per_fwhm=p["px_per_fwhm"], dark_e_per_px_s=p.get("dark_e_per_px_s", 0.0))
    return sc, cam, geometry(sc, cam)


# ---- R1 ---------------------------------------------------------------------------------------------
def r1_production(rp: ReservedPoint) -> list[dict]:
    """Production side of R1: per time, the offset-averaged rates and per-node records."""
    from ..production import eps_at
    p, _ = _point_dict(rp)
    sc, cam, geo = production_state(p)
    out = []
    for t in R1_TIMES_S:
        det = []
        r = eps_at(sc, cam, geo, t, "matched", details=det)
        out.append({"t": t, "fp": r.fp, "fn": r.fn, "nodes": det, "cam": cam})
    return out


def r1_band(nodes):
    """Offset-averaged lattice rates over +/- one lattice step around each node's threshold
    (the band rule of tests/mc_compare.py::mc_agrees for rebinned statistics)."""
    fp_lo = sum(n["weight"] * n["p0"].sf(n["threshold"] + n["p0"].h) for n in nodes)
    fp_hi = sum(n["weight"] * n["p0"].sf(n["threshold"] - n["p0"].h) for n in nodes)
    fn_lo = sum(n["weight"] * n["p1"].cdf_lt(n["threshold"] - n["p1"].h) for n in nodes)
    fn_hi = sum(n["weight"] * n["p1"].cdf_lt(n["threshold"] + n["p1"].h) for n in nodes)
    return (fp_lo, fp_hi), (fn_lo, fn_hi)


def r1_mc(nodes, cam, seed_key, shift_sign: int = 0) -> dict:
    """Pooled MC counts over the offset nodes. shift_sign != 0 applies the R5(a) wrong threshold."""
    from ..counts.mc import mc_error_rates
    from ..sensor.cmos import PixelReadout
    ws = np.array([n["weight"] for n in nodes])
    if not np.allclose(ws, ws[0]):
        raise ValueError("R1 pooling assumes equally weighted offset nodes (midpoint rule)")
    ro = PixelReadout(classes=cam.read_classes, gain_e_per_adu=None)
    k_fp = k_fn = 0
    for j, n in enumerate(nodes):
        T = n["threshold"]
        if shift_sign:
            T = T + shift_sign * R5_DELTA_FRAC * (n["p1"].mean() - n["p0"].mean())
        m = mc_error_rates(n["site"], ro, T, R1_N_PER_NODE, stream(MC_SEED, *seed_key, j), kind="matched",
                           weights=n["w"])
        k_fp += m["fp"]
        k_fn += m["fn"]
    return {"k_fp": k_fp, "k_fn": k_fn, "n": R1_N_PER_NODE * len(nodes)}


def testable(rate: float, n_total: int) -> bool:
    return rate * n_total >= R1_MIN_EXPECTED


def r1_agrees(band, k, n, cl) -> bool:
    """Band (lo, hi) of the production rate overlaps the MC Clopper-Pearson interval at level cl."""
    from ..fidelity import clopper_pearson
    lo, hi = clopper_pearson(k, n, cl)
    return band[0] <= hi and band[1] >= lo


# ---- R2 ---------------------------------------------------------------------------------------------
def same_result(a: dict, b: dict) -> bool:
    """Exact equality of two evaluator outputs, with NaN (unreached target) equal to NaN."""
    return a.keys() == b.keys() and all(
        (isinstance(a[k], float) and isinstance(b[k], float) and np.isnan(a[k]) and np.isnan(b[k])) or a[k] == b[k]
        for k in a)


def r2_eval(rp: ReservedPoint) -> dict:
    from .evaluate import sobol_point
    p, cfg = _point_dict(rp)
    seed = 0
    base = sobol_point(p, seed)
    # plumbing null: the explicit merge used for the refinement must reproduce the production path exactly
    base_explicit = sobol_point(_explicit(p), seed)
    if not same_result(base_explicit, base):
        raise AssertionError(f"explicit scenario merge differs from the production path at {rp.key}")
    q = _explicit(p, {"n_tau": R2_N_TAU})
    q["ds"], q["n_offsets"] = R2_DS, R2_N_OFFSETS
    ref = sobol_point(q, seed)
    return {"targets": cfg["fixed"]["targets"], "base": base, "refined": ref}


def r2_verdict(res: dict, add_to_refined: float = 0.0) -> dict:
    """Per target: reached flags identical, and |delta log10 t| <= R2_TOL_LOG10 where reached.
    add_to_refined perturbs the refined log10 t (negative control R5(b))."""
    out = {}
    for tg in res["targets"]:
        k = f"{tg:.0e}"
        rb, rr = res["base"][f"reached_{k}"], res["refined"][f"reached_{k}"]
        if rb != rr:
            out[k] = {"pass": False, "reason": "reached flag differs", "reached": (rb, rr)}
            continue
        if not rb:
            out[k] = {"pass": True, "reason": "unreached in both", "delta": None}
            continue
        d = res["refined"][f"log10t_{k}"] + add_to_refined - res["base"][f"log10t_{k}"]
        out[k] = {"pass": abs(d) <= R2_TOL_LOG10, "delta": d}
    return out


# ---- R3 ---------------------------------------------------------------------------------------------
def r3_verdict(rows_1024: list[dict], rows_512: list[dict]) -> dict:
    """rows_*: output of sweep.sobol.analyse (param, ST, ST_conf, ...). Criterion R3."""
    a = {r["param"]: r for r in rows_1024}
    b = {r["param"]: r for r in rows_512}
    checks = {}
    for name, r in a.items():
        if r["ST"] >= R3_MIN_ST:
            d = r["ST"] - b[name]["ST"]
            checks[name] = {"pass": abs(d) <= r["ST_conf"], "delta": d, "half_width": r["ST_conf"]}
    top = lambda rows: {x["param"] for x in sorted(rows, key=lambda x: -x["ST"])[:2]}
    top_ok = top(rows_1024) == top(rows_512)
    return {"indices": checks, "top2_same": top_ok,
            "pass": top_ok and all(c["pass"] for c in checks.values()), "n_checked": len(checks)}


# ---- driver -----------------------------------------------------------------------------------------
def _r1_job(args):
    rp, shift = args
    t0 = time.process_time()
    prod = r1_production(rp)
    res = []
    for e in prod:
        fpb, fnb = r1_band(e["nodes"])
        n_tot = R1_N_PER_NODE * len(e["nodes"])
        rec = {"point": rp.key, "t": e["t"], "fp": e["fp"], "fn": e["fn"], "fp_band": fpb, "fn_band": fnb,
               "fp_testable": testable(e["fp"], n_tot), "fn_testable": testable(e["fn"], n_tot), "n": n_tot}
        if rec["fp_testable"] or rec["fn_testable"]:
            rec.update(r1_mc(e["nodes"], e["cam"], (rp.design, rp.label, int(round(e["t"] * 1000))), shift))
        res.append(rec)
    return rp.key, res, time.process_time() - t0


def _r2_job(rp):
    t0 = time.process_time()
    r = r2_eval(rp)
    return rp.key, r, time.process_time() - t0


def run_r1(points, workers=1, shift=0):
    jobs = [(rp, shift) for rp in points]
    if workers <= 1:
        out = [_r1_job(j) for j in jobs]
    else:
        with get_context("spawn").Pool(workers) as pool:
            out = list(pool.imap_unordered(_r1_job, jobs, chunksize=1))
    recs = [r for _, rs, _ in out for r in rs]
    cpu = {k: c for k, _, c in out}
    M = sum(r["fp_testable"] for r in recs) + sum(r["fn_testable"] for r in recs)
    cl = 1 - R1_FAMILY_ALPHA / max(M, 1)
    for r in recs:
        for q in ("fp", "fn"):
            if r[f"{q}_testable"]:
                r[f"{q}_pass"] = r1_agrees(r[f"{q}_band"], r[f"k_{q}"], r["n"], cl)
    return {"M": M, "cl": cl, "records": recs, "cpu_s": cpu}


def run_r2(points, workers=1):
    if workers <= 1:
        out = [_r2_job(rp) for rp in points]
    else:
        with get_context("spawn").Pool(workers) as pool:
            out = list(pool.imap_unordered(_r2_job, points, chunksize=1))
    res = {k: {"raw": r, "verdict": r2_verdict(r)} for k, r, _ in out}
    return {"results": res, "cpu_s": {k: c for k, _, c in out}}


def summarise(points, r1=None, r2=None) -> dict:
    """R4: denominators. A point missing from a criterion that was run is a failure."""
    keys = [rp.key for rp in points]
    s = {"n_points": len(keys)}
    if r1 is not None:
        seen = {r["point"] for r in r1["records"]}
        tested = [(r, q) for r in r1["records"] for q in ("fp", "fn") if r[f"{q}_testable"]]
        s["R1"] = {"points_evaluated": len(seen), "comparisons_total": 2 * len(r1["records"]),
                   "testable": len(tested), "passed": sum(r[f"{q}_pass"] for r, q in tested),
                   "untestable": 2 * len(r1["records"]) - len(tested), "missing_points": sorted(set(keys) - seen),
                   "cl": r1["cl"]}
        s["R1"]["pass"] = (s["R1"]["passed"] == s["R1"]["testable"] and not s["R1"]["missing_points"])
    if r2 is not None:
        seen = set(r2["results"])
        flat = [v for k in seen for v in r2["results"][k]["verdict"].values()]
        s["R2"] = {"points_evaluated": len(seen), "comparisons": len(flat), "passed": sum(v["pass"] for v in flat),
                   "missing_points": sorted(set(keys) - seen)}
        s["R2"]["pass"] = (s["R2"]["passed"] == s["R2"]["comparisons"] and not s["R2"]["missing_points"])
    return s


def _select(points, spec):
    if not spec:
        return points
    want = set(spec.split(","))
    sel = [rp for rp in points if rp.key in want or f"{rp.design}:{points.index(rp)}" in want]
    if len(sel) != len(want):
        raise SystemExit(f"unknown point(s) in {spec}")
    return sel


def _json(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, tuple):
        return list(o)
    raise TypeError(type(o))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--reserved", default=str(ROOT / "configs/production/reserved_validation_2026-11-19.yaml"))
    ap.add_argument("--criteria", default="R1,R2")
    ap.add_argument("--points", default="", help="comma list of design:label (default: all 20)")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--perturb", choices=["none", "r1"], default="none",
                    help="r1: MC applies T +/- 0.10 x separation (negative control R5(a)); prints verdicts only")
    ap.add_argument("--out", default=str(ROOT / "runs/production/reserved_check"))
    a = ap.parse_args(argv)
    pts = _select(load_reserved(a.reserved), a.points)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    if a.perturb == "r1":
        rows = []
        for sign in (-1, +1):
            r1 = run_r1(pts, a.workers, shift=sign)
            rows += [{"sign": sign, "point": r["point"], "t": r["t"], **{f"{q}_pass": r.get(f"{q}_pass")
                      for q in ("fp", "fn")}} for r in r1["records"]]
        failed = [r for r in rows if r["fp_pass"] is False or r["fn_pass"] is False]
        print(json.dumps({"R5a_rows": rows, "R5a_pass": bool(failed)}, indent=1, default=_json))
        return
    crit = set(a.criteria.split(","))
    r1 = run_r1(pts, a.workers) if "R1" in crit else None
    r2 = run_r2(pts, a.workers) if "R2" in crit else None
    summ = summarise(pts, r1, r2)
    payload = {"reserved_file": a.reserved, "summary": summ}
    if r1 is not None:
        payload["R1"] = r1
    if r2 is not None:
        payload["R2"] = r2
    tag = "_".join(sorted(crit)) + ("_" + a.points.replace(":", "-").replace(",", "+") if a.points else "")
    (out / f"reserved_check_{tag}.json").write_text(json.dumps(payload, indent=1, default=_json))
    print(json.dumps(summ, indent=1, default=_json))


if __name__ == "__main__":
    main()

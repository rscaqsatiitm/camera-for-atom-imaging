"""Build reports/p1_production_note.md (production check-in) from the collected sweep HDF5 files.

    python -m p1_atomdet.production_report
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from camsim_common.io import load_h5  # noqa: E402
from camsim_common.placeholders import open_placeholders, placeholder_markdown  # noqa: E402
from camsim_common.schema import SensorModel, read_noise_median, read_noise_rms  # noqa: E402

from .validation.common import ROOT, fmt, md_table  # noqa: E402
from .validation.figures import C, INK, INK2, _fig  # noqa: E402

RUNS = ROOT / "runs" / "production"
FIG = ROOT / "reports" / "figures"


def _load(name):
    d, a = load_h5(RUNS / name / f"{name}.h5")
    P = {k.split("/", 1)[1]: v for k, v in d.items() if k.startswith("params/")}
    R = {k.split("/", 1)[1]: v for k, v in d.items() if k.startswith("results/")}
    Cst = {k.split("/", 1)[1]: v for k, v in d.items() if k.startswith("cost/")}
    for k, v in P.items():
        if v.dtype.kind == "S":
            P[k] = v.astype(str)
    return P, R, Cst, a


def _ms(x):
    return "n/r" if not np.isfinite(x) else f"{1000 * x:.1f}"


def preset_table(P, R):
    rows, best = [], {}
    names = list(dict.fromkeys(P["preset"]))
    for nm in names:
        sel = np.nonzero(P["preset"] == nm)[0]
        clf = "matched" if np.isfinite(R.get("matched_t_1e-04", np.full(len(P["preset"]), np.nan))[sel]).any() else "box"
        t4 = R[f"{clf}_t_1e-04"][sel]
        t4f = np.where(np.isfinite(t4), t4, np.inf)
        j = sel[int(np.argmin(t4f))]
        best[nm] = (j, clf)
        sm = SensorModel.from_yaml(ROOT / "camsim_common" / "presets" / f"{nm}.yaml")
        rows.append([sm.camera.replace("Teledyne Photometrics ", "").replace("Hamamatsu ", "").replace("Andor ", ""),
                     sm.mode.split(",")[0], f"{float(sm.qe_852.value):.3f}",
                     f"{read_noise_median(sm.read_noise):.2f}/{read_noise_rms(sm.read_noise):.2f}" if sm.emccd is None else "EM x1000",
                     f"{P['px_per_fwhm'][j]:g}", f"{R['magnification'][j]:.1f}", clf if sm.emccd is None or sm.emccd.mode == "analog" else "PC box",
                     _ms(R[f"{clf}_t_1e-03"][j]), _ms(R[f"{clf}_t_1e-04"][j]), fmt(R[f"{clf}_eps_40ms"][j]),
                     f"{fmt(R[f'{clf}_fp_40ms'][j])} / {fmt(R[f'{clf}_fn_40ms'][j])}"])
    order = np.argsort([float(r[8]) if r[8] != "n/r" else 1e9 for r in rows])
    rows = [rows[i] for i in order]
    tab = md_table(["camera", "mode", "QE 852", "read median/rms (e-)", "best px/FWHM", "magnification", "classifier",
                    "t(1e-3) ms", "t(1e-4) ms", "eps(40 ms)", "FP / FN (40 ms)"], rows)
    return tab, best, names


def sampling_table(P, R, names):
    pp = sorted(set(P["px_per_fwhm"]))
    rows = []
    for nm in names:
        r = [nm]
        for q in pp:
            sel = np.nonzero((P["preset"] == nm) & (P["px_per_fwhm"] == q))[0]
            if sel.size == 0:
                r.append("-")
                continue
            j = sel[0]
            v = R.get("matched_t_1e-04", np.full(len(P["preset"]), np.nan))[j]
            if not np.isfinite(v):
                v = R["box_t_1e-04"][j]
            r.append(_ms(v))
        rows.append(r)
    return md_table(["preset \\ px per FWHM"] + [f"{q:g}" for q in pp], rows)


def fig_presets(P, R, best, names, out):
    labs, t3, t4 = [], [], []
    for nm in names:
        j, clf = best[nm]
        labs.append(nm)
        t3.append(1000 * R[f"{clf}_t_1e-03"][j])
        t4.append(1000 * R[f"{clf}_t_1e-04"][j])
    o = np.argsort(np.where(np.isfinite(t4), t4, 1e9))
    fig, (ax,) = _fig(w=7.2, h=4.2)
    y = np.arange(len(o))
    ax.scatter(np.array(t3)[o], y, color=C[0], s=36, marker="o", label="t to eps = 1e-3", zorder=3)
    ax.scatter(np.array(t4)[o], y, color=C[1], s=36, marker="s", label="t to eps = 1e-4", zorder=3)
    for yy, a, b in zip(y, np.array(t3)[o], np.array(t4)[o]):
        ax.plot([a, b], [yy, yy], color=INK2, lw=0.8, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels([labs[i] for i in o], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("imaging time (ms), best pixel sampling", color=INK, fontsize=9)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def generic_views(P, R, out_png):
    """Heat maps of t(1e-4) (matched) in the (QE, read-noise median) plane; ppf = 1.0, background x1."""
    qe = np.array(sorted(set(P["qe"])))
    rn = np.array(sorted(set(P["read_median_e"])))
    fig, axs = _fig(w=8.4, h=3.6, n=2)
    tabs = {}
    for ax, tail in zip(axs, (0.0, 0.05)):
        Z = np.full((rn.size, qe.size), np.nan)
        for i, r in enumerate(rn):
            for j, q in enumerate(qe):
                sel = np.nonzero((P["qe"] == q) & (P["read_median_e"] == r) & (P["tail_fraction"] == tail)
                                 & (P["px_per_fwhm"] == 1.0) & (P["bg_scale"] == 1.0))[0]
                if sel.size:
                    Z[i, j] = 1000 * R["matched_t_1e-04"][sel[0]]
        tabs[tail] = Z
        Zp = np.where(np.isfinite(Z), Z, np.nanmax(Z) * 1.5 if np.isfinite(Z).any() else 1)
        ax.imshow(np.log10(Zp), origin="lower", cmap="Blues", aspect="auto")
        for i in range(rn.size):
            for j in range(qe.size):
                txt = "n/r" if not np.isfinite(Z[i, j]) else f"{Z[i, j]:.0f}"
                ax.text(j, i, txt, ha="center", va="center", fontsize=7,
                        color="#fcfcfb" if np.isfinite(Z[i, j]) and Z[i, j] > np.nanmedian(Z) else INK)
        ax.set_xticks(range(qe.size))
        ax.set_xticklabels([f"{q:g}" for q in qe], fontsize=8)
        ax.set_yticks(range(rn.size))
        ax.set_yticklabels([f"{r:g}" for r in rn], fontsize=8)
        ax.set_xlabel("QE at 852 nm", color=INK, fontsize=9)
        ax.set_ylabel("read noise, median (e- rms)", color=INK, fontsize=9)
        ax.set_title(f"tail fraction {tail:g}", fontsize=9, color=INK)
        ax.grid(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)
    return qe, rn, tabs


def one_at_a_time(P, R, ref, key="matched_t_1e-04"):
    """Range of t over each axis with the others at the reference point (pre-Sobol ranking)."""
    rows = []
    for ax, vals in ((k, sorted(set(P[k]))) for k in ref):
        ts = []
        for v in vals:
            m = np.ones(len(P["qe"]), bool)
            for k, rv in ref.items():
                m &= (P[k] == (v if k == ax else rv))
            sel = np.nonzero(m)[0]
            ts.append(1000 * R[key][sel[0]] if sel.size else np.nan)
        ts = np.array(ts)
        fin = ts[np.isfinite(ts)]
        rows.append((ax, vals, ts, (fin.max() / fin.min()) if fin.size > 1 else np.nan))
    rows.sort(key=lambda r: -np.nan_to_num(r[3], nan=1e9))
    return md_table(["parameter", "values", "t(1e-4) (ms), others at reference", "max/min"],
                    [[a, ", ".join(f"{v:g}" for v in vals), ", ".join(_ms(t / 1000) for t in ts), f"{r:.2f}" if np.isfinite(r) else "n/r (target not reached)"]
                     for a, vals, ts, r in rows])


def sensitivity_table(P, R):
    from .sweep.evaluate import SENSITIVITY_CASES
    rows = []
    for case in dict.fromkeys(P["case"]):
        for q in sorted(set(P["px_per_fwhm"])):
            sel = np.nonzero((P["case"] == case) & (P["px_per_fwhm"] == q))[0]
            if not sel.size:
                continue
            j = sel[0]
            clf = "matched" if np.isfinite(R.get("matched_t_1e-04", np.full(len(P["case"]), np.nan))[j]) else "box"
            preset, change = SENSITIVITY_CASES[case]
            rows.append([case, preset, ", ".join(f"{k}={v}" for k, v in change.items()) or "baseline", f"{q:g}", clf,
                         _ms(R[f"{clf}_t_1e-03"][j]), _ms(R[f"{clf}_t_1e-04"][j]), fmt(R[f"{clf}_eps_40ms"][j])])
    return md_table(["case", "preset", "change", "px/FWHM", "classifier", "t(1e-3) ms", "t(1e-4) ms", "eps(40 ms)"], rows)


def cost_summary(*costs):
    rows = []
    for name, Cst in costs:
        rows.append([name, len(Cst["cpu_h"]), f"{Cst['cpu_h'].sum():.2f}", f"{3600 * np.median(Cst['cpu_h']):.0f}",
                     f"{3600 * Cst['cpu_h'].max():.0f}", f"{Cst['peak_rss_mb'].max():.0f}"])
    return md_table(["study", "points", "CPU-h total", "median CPU-s/point", "max CPU-s/point", "peak RSS (MB)"], rows)


def build(out_md: Path = ROOT / "reports" / "p1_production_note.md") -> Path:
    FIG.mkdir(parents=True, exist_ok=True)
    Pp, Rp, Cp, _ = _load("prod_presets")
    Pg, Rg, Cg, _ = _load("prod_generic_reduced")
    Ps, Rs, Cs, _ = _load("prod_sensitivity")
    ptab, best, names = preset_table(Pp, Rp)
    stab = sampling_table(Pp, Rp, names)
    fig_presets(Pp, Rp, best, names, FIG / "prod_presets.png")
    qe, rn, maps = generic_views(Pg, Rg, FIG / "prod_generic_map.png")
    ref = {"qe": 0.6, "read_median_e": 1.0, "tail_fraction": 0.0, "px_per_fwhm": 0.7, "bg_scale": 1.0}
    oat4 = one_at_a_time(Pg, Rg, ref, "matched_t_1e-04")
    oat3 = one_at_a_time(Pg, Rg, ref, "matched_t_1e-03")
    sens = sensitivity_table(Ps, Rs)
    cost = cost_summary(("presets", Cp), ("generic (reduced)", Cg), ("sensitivity", Cs))

    # best sampling per (qe, read noise) in the generic grid (tail 0, bg x1)
    rows = []
    for r in sorted(set(Pg["read_median_e"])):
        row = [f"{r:g}"]
        for q in sorted(set(Pg["qe"])):
            m = (Pg["qe"] == q) & (Pg["read_median_e"] == r) & (Pg["tail_fraction"] == 0) & (Pg["bg_scale"] == 1.0)
            sel = np.nonzero(m)[0]
            t = Rg["matched_t_1e-04"][sel]
            if np.isfinite(t).any():
                k = sel[int(np.nanargmin(t))]
                row.append(f"{Pg['px_per_fwhm'][k]:g} ({1000 * Rg['matched_t_1e-04'][k]:.0f} ms)")
            else:
                row.append("n/r")
        rows.append(row)
    best_ppf = md_table(["read median (e-) \\ QE"] + [f"{q:g}" for q in sorted(set(Pg["qe"]))], rows)
    text = NOTE_TEMPLATE.format(ptab=ptab, stab=stab, best_ppf=best_ppf, oat4=oat4, oat3=oat3, sens=sens, cost=cost,
                                placeholders=placeholder_markdown(open_placeholders()),
                                n_open=len(open_placeholders()))
    out_md.write_text(text)
    return out_md


NOTE_TEMPLATE = ""  # filled below after the results are read (see bottom of file)


if __name__ == "__main__":
    print(build())

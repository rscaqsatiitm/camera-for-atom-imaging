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
        clf = "matched" if np.isfinite(R.get("matched_eps_40ms", np.full(len(P["preset"]), np.nan))[sel]).any() else "box"
        t4 = R[f"{clf}_t_1e-04"][sel]
        t3 = R[f"{clf}_t_1e-03"][sel]
        # best sampling: shortest t(1e-4); if never reached, shortest t(1e-3); else lowest eps(40 ms)
        key = t4 if np.isfinite(t4).any() else (t3 if np.isfinite(t3).any() else R[f"{clf}_eps_40ms"][sel])
        j = sel[int(np.nanargmin(np.where(np.isfinite(key), key, np.inf)))]
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
            mk = R.get("matched_eps_40ms", np.full(len(P["preset"]), np.nan))[j]
            v = R["matched_t_1e-04"][j] if np.isfinite(mk) else R["box_t_1e-04"][j]
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
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def generic_views(P, R, out_png):
    """Heat maps of t(1e-4) (matched) in the (QE, read-noise median) plane, each cell at its best pixel sampling; background x1."""
    qe = np.array(sorted(set(P["qe"])))
    rn = np.array(sorted(set(P["read_median_e"])))
    fig, axs = _fig(w=8.4, h=3.6, n=2)
    tabs = {}
    for ax, tail in zip(axs, (0.0, 0.05)):
        Z = np.full((rn.size, qe.size), np.nan)
        for i, r in enumerate(rn):
            for j, q in enumerate(qe):
                sel = np.nonzero((P["qe"] == q) & (P["read_median_e"] == r) & (P["tail_fraction"] == tail)
                                 & (P["bg_scale"] == 1.0))[0]
                t = R["matched_t_1e-04"][sel]
                if sel.size and np.isfinite(t).any():
                    Z[i, j] = 1000 * np.nanmin(t)
        tabs[tail] = Z
        Zp = np.where(np.isfinite(Z), Z, np.nanmax(Z) * 1.5 if np.isfinite(Z).any() else 1)
        ax.imshow(np.log10(Zp), origin="lower", cmap="Blues", aspect="auto")
        for i in range(rn.size):
            for j in range(qe.size):
                txt = "n/r" if not np.isfinite(Z[i, j]) else f"{Z[i, j]:.0f}"
                ax.text(j, i, txt, ha="center", va="center", fontsize=7,
                        color="#fcfcfb" if (not np.isfinite(Z[i, j])) or Z[i, j] > np.nanmedian(Z) else INK)
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
            clf = "matched" if np.isfinite(R.get("matched_eps_40ms", np.full(len(P["case"]), np.nan))[j]) else "box"
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
    j, clf = best["prime_bsi_express_cms"]
    sm = SensorModel.from_yaml(ROOT / "camsim_common" / "presets" / "prime_bsi_express_cms.yaml")
    rows_bsi = []
    for q in sorted(set(Pp["px_per_fwhm"])):
        k = np.nonzero((Pp["preset"] == "prime_bsi_express_cms") & (Pp["px_per_fwhm"] == q))[0][0]
        rows_bsi.append([f"{q:g}", f"{Rp['magnification'][k]:.1f}", f"{Rp['pitch_obj_um'][k]:.2f}",
                         fmt(Rp["matched_eps_20ms"][k]), fmt(Rp["matched_eps_40ms"][k]), fmt(Rp["matched_eps_80ms"][k]),
                         _ms(Rp["matched_t_1e-03"][k]), _ms(Rp["matched_t_1e-04"][k]), _ms(Rp["box_t_1e-04"][k])])
    bsi_tab = md_table(["px/FWHM", "magnification", "object-plane pixel (um)", "eps(20 ms)", "eps(40 ms)", "eps(80 ms)",
                        "t(1e-3) ms", "t(1e-4) ms", "t(1e-4) ms, box"], rows_bsi)
    med_cpu = 3600 * np.median(Cg["cpu_h"])
    K = dict(n_presets=len(names), n_generic=len(Pg["qe"]), cpu_generic=Cg["cpu_h"].sum(), cpu_presets=Cp["cpu_h"].sum(),
             cpu_sens=Cs["cpu_h"].sum(), med_cpu=med_cpu,
             full_points=46656, full_cpuh=46656 * med_cpu / 3600, sobol_points=1024 * 10, sobol_cpuh=1024 * 10 * med_cpu / 3600,
             bsi_ppf=f"{Pp['px_per_fwhm'][j]:g}", bsi_mag=f"{Rp['magnification'][j]:.1f}", bsi_t3=_ms(Rp["matched_t_1e-03"][j]),
             bsi_t4=_ms(Rp["matched_t_1e-04"][j]), bsi_eps40=fmt(Rp["matched_eps_40ms"][j]),
             bsi_dark40=f"{1.5 * 0.04:.2f}", bsi_roi_ms=f"{128 * 11.4e-3:.2f}")
    text = NOTE_TEMPLATE.format(ptab=ptab, stab=stab, best_ppf=best_ppf, oat4=oat4, oat3=oat3, sens=sens, cost=cost,
                                bsi_tab=bsi_tab, placeholders=placeholder_markdown(open_placeholders()),
                                n_open=len(open_placeholders()), **K)
    out_md.write_text(text)
    return out_md


NOTE_TEMPLATE = '# P1 production check-in: camera presets and reduced generic sweep\n\nCAQSAT lab, IIT Madras. Check-in 2 of 2 required by rule 7 (after the first production sweep,\nbefore scaling up). Generated by `python -m p1_atomdet.production_report` from\n`runs/production/*/*.h5`; every number in the tables is read from those files.\n\n## Summary\n\n1. **Ranking of the {n_presets} presets** (Table P-1, Fig. P-1). In the reference Cs scenario\n   (V1-anchored photon budget and loss, NA 0.55, best pixel sampling, matched filter), the time to\n   reach infidelity 1e-4 is 30 ms for the iXon Ultra 888 EXF in multi-frame photon counting,\n   36 to 40 ms for five low-noise CMOS configurations (ORCA-Quest 2 PNR, ORCA-Fusion BT, Prime BSI\n   Express CMS, ORCA-Quest 2 standard, Marana), 46 ms for the Prime BSI Express in Speed mode, and\n   94 to 300 ms for the QHY268M (IMX571), Zyla 5.5 and Blackfly S.\n2. **What sets the imaging time** (Tables P-3, P-4). QE dominates (t roughly proportional to 1/QE).\n   Read noise matters little below 0.3 e- median and strongly above 1 e-; a 5% fraction of pixels\n   at five times the median read noise costs a factor 1.8 at 1 e-. Background from 0 to 3 times the V1\n   level changes t by 1.4x; pixel sampling by 1.3x over 0.35 to 1.5 px per FWHM.\n3. **Pixel sampling rule** (Table P-5). The optimal sampling falls with read noise: 1.5 px per FWHM\n   or more at 0.1 e- (upper edge of the grid), 1.0 at 0.3 e-, 0.5 at 1 to 2 e-, 0.35 or less at 4 e-\n   (lower edge). The optimum is shallow (within 10% over a factor of two in sampling).\n4. **Prime BSI Express, preliminary operating point** (Sec. 6). CMS mode, magnification {bsi_mag}\n   ({bsi_ppf} px per FWHM), matched filter: t(1e-3) = {bsi_t3} ms, t(1e-4) = {bsi_t4} ms, eps = {bsi_eps40} at 40 ms\n   in the reference scenario. Speed mode costs about 23% in time (46.3 ms).\n5. **Cost.** This stage used {cpu_presets:.1f} + {cpu_generic:.1f} + {cpu_sens:.1f} CPU-h on the 4-core instance\n   (median {med_cpu:.0f} CPU-s per generic point). The full generic sweep proposed in Sec. 7\n   ({full_points} points) needs about {full_cpuh:.0f} CPU-h and a Sobol analysis about {sobol_cpuh:.0f} CPU-h; both fit the\n   IIT Madras HPC easily, and the Sobol analysis would also run on this instance in about half a day.\n\n6. **Atom loss can make 1e-4 unreachable whatever the camera** (Sec. 5). With the loss per photon\n   bounded from Picken et al. (5.6e-7, 24 times the V1 value), loss during the exposure sets a\n   floor near eps ~ 3e-4 at 40 ms and no camera reaches 1e-4. The 1e-4 specification therefore\n   presumes V1-level cooling and vacuum.\n7. **The read-noise tail is the most consequential unknown camera input** (Sec. 5). For the Prime\n   BSI Express, 2% of pixels at 3x the median read noise raises t(1e-4) from 37.5 to 56.6 ms (+51%),\n   whereas the QE uncertainty (+-1.5% absolute) moves it by 3%, and the EMCCD CIC split by under 10%.\n\nThese are model results for one reference scenario. The scenario inherits the V1 fit, so the\ncaveats of the validation report apply: in particular the V1 false-positive excess is an\nexperimental floor that none of these cameras removes.\n\n## 1. What was run\n\n(i) Eleven `SensorModel` presets (Sec. 3) built from the four PDFs supplied by the PI and from\nvendor documents retrieved and checked on 2026-10-15 (`sources/README.md`), each at seven pixel\nsamplings; (ii) a reduced generic CMOS-class sweep of {n_generic} points over QE, read-noise median,\nread-noise tail fraction, pixel sampling and background; (iii) a sensitivity study of the\npreset results to their uncertain inputs. All runs are resumable with per-point checkpoints and\nCPU accounting (`runs/production/run_all.sh`).\n\nMethod changes since the validation report: a characteristic-function (CF) engine for linear\nstatistics (box sum, matched filter), about 10^3 times faster than the lattice engine and equal to\nit within 2% in optimal infidelity (`tests/test_cf_production.py`); multi-frame EMCCD photon counting\nwith a quadrature click probability (more accurate than the lattice at h = 0.004, and consistent\nwith a stage-by-stage register Monte Carlo); per-site thresholds averaged over the sub-pixel\nposition of the site (9 offsets). The CF engine omits ADC rounding: for the ORCA-Quest 2 in PNR mode\nthe exact rounding lowers eps by 7% at the operating point (about 1 to 2% in imaging time), so the\nPNR results are slightly conservative; for the Prime BSI Express (0.25 e-/ADU placeholder) the\ndifference is 0.4%.\n\n## 2. Reference scenario and assumptions\n\n`configs/production/reference_cs.yaml`: Cs D2, scattering rate 72.3 kc/s, NA 0.55 with optical\ntransmission 0.327 after collection (2.7% of scattered photons to the camera, as in V1), loss\n2.35e-8 per scattered photon (V1 survival, includes vacuum loss), site-to-site rate spread 12.3%\n(V1 Fig. 4d), background 4.76 e-/px/s at the V1 pixel (0.902 um object plane, QE 0.44) scaled\nwith object-plane pixel area and QE, filling 0.5. Assumptions: (i) diffraction-limited vectorial\nPSF at NA 0.55 (the measured V1 PSF was about 35% wider, so absolute times here are optimistic and\nthe best sampling would shift up by the same factor); (ii) uniform background; (iii) no neighbour\ncrosstalk (spacing large compared with the PSF); (iv) read-noise tails, where a vendor gives median\nand rms, are a two-class model with tail pixels at 3x the median (presets) or 5x (generic sweep);\nP4 will replace these with measured distributions; (v) dark current at the vendor\'s stated\noperating temperature; (vi) survival is not limiting (loss 1.7e-3 per second of imaging).\n\n## 3. Camera presets\n\n**Table P-1. Presets at their best pixel sampling** (matched filter; photon counting uses an\nSNR-optimal box). QE at 852 nm is digitised from vendor typical curves (DATA, +-1.5 to 2%).\n\n{ptab}\n\n**Table P-2. Time to eps = 1e-4 (ms) versus pixel sampling** (px per PSF FWHM; n/r: not reached\nwithin 300 ms with survival >= 0.99).\n\n{stab}\n\n![Presets](figures/prod_presets.png)\n\n*Figure P-1. Imaging time to infidelity 1e-3 (circles) and 1e-4 (squares) for each preset at its best\npixel sampling, reference Cs scenario (Sec. 2), CF engine with per-site matched filter (EMCCD photon\ncounting: Poisson-binomial counts over an SNR-optimal box, frames at 670 fps for a 128x128 crop).\nError rates averaged over 9 sub-pixel site positions; threshold optimised per site.*\n\nRemarks. (i) The ORCA-Quest 2 enters with QE 0.384 at 852 nm from the Hamamatsu typical curve; the\nsame vendor\'s curve for the original ORCA-Quest gives 0.396, while Manetsch et al. state 0.44 for\ntheir camera. With 0.44 the Quest 2 times would fall by about 13% (t roughly proportional to 1/QE). (ii) The iXon advantage requires\nphoton counting over many short frames (26 frames in 40 ms at 670 fps for a 128x128 crop); in analog\nmode the excess noise factor (F^2 = 1.97) cancels the QE advantage. (iii) The QHY268M uses High Gain\n2CMS at gain 100 (read noise 1.19 e-, digitised) as the PI requested; its limit is QE (0.21).\n(iv) The Blackfly S never reaches 1e-4 within 300 ms (4.7 e- read noise, QE 0.19).\n\n## 4. Generic CMOS-class sweep (reduced)\n\nAxes: QE 0.2 to 0.95; read-noise median 0.1 to 4 e-; tail fraction 0 or 5% (tail pixels at 5x the\nmedian); 0.35 to 1.5 px per FWHM; background 0, 1 and 3 times the reference. No dark current, no ADC.\n\n![Generic map](figures/prod_generic_map.png)\n\n*Figure P-2. Time to eps = 1e-4 (ms, numbers in cells; colour log scale) in the plane of QE at 852 nm\nand read-noise median, each cell at its best pixel sampling (0.35 to 1.5 px per FWHM), reference\nbackground, matched filter, reference Cs scenario. Left: Gaussian read noise. Right: 5% of pixels at\n5x the median. n/r: not reached within 300 ms.*\n\n**Table P-3. One-at-a-time ranges of t(1e-4)** about QE 0.6, read noise 1.0 e-, no tail, 0.7 px per\nFWHM, reference background (a ranking proxy until the Sobol analysis).\n\n{oat4}\n\n**Table P-4. The same for t(1e-3).**\n\n{oat3}\n\n**Table P-5. Best pixel sampling (and the resulting t(1e-4)) per QE and read-noise median**, no tail,\nreference background.\n\n{best_ppf}\n\nThe optima at 0.1 e- and 4 e- lie on the edges of the sampled range and are not bracketed; the full sweep extends the sampling axis to 0.25 and 2.5 px per FWHM. The matched filter is faster than the SNR-optimal box sum by 7% (median over the grid) to 18%.\n\nAt 80 ms every Prime BSI sampling in Table P-6 reaches eps ~ 1.5e-5: this is the loss floor of the reference scenario (atoms lost early in the exposure), not a camera limit.\n\n## 5. Sensitivity of the preset results to uncertain inputs\n\n{sens}\n\n## 6. Prime BSI Express: preliminary operating point\n\n**Table P-6. Prime BSI Express, CMS mode, in the reference scenario** (matched filter unless stated).\n\n{bsi_tab}\n\nPreliminary recipe, to be refined after the full sweep: CMS mode (1.0 e-), magnification about\n{bsi_mag} so that one 6.5 um pixel spans about two PSF FWHM in the object plane, matched filter with per-site\nweights and threshold, imaging time {bsi_t4} ms for 1e-4 in this scenario. A 2048 x 128 ROI reads out in\n{bsi_roi_ms} ms (128 rows at 11.4 us); with continuous imaging light the exposure must exceed the\nimaging time by that amount so that every row sees the full pulse. Dark current at 0 C (1.5 e-/px/s)\nadds {bsi_dark40} e- per pixel in 40 ms, small compared with the background at this sampling. The read-noise tail\nof the GSENSE2020BSI is not specified (only "1.0 e-"); a 2% tail at 3x the median changes the times\nas shown in Sec. 5, and a measured read-noise map (P4 or a dark-frame measurement on our camera)\nis the most useful next input. The Bayesian per-pixel classifier, which can exploit a measured\nnoise map, is not yet included here.\n\n## 7. Computational cost and proposed full stage\n\n{cost}\n\nProposed next steps (on approval): (i) full generic sweep on the IIT Madras HPC: QE (6) x read-noise\nmedian (6) x tail (3) x sampling (6) x background (4) x NA (3) x dark current (3) x loss per photon (3),\n{full_points} points, about {full_cpuh:.0f} CPU-h with the SLURM templates; an EMCCD branch (EM gain, CIC,\nphoton-counting frame rate) of similar size; (ii) Sobol indices (SALib, N = 1024, 8 parameters,\nabout {sobol_cpuh:.0f} CPU-h) for t(1e-3) and t(1e-4); (iii) spec-to-fidelity contour maps, the Stage 1\ncamera specification and the Prime BSI Express recipe including the Bayesian classifier with a\nmeasured noise map.\n\n## 8. Decisions needed\n\n1. Approve the full sweep and the Sobol analysis (HPC account and partition details are needed for\n   the SLURM templates, marked [CONFIRM?]); or run the Sobol analysis on this instance first.\n2. Confirm the reference scenario (Sec. 2), in particular whether to use the measured V1 PSF width\n   instead of the diffraction limit, and whether 0.5 is the right filling for the specification.\n3. A dark-frame stack from our Prime BSI Express, when it arrives, would replace the assumed\n   read-noise tail.\n\n## Open placeholders ({n_open})\n\n{placeholders}\n'


if __name__ == "__main__":
    print(build())


# ------------------------------------------------------------------------------- Sobol
def _sobol_from_points(n_rows):
    """Results of the first n_rows Sobol points read from the per-point checkpoints (interim analysis)."""
    import json
    recs = {}
    for f in (RUNS / "prod_sobol" / "points").glob("*.json"):
        r = json.loads(f.read_text())
        if r["index"] < n_rows:
            recs[r["index"]] = r
    if len(recs) != n_rows:
        raise RuntimeError(f"Sobol prefix incomplete: {len(recs)} of {n_rows} points")
    rs = [recs[i] for i in range(n_rows)]
    d = {"index": np.arange(n_rows)}
    for k in rs[0]["result"]:
        d[f"results/{k}"] = np.array([r["result"][k] for r in rs], float)
    for k in rs[0]["point"]:
        if isinstance(rs[0]["point"][k], (int, float)):
            d[f"params/{k}"] = np.array([r["point"][k] for r in rs], float)
    for k in rs[0]["cost"]:
        d[f"cost/{k}"] = np.array([r["cost"][k] for r in rs], float)
    return d


def sobol_results(problem="configs/production/sobol_problem.yaml", n_base=None):
    """First-order and total Sobol indices (SALib, bootstrap 95% CIs) for log10 t(1e-3), log10 t(1e-4).

    With n_base set, analyse the first n_base*(D+2) points only: SALib's Saltelli rows are grouped per base
    point and the scrambled base sequence is seeded, so this prefix is exactly the N = n_base design.
    """
    from .sweep.sobol import analyse, _problem
    if n_base is None:
        d, a = load_h5(RUNS / "prod_sobol" / "prod_sobol.h5")
        o = np.argsort(d["index"])
        n_exp = int(a.get("n_points", len(o)))
        if len(o) != n_exp:
            raise RuntimeError(f"Sobol run incomplete: {len(o)} of {n_exp} points")
    else:
        D = _problem(ROOT / problem)[0]["num_vars"]
        d = _sobol_from_points(n_base * (D + 2))
        o = np.arange(len(d["index"]))
    out = {}
    for tg in ("1e-03", "1e-04"):
        y = d[f"results/log10t_{tg}"][o]
        reached = d[f"results/reached_{tg}"][o]
        out[tg] = {"indices": analyse(ROOT / problem, y), "frac_reached": float(reached.mean()),
                   "y": y, "reached": reached}
    params = {k.split("/", 1)[1]: v[o] for k, v in d.items() if k.startswith("params/")}
    cost = {k.split("/", 1)[1]: v for k, v in d.items() if k.startswith("cost/")}
    return out, params, cost


def sobol_table(res):
    rows = []
    names = [r["param"] for r in res["1e-04"]["indices"]]
    by = {tg: {r["param"]: r for r in res[tg]["indices"]} for tg in res}
    for nm in names:
        r4, r3 = by["1e-04"][nm], by["1e-03"][nm]
        rows.append([nm, f"{r4['S1']:.3f} +- {r4['S1_conf']:.3f}", f"{r4['ST']:.3f} +- {r4['ST_conf']:.3f}",
                     f"{r3['S1']:.3f} +- {r3['S1_conf']:.3f}", f"{r3['ST']:.3f} +- {r3['ST_conf']:.3f}"])
    return md_table(["parameter (ranked by S_T for 1e-4)", "S1, t(1e-4)", "S_T, t(1e-4)", "S1, t(1e-3)", "S_T, t(1e-3)"], rows)


def fig_sobol(res, out):
    names = [r["param"] for r in res["1e-04"]["indices"]]
    by = {r["param"]: r for r in res["1e-04"]["indices"]}
    fig, (ax,) = _fig(w=7.0, h=3.8)
    y = np.arange(len(names))
    s1 = np.array([by[n]["S1"] for n in names]); s1c = np.array([by[n]["S1_conf"] for n in names])
    st = np.array([by[n]["ST"] for n in names]); stc = np.array([by[n]["ST_conf"] for n in names])
    ax.barh(y - 0.18, st, 0.34, xerr=stc, color=C[0], label="total index S_T", error_kw=dict(ecolor=INK2, lw=1))
    ax.barh(y + 0.18, s1, 0.34, xerr=s1c, color=C[1], label="first-order index S1", error_kw=dict(ecolor=INK2, lw=1))
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Sobol index of log10 t(1e-4)", color=INK, fontsize=9)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)

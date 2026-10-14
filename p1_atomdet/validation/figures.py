"""Figures for the P1 validation report (static PNG, matplotlib).

Palette: categorical slots 1-4 of the reference palette (blue, orange, aqua, yellow) in fixed
order, validated for colour-vision deficiency on a light surface; every series also differs in
line style or marker and is named in a legend, and all plotted numbers appear in report tables.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from camsim_common.io import load_h5  # noqa: E402

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def _style(ax):
    ax.set_facecolor("#fcfcfb")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def _fig(w=6.4, h=3.6, n=1):
    fig, axs = plt.subplots(1, n, figsize=(w, h), facecolor="#fcfcfb")
    axs = np.atleast_1d(axs)
    for a in axs:
        _style(a)
    return fig, axs


def fig_v1_histogram(run_dir: Path, out: Path):
    d, _ = load_h5(run_dir / "v1_manetsch.h5")
    fig, (ax,) = _fig()
    ax.plot(d["inset_x"], d["inset_logd_data"], color=C[0], lw=2, label="Manetsch et al. Ext. Data Fig. 4b inset (digitised)")
    ax.plot(d["inset_x"], d["inset_logd_model"], color=C[1], lw=2, ls="--", label="model, box sum (baseline)")
    ax.set_xlabel("photons in 7x7 box (unweighted)", color=INK, fontsize=9)
    ax.set_ylabel("log10 probability per photon", color=INK, fontsize=9)
    ax.set_ylim(-6.2, -1.0)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def fig_v1_errors(run_dir: Path, out: Path):
    d, _ = load_h5(run_dir / "v1_manetsch.h5")
    res = json.loads((run_dir / "v1_manetsch.json").read_text())
    x0, p0, x1, p1 = d["p0_x"], d["p0_p"], d["p1_x"], d["p1_p"]
    T = np.arange(25, 60, 0.107)
    fp = np.array([p0[x0 >= t - 1e-9].sum() for t in T])
    fn = np.array([p1[x1 < t - 1e-9].sum() for t in T])
    F = res["fixed"]["filling"]
    fig, (ax,) = _fig()
    ax.semilogy(T, fp, color=C[0], lw=2, label="FP(T), model")
    ax.semilogy(T, fn, color=C[1], lw=2, ls="--", label="FN(T), model")
    ax.semilogy(T, (1 - F) * fp + F * fn, color=C[2], lw=2, ls=":", label="infidelity, model")
    t = res["targets"]
    ax.axhline(t["fp"], color=C[0], lw=1, alpha=0.6)
    ax.axhline(t["fn"], color=C[1], lw=1, alpha=0.6, ls="--")
    ax.text(26, t["fp"] * 1.25, "measured FP 7.01e-5", fontsize=7, color=INK2)
    ax.text(26, t["fn"] * 0.62, "measured FN 5.5e-5", fontsize=7, color=INK2)
    ax.set_ylim(1e-7, 1e-2)
    ax.set_xlabel("threshold T (photoelectrons, 7x7 box)", color=INK, fontsize=9)
    ax.set_ylabel("error probability", color=INK, fontsize=9)
    ax.legend(frameon=False, fontsize=8, loc="upper center")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def fig_v1_loss_shapes(run_dir: Path, out: Path):
    d, _ = load_h5(run_dir / "v1_manetsch.h5")
    res = json.loads((run_dir / "v1_manetsch.json").read_text())
    L = res["fixed"]["L"]
    fig, (ax,) = _fig()
    tt = np.linspace(0, 1, 200)
    ax.plot(tt, L * np.exp(-L * tt) / (1 - np.exp(-L)), color=C[0], lw=2, label="model A (constant hazard)")
    for i, (k, ls, mk) in enumerate((("thermostat", "--", "o"), ("heating", ":", "s"))):
        c, h = d[f"modelB_{k}_centres"], d[f"modelB_{k}_density"]
        n = res["model_b"][k]["n_lost"]
        ax.plot(c, h, color=C[i + 1], lw=2, ls=ls, marker=mk, ms=4, label=f"model B, {k} ({n} lost atoms)")
    ax.set_xlabel("loss time / exposure", color=INK, fontsize=9)
    ax.set_ylabel("probability density of loss time", color=INK, fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def fig_v2_tails(run_dir: Path, out: Path):
    d, _ = load_h5(run_dir / "v2_picken.h5")
    res = json.loads((run_dir / "v2_picken.json").read_text())
    from scipy import stats
    fig, axs = _fig(w=7.6, h=3.4, n=2)
    for ax, t in zip(axs, ("trap1", "trap2")):
        tr = res["traps"][t]
        x, p = d[f"{t}_p0_adu"], d[f"{t}_p0"]
        g = tr["g_e_per_adu"]
        # survival function of the background statistic vs the Gaussian of Table II
        sf = np.cumsum(p[::-1])[::-1]
        mu0 = tr["model_mu0_adu"]
        s0 = tr["model_sigma0_adu"]
        ax.semilogy(x, sf, color=C[0], lw=2, label="Poisson x Gaussian model")
        ax.semilogy(x, stats.norm.sf(x, mu0, s0), color=C[1], lw=2, ls="--", label="Gaussian (Table II)")
        ax.axvline(tr["cmin_adu"], color=INK2, lw=1)
        ax.text(tr["cmin_adu"] + 3, 1e-2, "c_min", fontsize=7, color=INK2)
        ax.set_xlim(150, 420)
        ax.set_ylim(1e-10, 1)
        ax.set_xlabel(f"{t.replace('trap', 'trap ')}: counts in 3x3 box (ADU)", color=INK, fontsize=9)
        ax.set_ylabel("P(S >= c | no atom)", color=INK, fontsize=9)
    axs[0].legend(frameon=False, fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def make_all(run_dir: Path, fig_dir: Path) -> dict:
    fig_dir.mkdir(parents=True, exist_ok=True)
    out = {}
    for name, fn in (("v1_histogram", fig_v1_histogram), ("v1_errors", fig_v1_errors),
                     ("v1_loss_shapes", fig_v1_loss_shapes), ("v2_tails", fig_v2_tails)):
        p = fig_dir / f"{name}.png"
        fn(run_dir, p)
        out[name] = p
    return out

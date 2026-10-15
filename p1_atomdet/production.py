"""Production-stage evaluation: imaging time to a target infidelity for a camera in a reference scenario.

Reference scenario (configs/production/reference_cs.yaml): Cs D2 imaging anchored on V1, with
scattering rate, optical transmission, loss per photon and rate spread taken from the V1 fit
(DERIVED), and the objective NA as a parameter. A camera enters through QE, read-noise
distribution, gain stage (EMCCD), dark current and pixel sampling (pixels per PSF FWHM).

Error rates are averaged over the sub-pixel position of the site (uniform), with per-site
classifier weights and thresholds; FP and FN are reported at the optimal threshold of the
averaged error. The CF engine (counts/cf.py) is used throughout; it omits ADC rounding, PRNU and
serial CIC, which the lattice engine covers for the camera presets.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .atoms.loss_phenom import TauDistribution
from .constants import SPECIES
from .counts.cf import CFNoise, box_mask_for_snr, cf_pair, matched_weights_analytic
from .counts.pmf import SiteModel
from .fidelity import ErrorRates, error_curves, time_to_target
from .optics.array import PSFSpec, pitch_for_pixels_per_fwhm, pixel_fractions
from .optics.collection import collection_fraction


@dataclass
class Scenario:
    species: str = "Cs133"
    rate_s: float = 72300.0           # scattering rate (V1, DERIVED)
    na: float = 0.55
    eta_opt: float = 0.327            # transmission after collection, V1: 0.027 / eta_coll(0.55) (DERIVED)
    dipole_pattern: str = "isotropic"
    p_loss_per_photon: float = 2.35e-8   # V1: -ln(0.999864)/(R t) (DERIVED; includes vacuum)
    tau_vac_s: float = np.inf
    rate_rel_sd: float = 0.1226       # V1 Fig. 4d (DATA)
    bg_rate_e_per_px_s: float = 4.76  # V1: 18.66 e- / 49 px / 80 ms at the V1 sampling (DATA); scale with sampling
    filling: float = 0.5
    psf_blur_px: float = 0.0
    psf_kind: str = "vectorial"       # vectorial (diffraction limit at na) | gaussian (measured waist)
    psf_waist_m: float | None = None  # object-plane 1/e^2 radius for psf_kind == gaussian, quoted at na_ref
    psf_na_ref: float = 0.55          # gaussian waist scales as na_ref/na when na differs (assumption)
    n_offsets: int = 3                # offsets per axis on [0, 0.5] (symmetry-reduced)
    offset_rule: str = "midpoint"     # midpoint (cell centres; default) | trapezoid (includes 0 and 0.5)
    ds: float = 0.02
    n_tau: int = 24


@dataclass
class Camera:
    qe: float
    read_classes: list                # [(class weight, [(w, mean, sigma)])]
    pixels_per_fwhm: float
    dark_e_per_px_s: float = 0.0
    gain_kind: str = "none"           # none | gamma
    f2: float = 2.0
    cic_per_px: float = 0.0           # EMCCD clock-induced charge (treated as extra input electrons)
    em_gain: float = 1.0
    label: str = ""


@dataclass
class Geometry:
    fracs: list
    weights_off: np.ndarray
    box: int
    pitch_m: float


def geometry(sc: Scenario, cam: Camera, supersample: int = 6) -> Geometry:
    sp = SPECIES[sc.species]
    if sc.psf_kind == "gaussian":
        psf = PSFSpec("gaussian", wavelength_m=float(sp.wavelength_m), na=sc.na,
                      waist_m=sc.psf_waist_m * sc.psf_na_ref / sc.na)
    else:
        psf = PSFSpec("vectorial", wavelength_m=float(sp.wavelength_m), na=sc.na, n_quad=48)
    pitch = pitch_for_pixels_per_fwhm(psf, cam.pixels_per_fwhm)
    psf.blur_sigma_m = sc.psf_blur_px * pitch
    box = 2 * int(np.ceil(2.5 * cam.pixels_per_fwhm)) + 1
    box = max(box, 3)
    if sc.offset_rule == "midpoint":
        g = (np.arange(sc.n_offsets) + 0.5) * 0.5 / sc.n_offsets
        wt1 = np.ones(sc.n_offsets)
    else:
        g = np.linspace(0.0, 0.5, sc.n_offsets)
        wt1 = np.ones(sc.n_offsets)
        if sc.n_offsets > 1:
            wt1[0] = wt1[-1] = 0.5
    wt = np.outer(wt1, wt1).ravel()
    wt /= wt.sum()
    fracs = [pixel_fractions(psf, pitch, box, (dx, dy), supersample) for dx in g for dy in g]
    return Geometry(fracs, wt, box, pitch)


def mu_per_s(sc: Scenario, cam: Camera) -> float:
    return sc.rate_s * float(collection_fraction(sc.na, sc.dipole_pattern)) * sc.eta_opt * cam.qe


def _bg_rate(sc: Scenario, cam: Camera, geo: Geometry) -> float:
    """Background per pixel scales with pixel area in the object plane; reference = V1 pixel (0.902 um)."""
    return sc.bg_rate_e_per_px_s * (geo.pitch_m / 0.902e-6) ** 2 * cam.qe / 0.44


def eps_at(sc: Scenario, cam: Camera, geo: Geometry, t: float, classifier: str = "matched") -> ErrorRates:
    noise = CFNoise.from_classes(cam.read_classes, gain_kind=cam.gain_kind, f2=cam.f2)
    mu = mu_per_s(sc, cam) * t
    bg = _bg_rate(sc, cam, geo) * t + cam.dark_e_per_px_s * t + cam.cic_per_px
    L = sc.p_loss_per_photon * sc.rate_s * t + (t / sc.tau_vac_s if np.isfinite(sc.tau_vac_s) else 0.0)
    tau = TauDistribution.from_hazard(L, sc.n_tau)
    fp_curves, fn_curves, T_ref = [], [], None
    for fr, a in zip(geo.fracs, geo.weights_off):
        site = SiteModel(frac=fr, mu=mu, bg=np.full(fr.shape, bg), tau=tau, rate_rel_sd=sc.rate_rel_sd, ds=sc.ds)
        w = matched_weights_analytic(site, noise) if classifier == "matched" else box_mask_for_snr(site, noise)
        p0, p1 = cf_pair(site, noise, w)
        # per-site optimal threshold (sites may differ in sub-pixel offset)
        T, fp, fn = error_curves(p0, p1)
        e = (1 - sc.filling) * fp + sc.filling * fn
        i = int(np.argmin(e))
        fp_curves.append(a * fp[i])
        fn_curves.append(a * fn[i])
    fp, fn = float(sum(fp_curves)), float(sum(fn_curves))
    return ErrorRates(float("nan"), fp, fn, (1 - sc.filling) * fp + sc.filling * fn, sc.filling)


def survival(sc: Scenario, t: float) -> float:
    return float(np.exp(-sc.p_loss_per_photon * sc.rate_s * t - (t / sc.tau_vac_s if np.isfinite(sc.tau_vac_s) else 0.0)))


def time_to(sc: Scenario, cam: Camera, target: float, classifier: str = "matched", t_lo=1e-3, t_hi=0.2,
            min_survival: float = 0.99, geo: Geometry | None = None):
    geo = geo or geometry(sc, cam)
    t, info = time_to_target(lambda tt: max(eps_at(sc, cam, geo, tt, classifier).eps, 1e-300), target, t_lo, t_hi,
                             survival_of_t=lambda tt: survival(sc, tt), min_survival=min_survival, n_scan=12)
    return t, info


def read_classes(median_e: float, tail_fraction: float = 0.0, tail_factor: float = 5.0):
    """Gaussian core plus a fraction of 'tail' pixels with rms tail_factor x median (P4 will replace this)."""
    cl = [(1.0 - tail_fraction, [(1.0, 0.0, median_e)])]
    if tail_fraction > 0:
        cl.append((tail_fraction, [(1.0, 0.0, tail_factor * median_e)]))
    return cl


# ----------------------------------------------------------------- EMCCD photon counting
def pc_click_prob(lam, em_gain: float, n_stages: int, f2: float, sigma_ro_e: float, thr_sigma: float,
                  cic_register: float = 0.0, n_gh: int = 40):
    """P(click) per pixel for Poisson(lam) input electrons (vectorised over lam).

    Register output (input-referred) for n electrons: Gamma(n a, theta), a = 1/(F^2 - 1), theta = F^2 - 1;
    Gaussian read noise sigma_in = sigma_ro/g; threshold tau = thr_sigma * sigma_in. Register-generated
    CIC clicks with probability c_reg * mean_j exp(-tau m^j) (exponential output of mean m^-j, stage j).
    """
    from scipy import special, stats
    lam = np.atleast_1d(np.asarray(lam, float))
    a, th = 1.0 / (f2 - 1.0), f2 - 1.0
    s_in = sigma_ro_e / em_gain
    tau = thr_sigma * s_in
    x, w = np.polynomial.hermite_e.hermegauss(n_gh)
    w = w / w.sum()
    nmax = int(lam.max() + 12 * np.sqrt(lam.max() + 1) + 20)
    n = np.arange(nmax + 1)
    thr = np.clip(tau - s_in * x, 0, None)  # (n_gh,)
    Pn = np.empty(n.size)
    Pn[0] = stats.norm.sf(thr_sigma)
    Pn[1:] = (special.gammaincc((n[1:] * a)[:, None], thr[None, :] / th) * w[None, :]).sum(1)
    p = stats.poisson.pmf(n[None, :], lam[:, None]) @ Pn
    if cic_register > 0:
        m = em_gain ** (1.0 / n_stages)
        j = np.arange(1, n_stages + 1)
        p_c = cic_register * np.mean(np.exp(-tau * m ** j.astype(float)))
        p = 1 - (1 - p) * (1 - p_c)
    return p


def pc_pair(site: SiteModel, p_of_lam, mask, n_frames: int = 1):
    """Box-mask click-count distributions (P0, P1) on the integer lattice.

    Multi-frame photon counting: the exposure is split into n_frames frames; each pixel's count is
    Binomial(n_frames, p(lam / n_frames)) given the signal multiplier s (loss timing across frames
    is not resolved). Pixel counts are convolved exactly; s is mixed as elsewhere.
    """
    from scipy import stats
    from .counts.lattice import Lattice
    f, b = site.frac.ravel(), site.bg.ravel()
    idx = np.nonzero(np.asarray(mask).ravel() > 0)[0]
    k = np.arange(n_frames + 1)

    def dist(lam):
        ps = p_of_lam(lam / n_frames)
        d = np.array([1.0])
        for q in ps:
            d = np.convolve(d, stats.binom.pmf(k, n_frames, q))
        return d

    p0 = dist(b[idx])
    sv, sw = site.signal_multiplier_nodes()
    p1 = np.zeros(p0.size)
    for s, a in zip(sv, sw):
        p1 += a * dist(b[idx] + site.mu * s * f[idx])
    return Lattice(p0, 0.0, 1.0), Lattice(p1, 0.0, 1.0)


def eps_at_pc(sc: Scenario, cam: "Camera", geo: Geometry, t: float, pc: dict) -> ErrorRates:
    pc = dict(pc)
    fps = pc.pop("frame_rate_hz", None)
    n_frames = max(1, int(np.floor(t * fps))) if fps else 1
    mu = mu_per_s(sc, cam) * t
    bg = _bg_rate(sc, cam, geo) * t + cam.dark_e_per_px_s * t + cam.cic_per_px * n_frames
    L = sc.p_loss_per_photon * sc.rate_s * t + (t / sc.tau_vac_s if np.isfinite(sc.tau_vac_s) else 0.0)
    tau = TauDistribution.from_hazard(L, sc.n_tau)
    fp_t, fn_t = 0.0, 0.0
    p_of_lam = lambda lam: pc_click_prob(lam, **pc)
    for fr, a in zip(geo.fracs, geo.weights_off):
        site = SiteModel(frac=fr, mu=mu, bg=np.full(fr.shape, bg), tau=tau, rate_rel_sd=sc.rate_rel_sd, ds=sc.ds)
        # mask: pixels by signal fraction, prefix maximising Gaussian d' with Bernoulli variances
        o = np.argsort(-fr.ravel())
        p0 = p_of_lam(np.full(o.size, bg / n_frames))
        p1 = p_of_lam((bg + mu * fr.ravel()[o]) / n_frames)
        dsum = np.cumsum(p1 - p0) * n_frames
        v = n_frames * np.cumsum(p0 * (1 - p0) + p1 * (1 - p1)) / 2
        n = int(np.argmax(dsum / np.sqrt(v))) + 1
        mask = np.zeros(fr.size)
        mask[o[:n]] = 1
        P0, P1 = pc_pair(site, p_of_lam, mask, n_frames)
        T, fp, fn = error_curves(P0, P1)
        e = (1 - sc.filling) * fp + sc.filling * fn
        i = int(np.argmin(e))
        fp_t += a * fp[i]
        fn_t += a * fn[i]
    return ErrorRates(float("nan"), fp_t, fn_t, (1 - sc.filling) * fp_t + sc.filling * fn_t, sc.filling)


# ----------------------------------------------------------------- presets -> Camera
def camera_from_sensor(sm, pixels_per_fwhm: float, species: str = "Cs133", cic_register_fraction=None,
                       pc_frame_rate_hz: float | None = None):
    """Map a SensorModel to the production Camera (plus PC parameters for photon-counting EMCCDs)."""
    from camsim_common.provenance import val
    from .sensor.emccd import enf2_finite
    qe = float(val(sm.qe_852) if species == "Cs133" else val(sm.qe_780))
    D = sm.dark()
    if sm.emccd is None:
        return Camera(qe=qe, read_classes=sm.read_noise.pixel_classes(), pixels_per_fwhm=pixels_per_fwhm,
                      dark_e_per_px_s=D, label=sm.name), None
    e = sm.emccd
    g, r = float(val(e.em_gain)), int(val(e.n_stages))
    fr = float(val(e.cic_register_fraction)) if cic_register_fraction is None else cic_register_fraction
    f2 = enf2_finite(g, r) if e.register_model == "stagewise" else 2.0
    s_ro = float(val(sm.read_noise.sigma_e))
    cic = float(val(e.cic))
    cam = Camera(qe=qe, read_classes=[(1.0, [(1.0, 0.0, s_ro / g)])], pixels_per_fwhm=pixels_per_fwhm,
                 dark_e_per_px_s=D, gain_kind="gamma", f2=f2, cic_per_px=cic * (1 - fr), em_gain=g, label=sm.name)
    pc = None
    if e.mode == "photon_counting":
        pc = dict(em_gain=g, n_stages=r, f2=f2, sigma_ro_e=s_ro, thr_sigma=e.pc_threshold_sigma, cic_register=cic * fr,
                  frame_rate_hz=pc_frame_rate_hz)
    return cam, pc


def eps_any(sc, cam, geo, t, pc=None, classifier="matched"):
    return eps_at_pc(sc, cam, geo, t, pc) if pc is not None else eps_at(sc, cam, geo, t, classifier)


def time_to_any(sc, cam, target, pc=None, classifier="matched", geo=None, t_lo=1e-3, t_hi=0.3, min_survival=0.99):
    geo = geo or geometry(sc, cam)
    return time_to_target(lambda tt: max(eps_any(sc, cam, geo, tt, pc, classifier).eps, 1e-300), target, t_lo, t_hi,
                          survival_of_t=lambda tt: survival(sc, tt), min_survival=min_survival, n_scan=12)


def times_to_targets(eps_of_t, targets, t_lo=1e-3, t_hi=1.0, n_scan=13, survival_of_t=None, min_survival=0.99):
    """Shortest t with eps(t) <= target for several targets from one shared log scan, refined by Brent.

    Returns {target: t or None} and the scan. eps may rise at long t (loss), so the first crossing is used.
    """
    from scipy import optimize
    ts = np.geomspace(t_lo, t_hi, n_scan)
    eps = np.array([eps_of_t(t) for t in ts])
    ok_s = np.ones(ts.size, bool) if survival_of_t is None else np.array([survival_of_t(t) >= min_survival for t in ts])
    out = {}
    for target in targets:
        ok = (eps <= target) & ok_s
        if not ok.any():
            out[target] = None
            continue
        j = int(np.argmax(ok))
        if j == 0:
            out[target] = float(ts[0])
            continue
        f = lambda lt: np.log(max(eps_of_t(np.exp(lt)), 1e-300)) - np.log(target)
        out[target] = float(np.exp(optimize.brentq(f, np.log(ts[j - 1]), np.log(ts[j]), xtol=2e-3)))
    return out, (ts, eps)

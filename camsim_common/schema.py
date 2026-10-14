"""SensorModel: the shared camera description used by P1 to P4.

One SensorModel describes one camera in one operating mode (e.g. Prime BSI Express CMS).
Every numeric field is a ProvValue, so provenance travels with the file.

Producers:  P2/P3 fill `qe` and `mtf`; P4 fills `read_noise`, `shutter.line_time` and `dark`.
Consumer:   P1 reads complete files and returns fidelity and imaging time.
"""
from __future__ import annotations

from pathlib import Path
from typing import Annotated, Literal, Optional, Union

import numpy as np
import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator

from .provenance import ProvValue, val

K_B_EV = 8.617333262e-5  # Boltzmann constant in eV/K (CODATA 2018 exact-derived; DERIVED from SI definition)


class _Base(BaseModel):
    model_config = ConfigDict(extra="forbid")


# --------------------------------------------------------------------------- QE
class QETable(_Base):
    """QE(lambda, T). values has shape (len(temperature_K), len(wavelength_nm))."""

    wavelength_nm: list[float]
    temperature_K: list[float] = Field(default_factory=lambda: [293.15])
    values: list[list[float]]
    prov: ProvValue  # provenance of the table as a whole (value field holds a short label)

    @model_validator(mode="after")
    def _shape(self):
        a = np.asarray(self.values, float)
        if a.shape != (len(self.temperature_K), len(self.wavelength_nm)):
            raise ValueError(f"QE table shape {a.shape} does not match axes")
        if np.any(a < 0) or np.any(a > 1):
            raise ValueError("QE values must lie in [0, 1]")
        return self

    def __call__(self, wavelength_nm: float, temperature_K: float | None = None) -> float:
        a = np.asarray(self.values, float)
        wl = np.asarray(self.wavelength_nm, float)
        per_T = np.array([np.interp(wavelength_nm, wl, row) for row in a])
        T = np.asarray(self.temperature_K, float)
        if temperature_K is None or len(T) == 1:
            return float(per_T[0])
        return float(np.interp(temperature_K, T, per_T))


# ---------------------------------------------------------------- read noise
class GaussianCore(_Base):
    """All pixels Gaussian with the same sigma (e- rms per read)."""

    kind: Literal["gaussian"] = "gaussian"
    sigma_e: ProvValue

    def pixel_classes(self):
        return [(1.0, [(1.0, 0.0, float(val(self.sigma_e)))])]


class RTSAmplitude(_Base):
    amplitude_e: ProvValue
    weight: float = 1.0


class GaussianRTSMixture(_Base):
    """Gaussian core plus a fraction f_tail of random-telegraph-signal (RTS) pixels.

    A tail pixel, in each read, shows the core Gaussian displaced by +A or -A with
    probability p_jump/2 each (CDS of a two-state trap), else by 0. The amplitude A of a
    given tail pixel is drawn from `amplitudes` with the listed weights.
    """

    kind: Literal["gaussian_rts"] = "gaussian_rts"
    sigma_core_e: ProvValue
    f_tail: ProvValue
    p_jump: ProvValue
    amplitudes: list[RTSAmplitude]

    def pixel_classes(self):
        s = float(val(self.sigma_core_e))
        f = float(val(self.f_tail))
        p = float(val(self.p_jump))
        wsum = sum(a.weight for a in self.amplitudes)
        out = [(1.0 - f, [(1.0, 0.0, s)])]
        for a in self.amplitudes:
            A = float(val(a.amplitude_e))
            out.append((f * a.weight / wsum, [(1 - p, 0.0, s), (p / 2, A, s), (p / 2, -A, s)]))
        return out


class EmpiricalSigmaCDF(_Base):
    """Per-pixel rms read noise given as an empirical CDF (e.g. digitised vendor histogram).
    Each pixel is Gaussian with its own sigma."""

    kind: Literal["empirical_sigma_cdf"] = "empirical_sigma_cdf"
    sigma_e: list[float]
    cdf: list[float]
    prov: ProvValue

    @model_validator(mode="after")
    def _mono(self):
        c = np.asarray(self.cdf)
        if len(c) != len(self.sigma_e) or np.any(np.diff(c) < 0) or abs(c[-1] - 1) > 1e-9:
            raise ValueError("cdf must be non-decreasing, end at 1, and match sigma_e")
        return self

    def pixel_classes(self):
        c = np.concatenate([[0.0], np.asarray(self.cdf)])
        w = np.diff(c)
        return [(float(wi), [(1.0, 0.0, float(si))]) for wi, si in zip(w, self.sigma_e) if wi > 0]


ReadNoise = Annotated[Union[GaussianCore, GaussianRTSMixture, EmpiricalSigmaCDF], Field(discriminator="kind")]


def per_pixel_rms(classes) -> tuple[np.ndarray, np.ndarray]:
    """Return (weights, rms) of the per-pixel rms distribution."""
    w, r = [], []
    for cw, comps in classes:
        m2 = sum(c * (mu**2 + s**2) for c, mu, s in comps)
        w.append(cw)
        r.append(np.sqrt(m2))
    return np.asarray(w), np.asarray(r)


def read_noise_median(rn) -> float:
    """Median over pixels of the per-pixel rms (the number vendors usually quote as 'median')."""
    w, r = per_pixel_rms(rn.pixel_classes())
    o = np.argsort(r)
    cw = np.cumsum(w[o])
    return float(r[o][np.searchsorted(cw, 0.5 * cw[-1])])


def read_noise_rms(rn) -> float:
    """Quadrature mean over pixels of the per-pixel rms (the vendor 'rms' figure)."""
    w, r = per_pixel_rms(rn.pixel_classes())
    return float(np.sqrt(np.sum(w * r**2) / np.sum(w)))


# --------------------------------------------------------------------- dark
class DarkArrhenius(_Base):
    """D(T) = D_ref * (T/T_ref)^1.5 * exp(-Ea/kB * (1/T - 1/T_ref))   [e-/pixel/s]."""

    D_ref: ProvValue
    T_ref_K: ProvValue
    Ea_eV: ProvValue
    T_op_K: ProvValue

    def __call__(self, T_K: float | None = None) -> float:
        T = float(val(self.T_op_K)) if T_K is None else T_K
        Tr = float(val(self.T_ref_K))
        return float(val(self.D_ref)) * (T / Tr) ** 1.5 * np.exp(-float(val(self.Ea_eV)) / K_B_EV * (1 / T - 1 / Tr))


# -------------------------------------------------------------------- EMCCD
class EMCCDBlock(_Base):
    em_gain: ProvValue
    cic: ProvValue  # events / pixel / frame
    n_stages: ProvValue
    cic_register_fraction: ProvValue  # fraction of CIC generated inside the multiplication register
    register_model: Literal["erlang", "stagewise"] = "erlang"
    mode: Literal["analog", "photon_counting"] = "analog"
    pc_threshold_sigma: float = 5.0  # photon-counting threshold in units of read-noise sigma (output referred)


# ---------------------------------------------------------------- geometry
class MTFBlock(_Base):
    """Either a tabulated MTF(f) or a Gaussian charge-spread kernel (sigma, um)."""

    freq_cyc_per_mm: Optional[list[float]] = None
    mtf: Optional[list[float]] = None
    kernel_sigma_um: Optional[ProvValue] = None
    prov: Optional[ProvValue] = None


class ShutterBlock(_Base):
    type: Literal["global", "rolling"]
    line_time_us: ProvValue
    frame_time_ms: ProvValue
    n_rows: int = 2048


class CosmicBlock(_Base):
    rate_per_cm2_s: ProvValue
    deposit_e_per_um: ProvValue
    thickness_um: ProvValue


class ADCBlock(_Base):
    gain_e_per_adu: ProvValue
    adc_bits: ProvValue
    full_well_e: ProvValue
    offset_adu: ProvValue


# ---------------------------------------------------------------- top level
class SensorModel(_Base):
    schema_version: Literal["1.0"] = "1.0"
    name: str
    camera: str
    mode: str
    sensor_type: Literal["sCMOS", "qCMOS", "EMCCD", "CCD", "CMOS"]
    notes: str = ""

    qe: QETable
    qe_780: ProvValue
    qe_852: ProvValue
    prnu: ProvValue  # fractional rms
    read_noise: ReadNoise
    read_noise_median_e: ProvValue  # stored, checked against the distribution
    read_noise_rms_e: ProvValue
    adc: ADCBlock
    dark: DarkArrhenius
    emccd: Optional[EMCCDBlock] = None
    pixel_pitch_um: ProvValue
    mtf: Optional[MTFBlock] = None
    shutter: ShutterBlock
    cosmic: CosmicBlock
    qcmos_pnr: bool = False  # photon-number-resolving output (rounded to integer electrons)

    @model_validator(mode="after")
    def _consistency(self):
        if self.sensor_type == "EMCCD" and self.emccd is None:
            raise ValueError("EMCCD sensor requires an emccd block")
        return self

    # --- derived checks
    def read_noise_check(self, rtol: float = 0.05) -> dict:
        """Compare the stored median/rms with those implied by the distribution."""
        med, rms = read_noise_median(self.read_noise), read_noise_rms(self.read_noise)
        sm, sr = float(val(self.read_noise_median_e)), float(val(self.read_noise_rms_e))
        return {
            "median_model": med, "median_stored": sm, "median_ok": abs(med - sm) <= rtol * max(sm, 1e-12),
            "rms_model": rms, "rms_stored": sr, "rms_ok": abs(rms - sr) <= rtol * max(sr, 1e-12),
        }

    # --- IO
    def to_yaml(self, path: str | Path | None = None) -> str:
        text = yaml.safe_dump(self.model_dump(mode="json", exclude_none=True), sort_keys=False, width=110)
        if path is not None:
            Path(path).write_text(text)
        return text

    @classmethod
    def from_yaml(cls, path_or_text: str | Path) -> "SensorModel":
        p = Path(str(path_or_text))
        text = p.read_text() if (len(str(path_or_text)) < 4096 and p.exists()) else str(path_or_text)
        return cls.model_validate(yaml.safe_load(text))

    def walk_values(self):
        """Yield (dotted_path, ProvValue) for every ProvValue in the model."""
        yield from _walk(self, self.name)


def _walk(obj, prefix):
    if isinstance(obj, ProvValue):
        yield prefix, obj
    elif isinstance(obj, BaseModel):
        for k in type(obj).model_fields:
            yield from _walk(getattr(obj, k), f"{prefix}.{k}")
    elif isinstance(obj, (list, tuple)):
        for i, x in enumerate(obj):
            yield from _walk(x, f"{prefix}[{i}]")

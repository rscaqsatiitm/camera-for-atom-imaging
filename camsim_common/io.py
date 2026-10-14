"""HDF5 output with the resolved configuration and git state embedded as attributes."""
from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import yaml


def git_state(repo: str | Path | None = None) -> dict:
    cwd = str(repo or Path(__file__).resolve().parents[1])
    def run(*a):
        try:
            return subprocess.run(["git", *a], cwd=cwd, capture_output=True, text=True, timeout=10).stdout.strip()
        except Exception:
            return ""
    return {"git_hash": run("rev-parse", "HEAD") or "uncommitted", "git_dirty": bool(run("status", "--porcelain"))}


def load_config(path: str | Path) -> dict:
    return yaml.safe_load(Path(path).read_text())


def _to_attr(v: Any):
    if isinstance(v, (dict, list, tuple)):
        return json.dumps(v, default=str)
    if v is None:
        return "None"
    return v


def save_h5(path: str | Path, datasets: dict[str, Any], config: dict, attrs: dict | None = None,
            include_git: bool = True) -> Path:
    """Write datasets (nested dicts become groups). The config is stored as YAML text.

    Written atomically (tmp file then rename), so an interrupted run never leaves a
    truncated file that a resumed sweep would mistake for a finished point. Timestamps are
    NOT written, so identical inputs give byte-comparable datasets.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with h5py.File(tmp, "w") as f:
        f.attrs["config_yaml"] = yaml.safe_dump(config, sort_keys=True)
        if include_git:
            for k, v in git_state().items():
                f.attrs[k] = v
        for k, v in (attrs or {}).items():
            f.attrs[k] = _to_attr(v)
        _write_group(f, datasets)
    tmp.replace(path)
    return path


def _write_group(g, d: dict):
    for k, v in d.items():
        if isinstance(v, dict):
            _write_group(g.create_group(k), v)
        else:
            g.create_dataset(k, data=np.asarray(v))


def load_h5(path: str | Path) -> tuple[dict, dict]:
    """Return (datasets, attrs). attrs['config'] holds the parsed config."""
    out: dict = {}
    with h5py.File(path, "r") as f:
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                out[name] = obj[()]
        f.visititems(visit)
        attrs = {k: f.attrs[k] for k in f.attrs}
    attrs["config"] = yaml.safe_load(attrs.get("config_yaml", "{}"))
    return out, attrs


def datasets_equal(a: str | Path, b: str | Path) -> bool:
    """True if two HDF5 files hold identical datasets and identical config (reproducibility check)."""
    da, aa = load_h5(a)
    db, ab = load_h5(b)
    if aa["config"] != ab["config"] or set(da) != set(db):
        return False
    return all(np.array_equal(da[k], db[k]) for k in da)

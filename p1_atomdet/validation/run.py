"""Run validation stages from YAML configs.

    python -m p1_atomdet.validation.run --config configs/validation/v1_manetsch.yaml [more.yaml] [--out runs/validation]

Each config produces <out>/<name>.h5 (arrays + embedded config + git hash) and <out>/<name>.json
(scalar results and markdown tables), and appends timing/memory to <out>/runlog.jsonl.
"""
from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path

import numpy as np

from camsim_common.io import load_config, save_h5
from camsim_common.runlog import runlog

from .common import ROOT

MODULES = {"v1_manetsch": "v1_manetsch", "v2_picken": "v2_picken", "v3_analytic": "v3_analytic",
           "convergence": "convergence"}


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return o


def run_one(cfg_path: Path, out: Path) -> dict:
    cfg = load_config(cfg_path)
    mod = importlib.import_module(f"p1_atomdet.validation.{MODULES[cfg['name']]}")
    with runlog(cfg["name"], out / "runlog.jsonl", {"config": str(cfg_path)}) as rec:
        res, arrays = mod.run(cfg, out)
        res["tables"] = mod.tables(res)
    res["runlog"] = {k: rec[k] for k in ("wall_s", "cpu_s", "peak_rss_mb", "n_cpu")}
    save_h5(out / f"{cfg['name']}.h5", arrays, cfg)
    (out / f"{cfg['name']}.json").write_text(json.dumps(_jsonable(res), indent=1))
    return res


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", nargs="+", required=True)
    ap.add_argument("--out", default=str(ROOT / "runs" / "validation"))
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for c in a.config:
        res = run_one(Path(c), out)
        print(f"== {c}: wall {res['runlog']['wall_s']:.1f} s, peak RSS {res['runlog']['peak_rss_mb']:.0f} MB")
        for name, tab in res["tables"].items():
            print(f"-- {name}\n{tab}\n")


if __name__ == "__main__":
    main()

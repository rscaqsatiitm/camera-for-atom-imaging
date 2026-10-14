"""Resumable parameter sweeps with per-point checkpoints and CPU-hour accounting.

    python -m p1_atomdet.sweep.grid --config configs/sweep/<name>.yaml [--workers 4] [--shard i/n]
    python -m p1_atomdet.sweep.grid --config ... --collect      # merge finished points into one HDF5

Config keys:
  name, seed, evaluator ("module:function"), fixed {..}, axes {param: [values]}  (full factorial)
  or points_file (a .npy/.csv of explicit points, e.g. a Saltelli sample) with columns = axes names.
Each point is written atomically to <out>/points/<index>.json together with wall time, CPU
seconds and peak RSS; points whose file exists are skipped, so a killed run resumes where it
stopped. --shard i/n evaluates every n-th point starting at i (SLURM array tasks).
"""
from __future__ import annotations

import argparse
import importlib
import itertools
import json
import os
from multiprocessing import get_context
from pathlib import Path

import numpy as np

from camsim_common.io import load_config, save_h5
from camsim_common.rng import numba_seed
from camsim_common.runlog import runlog

ROOT = Path(__file__).resolve().parents[2]


def grid_points(cfg: dict) -> list[dict]:
    if "points_file" in cfg:
        arr = np.loadtxt(ROOT / cfg["points_file"], delimiter=",", ndmin=2)
        names = list(cfg["axes"].keys())
        return [dict(zip(names, map(float, row))) for row in arr]
    names = list(cfg["axes"].keys())
    vals = [cfg["axes"][k] for k in names]
    return [dict(zip(names, combo)) for combo in itertools.product(*vals)]


def _resolve(path: str):
    mod, fn = path.split(":")
    return getattr(importlib.import_module(mod), fn)


def _point_file(out: Path, i: int) -> Path:
    return out / "points" / f"{i:07d}.json"


def _work(args):
    i, point, cfg, out = args
    f = _point_file(out, i)
    if f.exists():
        return i, "skipped"
    ev = _resolve(cfg["evaluator"])
    seed = numba_seed(cfg["seed"], cfg["name"], i)
    with runlog(f"{cfg['name']}[{i}]") as rec:
        result = ev({**cfg.get("fixed", {}), **point}, seed)
    payload = {"index": i, "point": point, "result": result, "seed": seed,
               "cost": {k: rec[k] for k in ("wall_s", "cpu_s", "cpu_h", "peak_rss_mb")}}
    tmp = f.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, default=float))
    os.replace(tmp, f)
    return i, "done"


def run(cfg: dict, out: Path, workers: int = 1, shard: tuple[int, int] = (0, 1), limit: int | None = None):
    (out / "points").mkdir(parents=True, exist_ok=True)
    (out / "config.json").write_text(json.dumps(cfg, default=str, indent=1))
    pts = grid_points(cfg)
    todo = [(i, p, cfg, out) for i, p in enumerate(pts) if i % shard[1] == shard[0] and not _point_file(out, i).exists()]
    if limit is not None:
        todo = todo[:limit]
    if workers <= 1:
        status = [_work(t) for t in todo]
    else:
        with get_context("spawn").Pool(workers) as pool:
            status = list(pool.imap_unordered(_work, todo, chunksize=1))
    return {"n_points": len(pts), "evaluated": sum(s == "done" for _, s in status),
            "remaining": sum(not _point_file(out, i).exists() for i in range(len(pts)))}


def collect(cfg: dict, out: Path) -> Path:
    pts = grid_points(cfg)
    recs = []
    for i in range(len(pts)):
        f = _point_file(out, i)
        if f.exists():
            recs.append(json.loads(f.read_text()))
    if not recs:
        raise RuntimeError("no finished points")
    keys = sorted(recs[0]["result"].keys())
    data = {"index": np.array([r["index"] for r in recs]),
            "params": {k: np.array([r["point"][k] for r in recs], float) for k in recs[0]["point"]},
            "results": {k: np.array([r["result"].get(k, np.nan) for r in recs], float) for k in keys},
            "cost": {k: np.array([r["cost"][k] for r in recs], float) for k in recs[0]["cost"]}}
    return save_h5(out / f"{cfg['name']}.h5", data, cfg, attrs={"n_points": len(pts), "n_done": len(recs),
                                                               "cpu_h_total": float(data["cost"]["cpu_h"].sum())})


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--limit", type=int, default=None, help="evaluate at most this many points (testing)")
    ap.add_argument("--collect", action="store_true")
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    out = Path(a.out or ROOT / "runs" / "sweeps" / cfg["name"])
    if a.collect:
        print(collect(cfg, out))
        return
    i, n = map(int, a.shard.split("/"))
    print(run(cfg, out, a.workers, (i, n), a.limit))


if __name__ == "__main__":
    main()

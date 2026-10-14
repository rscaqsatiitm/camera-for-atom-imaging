"""Global sensitivity (Sobol indices) with SALib's Saltelli sampling.

    python -m p1_atomdet.sweep.sobol sample  --problem configs/sweep/sobol_problem.yaml --n 1024 --out points.csv
    (evaluate points.csv with sweep.grid using points_file)
    python -m p1_atomdet.sweep.sobol analyse --problem ... --h5 runs/sweeps/<name>/<name>.h5 --key log10_t_1e3

Parameters with log-uniform ranges are sampled in log10 and converted before evaluation.
Executed only after the validation check-in (production stage).
"""
from __future__ import annotations

import argparse

import numpy as np
import yaml


def _problem(path):
    spec = yaml.safe_load(open(path))
    names = list(spec["parameters"].keys())
    bounds, logs = [], []
    for n in names:
        d = spec["parameters"][n]
        lo, hi = d["range"]
        if d.get("log", False):
            lo, hi = np.log10(lo), np.log10(hi)
        bounds.append([lo, hi])
        logs.append(bool(d.get("log", False)))
    return {"num_vars": len(names), "names": names, "bounds": bounds}, logs


def sample(problem_path, n, out, seed=0):
    from SALib.sample import sobol as sobol_sample
    prob, logs = _problem(problem_path)
    X = sobol_sample.sample(prob, n, calc_second_order=False, seed=seed)
    X[:, logs] = 10 ** X[:, logs]
    np.savetxt(out, X, delimiter=",")
    return X.shape


def analyse(problem_path, y):
    from SALib.analyze import sobol as sobol_an
    prob, _ = _problem(problem_path)
    Si = sobol_an.analyze(prob, np.asarray(y, float), calc_second_order=False)
    order = np.argsort(-Si["ST"])
    return [{"param": prob["names"][i], "S1": Si["S1"][i], "S1_conf": Si["S1_conf"][i], "ST": Si["ST"][i],
             "ST_conf": Si["ST_conf"][i]} for i in order]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["sample", "analyse"])
    ap.add_argument("--problem", required=True)
    ap.add_argument("--n", type=int, default=1024)
    ap.add_argument("--out", default="points.csv")
    ap.add_argument("--h5")
    ap.add_argument("--key")
    a = ap.parse_args(argv)
    if a.cmd == "sample":
        print(sample(a.problem, a.n, a.out))
    else:
        from camsim_common.io import load_h5
        d, _ = load_h5(a.h5)
        o = np.argsort(d["index"])
        for r in analyse(a.problem, d[f"results/{a.key}"][o]):
            print(r)


if __name__ == "__main__":
    main()

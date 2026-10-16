"""Resume test: a sweep interrupted midway restarts and skips completed points."""
import json

from p1_atomdet.sweep.grid import collect, grid_points, run

CFG = {"name": "toy", "seed": 5, "evaluator": "p1_atomdet.sweep.evaluate:toy_point", "fixed": {"b": 1.0},
       "axes": {"a": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]}}


def test_resume_skips_completed_points(tmp_path):
    first = run(CFG, tmp_path, limit=3)            # simulate a run killed after 3 points
    assert first["evaluated"] == 3 and first["remaining"] == 3
    done_before = {p.name: p.read_text() for p in (tmp_path / "points").glob("*.json")}
    second = run(CFG, tmp_path)
    assert second["evaluated"] == 3 and second["remaining"] == 0
    for name, text in done_before.items():            # completed points untouched
        assert (tmp_path / "points" / name).read_text() == text
    rec = json.loads((tmp_path / "points" / "0000004.json").read_text())
    assert rec["result"]["y"] == 17.0 and "cpu_h" in rec["cost"]
    h5 = collect(CFG, tmp_path)
    assert h5.exists()


def test_seeds_independent_of_execution_order(tmp_path):
    a = run(CFG, tmp_path / "a", shard=(1, 2))
    run(CFG, tmp_path / "a", shard=(0, 2))
    run(CFG, tmp_path / "b")
    for i in range(len(grid_points(CFG))):
        ra = json.loads((tmp_path / "a" / "points" / f"{i:07d}.json").read_text())
        rb = json.loads((tmp_path / "b" / "points" / f"{i:07d}.json").read_text())
        assert ra["result"] == rb["result"]


def test_times_to_targets_early_stop_matches_full_scan():
    """Stopping the scan after the last first-crossing leaves every crossing time unchanged."""
    import numpy as np
    from p1_atomdet.production import times_to_targets
    calls = []

    def eps(t):  # falls as exp(-t/tau), then rises (loss); minimum near 0.2 s
        calls.append(t)
        return 0.5 * np.exp(-t / 0.012) + 2e-5 * (t / 0.2) ** 4
    res, (ts, _) = times_to_targets(eps, [1e-3, 1e-4], t_lo=1e-3, t_hi=5.0)
    for tg, tau in ((1e-3, 0.012 * np.log(500)), (1e-4, 0.012 * np.log(5000))):
        assert abs(res[tg] / tau - 1) < 5e-3
    assert ts[-1] < 0.2 and max(calls) < 0.2           # never evaluated the expensive long-t points
    res_s, _ = times_to_targets(eps, [1e-5], t_lo=1e-3, t_hi=5.0,
                                survival_of_t=lambda t: np.exp(-t), min_survival=0.99)
    assert res_s[1e-5] is None                          # unreachable before survival drops below 0.99

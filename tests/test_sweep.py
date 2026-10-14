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

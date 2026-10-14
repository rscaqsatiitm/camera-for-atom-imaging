"""Wall-clock, CPU time and peak memory logging for runs and sweep points."""
from __future__ import annotations

import json
import os
import resource
import time
from contextlib import contextmanager
from pathlib import Path

import psutil


def _cpu_seconds() -> float:
    ru = resource.getrusage(resource.RUSAGE_SELF)
    rc = resource.getrusage(resource.RUSAGE_CHILDREN)
    return ru.ru_utime + ru.ru_stime + rc.ru_utime + rc.ru_stime


def _peak_rss_mb() -> float:
    # ru_maxrss is in kB on Linux
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


@contextmanager
def runlog(label: str, logfile: str | Path | None = None, extra: dict | None = None):
    """Context manager; yields a dict that is filled with timing on exit and appended to logfile (JSON lines)."""
    rec = {"label": label, "pid": os.getpid(), "n_cpu": psutil.cpu_count(logical=True)}
    if extra:
        rec.update(extra)
    t0, c0 = time.perf_counter(), _cpu_seconds()
    try:
        yield rec
    finally:
        rec["wall_s"] = time.perf_counter() - t0
        rec["cpu_s"] = _cpu_seconds() - c0
        rec["cpu_h"] = rec["cpu_s"] / 3600.0
        rec["peak_rss_mb"] = _peak_rss_mb()
        if logfile is not None:
            p = Path(logfile)
            p.parent.mkdir(parents=True, exist_ok=True)
            with p.open("a") as f:
                f.write(json.dumps(rec, default=str) + "\n")

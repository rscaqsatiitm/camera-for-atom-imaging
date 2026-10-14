"""Scan presets and configs for open placeholders and build the provenance ledger.

    python -m camsim_common.placeholders [paths...]        # list open [VALUE?] entries
    python -m camsim_common.placeholders --ledger [paths]  # full provenance ledger (markdown)

A placeholder is any ProvValue-like mapping with prov == PLACEHOLDER, or any string that
contains "[VALUE?". The same scan is used by the report generator, so the report list and
this command cannot drift apart.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PATHS = [ROOT / "camsim_common" / "presets", ROOT / "configs"]


def _iter_files(paths):
    for p in paths:
        p = Path(p)
        if p.is_dir():
            yield from sorted(x for x in p.rglob("*") if x.suffix in (".yaml", ".yml", ".json"))
        elif p.exists():
            yield p


def _walk(obj, prefix):
    if isinstance(obj, dict):
        if "prov" in obj and "value" in obj:
            yield prefix, obj
            return
        for k, v in obj.items():
            yield from _walk(v, f"{prefix}.{k}" if prefix else str(k))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from _walk(v, f"{prefix}[{i}]")
    elif isinstance(obj, str) and "[VALUE?" in obj:
        yield prefix, {"value": obj, "unit": "", "prov": "PLACEHOLDER", "pointer": "", "note": "inline string"}


def scan(paths=None) -> list[dict]:
    """Return a list of records {file, path, value, unit, prov, pointer, note} for every tagged value."""
    rows = []
    for f in _iter_files(paths or DEFAULT_PATHS):
        try:
            data = yaml.safe_load(f.read_text())
        except Exception as e:  # pragma: no cover
            rows.append({"file": str(f), "path": "<parse error>", "value": str(e), "unit": "", "prov": "ERROR",
                         "pointer": "", "note": ""})
            continue
        for path, rec in _walk(data, ""):
            rows.append({"file": str(Path(f).relative_to(ROOT)) if Path(f).is_relative_to(ROOT) else str(f),
                         "path": path, **{k: rec.get(k, "") for k in ("value", "unit", "prov", "pointer", "note")}})
    return rows


def open_placeholders(paths=None) -> list[dict]:
    return [r for r in scan(paths) if r["prov"] == "PLACEHOLDER"]


def fmt_value(v) -> str:
    if isinstance(v, list):
        return f"[{len(v)} values]" if len(v) > 4 else str(v)
    return str(v)


def ledger_markdown(rows) -> str:
    out = ["| value | location | type | pointer |", "|---|---|---|---|"]
    for r in rows:
        v = f"{fmt_value(r['value'])} {r['unit']}".strip()
        if r["prov"] == "PLACEHOLDER":
            v = f"[VALUE?: approx. {v}, verify]"
        out.append(f"| {v} | {r['file']}: `{r['path']}` | {r['prov']} | {r['pointer'] or r['note']} |")
    return "\n".join(out)


def placeholder_markdown(rows) -> str:
    out = []
    for r in rows:
        v = f"{fmt_value(r['value'])} {r['unit']}".strip()
        note = f" ({r['note']})" if r["note"] else ""
        out.append(f"- `{r['file']}`: `{r['path']}` = [VALUE?: approx. {v}, verify]{note}")
    return "\n".join(out) if out else "(none)"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="*")
    ap.add_argument("--ledger", action="store_true")
    a = ap.parse_args(argv)
    paths = a.paths or DEFAULT_PATHS
    if a.ledger:
        print(ledger_markdown(scan(paths)))
    else:
        rows = open_placeholders(paths)
        print(placeholder_markdown(rows))
        print(f"\n{len(rows)} open placeholders", file=sys.stderr)


if __name__ == "__main__":
    main()

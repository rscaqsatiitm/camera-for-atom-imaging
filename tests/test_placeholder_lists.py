"""The open [VALUE?] lists printed in the reports must equal the scanner output (exact match).

Criterion fixed on 2026-11-19 before the first run: line-for-line equality, order included.
  - reports/p1_validation.md lists the validation-stage files: the two validation presets and
    configs/validation/ (40 entries on 2026-11-19).
  - reports/p1_production_note.md lists the full default scan, presets/ and configs/ (116 entries).
"""
from pathlib import Path

from camsim_common.placeholders import DEFAULT_PATHS, ROOT, open_placeholders, placeholder_markdown

VALIDATION_SCOPE = [ROOT / "camsim_common" / "presets" / "orca_quest_c15550_20up_manetsch.yaml",
                    ROOT / "camsim_common" / "presets" / "zyla55_picken_globalshutter_200MHz.yaml",
                    ROOT / "configs" / "validation"]


def _report_list(path: Path, heading: str) -> list[str]:
    text = path.read_text().splitlines()
    i = next(k for k, line in enumerate(text) if line.startswith(heading))
    out = []
    for line in text[i + 1:]:
        if line.startswith("## "):
            break
        if line.startswith("- `"):
            out.append(line)
    return out


def _scan_list(paths) -> list[str]:
    return [line for line in placeholder_markdown(open_placeholders(paths)).splitlines() if line.startswith("- `")]


def test_validation_report_list_matches_scanner():
    rep = _report_list(ROOT / "reports" / "p1_validation.md", "## Open placeholders")
    assert rep, "no placeholder list found in the validation report"
    assert rep == _scan_list(VALIDATION_SCOPE)


def test_production_note_list_matches_full_scan():
    rep = _report_list(ROOT / "reports" / "p1_production_note.md", "## Open placeholders")
    assert rep, "no placeholder list found in the production note"
    assert rep == _scan_list(DEFAULT_PATHS)

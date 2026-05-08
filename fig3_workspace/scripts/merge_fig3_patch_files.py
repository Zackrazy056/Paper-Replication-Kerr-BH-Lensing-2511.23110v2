"""Merge Fig.3 patch JSON files by (a, omega, l).

Later input files win when quality ties. Quality is ranked by validity and by
the available angular-grid length, so a valid 281-point row replaces an older
short-grid row automatically.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def load_rows(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    rows = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(rows, list):
        raise ValueError(f"{path} is not a JSON list")
    return [row for row in rows if isinstance(row, dict)]


def row_key(row: dict[str, Any]) -> tuple[int, float, float] | None:
    if not {"l", "a", "omega"}.issubset(row):
        return None
    return int(row["l"]), round(float(row["a"]), 3), round(float(row["omega"]), 10)


def theta_len(row: dict[str, Any]) -> int:
    s_theta = row.get("S_theta", [])
    s_pi = row.get("S_pi_minus_theta", [])
    if not isinstance(s_theta, list) or not isinstance(s_pi, list):
        return -1
    return min(len(s_theta), len(s_pi))


def quality(row: dict[str, Any]) -> tuple[int, int]:
    return (1 if row.get("valid", True) is True else 0, theta_len(row))


def merge_rows(paths: list[Path]) -> list[dict[str, Any]]:
    merged: dict[tuple[int, float, float], dict[str, Any]] = {}
    passthrough: list[dict[str, Any]] = []

    for path in paths:
        for row in load_rows(path):
            key = row_key(row)
            if key is None:
                passthrough.append(row)
                continue
            previous = merged.get(key)
            if previous is None or quality(row) >= quality(previous):
                merged[key] = row

    ordered = sorted(merged.values(), key=lambda r: (round(float(r["a"]), 3), round(float(r["omega"]), 10), int(r["l"])))
    return [*passthrough, *ordered]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", help="Patch JSON files to merge, in priority order.")
    parser.add_argument("--base", help="Existing base JSON to read before inputs.")
    parser.add_argument("--out", required=True, help="Merged output JSON.")
    args = parser.parse_args()

    paths = []
    if args.base:
        paths.append(Path(args.base))
    paths.extend(Path(item) for item in args.inputs)

    rows = merge_rows(paths)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    print(f"Merged {len(paths)} files -> {out_path} ({len(rows)} rows)")


if __name__ == "__main__":
    main()

"""Audit Fig.56 spherical-series-reduction cache coverage.

The strict Fig.56/Table I path needs one spherical SR row per (a, omega)
combination present in the cached scattering data. This script reports the
missing combinations and writes a small machine-readable manifest.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"

INPUT_FILES = [
    "fig4_data.json",
    "fig56_patch_data.json",
    "fig56_patch_a_minus099_theta30_hi.json",
]
SR_FILE = "fig56_sr_spherical.json"


def _key(a: float, omega: float) -> tuple[float, float]:
    return round(float(a), 3), round(float(omega), 10)


def load_json(path: Path) -> list[dict]:
    if not path.exists():
        return []
    return json.loads(path.read_text(encoding="utf-8"))


def expected_combos() -> set[tuple[float, float]]:
    combos: set[tuple[float, float]] = set()
    for name in INPUT_FILES:
        for row in load_json(DATA_DIR / name):
            if not row.get("valid", True):
                continue
            required = {"a", "omega", "l", "Re_C", "B_ratio_re", "B_ratio_im", "S_gamma"}
            if not required.issubset(row):
                continue
            l_val = int(row["l"])
            if 2 <= l_val <= 60:
                combos.add(_key(row["a"], row["omega"]))
    return combos


def existing_combos() -> set[tuple[float, float]]:
    rows = load_json(DATA_DIR / SR_FILE)
    return {_key(row["a"], row["omega"]) for row in rows if "a" in row and "omega" in row}


def group_by_a(items: set[tuple[float, float]]) -> dict[str, list[float]]:
    grouped: dict[float, list[float]] = defaultdict(list)
    for a_val, omega in sorted(items):
        grouped[a_val].append(omega)
    return {f"{a_val:.3f}": values for a_val, values in sorted(grouped.items())}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strict", action="store_true", help="Exit nonzero when coverage is incomplete.")
    parser.add_argument(
        "--out",
        default=str(RESULTS_DIR / "fig56_sr_coverage.json"),
        help="Coverage manifest path.",
    )
    args = parser.parse_args()

    expected = expected_combos()
    existing = existing_combos()
    missing = expected - existing
    extra = existing - expected

    manifest = {
        "input_files": INPUT_FILES,
        "sr_file": SR_FILE,
        "expected_count": len(expected),
        "existing_count": len(existing),
        "missing_count": len(missing),
        "extra_count": len(extra),
        "missing_by_a": group_by_a(missing),
        "existing_by_a": group_by_a(existing),
        "extra_by_a": group_by_a(extra),
        "wolfram_batch_hint": (
            "$env:FIG56_SR_A_FILTER='0.8'; $env:FIG56_SR_BATCH='8'; "
            "wolframscript -file mathematica/generate_fig56_sr_spherical.wls"
        ),
    }

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(
        "Fig.56 SR coverage: "
        f"{len(existing)}/{len(expected)} existing, {len(missing)} missing, {len(extra)} extra"
    )
    for a_key, values in manifest["missing_by_a"].items():
        print(f"  missing a={a_key}: {len(values)} omega points")
    print(f"Coverage manifest: {out_path}")

    if args.strict and missing:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

"""Audit Fig.3 numerical data coverage and current curve stability.

This script is intentionally diagnostic: it does not try to decide the final
physics method. It records which cached modes are actually usable by
``main_fig3.py`` and how much of the plotted curve is hidden by the paper-style
0..80 y-axis clipping.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import mpmath as mp
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.physics_fig3 import Fig3Engine, M

DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"
MODE_FILES = [
    DATA_DIR / "fig3_data.json",
    DATA_DIR / "fig3_patch_w2_a0_a099.json",
]
SR_FILE = DATA_DIR / "fig3_sr_spherical.json"
THETA_DEG = np.linspace(40.0, 180.0, 281)
THETA_RAD_MP = [mp.mpf(str(theta)) * mp.pi / 180 for theta in THETA_DEG]
CASES = [(0.0, 0.6), (0.99, 0.6), (0.0, 2.0), (0.99, 2.0)]


def load_rows(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    rows = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(rows, list):
        return []
    return [row for row in rows if isinstance(row, dict)]


def combo_key(a_val: float, omega_val: float) -> tuple[str, str]:
    return f"{float(a_val):.3f}", f"{float(omega_val):.4f}"


def mode_key(row: dict[str, Any]) -> tuple[int, str, str]:
    return int(row["l"]), *combo_key(float(row["a"]), float(row["omega"]))


def theta_len(row: dict[str, Any], name: str = "S_theta") -> int:
    vals = row.get(name, [])
    return len(vals) if isinstance(vals, list) else -1


def raw_mode_coverage() -> dict[str, Any]:
    report: dict[str, Any] = {}
    for path in MODE_FILES:
        rows = load_rows(path)
        grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
        for row in rows:
            if not {"l", "a", "omega"}.issubset(row):
                continue
            grouped[combo_key(row["a"], row["omega"])].append(row)

        file_report: dict[str, Any] = {}
        for key, key_rows in sorted(grouped.items()):
            l_values = sorted({int(row["l"]) for row in key_rows})
            valid_rows = [row for row in key_rows if row.get("valid", True) is True]
            invalid_l = sorted(int(row["l"]) for row in key_rows if row.get("valid", True) is not True)
            theta_lengths = sorted({theta_len(row) for row in key_rows})
            combo_name = f"a={key[0]},omega={key[1]}"
            file_report[combo_name] = {
                "rows": len(key_rows),
                "valid_rows": len(valid_rows),
                "l_min": min(l_values) if l_values else None,
                "l_max": max(l_values) if l_values else None,
                "theta_lengths": theta_lengths,
                "invalid_l_head": invalid_l[:20],
                "invalid_count": len(invalid_l),
            }
        report[path.name] = file_report
    return report


def effective_mode_coverage() -> dict[str, Any]:
    best: dict[tuple[int, str, str], dict[str, Any]] = {}
    for path in MODE_FILES:
        for row in load_rows(path):
            if row.get("valid", True) is not True:
                continue
            if not {"l", "a", "omega", "S_theta", "S_pi_minus_theta"}.issubset(row):
                continue
            key = mode_key(row)
            row_len = min(theta_len(row, "S_theta"), theta_len(row, "S_pi_minus_theta"))
            previous = best.get(key)
            if previous is None or row_len >= int(previous["theta_len"]):
                best[key] = {"theta_len": row_len, "source": path.name}

    grouped: dict[tuple[str, str], list[tuple[int, dict[str, Any]]]] = defaultdict(list)
    for (l_val, a_str, omega_str), info in best.items():
        grouped[(a_str, omega_str)].append((l_val, info))

    report: dict[str, Any] = {}
    for key, values in sorted(grouped.items()):
        l_values = sorted(l_val for l_val, _ in values)
        full_l = sorted(l_val for l_val, info in values if int(info["theta_len"]) >= len(THETA_DEG))
        short_l = sorted(l_val for l_val, info in values if int(info["theta_len"]) < len(THETA_DEG))
        source_counts: dict[str, int] = defaultdict(int)
        theta_counts: dict[str, int] = defaultdict(int)
        for _, info in values:
            source_counts[str(info["source"])] += 1
            theta_counts[str(info["theta_len"])] += 1
        combo_name = f"a={key[0]},omega={key[1]}"
        report[combo_name] = {
            "available_modes": len(l_values),
            "l_min": min(l_values) if l_values else None,
            "l_max": max(l_values) if l_values else None,
            "usable_281_modes": len(full_l),
            "usable_281_l_head": full_l[:20],
            "usable_281_l_tail": full_l[-20:],
            "short_grid_modes": len(short_l),
            "short_grid_l_head": short_l[:20],
            "missing_l_2_to_88": [
                l_val for l_val in range(2, 89) if (l_val, key[0], key[1]) not in best
            ],
            "source_counts": dict(sorted(source_counts.items())),
            "theta_len_counts": dict(sorted(theta_counts.items(), key=lambda item: int(item[0]))),
        }
    return report


def sr_cache_summary() -> dict[str, Any]:
    rows = load_rows(SR_FILE)
    report: dict[str, Any] = {}
    for row in rows:
        if not {"a", "omega"}.issubset(row):
            continue
        key = combo_key(row["a"], row["omega"])
        combo_name = f"a={key[0]},omega={key[1]}"
        report[combo_name] = {
            "k_sr": row.get("k_sr"),
            "l_mode_max": row.get("l_mode_max"),
            "l_basis_max": row.get("l_basis_max"),
            "fit_residual": row.get("fit_residual"),
            "fit_condition": row.get("fit_condition"),
            "theta_points": theta_len(row, "theta_deg_list"),
        }
    return report


def stats(values: np.ndarray) -> dict[str, Any]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return {"finite": False}
    return {
        "finite": True,
        "min": float(np.min(finite)),
        "max": float(np.max(finite)),
        "mean": float(np.mean(finite)),
        "p95": float(np.percentile(finite, 95)),
        "points_gt_80": int(np.sum(finite > 80.0)),
        "points_gt_200": int(np.sum(finite > 200.0)),
        "theta_at_max_deg": float(THETA_DEG[int(np.argmax(values))]),
    }


def split_precomputed_sr(
    engine: Fig3Engine, a_val: float, omega_val: float, l_max: int, k_sr: int
) -> dict[str, Any] | None:
    a_str, omega_str = combo_key(a_val, omega_val)
    f_curve = engine._find_spherical_sr_curve(a_str, omega_str, THETA_RAD_MP, k=k_sr)
    if f_curve is None:
        return None

    available_l = engine._available_l_values(a_str, omega_str)
    if not available_l:
        return None
    use_l_max = min(int(l_max), max(available_l))
    omega_mp = mp.mpf(str(omega_val))
    pref = engine._scattering_prefactor(omega_mp, mode="eq12")

    coeff_g: dict[int, mp.mpc] = {}
    used_g: list[int] = []
    skipped_g: list[dict[str, Any]] = []
    for l_val in range(2, use_l_max + 1):
        mode = engine.data_dict.get((l_val, a_str, omega_str))
        if mode is None:
            skipped_g.append({"l": l_val, "reason": "missing"})
            continue
        if len(mode["S_pi_minus_theta"]) < len(THETA_RAD_MP):
            skipped_g.append({"l": l_val, "reason": "short_theta_grid"})
            continue

        term1 = mp.power(-1, l_val + 1)
        bracket_g = (
            term1
            * 3
            * mp.j
            * M
            * omega_mp
            / (4 * mp.power(omega_mp, 4))
            * mode["B_ratio"]
        )
        coeff_g[l_val] = pref * mp.power(-1, l_val) * mode["S_0"] * mp.conj(bracket_g)
        used_g.append(l_val)

    f2 = np.zeros(len(THETA_RAD_MP), dtype=float)
    g2 = np.zeros(len(THETA_RAD_MP), dtype=float)
    for idx in range(len(THETA_RAD_MP)):
        sum_g = mp.mpc("0.0")
        for l_val, coeff in coeff_g.items():
            sum_g += coeff * engine.data_dict[(l_val, a_str, omega_str)]["S_pi_minus_theta"][idx]
        f2[idx] = float(abs(f_curve[idx]) ** 2)
        g2[idx] = float(abs(sum_g) ** 2)

    return {
        "used_g_modes": len(used_g),
        "used_g_l_head": used_g[:20],
        "used_g_l_tail": used_g[-20:],
        "skipped_g_modes": len(skipped_g),
        "skipped_g_head": skipped_g[:20],
        "f_abs2": stats(f2),
        "g_abs2": stats(g2),
        "g_fraction_p95": float(np.percentile(g2 / np.maximum(f2 + g2, 1e-300), 95)),
    }


def curve_summary(engine: Fig3Engine, l_max: int, k_sr: int) -> dict[str, Any]:
    report: dict[str, Any] = {}
    for a_val, omega_val in CASES:
        a_str, omega_str = combo_key(a_val, omega_val)
        y_values = np.asarray(
            engine.calc_angular_distribution_series_reduction(
                a_val=a_val,
                omega_val=omega_val,
                theta_rad_list=THETA_RAD_MP,
                l_max=l_max,
                k=k_sr,
                m=2,
                prefactor_mode="eq12",
            ),
            dtype=float,
        )
        combo_name = f"a={a_str},omega={omega_str}"
        report[combo_name] = {
            "spherical_f_sr_available": engine._find_spherical_sr_curve(
                a_str, omega_str, THETA_RAD_MP, k=k_sr
            )
            is not None,
            "series_reduction_total": stats(y_values),
        }
        split = split_precomputed_sr(engine, a_val, omega_val, l_max=l_max, k_sr=k_sr)
        if split is not None:
            report[combo_name]["precomputed_sr_split"] = split
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--l-max", type=int, default=88, help="Fig.3 l_max used for the audit.")
    parser.add_argument("--k-sr", type=int, default=2, help="Series-reduction order.")
    parser.add_argument(
        "--out",
        default=str(RESULTS_DIR / "fig3_numeric_audit.json"),
        help="Output JSON path.",
    )
    args = parser.parse_args()

    data_files = [str(path) for path in [*MODE_FILES, SR_FILE]]
    engine = Fig3Engine(data_files)

    manifest = {
        "theta_grid": {"min_deg": 40.0, "max_deg": 180.0, "points": len(THETA_DEG)},
        "settings": {"l_max": args.l_max, "k_sr": args.k_sr, "y_clip_max": 80.0},
        "raw_mode_coverage": raw_mode_coverage(),
        "effective_mode_coverage": effective_mode_coverage(),
        "spherical_sr_cache": sr_cache_summary(),
        "curve_summary": curve_summary(engine, l_max=args.l_max, k_sr=args.k_sr),
    }

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"Fig.3 numeric audit written to: {out_path}")
    for combo, info in manifest["curve_summary"].items():
        total = info["series_reduction_total"]
        print(
            f"  {combo}: max={total.get('max', float('nan')):.3g}, "
            f"p95={total.get('p95', float('nan')):.3g}, "
            f"points>80={total.get('points_gt_80', 0)}, "
            f"spherical_f={info['spherical_f_sr_available']}"
        )


if __name__ == "__main__":
    main()

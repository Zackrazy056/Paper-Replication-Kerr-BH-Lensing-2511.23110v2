"""Plot Fig.3-style 2x2 method comparison grid.

Panels:
  top row: a = 0, bottom row: a = 0.99M
  left column: M omega = 0.6, right column: M omega = 2.0

Each panel compares series reduction with Cesaro alpha=2 and alpha=5 on the
same linear axes.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from src.physics_fig3 import Fig3Engine

ROOT = Path(__file__).resolve().parent
DATA_FILES = [
    ROOT / "data" / "fig3_data.json",
    ROOT / "data" / "fig3_patch_w2_a0_a099.json",
    ROOT / "data" / "fig3_sr_spherical.json",
]


def _stats(values: np.ndarray) -> dict[str, float | int | bool]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return {"finite": False}
    return {
        "finite": True,
        "min": float(np.min(finite)),
        "max": float(np.max(finite)),
        "p95": float(np.percentile(finite, 95)),
        "points_gt_80": int(np.sum(finite > 80.0)),
    }


def _compute_curves(
    engine: Fig3Engine,
    theta_rad: np.ndarray,
    a_val: float,
    omega: float,
    l_max: int,
) -> dict[str, np.ndarray]:
    curves: dict[str, np.ndarray] = {}
    curves["Series reduction, k=2"] = np.asarray(
        engine.calc_angular_distribution_series_reduction(
            a_val=a_val,
            omega_val=omega,
            theta_rad_list=theta_rad,
            l_max=l_max,
            k=2,
            m=2,
            prefactor_mode="eq12",
        ),
        dtype=float,
    )
    curves["Cesaro, alpha=2"] = np.asarray(
        engine.calc_angular_distribution_cesaro(
            a_val=a_val,
            omega_val=omega,
            num_angles=len(theta_rad),
            l_max=l_max,
            alpha=2,
            m=2,
            prefactor_mode="eq12",
        ),
        dtype=float,
    )
    curves["Cesaro, alpha=5"] = np.asarray(
        engine.calc_angular_distribution_cesaro(
            a_val=a_val,
            omega_val=omega,
            num_angles=len(theta_rad),
            l_max=l_max,
            alpha=5,
            m=2,
            prefactor_mode="eq12",
        ),
        dtype=float,
    )
    return curves


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--l-max", type=int, default=88, help="Mode truncation.")
    parser.add_argument(
        "--out",
        default=str(ROOT / "results" / "Fig3_paper_methods_grid.png"),
        help="Output PNG path.",
    )
    parser.add_argument(
        "--metrics",
        default=str(ROOT / "results" / "fig3_paper_methods_grid_metrics.json"),
        help="Output metrics JSON path.",
    )
    args = parser.parse_args()

    theta_deg = np.linspace(40.0, 180.0, 281)
    theta_rad = np.deg2rad(theta_deg)
    engine = Fig3Engine([str(path) for path in DATA_FILES])

    spins = [(0.0, r"$a=0$"), (0.99, r"$a=0.99M$")]
    omegas = [(0.6, r"$M\omega=0.6$"), (2.0, r"$M\omega=2.0$")]
    styles = {
        "Series reduction, k=2": {
            "label": r"Series reduction, $k=2$",
            "color": "red",
            "linestyle": "-",
            "linewidth": 2.0,
        },
        "Cesaro, alpha=2": {
            "label": r"Cesaro, $\alpha=2$",
            "color": "blue",
            "linestyle": "-.",
            "linewidth": 1.9,
        },
        "Cesaro, alpha=5": {
            "label": r"Cesaro, $\alpha=5$",
            "color": "black",
            "linestyle": "-.",
            "linewidth": 1.9,
        },
    }

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), sharex=True, sharey=True)
    metrics: dict[str, dict[str, object]] = {}

    for row_idx, (a_val, a_label) in enumerate(spins):
        for col_idx, (omega, omega_label) in enumerate(omegas):
            ax = axes[row_idx, col_idx]
            curves = _compute_curves(engine, theta_rad, a_val, omega, args.l_max)
            combo_key = f"Momega={omega:.1f},a={a_val:.2f}"
            metrics[combo_key] = {}

            for curve_name, y_values in curves.items():
                style = styles[curve_name]
                ax.plot(
                    theta_deg,
                    y_values,
                    label=style["label"],
                    color=style["color"],
                    linestyle=style["linestyle"],
                    linewidth=style["linewidth"],
                    clip_on=True,
                )
                metrics[combo_key][curve_name] = _stats(y_values)

            ax.set_title(f"{omega_label}, {a_label}", fontsize=14)
            ax.set_xlim(30.0, 180.0)
            ax.set_ylim(0.0, 80.0)
            ax.set_xticks(np.arange(40, 181, 20))
            ax.set_yticks(np.arange(0, 81, 20))
            ax.grid(True, linestyle=":", linewidth=0.7, alpha=0.55)
            ax.tick_params(axis="both", direction="in", length=6, width=1.1, labelsize=11)

            if row_idx == 1:
                ax.set_xlabel(r"$\theta$ (degrees)", fontsize=13)
            if col_idx == 0:
                ax.set_ylabel(r"$M^{-2}\,d\sigma/d\Omega$", fontsize=13)
            if row_idx == 0 and col_idx == 0:
                ax.legend(loc="upper right", fontsize=10, framealpha=0.92)

    fig.subplots_adjust(left=0.09, right=0.985, bottom=0.09, top=0.94, wspace=0.08, hspace=0.18)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)

    metrics_path = Path(args.metrics)
    metrics_path.parent.mkdir(parents=True, exist_ok=True)
    metrics_path.write_text(json.dumps(metrics, indent=2), encoding="utf-8")

    print(f"[+] Saved Fig.3 methods grid: {out_path}")
    print(f"[+] Saved metrics: {metrics_path}")
    for combo_key, combo_metrics in metrics.items():
        sr_stats = combo_metrics["Series reduction, k=2"]
        print(
            f"    {combo_key}: SR max={sr_stats.get('max', float('nan')):.3g}, "
            f"points>80={sr_stats.get('points_gt_80', 0)}"
        )


if __name__ == "__main__":
    main()

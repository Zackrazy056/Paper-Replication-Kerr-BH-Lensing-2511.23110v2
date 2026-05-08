"""Diagnose Fig.3 convergence at one scattering angle as l_max varies."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.physics_fig3 import Fig3Engine


DATA_FILES = [
    ROOT / "data" / "fig3_data.json",
    ROOT / "data" / "fig3_patch_w2_a0_a099.json",
    ROOT / "data" / "fig3_sr_spherical.json",
]


def compute_value(engine: Fig3Engine, method: str, a_val: float, omega: float, theta_rad: np.ndarray, l_max: int) -> float:
    if method == "sr":
        values = engine.calc_angular_distribution_series_reduction(
            a_val=a_val,
            omega_val=omega,
            theta_rad_list=theta_rad,
            l_max=l_max,
            k=2,
            m=2,
            prefactor_mode="eq12",
        )
    elif method == "cesaro2":
        values = engine.calc_angular_distribution_cesaro(
            a_val=a_val,
            omega_val=omega,
            num_angles=len(theta_rad),
            l_max=l_max,
            alpha=2,
            m=2,
            prefactor_mode="eq12",
        )
    elif method == "cesaro5":
        values = engine.calc_angular_distribution_cesaro(
            a_val=a_val,
            omega_val=omega,
            num_angles=len(theta_rad),
            l_max=l_max,
            alpha=5,
            m=2,
            prefactor_mode="eq12",
        )
    elif method == "plain":
        values = engine.calc_angular_distribution(
            a_val=a_val,
            omega_val=omega,
            num_angles=len(theta_rad),
            prefactor_mode="eq12",
            l_max=l_max,
        )
    else:
        raise ValueError(f"Unknown method: {method}")
    return float(values[0])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a", type=float, default=0.99, help="Spin parameter.")
    parser.add_argument("--omega", type=float, default=2.0, help="Dimensionless frequency Momega.")
    parser.add_argument("--theta-deg", type=float, default=120.0, help="Scattering angle in degrees.")
    parser.add_argument("--l-min", type=int, default=2)
    parser.add_argument("--l-max", type=int, default=88)
    parser.add_argument(
        "--methods",
        default="sr,cesaro2,cesaro5",
        help="Comma-separated methods: sr, cesaro2, cesaro5, plain.",
    )
    parser.add_argument(
        "--out-prefix",
        default=str(ROOT / "results" / "fig3_partial_sum_theta120"),
        help="Output prefix for CSV and PNG.",
    )
    args = parser.parse_args()

    methods = [item.strip() for item in args.methods.split(",") if item.strip()]
    theta_rad = np.asarray([np.deg2rad(args.theta_deg)], dtype=float)
    engine = Fig3Engine([str(path) for path in DATA_FILES])

    rows: list[dict[str, float | int | str]] = []
    for l_max in range(args.l_min, args.l_max + 1):
        for method in methods:
            rows.append(
                {
                    "a": args.a,
                    "omega": args.omega,
                    "theta_deg": args.theta_deg,
                    "method": method,
                    "l_max": l_max,
                    "dsdo": compute_value(engine, method, args.a, args.omega, theta_rad, l_max),
                }
            )

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = out_prefix.with_suffix(".csv")
    png_path = out_prefix.with_suffix(".png")

    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["a", "omega", "theta_deg", "method", "l_max", "dsdo"])
        writer.writeheader()
        writer.writerows(rows)

    fig, ax = plt.subplots(figsize=(8, 5))
    for method in methods:
        method_rows = [row for row in rows if row["method"] == method]
        ax.plot([row["l_max"] for row in method_rows], [row["dsdo"] for row in method_rows], label=method)
    ax.set_xlabel(r"$l_{\max}$")
    ax.set_ylabel(r"$M^{-2}d\sigma/d\Omega$")
    ax.set_title(rf"$a={args.a:g},\ M\omega={args.omega:g},\ \theta={args.theta_deg:g}^\circ$")
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend()
    fig.tight_layout()
    fig.savefig(png_path, dpi=200)
    plt.close(fig)

    print(f"Saved CSV: {csv_path}")
    print(f"Saved PNG: {png_path}")


if __name__ == "__main__":
    main()

import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = os.path.dirname(__file__)
REPO_DIR = os.path.dirname(SCRIPT_DIR)
if REPO_DIR not in sys.path:
    sys.path.insert(0, REPO_DIR)

from src.physics_fig3 import Fig3Engine


def roughness_metric(y):
    arr = np.asarray(y, dtype=float)
    if arr.size < 3:
        return 0.0
    d2 = np.diff(arr, n=2)
    denom = float(np.std(arr))
    if denom <= 0:
        return 0.0
    return float(np.std(d2) / denom)


def main():
    repo_dir = REPO_DIR
    data_files = [
        os.path.join(repo_dir, "data", "fig3_data.json"),
        os.path.join(repo_dir, "data", "fig3_patch_w2_a0_a099.json"),
    ]
    engine = Fig3Engine(data_files)

    theta_deg = np.linspace(40, 180, 281)
    theta_rad = np.deg2rad(theta_deg)
    omega = 2.0
    l_max = 74

    cases = [
        {"a": 0.0, "label": "a = 0"},
        {"a": 0.99, "label": "a = 0.99M"},
    ]
    methods = [
        {
            "name": r"Series reduction ($k=2$)",
            "kind": "sr",
            "color": "red",
            "ls": "-",
        },
        {
            "name": r"Cesaro ($\alpha=2$)",
            "kind": "cesaro",
            "alpha": 2,
            "color": "blue",
            "ls": "-.",
        },
        {
            "name": r"Cesaro ($\alpha=5$)",
            "kind": "cesaro",
            "alpha": 5,
            "color": "black",
            "ls": "-.",
        },
    ]

    fig, axes = plt.subplots(2, 1, figsize=(8, 10), sharex=True, sharey=True)
    fig.patch.set_facecolor("#d3d3d3")
    metrics = []

    for i, case in enumerate(cases):
        ax = axes[i]
        ax.set_facecolor("#d3d3d3")
        for method in methods:
            if method["kind"] == "sr":
                y = engine.calc_angular_distribution_series_reduction(
                    a_val=case["a"],
                    omega_val=omega,
                    theta_rad_list=theta_rad,
                    l_max=l_max,
                    k=2,
                    m=2,
                    prefactor_mode="eq12",
                )
            else:
                y = engine.calc_angular_distribution_cesaro(
                    a_val=case["a"],
                    omega_val=omega,
                    num_angles=len(theta_rad),
                    l_max=l_max,
                    alpha=method["alpha"],
                    m=2,
                    prefactor_mode="eq12",
                )

            y_arr = np.asarray(y, dtype=float)
            y_clip = np.clip(y_arr, 0.0, 80.0)
            ax.plot(
                theta_deg,
                y_clip,
                color=method["color"],
                linestyle=method["ls"],
                linewidth=2,
                label=method["name"],
            )

            metrics.append(
                {
                    "a": case["a"],
                    "omega": omega,
                    "method": method["name"],
                    "l_max": l_max,
                    "roughness": roughness_metric(y_arr),
                    "raw_min": float(np.min(y_arr)),
                    "raw_max": float(np.max(y_arr)),
                }
            )

        ax.set_xlim(40, 180)
        ax.set_ylim(0, 80)
        ax.set_yticks(np.arange(0, 81, 20))
        ax.grid(True, linestyle=":", alpha=0.7)
        ax.set_title(rf"$M\omega = 2.0, {case['label']}$", fontsize=30)
        if i == 0:
            ax.legend(loc="upper center", fontsize=20, framealpha=0.9)

    axes[-1].set_xlabel(r"$\theta$ (degrees)", fontsize=22)
    axes[0].tick_params(labelsize=18, width=1.6, length=8, direction="in")
    axes[1].tick_params(labelsize=18, width=1.6, length=8, direction="in")

    out_dir = os.path.join(repo_dir, "results")
    os.makedirs(out_dir, exist_ok=True)

    fig_path = os.path.join(out_dir, "Fig3_w2_three_methods_compare.png")
    plt.tight_layout()
    plt.savefig(fig_path, dpi=300)
    plt.close(fig)

    json_path = os.path.join(out_dir, "fig3_w2_method_audit_metrics.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2, ensure_ascii=False)

    print(f"[+] Saved figure: {fig_path}")
    print(f"[+] Saved metrics: {json_path}")


if __name__ == "__main__":
    main()

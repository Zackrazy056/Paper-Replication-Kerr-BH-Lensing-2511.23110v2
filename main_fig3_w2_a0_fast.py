import os

import matplotlib.pyplot as plt
import numpy as np

from src.physics_fig3 import Fig3Engine


def main():
    print("=== Fig.3 Fast Single Case: a=0, Momega=2.0 ===")

    # Must match mathematica/generate_fig3_w2_a0_fast.wls (1 degree spacing).
    theta_deg = np.arange(40, 181, 1)
    theta_rad = np.deg2rad(theta_deg)

    base_file = os.path.join(os.path.dirname(__file__), "data", "fig3_data_w2_a0_fast.json")
    if not os.path.exists(base_file):
        base_file = os.path.join(os.path.dirname(__file__), "data", "fig3_data.json")
        print("[*] Fast data not found, fallback to data/fig3_data.json")
    data_files = [
        base_file,
        os.path.join(os.path.dirname(__file__), "data", "fig3_patch_w2_a0_a099.json"),
        os.path.join(os.path.dirname(__file__), "data", "fig3_sr_spherical.json"),
    ]
    engine = Fig3Engine(data_files)

    dsdo = engine.calc_angular_distribution_series_reduction(
        a_val=0.0,
        omega_val=2.0,
        theta_rad_list=theta_rad,
        l_max=60,
        k=2,
        m=2,
        prefactor_mode="eq12",
    )

    if not dsdo:
        raise RuntimeError("No curve generated. Please check data file validity.")

    y = np.maximum(np.array(dsdo, dtype=float), 1e-30)

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.semilogy(theta_deg, y, color="red", linestyle="-", linewidth=2)
    ax.set_xlim(40, 180)
    ax.set_xlabel(r"$\theta$ (degrees)", fontsize=13)
    ax.set_ylabel(r"$M^{-2}\,d\sigma/d\Omega$", fontsize=13)
    ax.set_title(r"Fast trend: $M\omega=2.0,\ a=0$ (SR $k=2$, $l_{\max}=60$)", fontsize=13)
    ax.grid(True, linestyle="--", alpha=0.5)

    out_path = os.path.join(os.path.dirname(__file__), "results", "Fig3_fast_Mw2p0_a0.png")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    print(f"[+] Done. Figure saved to {out_path}")


if __name__ == "__main__":
    main()

import os

import matplotlib.pyplot as plt
import numpy as np

from src.physics_fig3 import Fig3Engine


def main():
    print("=== Fig.3 Series-Reduction Single Test ===")

    # Match Mathematica grid: 40..180 deg with step 2 deg
    theta_deg = np.arange(40, 182, 2)
    theta_rad = np.deg2rad(theta_deg)

    data_files = [
        os.path.join(os.path.dirname(__file__), "data", "fig3_data_single_case.json"),
        os.path.join(os.path.dirname(__file__), "data", "fig3_sr_spherical.json"),
    ]
    engine = Fig3Engine(data_files)

    a_val = 0.0
    omega_val = 0.6
    print(f"[*] Computing series-reduction curve for a={a_val}, Momega={omega_val}")
    dsdo = engine.calc_angular_distribution_series_reduction(
        a_val=a_val,
        omega_val=omega_val,
        theta_rad_list=theta_rad,
        l_max=70,
        k=2,
        m=2,
        prefactor_mode="eq12",
    )

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.plot(theta_deg, dsdo, color="red", linestyle="-", linewidth=2)
    ax.set_xlim(40, 180)
    ax.set_ylim(0, 80)
    ax.set_yticks(np.arange(0, 81, 20))
    ax.set_xlabel(r"$\theta$ (degrees)", fontsize=13)
    ax.set_ylabel(r"$M^{-2}\,d\sigma/d\Omega$", fontsize=13)
    ax.set_title(r"Series reduction: $M\omega=0.6,\ a=0$", fontsize=14)
    ax.grid(True, linestyle="--", alpha=0.5)

    out_path = os.path.join(os.path.dirname(__file__), "results", "Fig3_series_test_Mw0p6_a0.png")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    print(f"[+] Done. Figure saved to {out_path}")


if __name__ == "__main__":
    main()

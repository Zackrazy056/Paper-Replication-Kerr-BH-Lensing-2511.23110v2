import os

import matplotlib.pyplot as plt
import numpy as np

from src.physics_fig3 import Fig3Engine


def _combo_filename(omega, a_val):
    a_tag = f"{a_val:+.3f}".replace("+", "p").replace("-", "m")
    w_tag = f"{omega:.1f}".replace(".", "p")
    return f"Fig3_combo_Mw{w_tag}_a{a_tag}.png"


def main():
    print("=== Kerr GW Lensing Fig.3 Reproduction ===")

    # Data grid must match Mathematica output (0.5 degree spacing).
    theta_deg_data = np.linspace(40, 180, 281)
    theta_rad_data = np.deg2rad(theta_deg_data)
    theta_deg_plot = theta_deg_data

    data_files = [
        os.path.join(os.path.dirname(__file__), "data", "fig3_data.json"),
        os.path.join(os.path.dirname(__file__), "data", "fig3_patch_w2_a0_a099.json"),
        os.path.join(os.path.dirname(__file__), "data", "fig3_sr_spherical.json"),
    ]
    engine = Fig3Engine(data_files)

    # Numerical settings.
    l_max_sr = int(os.environ.get("FIG3_LMAX", "88"))
    k_low = 2
    k_high = 2
    y_min = 0.0
    y_max = 80.0
    method = os.environ.get("FIG3_METHOD", "sr").strip().lower()
    cesaro_alpha = int(os.environ.get("FIG3_CESARO_ALPHA", "5"))
    if method not in {"sr", "cesaro", "plain"}:
        raise ValueError("FIG3_METHOD must be one of: sr, cesaro, plain")

    omegas = [0.5, 2.0]
    a_configs = [
        {"val": 0.0, "label": "a = 0", "color": "navy", "ls": "-"},
        {"val": 0.99, "label": "a = 0.99M", "color": "crimson", "ls": "-."},
    ]

    combo_out_dir = os.path.join(os.path.dirname(__file__), "results", "fig3_combos")
    os.makedirs(combo_out_dir, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    print(f"[*] Method: {method}, l_max={l_max_sr}, cesaro_alpha={cesaro_alpha}")

    for i, omega in enumerate(omegas):
        ax = axes[i]
        k_sr = k_high if omega >= 2.0 else k_low
        print(f"[*] Computing angular distribution for Momega = {omega}")

        for config in a_configs:
            if method == "sr":
                dsdo = engine.calc_angular_distribution_series_reduction(
                    a_val=config["val"],
                    omega_val=omega,
                    theta_rad_list=theta_rad_data,
                    l_max=l_max_sr,
                    k=k_sr,
                    m=2,
                    prefactor_mode="eq12",
                )
            elif method == "cesaro":
                dsdo = engine.calc_angular_distribution_cesaro(
                    a_val=config["val"],
                    omega_val=omega,
                    num_angles=len(theta_rad_data),
                    l_max=l_max_sr,
                    alpha=cesaro_alpha,
                    m=2,
                    prefactor_mode="eq12",
                )
            else:
                dsdo = engine.calc_angular_distribution(
                    a_val=config["val"],
                    omega_val=omega,
                    num_angles=len(theta_rad_data),
                    prefactor_mode="eq12",
                    l_max=l_max_sr,
                )
            if not dsdo:
                print(f"    [!] Empty curve for Momega={omega}, a={config['val']:+.3f}")
                continue

            dsdo_plot_raw = np.array(dsdo, dtype=float)
            dsdo_plot = np.clip(dsdo_plot_raw, y_min, y_max)

            ax.plot(
                theta_deg_plot,
                dsdo_plot,
                label=f"${config['label']}$",
                color=config["color"],
                linestyle=config["ls"],
                linewidth=2,
            )

            fig_single, ax_single = plt.subplots(1, 1, figsize=(7, 5))
            ax_single.plot(
                theta_deg_plot, dsdo_plot, color=config["color"], linestyle="-", linewidth=2.2
            )
            ax_single.set_xlabel(r"$\theta$ (degrees)", fontsize=13)
            ax_single.set_ylabel(r"$M^{-2} d\sigma/d\Omega$", fontsize=13)
            if method == "sr":
                method_label = rf"SR, $k={k_sr}$"
            elif method == "cesaro":
                method_label = rf"Cesaro, $\alpha={cesaro_alpha}$"
            else:
                method_label = "plain sum"
            ax_single.set_title(rf"$M\omega={omega},\ {config['label']}$ ({method_label})", fontsize=13)
            ax_single.set_xlim([40, 180])
            ax_single.set_ylim([y_min, y_max])
            ax_single.set_yticks(np.arange(0, 81, 20))
            ax_single.grid(True, linestyle="--", alpha=0.5)
            single_path = os.path.join(combo_out_dir, _combo_filename(omega, config["val"]))
            fig_single.tight_layout()
            fig_single.savefig(single_path, dpi=300)
            plt.close(fig_single)

            print(
                f"    [+] Saved combo curve: Momega={omega:.1f}, a={config['val']:+.3f}, "
                f"min(raw)={float(np.min(dsdo_plot_raw)):.3e}, max(raw)={float(np.max(dsdo_plot_raw)):.3e}"
            )

        ax.set_xlabel(r"$\theta$ (degrees)", fontsize=14)
        if i == 0:
            ax.set_ylabel(r"$M^{-2} d\sigma/d\Omega$", fontsize=14)
        if method == "sr":
            panel_label = r"SR $k=2$"
        elif method == "cesaro":
            panel_label = rf"Cesaro $\alpha={cesaro_alpha}$"
        else:
            panel_label = "plain sum"
        ax.set_title(r"$M\omega = {}$ ({})".format(omega, panel_label), fontsize=16)
        ax.set_xlim([40, 180])
        ax.set_ylim([y_min, y_max])
        ax.set_yticks(np.arange(0, 81, 20))
        ax.grid(True, linestyle="--", alpha=0.5)
        if i == 0:
            ax.legend(loc="upper right", fontsize=12)

    plt.tight_layout()

    if method == "sr":
        out_name = "Fig3_reproduction_SR.png"
    elif method == "cesaro":
        out_name = f"Fig3_reproduction_Cesaro_a{cesaro_alpha}.png"
    else:
        out_name = "Fig3_reproduction_plain.png"
    out_path = os.path.join(os.path.dirname(__file__), "results", out_name)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.savefig(out_path, dpi=300)
    print(f"[+] Done. Figure saved to {out_path}")
    print(f"[+] Individual combinations saved under: {combo_out_dir}")


if __name__ == "__main__":
    main()

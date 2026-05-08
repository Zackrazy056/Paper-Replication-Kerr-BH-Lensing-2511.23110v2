import os

import matplotlib.pyplot as plt

from src.physics_fig4 import Fig4A0Engine


def main():
    print("=== Fig.4 Reproduction (a=0) ===")

    data_file = os.path.join(os.path.dirname(__file__), "data", "fig4_a0_data.json")
    engine = Fig4A0Engine(data_file)

    # First-pass settings: prioritize trend and speed.
    l_max = 60
    k_sr = 2
    r_sl = 100
    x, y_sr, y_c2, y_c5 = engine.build_curves(a_val=0.0, l_max=l_max, k_sr=k_sr, r_sl=r_sl)
    if not x:
        raise RuntimeError("No valid points generated for Fig.4 a=0.")

    fig, ax = plt.subplots(1, 1, figsize=(7.2, 5.2))
    ax.loglog(x, y_sr, color="black", linewidth=2.0, label="Series reduction (k=2)")
    ax.loglog(x, y_c2, color="blue", linewidth=1.8, label="Cesaro (α = 2)")
    ax.loglog(x, y_c5, color="red", linewidth=1.8, label="Cesaro (α = 5)")

    ax.set_xlabel(r"$M\omega$", fontsize=14)
    ax.set_ylabel(r"$|f(\gamma=0,\theta_o=\phi_o=\pi/6)|/r_{\rm SL}$", fontsize=14)
    ax.set_title(r"$a=0$", fontsize=18)
    ax.grid(True, which="both", linestyle=":", alpha=0.5)
    ax.legend(loc="lower left", fontsize=12)

    out_path = os.path.join(os.path.dirname(__file__), "results", "Fig4_a0_reproduction.png")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    print(f"[+] Done. Figure saved to {out_path}")


if __name__ == "__main__":
    main()

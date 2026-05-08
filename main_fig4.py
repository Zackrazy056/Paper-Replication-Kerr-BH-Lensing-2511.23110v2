import os

import matplotlib.pyplot as plt
import json

from src.physics_fig4 import Fig4Engine


def _next_iter_path(results_dir, stem, ext):
    i = 1
    while True:
        cand = os.path.join(results_dir, f"{stem}_{i:02d}.{ext}")
        if not os.path.exists(cand):
            return cand
        i += 1


def main():
    print("=== Fig.4 Reproduction (a = 0, 0.8, 0.99, 0.999) ===")

    data_file = os.path.join(os.path.dirname(__file__), "data", "fig4_data.json")
    engine = Fig4Engine(data_file)

    l_max = 60
    r_sl = 100
    a_list = [0.0, 0.8, 0.99, 0.999]
    titles = [r"$a=0$", r"$a=0.8M$", r"$a=0.99M$", r"$a=0.999M$"]
    sr_dump = []

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True, sharey=True)
    axes = axes.ravel()

    for i, (a_val, title) in enumerate(zip(a_list, titles)):
        ax = axes[i]
        x, y_sr, y_c2, y_c5 = engine.build_curves(a_val=a_val, l_max=l_max, r_sl=r_sl)
        if not x:
            print(f"[!] No valid curve points for a={a_val}")
            continue
        sr_dump.append(
            {
                "a": a_val,
                "omega": x,
                "sr_k2": y_sr,
            }
        )

        ax.loglog(x, y_sr, color="black", linewidth=2.0, label="Series reduction (k=2)")
        ax.loglog(x, y_c2, color="blue", linewidth=1.8, label="Cesaro (alpha = 2)")
        ax.loglog(x, y_c5, color="red", linewidth=1.8, label="Cesaro (alpha = 5)")

        ax.set_title(title, fontsize=16)
        ax.grid(True, which="both", linestyle=":", alpha=0.5)
        if i in (2, 3):
            ax.set_xlabel(r"$M\omega$", fontsize=14)
        if i in (0, 2):
            ax.set_ylabel(r"$|f(\gamma=0,\theta_o=\phi_o=\pi/6)|/r_{\rm SL}$", fontsize=13)
        if i == 0:
            ax.legend(loc="lower left", fontsize=11)

    plt.tight_layout()
    results_dir = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(results_dir, exist_ok=True)
    out_path = _next_iter_path(results_dir, "Fig4_reproduction", "png")
    plt.savefig(out_path, dpi=300)
    sr_out = _next_iter_path(results_dir, "fig4_sr_curves", "json")
    with open(sr_out, "w", encoding="utf-8") as f:
        json.dump(sr_dump, f, ensure_ascii=False)
    print(f"[+] Done. Figure saved to {out_path}")
    print(f"[+] SR-only curves saved to {sr_out}")


if __name__ == "__main__":
    main()

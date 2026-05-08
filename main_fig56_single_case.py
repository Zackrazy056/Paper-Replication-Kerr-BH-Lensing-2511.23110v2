import math
import os

import matplotlib.pyplot as plt
import numpy as np

from src.physics_fig56 import Fig56Engine


def main() -> None:
    root = os.path.dirname(__file__)
    data_files = [
        os.path.join(root, "data", "fig4_data.json"),
        os.path.join(root, "data", "fig56_patch_data.json"),
        os.path.join(root, "data", "fig56_patch_a_minus099_theta30_hi.json"),
        os.path.join(root, "data", "fig56_sr_spherical_a0_theta30.json"),
    ]
    engine = Fig56Engine(data_files)

    theta_deg = 30.0
    a_val = 0.0
    freq = np.linspace(180.0, 640.0, 1800)

    wf = engine.build_waveforms(
        a_val=a_val,
        theta_deg=theta_deg,
        freq_hz=freq,
        l_max=60,
        k_sr=2,
        use_g_if_available=True,
        stabilize_f=False,
        require_sr_f=True,
    )

    psd = engine._aLIGO_psd_approx(freq)
    mm_plus = engine.mismatch(wf["h_dir_plus"], wf["h_obs_plus"], freq, psd)
    mm_cross = engine.mismatch(wf["h_dir_cross"], wf["h_obs_cross"], freq, psd)

    out_dir = os.path.join(root, "results")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "Fig6_theta30_a0_only.png")

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    ax.semilogy(freq, np.abs(wf["h_dir_plus"]), color="black", lw=1.9, label=r"$|h^{\rm dir}_{+}|$")
    ax.semilogy(freq, np.abs(wf["h_dir_cross"]), color="black", lw=1.1, ls="--", label=r"$|h^{\rm dir}_{\times}|$")
    ax.semilogy(freq, np.abs(wf["h_obs_plus"]), color="tab:blue", lw=1.9, label=r"$|h^{\rm obs}_{+}|$")
    ax.semilogy(freq, np.abs(wf["h_obs_cross"]), color="tab:red", lw=1.6, label=r"$|h^{\rm obs}_{\times}|$")

    ax.set_xlim(180, 640)
    ax.set_ylim(1e-8, 1e-2)
    ax.grid(True, which="both", ls=":", alpha=0.45)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel(r"$|h_{+,\times}(f)|$")
    ax.set_title(r"$\theta_o=30^\circ,\ a=0$")
    ax.legend(loc="lower left", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)

    print(f"[+] Saved: {out_path}")
    print(f"[*] mismatch_plus  = {mm_plus:.6f}")
    print(f"[*] mismatch_cross = {mm_cross:.6f}")


if __name__ == "__main__":
    main()

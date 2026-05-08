import csv
import json
import math
import os

import matplotlib.pyplot as plt
import numpy as np

from src.physics_fig56 import Fig56Engine


def _env_bool(name: str, default: bool) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}


def _env_int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    if raw is None or raw.strip() == "":
        return default
    return int(raw)


def _env_float(name: str, default: float) -> float:
    raw = os.environ.get(name)
    if raw is None or raw.strip() == "":
        return default
    return float(raw)


def draw_fig5(out_path: str) -> None:
    fig, ax = plt.subplots(1, 1, figsize=(8, 4.8))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")

    src = np.array([1.2, 3.0])
    lens = np.array([4.8, 3.0])
    obs = np.array([9.0, 4.6])

    ax.scatter(*src, s=180, c="royalblue", zorder=3)
    ax.scatter(*lens, s=260, c="black", zorder=3)
    ax.scatter(*obs, s=180, c="darkorange", zorder=3)
    ax.text(src[0] - 0.25, src[1] + 0.35, "S", fontsize=14, weight="bold")
    ax.text(lens[0] - 0.15, lens[1] + 0.4, "L", fontsize=14, weight="bold", color="white")
    ax.text(obs[0] + 0.1, obs[1] + 0.2, "O", fontsize=14, weight="bold")

    ax.plot([src[0], obs[0]], [src[1], obs[1]], color="gray", lw=1.8, ls="--", label="direct")
    ax.plot([src[0], lens[0]], [src[1], lens[1]], color="crimson", lw=2.2)
    ax.plot([lens[0], obs[0]], [lens[1], obs[1]], color="crimson", lw=2.2, label="scattered")

    ax.annotate(
        "",
        xy=lens + np.array([0, 1.1]),
        xytext=lens + np.array([0, 0.2]),
        arrowprops=dict(arrowstyle="-|>", lw=2.0, color="black"),
    )
    ax.text(lens[0] + 0.18, lens[1] + 0.95, r"$\vec{S}_{\rm BH}$", fontsize=12)

    theta = math.degrees(math.atan2(obs[1] - lens[1], obs[0] - lens[0]))
    arc = np.linspace(0, math.radians(theta), 80)
    r_arc = 1.2
    ax.plot(lens[0] + r_arc * np.cos(arc), lens[1] + r_arc * np.sin(arc), c="green", lw=1.5)
    ax.text(lens[0] + 1.05, lens[1] + 0.45, r"$\theta_o$", color="green", fontsize=12)

    ax.text(2.6, 2.55, r"$r_{\rm SL}$", color="crimson", fontsize=11)
    ax.text(6.6, 3.8, r"$r_{\rm LO}$", color="crimson", fontsize=11)
    ax.text(4.5, 4.95, r"$r_{\rm SO}$ (direct)", color="gray", fontsize=10)
    ax.set_title("Fig.5 Geometry (S-L-O setup)", fontsize=14)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def draw_fig6(
    engine: Fig56Engine,
    out_path: str,
    l_max: int,
    k_sr: int,
    use_g_if_available: bool,
    prefer_sr_f: bool,
    require_sr_f: bool,
) -> None:
    freq = np.linspace(20.0, 620.0, 2400)
    strain_scale = _env_float("FIG6_STRAIN_SCALE", 1e-20)
    configs = [
        (30.0, 0.0, r"$\theta_o=30^\circ,\ a=0$"),
        (30.0, 0.999, r"$\theta_o=30^\circ,\ a=0.999M$"),
        (60.0, 0.0, r"$\theta_o=60^\circ,\ a=0$"),
        (60.0, 0.999, r"$\theta_o=60^\circ,\ a=0.999M$"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    axes = axes.ravel()

    for ax, (theta_deg, a_val, title) in zip(axes, configs):
        wf = engine.build_waveforms(
            a_val=a_val,
            theta_deg=theta_deg,
            freq_hz=freq,
            l_max=l_max,
            k_sr=k_sr,
            use_g_if_available=use_g_if_available,
            prefer_sr_f=prefer_sr_f,
            require_sr_f=require_sr_f,
        )
        ax.semilogy(
            freq,
            strain_scale * np.abs(wf["h_dir_plus"]),
            color="tab:blue",
            lw=1.35,
            label=r"$|h^{\rm dir}_{+}|$",
        )
        ax.semilogy(
            freq,
            strain_scale * np.abs(wf["h_dir_cross"]),
            color="tab:orange",
            lw=1.35,
            ls="--",
            label=r"$|h^{\rm dir}_{\times}|$",
        )
        ax.semilogy(
            freq,
            strain_scale * np.abs(wf["h_obs_plus"]),
            color="tab:green",
            lw=1.35,
            label=r"$|h^{\rm obs}_{+}|$",
        )
        ax.semilogy(
            freq,
            strain_scale * np.abs(wf["h_obs_cross"]),
            color="tab:red",
            lw=1.35,
            ls="--",
            label=r"$|h^{\rm obs}_{\times}|$",
        )

        ax.set_title(title, fontsize=12)
        ax.set_xlim(20, 620)
        ax.set_ylim(1e-28, 1e-22)
        ax.set_xticks([200, 400, 600])
        ax.grid(True, which="both", ls=":", lw=0.45, alpha=0.55)

    axes[2].set_xlabel("Frequency [Hz]", fontsize=12)
    axes[3].set_xlabel("Frequency [Hz]", fontsize=12)
    axes[0].set_ylabel(r"$|h_{+,\times}(\nu)|$", fontsize=12)
    axes[2].set_ylabel(r"$|h_{+,\times}(\nu)|$", fontsize=12)
    axes[0].legend(loc="upper right", fontsize=9, frameon=True)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def export_table1(
    engine: Fig56Engine,
    out_csv: str,
    out_json: str,
    l_max: int,
    k_sr: int,
    n_freq: int,
    use_g_if_available: bool,
    prefer_sr_f: bool,
    require_sr_f: bool,
) -> None:
    a_list = [0.0, 0.8, 0.99, 0.999, -0.99]
    theta_list = [30.0, 60.0]
    rows = engine.build_table1(
        a_list=a_list,
        theta_list_deg=theta_list,
        f_low=20.0,
        f_high=640.0,
        n_freq=n_freq,
        l_max=l_max,
        k_sr=k_sr,
        use_g_if_available=use_g_if_available,
        prefer_sr_f=prefer_sr_f,
        require_sr_f=require_sr_f,
    )

    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(rows, f, ensure_ascii=False, indent=2)

    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=["theta_deg", "a", "mismatch_plus", "mismatch_cross"]
        )
        writer.writeheader()
        writer.writerows(rows)

    print("[*] Table 1 reproduction rows:")
    grouped = {}
    for row in rows:
        grouped.setdefault(int(row["theta_deg"]), []).append(row)
    for theta in sorted(grouped):
        print(f"    theta_o = {theta} deg")
        for row in sorted(grouped[theta], key=lambda r: r["a"]):
            print(
                "      a={:+.3f}  M_plus={:.4f}  M_cross={:.4f}".format(
                    row["a"], row["mismatch_plus"], row["mismatch_cross"]
                )
            )


def main() -> None:
    root = os.path.dirname(__file__)
    data_files = [
        os.path.join(root, "data", "fig4_data.json"),
        os.path.join(root, "data", "fig56_patch_data.json"),
        os.path.join(root, "data", "fig56_patch_a_minus099_theta30_hi.json"),
        os.path.join(root, "data", "fig56_sr_spherical.json"),
    ]
    engine = Fig56Engine(data_files)

    l_max = _env_int("FIG56_LMAX", 60)
    k_sr = _env_int("FIG56_K_SR", 2)
    table_n_freq = _env_int("FIG56_TABLE_N_FREQ", 2400)
    use_g_if_available = _env_bool("FIG56_USE_G", True)
    require_sr_f = _env_bool("FIG56_REQUIRE_SR_F", False)
    prefer_sr_f = _env_bool("FIG56_PREFER_SR_F", False) or require_sr_f

    print("=== Fig.5/Fig.6/Table I Reproduction ===")
    print(
        "[*] Settings: "
        f"l_max={l_max}, k_sr={k_sr}, table_n_freq={table_n_freq}, "
        f"use_g_if_available={use_g_if_available}, prefer_sr_f={prefer_sr_f}, "
        f"require_sr_f={require_sr_f}"
    )
    if not prefer_sr_f:
        print(
            "[*] FIG56_PREFER_SR_F=0: using the legacy modal path for f across all "
            "spin cases so cached configurations rerun consistently."
        )
    elif not require_sr_f:
        print(
            "[*] FIG56_PREFER_SR_F=1 and FIG56_REQUIRE_SR_F=0: spherical SR rows "
            "will be used when present, with legacy fallback for missing rows."
        )

    required = [
        (0.0, math.pi / 6),
        (0.999, math.pi / 6),
        (0.0, math.pi / 3),
        (0.999, math.pi / 3),
        (0.8, math.pi / 3),
        (0.99, math.pi / 3),
        (-0.99, math.pi / 6),
        (-0.99, math.pi / 3),
    ]
    missing = [(a, th) for (a, th) in required if not engine.has_combo(a, th)]
    if missing:
        print("[!] Missing required scattering data combinations:")
        for a_val, th in missing:
            print(f"    a={a_val:+.3f}, theta={math.degrees(th):.1f} deg")
        print("[!] Run: wolframscript -file mathematica/generate_fig56_patch_data.wls")
        return

    out_dir = os.path.join(root, "results")
    os.makedirs(out_dir, exist_ok=True)

    fig5_path = os.path.join(out_dir, "Fig5_geometry_reproduction.png")
    fig6_path = os.path.join(out_dir, "Fig6_reproduction.png")
    table_csv = os.path.join(out_dir, "Table1_reproduction.csv")
    table_json = os.path.join(out_dir, "Table1_reproduction.json")

    try:
        draw_fig5(fig5_path)
        draw_fig6(
            engine,
            fig6_path,
            l_max=l_max,
            k_sr=k_sr,
            use_g_if_available=use_g_if_available,
            prefer_sr_f=prefer_sr_f,
            require_sr_f=require_sr_f,
        )
        export_table1(
            engine,
            table_csv,
            table_json,
            l_max=l_max,
            k_sr=k_sr,
            n_freq=table_n_freq,
            use_g_if_available=use_g_if_available,
            prefer_sr_f=prefer_sr_f,
            require_sr_f=require_sr_f,
        )
    except ValueError as exc:
        if require_sr_f and "Missing SR spherical f" in str(exc):
            print(f"[!] {exc}")
            print("[!] Run `python reproduce.py sr-audit` to inspect missing Fig.56 SR rows.")
            print("[!] After activating Wolfram Engine, rerun mathematica/generate_fig56_sr_spherical.wls.")
            raise SystemExit(1) from None
        raise

    print(f"[+] Saved: {fig5_path}")
    print(f"[+] Saved: {fig6_path}")
    print(f"[+] Saved: {table_csv}")
    print(f"[+] Saved: {table_json}")


if __name__ == "__main__":
    main()

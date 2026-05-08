"""Paper-standard Fig.6 waveform plot.

This driver uses PyCBC IMRPhenomD waveforms for the direct component and the
incident lens-side waveform, then applies the cached Teukolsky scattering
factors. It is intentionally separate from ``main_fig56.py``, whose Fig.6 path
uses a lightweight analytic waveform for local reproduction.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from main_table1_paper_standard import (
    M1_REF,
    M2_REF,
    apply_scattering,
    build_engine,
    require_optional_deps,
    waveform_fd,
)


ROOT = Path(__file__).resolve().parent


def _plot_abs(values: np.ndarray) -> np.ndarray:
    amp = np.abs(values)
    return np.where(amp > 0.0, amp, np.nan)


def build_panel_waveforms(args: argparse.Namespace) -> list[dict[str, object]]:
    deps = require_optional_deps()
    engine = build_engine()
    freq = np.linspace(float(args.f_low), float(args.f_high), int(args.n_freq))
    panels: list[dict[str, object]] = []

    for theta_deg, a_val in [(30.0, 0.0), (30.0, 0.999), (60.0, 0.0), (60.0, 0.999)]:
        iota_obs = math.pi - math.radians(theta_deg)
        h_dir_plus, h_dir_cross = waveform_fd(
            deps,
            freq,
            M1_REF,
            M2_REF,
            0.0,
            0.0,
            iota_obs,
            distance_mpc=float(args.distance_mpc),
        )
        h_lens_plus, _ = waveform_fd(
            deps,
            freq,
            M1_REF,
            M2_REF,
            0.0,
            0.0,
            0.0,
            distance_mpc=float(args.distance_mpc),
        )
        # Match the +2-helicity convention used by the scattering model.
        h_lens_cross = 1j * h_lens_plus

        scattered = apply_scattering(
            engine,
            freq,
            h_dir_plus,
            h_dir_cross,
            h_lens_plus,
            h_lens_cross,
            a_val=a_val,
            theta_deg=theta_deg,
            l_max=int(args.l_max),
            require_sr_f=bool(args.strict_sr),
        )
        panel = {
            "theta_deg": theta_deg,
            "a": a_val,
            "freq_hz": freq,
            "h_dir_plus": h_dir_plus,
            "h_dir_cross": h_dir_cross,
            "h_obs_plus": scattered["h_obs_plus"],
            "h_obs_cross": scattered["h_obs_cross"],
        }
        panels.append(panel)
    return panels


def write_summary(panels: list[dict[str, object]], out_path: Path) -> None:
    rows: list[dict[str, float]] = []
    for panel in panels:
        row: dict[str, float] = {
            "theta_deg": float(panel["theta_deg"]),
            "a": float(panel["a"]),
        }
        for key in ["h_dir_plus", "h_dir_cross", "h_obs_plus", "h_obs_cross"]:
            amp = _plot_abs(np.asarray(panel[key]))
            finite = amp[np.isfinite(amp)]
            row[f"{key}_min"] = float(np.min(finite)) if finite.size else float("nan")
            row[f"{key}_max"] = float(np.max(finite)) if finite.size else float("nan")
        rows.append(row)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(rows, indent=2), encoding="utf-8")


def draw_fig6_paper(panels: list[dict[str, object]], out_path: Path, args: argparse.Namespace) -> None:
    titles = [
        r"$\theta_o=30^\circ,\ a=0$",
        r"$\theta_o=30^\circ,\ a=0.999M$",
        r"$\theta_o=60^\circ,\ a=0$",
        r"$\theta_o=60^\circ,\ a=0.999M$",
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), sharex=True, sharey=True)
    axes = axes.ravel()

    for ax, panel, title in zip(axes, panels, titles):
        freq = np.asarray(panel["freq_hz"])
        ax.semilogy(freq, _plot_abs(np.asarray(panel["h_dir_plus"])), color="tab:blue", lw=1.25, label=r"$|h^{\rm dir}_{+}|$")
        ax.semilogy(
            freq,
            _plot_abs(np.asarray(panel["h_dir_cross"])),
            color="tab:orange",
            lw=1.25,
            ls="--",
            label=r"$|h^{\rm dir}_{\times}|$",
        )
        ax.semilogy(freq, _plot_abs(np.asarray(panel["h_obs_plus"])), color="tab:green", lw=1.25, label=r"$|h^{\rm obs}_{+}|$")
        ax.semilogy(
            freq,
            _plot_abs(np.asarray(panel["h_obs_cross"])),
            color="tab:red",
            lw=1.25,
            ls="--",
            label=r"$|h^{\rm obs}_{\times}|$",
        )
        ax.set_title(title, fontsize=12)
        ax.set_xlim(float(args.f_low), float(args.f_high))
        ax.set_ylim(float(args.y_min), float(args.y_max))
        ax.set_xticks([200, 400, 600])
        ax.tick_params(axis="x", which="both", labelbottom=True)
        ax.grid(True, which="both", ls=":", lw=0.45, alpha=0.55)

    axes[2].set_xlabel("Frequency [Hz]", fontsize=12)
    axes[3].set_xlabel("Frequency [Hz]", fontsize=12)
    axes[0].set_ylabel(r"$|h_{+,\times}(\nu)|$", fontsize=12)
    axes[2].set_ylabel(r"$|h_{+,\times}(\nu)|$", fontsize=12)
    axes[0].legend(loc="upper right", fontsize=9, frameon=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-deps", action="store_true", help="Only check PyCBC/LALSuite/SciPy imports.")
    parser.add_argument("--strict-sr", action="store_true", help="Require strict spherical-SR rows.")
    parser.add_argument("--n-freq", type=int, default=2400)
    parser.add_argument("--f-low", type=float, default=20.0)
    parser.add_argument("--f-high", type=float, default=620.0)
    parser.add_argument("--y-min", type=float, default=1e-28)
    parser.add_argument("--y-max", type=float, default=1e-22)
    parser.add_argument("--l-max", type=int, default=60)
    parser.add_argument("--distance-mpc", type=float, default=402.6)
    parser.add_argument("--out", default=str(ROOT / "results" / "Fig6_paper_standard.png"))
    parser.add_argument("--summary", default=str(ROOT / "results" / "Fig6_paper_standard_summary.json"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        if args.check_deps:
            require_optional_deps()
            print("Paper-standard dependencies are available.")
            return
        panels = build_panel_waveforms(args)
    except RuntimeError as exc:
        print(f"[!] {exc}")
        raise SystemExit(2) from None

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    draw_fig6_paper(panels, out_path, args)
    summary_path = Path(args.summary)
    write_summary(panels, summary_path)
    print(f"[+] Saved: {out_path}")
    print(f"[+] Saved: {summary_path}")


if __name__ == "__main__":
    main()

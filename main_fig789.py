"""Paper-standard Fig.7-Fig.9 contour sweep.

The target quantities are:

- Fig.7: optimized mismatch, 1 - FF, over theta_o x lens spin.
- Fig.8: optimal-SNR ratio, ||h_obs|| / ||h_dir||.
- Fig.9: required SNR from ln B ~= rho^2 / 2 * (1 - FF^2).

The current cached spherical-SR scattering data is exact at theta_o=30,60 deg.
By default this driver uses those strict-SR anchor curves and linearly
interpolates/extrapolates the complex scattering amplitudes in theta. Use
``--theta-mode exact`` only when the scattering cache has all requested theta
values.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from main_table1_paper_standard import (
    LENS_MASS_MSUN,
    M1_REF,
    M2_REF,
    apply_scattering,
    build_engine,
    chirp_q_to_masses,
    masses_to_chirp_q,
    require_optional_deps,
    waveform_fd,
)
from src.physics_fig56 import M_SUN_SEC, Fig56Engine


ROOT = Path(__file__).resolve().parent
LN_B_THRESHOLD = 10.0


@dataclass
class GridResult:
    theta_deg: float
    a: float
    ff_plus: float
    ff_cross: float
    mismatch_plus: float
    mismatch_cross: float
    snr_ratio_plus: float
    snr_ratio_cross: float
    required_snr_plus: float
    required_snr_cross: float
    rho_dir_plus: float
    rho_dir_cross: float
    rho_obs_plus: float
    rho_obs_cross: float
    optimized: bool
    theta_mode: str


def _parse_float_list(raw: str | None, default: list[float]) -> list[float]:
    if raw is None or raw.strip() == "":
        return default
    return [float(part.strip()) for part in raw.split(",") if part.strip()]


def _theta_grid(args: argparse.Namespace) -> np.ndarray:
    if args.theta:
        return np.array(args.theta, dtype=float)
    return np.linspace(float(args.theta_min), float(args.theta_max), int(args.theta_count))


def _inner_norm(engine: Fig56Engine, h: np.ndarray, freq_hz: np.ndarray, psd: np.ndarray) -> float:
    inv_psd = 1.0 / np.maximum(psd, 1e-60)
    df = np.gradient(freq_hz)
    return float(math.sqrt(max(4.0 * np.sum((np.abs(h) ** 2) * inv_psd * df), 0.0)))


def _required_snr(ff: float, ln_b_threshold: float = LN_B_THRESHOLD) -> float:
    ff = float(np.clip(ff, 0.0, 1.0))
    denom = 1.0 - ff * ff
    if denom <= 1e-15:
        return float("inf")
    return float(math.sqrt((2.0 * float(ln_b_threshold)) / denom))


def _complex_theta_interp(
    theta_deg: float,
    anchors_deg: np.ndarray,
    values_by_anchor: list[np.ndarray],
) -> np.ndarray:
    if len(anchors_deg) == 1:
        return values_by_anchor[0].copy()
    order = np.argsort(anchors_deg)
    x = anchors_deg[order]
    vals = [values_by_anchor[int(i)] for i in order]
    x0 = float(x[0])
    x1 = float(x[-1])
    if x1 == x0:
        return vals[0].copy()
    weight = (float(theta_deg) - x0) / (x1 - x0)
    return vals[0] + weight * (vals[-1] - vals[0])


def anchor_scattering(
    engine: Fig56Engine,
    freq_hz: np.ndarray,
    a_val: float,
    theta_deg: float,
    l_max: int,
    strict_sr: bool,
    anchor_thetas: list[float],
) -> tuple[np.ndarray, np.ndarray]:
    f_anchors: list[np.ndarray] = []
    g_anchors: list[np.ndarray] = []
    for anchor in anchor_thetas:
        omega_grid, f_grid, g_grid = engine.build_scatter_grid(
            a_val=a_val,
            theta_obs=math.radians(anchor),
            l_max=l_max,
            k_sr=2,
            use_g_if_available=True,
            stabilize_f=False,
            prefer_sr_f=True,
            require_sr_f=strict_sr,
        )
        if len(omega_grid) == 0:
            raise ValueError(f"No anchor scattering data for a={a_val:+.3f}, theta={anchor:.3f} deg.")
        mw = engine.mw_from_hz(freq_hz, lens_mass_msun=LENS_MASS_MSUN)
        f_anchors.append(engine._interp_complex(omega_grid, f_grid, mw))
        g_anchors.append(engine._interp_complex(omega_grid, g_grid, mw))

    anchors = np.array(anchor_thetas, dtype=float)
    f_interp = _complex_theta_interp(theta_deg, anchors, f_anchors)
    g_interp = _complex_theta_interp(theta_deg, anchors, g_anchors)
    return f_interp, g_interp


def build_observed_waveforms(
    deps,
    engine: Fig56Engine,
    freq_hz: np.ndarray,
    a_val: float,
    theta_deg: float,
    l_max: int,
    strict_sr: bool,
    theta_mode: str,
    anchor_thetas: list[float],
) -> dict[str, np.ndarray]:
    iota_obs = math.pi - math.radians(theta_deg)
    h_dir_plus, h_dir_cross = waveform_fd(deps, freq_hz, M1_REF, M2_REF, 0.0, 0.0, iota_obs)
    h_lens_plus, _ = waveform_fd(deps, freq_hz, M1_REF, M2_REF, 0.0, 0.0, 0.0)
    h_lens_cross = 1j * h_lens_plus

    if theta_mode == "exact":
        scattered = apply_scattering(
            engine,
            freq_hz,
            h_dir_plus,
            h_dir_cross,
            h_lens_plus,
            h_lens_cross,
            a_val=a_val,
            theta_deg=theta_deg,
            l_max=l_max,
            require_sr_f=strict_sr,
        )
        return {
            "h_dir_plus": h_dir_plus,
            "h_dir_cross": h_dir_cross,
            **scattered,
        }

    f_interp, g_interp = anchor_scattering(
        engine,
        freq_hz,
        a_val=a_val,
        theta_deg=theta_deg,
        l_max=l_max,
        strict_sr=strict_sr,
        anchor_thetas=anchor_thetas,
    )
    theta = math.radians(theta_deg)
    dt = 2.0 * 100.0 * (LENS_MASS_MSUN * M_SUN_SEC) * math.sin(theta / 2.0) ** 2
    phase_delay = np.exp(2j * math.pi * freq_hz * dt)
    h_sw_plus = ((f_interp + g_interp) / 100.0) * h_lens_plus
    h_sw_cross = ((f_interp - g_interp) / 100.0) * h_lens_cross
    return {
        "h_dir_plus": h_dir_plus,
        "h_dir_cross": h_dir_cross,
        "h_obs_plus": h_dir_plus + phase_delay * h_sw_plus,
        "h_obs_cross": h_dir_cross + phase_delay * h_sw_cross,
        "f_scatter": f_interp,
        "g_scatter": g_interp,
    }


def fitting_factor_with_params(
    deps,
    engine: Fig56Engine,
    target: np.ndarray,
    freq_hz: np.ndarray,
    psd: np.ndarray,
    theta_deg: float,
    pol: str,
    optimize: bool,
    maxiter: int,
    x0: np.ndarray | None,
    multistart: bool,
) -> tuple[float, np.ndarray]:
    chirp0, q0 = masses_to_chirp_q(M1_REF, M2_REF)
    baseline = np.array([chirp0, q0, 0.0, 0.0], dtype=float)
    iota = math.pi - math.radians(theta_deg)

    def template_from_params(params: np.ndarray) -> np.ndarray | None:
        chirp, q, spin1z, spin2z = params
        if not (20.0 <= chirp <= 40.0 and 0.4 <= q <= 1.0 and -0.99 <= spin1z <= 0.99 and -0.99 <= spin2z <= 0.99):
            return None
        mass1, mass2 = chirp_q_to_masses(chirp, q)
        hp, hc = waveform_fd(deps, freq_hz, mass1, mass2, spin1z, spin2z, iota)
        return hp if pol == "plus" else hc

    def match_for(params: np.ndarray) -> float:
        h_template = template_from_params(params)
        if h_template is None:
            return 0.0
        return 1.0 - engine.mismatch(h_template, target, freq_hz, psd)

    if not optimize:
        baseline_match = match_for(baseline)
        return float(np.clip(baseline_match, 0.0, 1.0)), baseline

    starts = [baseline]
    if x0 is not None:
        starts.append(np.array(x0, dtype=float))
    if multistart:
        starts.extend(
            [
                np.array([chirp0 * 0.98, q0, 0.0, 0.0], dtype=float),
                np.array([chirp0 * 1.02, q0, 0.0, 0.0], dtype=float),
                np.array([chirp0, q0 * 0.95, 0.0, 0.0], dtype=float),
                np.array([chirp0, min(q0 * 1.05, 1.0), 0.0, 0.0], dtype=float),
                np.array([chirp0, q0, 0.15, 0.15], dtype=float),
                np.array([chirp0, q0, -0.15, -0.15], dtype=float),
            ]
        )

    unique_starts: list[np.ndarray] = []
    for start in starts:
        if any(np.allclose(start, prev, rtol=0.0, atol=1e-12) for prev in unique_starts):
            continue
        unique_starts.append(start)

    best = -math.inf
    params = baseline.copy()
    for start in unique_starts:
        start_match = match_for(start)
        if start_match > best:
            best = start_match
            params = start.copy()

        def objective(candidate: np.ndarray) -> float:
            return -match_for(candidate)

        result = deps.minimize(
            objective,
            start,
            method="Nelder-Mead",
            options={"maxiter": int(maxiter), "xatol": 1e-3, "fatol": 1e-4},
        )
        result_match = -float(result.fun) if np.isfinite(result.fun) else 0.0
        if result_match > best:
            best = result_match
            params = np.array(result.x, dtype=float)
    return float(np.clip(best, 0.0, 1.0)), params


def run_grid(args: argparse.Namespace) -> list[GridResult]:
    deps = require_optional_deps()
    engine = build_engine()
    freq = np.linspace(float(args.f_low), float(args.f_high), int(args.n_freq))
    psd = engine._aLIGO_psd_approx(freq)
    theta_vals = _theta_grid(args)
    spin_vals = np.array(args.spin, dtype=float)
    anchor_thetas = _parse_float_list(args.anchor_theta, [30.0, 60.0])

    rows: list[GridResult] = []
    last_params: dict[str, np.ndarray | None] = {"plus": None, "cross": None}
    n_total = len(theta_vals) * len(spin_vals)
    idx = 0
    for theta_deg in theta_vals:
        for a_val in spin_vals:
            idx += 1
            print(f"[*] ({idx}/{n_total}) theta={theta_deg:.2f} deg, a={a_val:+.3f}", flush=True)
            wf = build_observed_waveforms(
                deps,
                engine,
                freq,
                a_val=float(a_val),
                theta_deg=float(theta_deg),
                l_max=int(args.l_max),
                strict_sr=bool(args.strict_sr),
                theta_mode=args.theta_mode,
                anchor_thetas=anchor_thetas,
            )
            ff_plus, last_params["plus"] = fitting_factor_with_params(
                deps,
                engine,
                wf["h_obs_plus"],
                freq,
                psd,
                float(theta_deg),
                "plus",
                bool(args.optimize),
                int(args.maxiter),
                last_params["plus"],
                bool(args.multistart),
            )
            ff_cross, last_params["cross"] = fitting_factor_with_params(
                deps,
                engine,
                wf["h_obs_cross"],
                freq,
                psd,
                float(theta_deg),
                "cross",
                bool(args.optimize),
                int(args.maxiter),
                last_params["cross"],
                bool(args.multistart),
            )

            rho_dir_plus = _inner_norm(engine, wf["h_dir_plus"], freq, psd)
            rho_dir_cross = _inner_norm(engine, wf["h_dir_cross"], freq, psd)
            rho_obs_plus = _inner_norm(engine, wf["h_obs_plus"], freq, psd)
            rho_obs_cross = _inner_norm(engine, wf["h_obs_cross"], freq, psd)
            ratio_plus = rho_obs_plus / max(rho_dir_plus, 1e-30)
            ratio_cross = rho_obs_cross / max(rho_dir_cross, 1e-30)

            row = GridResult(
                theta_deg=float(theta_deg),
                a=float(a_val),
                ff_plus=float(ff_plus),
                ff_cross=float(ff_cross),
                mismatch_plus=float(1.0 - ff_plus),
                mismatch_cross=float(1.0 - ff_cross),
                snr_ratio_plus=float(ratio_plus),
                snr_ratio_cross=float(ratio_cross),
                required_snr_plus=_required_snr(ff_plus, args.ln_b_threshold),
                required_snr_cross=_required_snr(ff_cross, args.ln_b_threshold),
                rho_dir_plus=float(rho_dir_plus),
                rho_dir_cross=float(rho_dir_cross),
                rho_obs_plus=float(rho_obs_plus),
                rho_obs_cross=float(rho_obs_cross),
                optimized=bool(args.optimize),
                theta_mode=str(args.theta_mode),
            )
            rows.append(row)
            print(
                "    M+=(%.5g), Mx=(%.5g), rho_ratio+=(%.5g), rho_ratio_x=(%.5g), rho_req+=(%.4g), rho_req_x=(%.4g)"
                % (
                    row.mismatch_plus,
                    row.mismatch_cross,
                    row.snr_ratio_plus,
                    row.snr_ratio_cross,
                    row.required_snr_plus,
                    row.required_snr_cross,
                ),
                flush=True,
            )
    return rows


def _rows_to_dicts(rows: list[GridResult]) -> list[dict[str, float | bool | str]]:
    return [row.__dict__.copy() for row in rows]


def _grid_array(rows: list[GridResult], field: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    theta = np.array(sorted({row.theta_deg for row in rows}), dtype=float)
    spin = np.array(sorted({row.a for row in rows}), dtype=float)
    z = np.full((len(theta), len(spin)), np.nan, dtype=float)
    theta_index = {v: i for i, v in enumerate(theta)}
    spin_index = {v: i for i, v in enumerate(spin)}
    for row in rows:
        z[theta_index[row.theta_deg], spin_index[row.a]] = float(getattr(row, field))
    return spin, theta, z


def _plot_two_panel(
    rows: list[GridResult],
    field_plus: str,
    field_cross: str,
    title: str,
    cbar_label: str,
    out_path: Path,
    cap: float | None = None,
) -> None:
    spin, theta, z_plus = _grid_array(rows, field_plus)
    _, _, z_cross = _grid_array(rows, field_cross)
    if cap is not None:
        z_plus = np.clip(z_plus, 0.0, cap)
        z_cross = np.clip(z_cross, 0.0, cap)

    vals = np.concatenate([z_plus[np.isfinite(z_plus)], z_cross[np.isfinite(z_cross)]])
    if vals.size == 0:
        raise ValueError(f"No finite values for {title}.")
    vmin = float(np.nanmin(vals))
    vmax = float(np.nanmax(vals))
    if math.isclose(vmin, vmax):
        vmax = vmin + 1e-6
    levels = np.linspace(vmin, vmax, 18)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), sharex=True, sharey=True)
    spin_mesh, theta_mesh = np.meshgrid(spin, theta)
    use_scatter = len(spin) < 2 or len(theta) < 2
    mappable = None
    for ax, z, panel_title in [
        (axes[0], z_plus, r"$+$ polarization"),
        (axes[1], z_cross, r"$\times$ polarization"),
    ]:
        if use_scatter:
            finite = np.isfinite(z)
            mappable = ax.scatter(
                spin_mesh[finite],
                theta_mesh[finite],
                c=z[finite],
                s=95,
                marker="s",
                cmap="viridis",
                vmin=vmin,
                vmax=vmax,
            )
            x_pad = 0.05 if len(spin) == 1 else 0.0
            y_pad = 1.0 if len(theta) == 1 else 0.0
            ax.set_xlim(float(spin[0]) - x_pad, float(spin[-1]) + x_pad)
            ax.set_ylim(float(theta[0]) - y_pad, float(theta[-1]) + y_pad)
        else:
            mappable = ax.contourf(spin, theta, z, levels=levels, cmap="viridis", extend="both")
            ax.contour(spin, theta, z, levels=levels[::3], colors="white", linewidths=0.45, alpha=0.65)
            ax.scatter(spin_mesh, theta_mesh, s=10, c="black", alpha=0.45)
        ax.set_title(panel_title)
        ax.set_xlabel(r"Lens spin $a/M$")
        ax.grid(True, ls=":", lw=0.4, alpha=0.5)
    axes[0].set_ylabel(r"$\theta_o$ [deg]")
    fig.suptitle(title)
    cbar = fig.colorbar(mappable, ax=axes.ravel().tolist(), shrink=0.94)
    cbar.set_label(cbar_label)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_outputs(rows: list[GridResult], args: argparse.Namespace) -> None:
    out_base = Path(args.out)
    out_base.parent.mkdir(parents=True, exist_ok=True)
    data = _rows_to_dicts(rows)
    out_base.write_text(json.dumps(data, indent=2), encoding="utf-8")
    csv_path = out_base.with_suffix(".csv")
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(data[0].keys()))
        writer.writeheader()
        writer.writerows(data)

    prefix = out_base.with_suffix("")
    fig7 = prefix.with_name(prefix.name + "_Fig7_mismatch_contour.png")
    fig8 = prefix.with_name(prefix.name + "_Fig8_snr_ratio_contour.png")
    fig9 = prefix.with_name(prefix.name + "_Fig9_required_snr_contour.png")
    _plot_two_panel(
        rows,
        "mismatch_plus",
        "mismatch_cross",
        "Fig.7 reproduction: optimized mismatch",
        r"$1-\mathrm{FF}$",
        fig7,
    )
    _plot_two_panel(
        rows,
        "snr_ratio_plus",
        "snr_ratio_cross",
        "Fig.8 reproduction: optimal-SNR ratio",
        r"$\rho_{\rm obs}/\rho_{\rm dir}$",
        fig8,
        cap=args.snr_ratio_cap,
    )
    _plot_two_panel(
        rows,
        "required_snr_plus",
        "required_snr_cross",
        "Fig.9 reproduction: required SNR",
        r"$\rho_{\rm req}$",
        fig9,
        cap=args.required_snr_cap,
    )

    print(f"[+] Saved: {out_base}")
    print(f"[+] Saved: {csv_path}")
    print(f"[+] Saved: {fig7}")
    print(f"[+] Saved: {fig8}")
    print(f"[+] Saved: {fig9}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-deps", action="store_true", help="Only check PyCBC/LALSuite/SciPy imports.")
    parser.add_argument("--strict-sr", action="store_true", help="Require strict spherical-SR anchor rows.")
    parser.add_argument("--optimize", action="store_true", help="Run Nelder-Mead FF optimization.")
    parser.add_argument(
        "--multistart",
        action="store_true",
        help="Try deterministic extra optimizer starts to reduce local-minimum and grid-order sensitivity.",
    )
    parser.add_argument("--theta-mode", choices=["anchor", "exact"], default="anchor")
    parser.add_argument("--anchor-theta", default="30,60", help="Comma-separated anchor theta values in degrees.")
    parser.add_argument("--theta", type=float, nargs="+", help="Explicit theta grid in degrees.")
    parser.add_argument("--theta-min", type=float, default=30.0)
    parser.add_argument("--theta-max", type=float, default=89.0)
    parser.add_argument("--theta-count", type=int, default=9)
    parser.add_argument("--spin", type=float, nargs="+", default=[-0.99, 0.0, 0.8, 0.99, 0.999])
    parser.add_argument("--n-freq", type=int, default=1000)
    parser.add_argument("--f-low", type=float, default=20.0)
    parser.add_argument("--f-high", type=float, default=640.0)
    parser.add_argument("--l-max", type=int, default=60)
    parser.add_argument("--maxiter", type=int, default=100)
    parser.add_argument("--ln-b-threshold", type=float, default=LN_B_THRESHOLD)
    parser.add_argument("--snr-ratio-cap", type=float, default=2.0)
    parser.add_argument("--required-snr-cap", type=float, default=300.0)
    parser.add_argument("--out", default=str(ROOT / "results" / "Fig789_grid.json"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        if args.check_deps:
            require_optional_deps()
            print("Paper-standard dependencies are available.")
            return
        rows = run_grid(args)
        save_outputs(rows, args)
    except RuntimeError as exc:
        print(f"[!] {exc}")
        raise SystemExit(2) from None


if __name__ == "__main__":
    main()

"""Paper-standard Table I driver scaffold.

This script is intentionally stricter than main_fig56.py. It requires PyCBC,
LALSuite/lalsimulation, and SciPy so that waveform generation can use
IMRPhenomD and fitting-factor optimization can use scipy.optimize, matching the
paper methodology more closely.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from src.physics_fig56 import Fig56Engine


ROOT = Path(__file__).resolve().parent
M1_REF = 36.9
M2_REF = 32.8
LENS_MASS_MSUN = 100.0


@dataclass
class OptionalDeps:
    scipy: object
    minimize: object
    get_fd_waveform: object


def require_optional_deps() -> OptionalDeps:
    missing = [
        name
        for name in ("scipy", "pycbc", "lal", "lalsimulation")
        if importlib.util.find_spec(name) is None
    ]
    if missing:
        raise RuntimeError(
            "Missing paper-standard dependencies: "
            + ", ".join(missing)
            + ". On Windows, LALSuite wheels may be unavailable; use a Linux/WSL/cluster "
            "environment if PyCBC cannot be installed locally."
        )

    from scipy.optimize import minimize
    from pycbc.waveform import get_fd_waveform

    import scipy

    return OptionalDeps(scipy=scipy, minimize=minimize, get_fd_waveform=get_fd_waveform)


def build_engine() -> Fig56Engine:
    data_files = [
        ROOT / "data" / "fig4_data.json",
        ROOT / "data" / "fig56_patch_data.json",
        ROOT / "data" / "fig56_patch_a_minus099_theta30_hi.json",
        ROOT / "data" / "fig56_sr_spherical.json",
        ROOT / "data" / "fig789_patch_data.json",
        ROOT / "data" / "fig789_sr_spherical.json",
    ]
    return Fig56Engine([str(path) for path in data_files])


def chirp_q_to_masses(chirp_mass: float, q: float) -> tuple[float, float]:
    q = float(q)
    eta = q / (1.0 + q) ** 2
    total = chirp_mass / (eta ** (3.0 / 5.0))
    m1 = total / (1.0 + q)
    m2 = q * m1
    return float(max(m1, m2)), float(min(m1, m2))


def masses_to_chirp_q(m1: float, m2: float) -> tuple[float, float]:
    m_high = max(m1, m2)
    m_low = min(m1, m2)
    total = m_high + m_low
    eta = (m_high * m_low) / total**2
    chirp = total * eta ** (3.0 / 5.0)
    q = m_low / m_high
    return float(chirp), float(q)


def waveform_fd(
    deps: OptionalDeps,
    freq_hz: np.ndarray,
    mass1: float,
    mass2: float,
    spin1z: float,
    spin2z: float,
    inclination: float,
    distance_mpc: float = 402.6,
) -> tuple[np.ndarray, np.ndarray]:
    delta_f = float(freq_hz[1] - freq_hz[0])
    f_final = float(freq_hz[-1])
    hp, hc = deps.get_fd_waveform(
        approximant="IMRPhenomD",
        mass1=float(mass1),
        mass2=float(mass2),
        spin1z=float(spin1z),
        spin2z=float(spin2z),
        inclination=float(inclination),
        distance=float(distance_mpc),
        delta_f=delta_f,
        f_lower=float(freq_hz[0]),
        f_final=f_final,
    )
    src_freq = np.asarray(hp.sample_frequencies, dtype=float)
    hp_arr = np.asarray(hp, dtype=complex)
    hc_arr = np.asarray(hc, dtype=complex)
    hp_interp = np.interp(freq_hz, src_freq, hp_arr.real) + 1j * np.interp(freq_hz, src_freq, hp_arr.imag)
    hc_interp = np.interp(freq_hz, src_freq, hc_arr.real) + 1j * np.interp(freq_hz, src_freq, hc_arr.imag)
    return hp_interp, hc_interp


def apply_scattering(
    engine: Fig56Engine,
    freq_hz: np.ndarray,
    h_dir_plus: np.ndarray,
    h_dir_cross: np.ndarray,
    h_lens_plus: np.ndarray,
    h_lens_cross: np.ndarray,
    a_val: float,
    theta_deg: float,
    l_max: int,
    require_sr_f: bool,
) -> dict[str, np.ndarray]:
    theta = math.radians(theta_deg)
    omega_grid, f_grid, g_grid = engine.build_scatter_grid(
        a_val=a_val,
        theta_obs=theta,
        l_max=l_max,
        k_sr=2,
        use_g_if_available=True,
        stabilize_f=False,
        prefer_sr_f=True,
        require_sr_f=require_sr_f,
    )
    mw = engine.mw_from_hz(freq_hz, lens_mass_msun=LENS_MASS_MSUN)
    f_interp = engine._interp_complex(omega_grid, f_grid, mw)
    g_interp = engine._interp_complex(omega_grid, g_grid, mw)
    dt = 2.0 * 100.0 * (LENS_MASS_MSUN * 4.925490947e-6) * math.sin(theta / 2.0) ** 2
    phase_delay = np.exp(2j * math.pi * freq_hz * dt)

    h_sw_plus = ((f_interp + g_interp) / 100.0) * h_lens_plus
    h_sw_cross = ((f_interp - g_interp) / 100.0) * h_lens_cross
    return {
        "h_obs_plus": h_dir_plus + phase_delay * h_sw_plus,
        "h_obs_cross": h_dir_cross + phase_delay * h_sw_cross,
        "f_scatter": f_interp,
        "g_scatter": g_interp,
    }


def fitting_factor(
    deps: OptionalDeps,
    engine: Fig56Engine,
    target: np.ndarray,
    freq_hz: np.ndarray,
    psd: np.ndarray,
    theta_deg: float,
    pol: str,
    optimize: bool,
) -> float:
    chirp0, q0 = masses_to_chirp_q(M1_REF, M2_REF)
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

    x0 = np.array([chirp0, q0, 0.0, 0.0], dtype=float)
    if not optimize:
        return match_for(x0)

    def objective(params: np.ndarray) -> float:
        return -match_for(params)

    result = deps.minimize(
        objective,
        x0,
        method="Nelder-Mead",
        options={"maxiter": 120, "xatol": 1e-3, "fatol": 1e-4},
    )
    best = max(match_for(x0), -float(result.fun) if result.success or np.isfinite(result.fun) else 0.0)
    return float(min(max(best, 0.0), 1.0))


def run_table(args: argparse.Namespace) -> list[dict[str, float]]:
    deps = require_optional_deps()
    engine = build_engine()
    freq = np.linspace(args.f_low, args.f_high, args.n_freq)
    psd = engine._aLIGO_psd_approx(freq)

    rows: list[dict[str, float]] = []
    for theta_deg in args.theta:
        iota_obs = math.pi - math.radians(theta_deg)
        h_dir_plus, h_dir_cross = waveform_fd(deps, freq, M1_REF, M2_REF, 0.0, 0.0, iota_obs)
        h_lens_plus, h_lens_cross = waveform_fd(deps, freq, M1_REF, M2_REF, 0.0, 0.0, 0.0)
        # Enforce the positive-helicity incident-wave convention used in the scattering model.
        h_lens_cross = 1j * h_lens_plus

        for a_val in args.spin:
            scattered = apply_scattering(
                engine,
                freq,
                h_dir_plus,
                h_dir_cross,
                h_lens_plus,
                h_lens_cross,
                a_val=a_val,
                theta_deg=theta_deg,
                l_max=args.l_max,
                require_sr_f=args.strict_sr,
            )
            ff_plus = fitting_factor(
                deps,
                engine,
                scattered["h_obs_plus"],
                freq,
                psd,
                theta_deg,
                "plus",
                optimize=args.optimize,
            )
            ff_cross = fitting_factor(
                deps,
                engine,
                scattered["h_obs_cross"],
                freq,
                psd,
                theta_deg,
                "cross",
                optimize=args.optimize,
            )
            row = {
                "theta_deg": float(theta_deg),
                "a": float(a_val),
                "mismatch_plus": float(1.0 - ff_plus),
                "mismatch_cross": float(1.0 - ff_cross),
                "optimized": bool(args.optimize),
            }
            rows.append(row)
            print(
                "theta={theta_deg:.1f} a={a:+.3f} M_plus={mp:.5g} M_cross={mc:.5g}".format(
                    theta_deg=theta_deg,
                    a=a_val,
                    mp=row["mismatch_plus"],
                    mc=row["mismatch_cross"],
                )
            )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-deps", action="store_true", help="Only check optional dependencies.")
    parser.add_argument("--strict-sr", action="store_true", help="Require complete spherical SR data.")
    parser.add_argument("--optimize", action="store_true", help="Run Nelder-Mead fitting-factor optimization.")
    parser.add_argument("--n-freq", type=int, default=1600)
    parser.add_argument("--f-low", type=float, default=20.0)
    parser.add_argument("--f-high", type=float, default=640.0)
    parser.add_argument("--l-max", type=int, default=60)
    parser.add_argument("--theta", type=float, nargs="+", default=[30.0, 60.0])
    parser.add_argument("--spin", type=float, nargs="+", default=[0.0, 0.8, 0.99, 0.999, -0.99])
    parser.add_argument(
        "--out",
        default=str(ROOT / "results" / "Table1_paper_standard.json"),
        help="Output JSON path.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        if args.check_deps:
            require_optional_deps()
            print("Paper-standard dependencies are available.")
            return

        rows = run_table(args)
    except RuntimeError as exc:
        print(f"[!] {exc}")
        raise SystemExit(2) from None

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    csv_path = out_path.with_suffix(".csv")
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["theta_deg", "a", "mismatch_plus", "mismatch_cross", "optimized"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"[+] Saved: {out_path}")
    print(f"[+] Saved: {csv_path}")


if __name__ == "__main__":
    main()

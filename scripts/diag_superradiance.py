import argparse
import sys
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from interface import nearest_available_omega
from physics import (
    calc_absorption_cross_section,
    get_absorption_mode_rows,
    omega_tilde,
    superradiant_limit,
)


def _fmt(x, digits=10):
    return mp.nstr(x, digits)


def _target_frequencies(m_omega_h):
    return [
        ("0.5*mOmegaH", mp.mpf("0.5") * m_omega_h),
        ("0.95*mOmegaH", mp.mpf("0.95") * m_omega_h),
        ("1.2*mOmegaH", mp.mpf("1.2") * m_omega_h),
    ]


def run_diag(a, m, l_max):
    a_mp = mp.mpf(str(a))
    m_omega_h = superradiant_limit(a_mp, m=m)

    print(f"\n=== a={a_mp}, m={m}, mOmega_H={_fmt(m_omega_h, 14)} ===")
    for label, omega_target in _target_frequencies(m_omega_h):
        omega_used = nearest_available_omega(a_mp, omega_target, m=m)
        tilde = omega_tilde(a_mp, omega_used, m=m)
        mode_data = get_absorption_mode_rows(
            a_mp,
            omega_used,
            m=m,
            l_min=2,
            l_max=l_max,
            parities=(-1, 1),
            allow_nearest_omega=False,
        )
        rows = sorted(mode_data["rows"], key=lambda r: (r["l"], r["P"]))

        sigma_raw = calc_absorption_cross_section(
            a_mp,
            omega_used,
            P=None,
            m=m,
            l_max=l_max,
            clip_non_superradiant_negative=False,
        )
        sigma_clip = calc_absorption_cross_section(
            a_mp,
            omega_used,
            P=None,
            m=m,
            l_max=l_max,
            clip_non_superradiant_negative=True,
        )

        gt1_rows = [r for r in rows if r["sfac_modsq"] > 1]
        sr_gt1_rows = [r for r in rows if r["omega_tilde"] < 0 and r["sfac_modsq"] > 1]

        print(
            f"\n[{label}] target={_fmt(omega_target, 12)}, "
            f"used={_fmt(omega_used, 12)}, omega_tilde={_fmt(tilde, 12)}"
        )
        print(
            f"rows={len(rows)}, |S|^2>1: total={len(gt1_rows)}, "
            f"superradiant={len(sr_gt1_rows)}, "
            f"sigma_raw={sigma_raw:.8g}, sigma_clip={sigma_clip:.8g}"
        )
        print("l  P  |S|^2         1-|S|^2       K_l^P")
        for row in rows:
            print(
                f"{row['l']:2d} {row['P']:2d} "
                f"{_fmt(row['sfac_modsq'], 10):>12} "
                f"{_fmt(row['one_minus_modsq'], 10):>12} "
                f"{_fmt(row['K_lP'], 10):>12}"
            )


def main():
    parser = argparse.ArgumentParser(
        description="Diagnose superradiance mode-by-mode for Fig.2 left panel."
    )
    parser.add_argument(
        "--a-list",
        type=float,
        nargs="+",
        default=[0.8, 0.99, 0.999],
        help="Spin values to diagnose.",
    )
    parser.add_argument("--m", type=int, default=2, help="Azimuthal mode m.")
    parser.add_argument("--l-max", type=int, default=10, help="Max l in diagnostic table.")
    parser.add_argument("--dps", type=int, default=80, help="mpmath decimal precision.")
    args = parser.parse_args()

    mp.mp.dps = max(mp.mp.dps, int(args.dps))
    for a in args.a_list:
        run_diag(a, args.m, args.l_max)


if __name__ == "__main__":
    main()

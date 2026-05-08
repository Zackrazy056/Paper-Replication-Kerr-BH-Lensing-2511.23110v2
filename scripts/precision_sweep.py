import argparse
import importlib
import sys
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))


def _fmt(x, digits=12):
    return mp.nstr(x, digits)


def _load_modules_at_dps(dps):
    mp.mp.dps = int(dps)
    interface = importlib.import_module("interface")
    physics = importlib.import_module("physics")
    interface = importlib.reload(interface)
    physics = importlib.reload(physics)
    return interface, physics


def main():
    parser = argparse.ArgumentParser(
        description="Precision sweep for |Sfac|^2 stability in Fig.2 left-panel chain."
    )
    parser.add_argument("--a", type=float, required=True, help="Spin parameter a/M.")
    parser.add_argument(
        "--omega",
        type=float,
        default=None,
        help="Target omega. If omitted, uses 0.5*m*Omega_H.",
    )
    parser.add_argument("--m", type=int, default=2, help="Azimuthal mode.")
    parser.add_argument("--l-min", type=int, default=2, help="Min l.")
    parser.add_argument("--l-max", type=int, default=6, help="Max l.")
    parser.add_argument(
        "--dps-list",
        type=int,
        nargs="+",
        default=[50, 80, 120, 160],
        help="Decimal precision list.",
    )
    parser.add_argument(
        "--nearest",
        action="store_true",
        help="Project omega to nearest tabulated grid point.",
    )
    args = parser.parse_args()

    mode_history = {}
    omega_history = {}

    for dps in args.dps_list:
        interface, physics = _load_modules_at_dps(dps)
        a_mp = mp.mpf(str(args.a))
        m_omega_h = physics.superradiant_limit(a_mp, m=args.m)
        omega_target = (
            mp.mpf(str(args.omega))
            if args.omega is not None
            else mp.mpf("0.5") * m_omega_h
        )
        omega_used = (
            interface.nearest_available_omega(a_mp, omega_target, m=args.m)
            if args.nearest
            else omega_target
        )
        omega_history[dps] = omega_used

        mode_data = physics.get_absorption_mode_rows(
            a_mp,
            omega_used,
            m=args.m,
            l_min=args.l_min,
            l_max=args.l_max,
            parities=(-1, 1),
            allow_nearest_omega=False,
        )
        rows = mode_data["rows"]

        print(
            f"\n[dps={dps}] a={a_mp}, target_omega={_fmt(omega_target)}, "
            f"used_omega={_fmt(omega_used)}, mOmega_H={_fmt(m_omega_h)}"
        )
        for row in sorted(rows, key=lambda r: (r["l"], r["P"])):
            key = (row["l"], row["P"])
            mode_history.setdefault(key, []).append((dps, row["sfac_modsq"]))
            print(
                f"l={row['l']:2d}, P={row['P']:2d}, "
                f"|S|^2={_fmt(row['sfac_modsq'], 14)}"
            )

    if not mode_history:
        print("\nNo mode rows were found for this (a, omega) setup.")
        return

    print("\n=== Relative Drift vs Highest Precision ===")
    for key in sorted(mode_history):
        seq = mode_history[key]
        seq_sorted = sorted(seq, key=lambda t: t[0])
        dps_ref, val_ref = seq_sorted[-1]
        print(f"mode l={key[0]}, P={key[1]}, ref_dps={dps_ref}, ref={_fmt(val_ref, 14)}")
        for dps, val in seq_sorted:
            if val_ref == 0:
                rel = mp.fabs(val - val_ref)
            else:
                rel = mp.fabs((val - val_ref) / val_ref)
            print(f"  dps={dps:3d}, rel_diff={_fmt(rel, 8)}")


if __name__ == "__main__":
    main()

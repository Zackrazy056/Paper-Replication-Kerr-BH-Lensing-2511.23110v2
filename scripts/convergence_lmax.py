import argparse
import sys
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from interface import get_available_omegas
from physics import calc_backward_scattering


def _fmt(x, d=8):
    return f"{float(x):.{d}f}"


def main():
    parser = argparse.ArgumentParser(
        description="Convergence test of backward scattering peak vs l_max."
    )
    parser.add_argument("--a", type=float, default=0.999)
    parser.add_argument("--w-min", type=float, default=0.7)
    parser.add_argument("--w-max", type=float, default=1.05)
    parser.add_argument(
        "--lmax-list", type=int, nargs="+", default=[40, 60, 80, 100, 120]
    )
    parser.add_argument("--dps", type=int, default=80)
    args = parser.parse_args()

    mp.mp.dps = max(mp.mp.dps, args.dps)
    a = mp.mpf(str(args.a))
    w_all = [mp.mpf(str(w)) for w in get_available_omegas(a, m=2)]
    w_grid = [w for w in w_all if args.w_min <= float(w) <= args.w_max]
    if not w_grid:
        print("No omega points found in the requested range.")
        return

    print(
        f"a={args.a}, omega_range=[{args.w_min}, {args.w_max}], "
        f"points={len(w_grid)}, data_lmax_available=25 (from table)"
    )

    for lmax in args.lmax_list:
        vals = [
            calc_backward_scattering(a, w, l_min=2, l_max=lmax, early_stop=False)
            for w in w_grid
        ]
        vmax = max(vals)
        w_peak = w_grid[vals.index(vmax)]
        truncated = " (saturated by data l<=25)" if lmax > 25 else ""
        print(
            f"lmax={lmax:3d} peak={_fmt(vmax, 10)} at Momega={_fmt(w_peak, 4)}{truncated}"
        )


if __name__ == "__main__":
    main()

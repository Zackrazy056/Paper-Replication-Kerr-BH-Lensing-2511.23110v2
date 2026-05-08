import argparse
import sys
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from interface import get_available_omegas, get_teukolsky_data_mock
from physics import calc_phase_shift, superradiant_limit


def sigma_variant(a, omega, l_max=25, weight_mode="re_s2"):
    s = mp.mpf("0")
    for l in range(2, l_max + 1):
        d = get_teukolsky_data_mock(l, 2, -2, a, omega)
        if mp.fabs(d["Re_C"]) < mp.mpf("1e-20"):
            continue

        phase = calc_phase_shift(l, a, omega, 1, d["Re_C"], d["B_ratio"])
        kernel = mp.re(1 - mp.fabs(phase) ** 2)

        if weight_mode == "re_s2":
            weight = mp.re(d["S_0"] ** 2)
        elif weight_mode == "abs_s2":
            weight = mp.fabs(d["S_0"]) ** 2
        else:
            raise ValueError("Unknown weight mode")

        s += weight * kernel
    return (4 * mp.pi**2 / omega**2) * s / mp.pi


def main():
    parser = argparse.ArgumentParser(
        description="Low-frequency sensitivity test: omega_min and S^2 vs |S|^2."
    )
    parser.add_argument("--a-list", type=float, nargs="+", default=[0.8, 0.99, 0.999])
    parser.add_argument("--w-min-target", type=float, default=0.02)
    parser.add_argument("--w-max", type=float, default=0.5)
    parser.add_argument("--l-max", type=int, default=25)
    parser.add_argument("--dps", type=int, default=80)
    args = parser.parse_args()

    mp.mp.dps = max(mp.mp.dps, args.dps)
    print(
        f"target omega_min={args.w_min_target}, omega_max={args.w_max}, "
        f"l_max={args.l_max}"
    )

    for a_val in args.a_list:
        a = mp.mpf(str(a_val))
        m_omega_h = superradiant_limit(a, m=2)
        omega_all = [mp.mpf(str(w)) for w in get_available_omegas(a, m=2)]
        omega_used = [w for w in omega_all if float(w) <= args.w_max]
        if not omega_used:
            print(f"\na={a_val}: no omega points <= {args.w_max}")
            continue

        w_min_avail = float(min(omega_used))
        coverage_note = ""
        if w_min_avail > args.w_min_target:
            coverage_note = (
                f" [coverage gap: requested {args.w_min_target}, available {w_min_avail}]"
            )

        diffs = []
        vals_re = []
        vals_abs = []
        for w in omega_used:
            v_re = sigma_variant(a, w, l_max=args.l_max, weight_mode="re_s2")
            v_abs = sigma_variant(a, w, l_max=args.l_max, weight_mode="abs_s2")
            vals_re.append(v_re)
            vals_abs.append(v_abs)
            diffs.append(mp.fabs(v_re - v_abs))

        v0_re = vals_re[0]
        v0_abs = vals_abs[0]
        print(
            f"\na={a_val}, mOmega_H={float(m_omega_h):.8f}, "
            f"omega_points={len(omega_used)}{coverage_note}"
        )
        print(
            f"  at omega_min_available={w_min_avail:.4f}: "
            f"re[S^2]={float(v0_re):.8f}, |S|^2={float(v0_abs):.8f}"
        )
        print(
            f"  max |re[S^2]-|S|^2| over range: {float(max(diffs)):.3e}, "
            f"min sigma(re[S^2])={float(min(vals_re)):.8f}"
        )


if __name__ == "__main__":
    main()

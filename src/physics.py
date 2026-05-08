import mpmath as mp

try:
    from .interface import get_teukolsky_data, get_teukolsky_data_mock, nearest_available_omega
except ImportError:
    from interface import get_teukolsky_data, get_teukolsky_data_mock, nearest_available_omega

mp.mp.dps = max(mp.mp.dps, 50)
M = mp.mpf("1.0")
_RE_C_ZERO = mp.mpf("1e-20")
_ABS_TOL = mp.mpf("1e-8")
_CLIP_TOL = mp.mpf("1e-10")


def kerr_r_plus(a):
    a_mp = mp.mpf(a)
    disc = M**2 - a_mp**2
    if disc < 0:
        disc = mp.mpf("0.0")
    return M + mp.sqrt(disc)


def kerr_omega_h(a):
    a_mp = mp.mpf(a)
    if a_mp == 0:
        return mp.mpf("0.0")
    r_plus = kerr_r_plus(a_mp)
    if r_plus == 0:
        return mp.mpf("0.0")
    return a_mp / (2 * M * r_plus)


def superradiant_limit(a, m=2):
    return mp.mpf(m) * kerr_omega_h(a)


def omega_tilde(a, omega, m=2):
    return mp.mpf(omega) - superradiant_limit(a, m=m)


def calc_phase_shift(l, a, omega, P, re_c, b_ratio):
    """Complex S-factor eta*exp(2 i delta) under current convention."""
    term1 = mp.power(-1, l + 1)
    term2 = re_c + 12 * mp.j * M * omega * P * mp.power(-1, l)
    term3 = 16 * mp.power(omega, 4)
    phase_shift = term1 * (term2 / term3) * b_ratio
    return phase_shift


def _build_absorption_mode_row(l, a, omega, m, parity, data):
    phase = calc_phase_shift(l, a, omega, parity, data["Re_C"], data["B_ratio"])
    sfac_modsq = mp.power(mp.fabs(phase), 2)
    one_minus_modsq = mp.re(mp.mpf("1.0") - sfac_modsq)
    s0_sq = mp.re(mp.power(data["S_0"], 2))
    term = s0_sq * one_minus_modsq

    return {
        "l": int(l),
        "m": int(m),
        "P": int(parity),
        "omega": mp.mpf(omega),
        "omega_tilde": omega_tilde(a, omega, m=m),
        "S0_sq": s0_sq,
        "sfac_modsq": sfac_modsq,
        "one_minus_modsq": one_minus_modsq,
        "term": term,
        "K_lP": -term,
    }


def get_absorption_mode_rows(
    a,
    omega,
    m=2,
    l_min=2,
    l_max=10,
    parities=(-1, 1),
    allow_nearest_omega=False,
):
    """
    Return per-(l,P) absorption diagnostics for a given (a, omega).
    """
    omega_mp = mp.mpf(omega)
    if allow_nearest_omega:
        omega_mp = nearest_available_omega(a, omega_mp, m=m)

    rows = []
    for l in range(int(l_min), int(l_max) + 1):
        data = get_teukolsky_data(
            l, m, -2, a, omega_mp, allow_nearest_omega=allow_nearest_omega
        )
        if mp.fabs(data["Re_C"]) < _RE_C_ZERO:
            continue

        for parity in parities:
            rows.append(_build_absorption_mode_row(l, a, omega_mp, m, parity, data))

    return {"omega_used": omega_mp, "rows": rows}


def calc_absorption_cross_section(
    a,
    omega,
    P=1,
    m=2,
    l_min=2,
    l_max=100,
    clip_non_superradiant_negative=True,
    return_breakdown=False,
):
    """
    Absorption cross-section for the Fig.2 left panel chain.
    - P=+1/-1 computes a single parity branch.
    - P=None computes both parities and sums them.
    """
    omega_mp = mp.mpf(omega)
    if omega_mp == 0:
        if return_breakdown:
            return {"sigma_a": 0.0, "rows": [], "omega_used": omega_mp}
        return 0.0

    if P is None:
        parities = (-1, 1)
    else:
        parities = (int(P),)

    sum_sigma = mp.mpf("0.0")
    consecutive_small = 0
    rows = []

    for l in range(int(l_min), int(l_max) + 1):
        data = get_teukolsky_data_mock(l, m, -2, a, omega_mp)
        if mp.fabs(data["Re_C"]) < _RE_C_ZERO:
            break

        l_term_total = mp.mpf("0.0")
        for parity in parities:
            row = _build_absorption_mode_row(l, a, omega_mp, m, parity, data)
            one_minus_modsq = row["one_minus_modsq"]
            clipped = False

            if (
                clip_non_superradiant_negative
                and row["omega_tilde"] >= 0
                and one_minus_modsq < -_CLIP_TOL
            ):
                one_minus_modsq = mp.mpf("0.0")
                clipped = True
            elif mp.fabs(one_minus_modsq) < _CLIP_TOL:
                one_minus_modsq = mp.mpf("0.0")

            term_used = row["S0_sq"] * one_minus_modsq
            l_term_total += term_used
            sum_sigma += term_used

            if return_breakdown:
                row["one_minus_modsq_used"] = one_minus_modsq
                row["term_used"] = term_used
                row["clipped_non_sr"] = clipped
                rows.append(row)

        if mp.fabs(l_term_total) < _ABS_TOL:
            consecutive_small += 1
            if consecutive_small >= 3:
                break
        else:
            consecutive_small = 0

    sigma_a = (4 * mp.power(mp.pi, 2) / mp.power(omega_mp, 2)) * sum_sigma
    sigma_a_f = float(mp.re(sigma_a))

    if return_breakdown:
        return {
            "sigma_a": sigma_a_f,
            "omega_used": omega_mp,
            "m_omega_h": superradiant_limit(a, m=m),
            "rows": rows,
        }
    return sigma_a_f


def calc_backward_scattering(a, omega, l_min=2, l_max=100, early_stop=True):
    """Backward-scattering branch (m=+2) with configurable l truncation."""
    sum_f = mp.mpc("0.0")
    sum_g = mp.mpc("0.0")
    for l in range(int(l_min), int(l_max) + 1):
        data = get_teukolsky_data_mock(l, 2, -2, a, omega)
        if mp.fabs(data["Re_C"]) < _RE_C_ZERO:
            continue

        bracket_f = (
            mp.power(-1, l + 1) * data["Re_C"] / (16 * mp.power(omega, 4)) * data["B_ratio"]
            - 1
        )
        term_f = data["S_0"] * data["S_pi"] * bracket_f

        bracket_g = (
            mp.power(-1, l + 1)
            * 3
            * mp.j
            * M
            * omega
            / (4 * mp.power(omega, 4))
            * data["B_ratio"]
        )
        term_g = mp.power(-1, l) * data["S_0"] * data["S_0"] * mp.conj(bracket_g)

        sum_f += term_f
        sum_g += term_g

        if early_stop and mp.fabs(term_f) < _ABS_TOL and mp.fabs(term_g) < _ABS_TOL:
            break

    coeff = 2 * mp.pi / (mp.j * omega)
    dsigma_domega = mp.power(mp.fabs(coeff * sum_f), 2) + mp.power(mp.fabs(coeff * sum_g), 2)
    return float(dsigma_domega)


if __name__ == "__main__":
    pass

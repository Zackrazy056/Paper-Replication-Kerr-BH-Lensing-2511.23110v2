import json
import os
import re

import mpmath as mp

mp.mp.dps = max(mp.mp.dps, 50)

DATA_FILE = os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "data", "teukolsky_data_real.json"
)

SILENT = os.environ.get("KERR_GW_INTERFACE_SILENT", "0") == "1"

ZERO_RECORD = {
    "Re_C": mp.mpf("0.0"),
    "B_ratio": mp.mpc("0.0", "0.0"),
    "S_0": mp.mpf("0.0"),
    "S_pi": mp.mpf("0.0"),
}

real_data_dict = {}


def _log(message):
    if not SILENT:
        print(message)


def _is_nonfinite_token(value):
    s = str(value)
    return any(
        token in s
        for token in ("ComplexInfinity", "DirectedInfinity", "Indeterminate", "Infinity")
    )


def _sanitize_mathematica_string(value):
    s = str(value).strip().replace(" ", "")
    if _is_nonfinite_token(s):
        return None
    s = re.sub(r"`[-+0-9.]*", "", s)
    s = s.replace("*^", "e")
    return s


def parse_mathematica_number(num_str):
    if isinstance(num_str, mp.mpf):
        return num_str
    if isinstance(num_str, (int, float)):
        return mp.mpf(num_str)

    cleaned = _sanitize_mathematica_string(num_str)
    if cleaned is None or cleaned == "":
        return mp.mpf("0.0")
    if "I" in cleaned or "i" in cleaned:
        return mp.mpf("0.0")
    return mp.mpf(cleaned)


def _load_real_data():
    global real_data_dict
    try:
        with open(DATA_FILE, "r", encoding="utf-8") as f:
            raw = json.load(f)
    except FileNotFoundError:
        _log(f"[!] data file not found: {DATA_FILE}")
        return

    parse_fail_count = 0
    nonfinite_rows = 0
    invalid_flag_rows = 0
    rows_without_m = 0

    nonfinite_fields = {
        "Re_C": 0,
        "B_ratio_re": 0,
        "B_ratio_im": 0,
        "S_0": 0,
        "S_pi": 0,
    }

    invalid_by_m = {}
    loaded_by_m = {}

    for item in raw:
        m_val = item.get("m", None)
        if m_val is None:
            rows_without_m += 1
            m_int = 2
        else:
            m_int = int(m_val)

        key = (
            int(item["l"]),
            m_int,
            f"{float(item['a']):.3f}",
            f"{float(item['omega']):.4f}",
        )

        if "valid" in item and not bool(item["valid"]):
            invalid_flag_rows += 1
            invalid_by_m[m_int] = invalid_by_m.get(m_int, 0) + 1

        row_has_nonfinite = False
        for field in nonfinite_fields:
            if _is_nonfinite_token(item.get(field, "")):
                nonfinite_fields[field] += 1
                row_has_nonfinite = True

        if row_has_nonfinite:
            nonfinite_rows += 1
            real_data_dict[key] = ZERO_RECORD.copy()
            loaded_by_m[m_int] = loaded_by_m.get(m_int, 0) + 1
            continue

        try:
            b_re = parse_mathematica_number(item["B_ratio_re"])
            b_im = parse_mathematica_number(item["B_ratio_im"])
            real_data_dict[key] = {
                "Re_C": parse_mathematica_number(item["Re_C"]),
                "B_ratio": mp.mpc(b_re, b_im),
                "S_0": parse_mathematica_number(item["S_0"]),
                "S_pi": parse_mathematica_number(item["S_pi"]),
            }
        except Exception:
            parse_fail_count += 1
            real_data_dict[key] = ZERO_RECORD.copy()

        loaded_by_m[m_int] = loaded_by_m.get(m_int, 0) + 1

    _log(
        "[*] loaded Teukolsky table: "
        f"{len(real_data_dict)} records, parse-fallback {parse_fail_count}."
    )
    _log(
        "[*] data quality: "
        f"invalid-flag rows={invalid_flag_rows}, nonfinite rows={nonfinite_rows}, "
        f"rows-without-m={rows_without_m}."
    )
    _log(f"[*] rows by m: loaded={loaded_by_m}, invalid={invalid_by_m}")

    if nonfinite_rows > 0:
        _log(
            "[WARN] nonfinite Mathematica rows detected and zeroed. "
            f"field counts={nonfinite_fields}"
        )


_load_real_data()


def get_available_omegas(a, m=2):
    """
    Return sorted omega grid (as Python floats) available for a given (a, m).
    """
    a_key = f"{float(a):.3f}"
    m_int = int(m)

    omega_set = {
        float(key[3])
        for key in real_data_dict
        if key[2] == a_key and key[1] == m_int
    }
    return sorted(omega_set)


def nearest_available_omega(a, omega, m=2):
    """
    Return nearest tabulated omega as mpmath number.
    If no grid exists, return the input omega unchanged.
    """
    omega_grid = get_available_omegas(a, m=m)
    if not omega_grid:
        return mp.mpf(omega)

    omega_f = float(omega)
    nearest = min(omega_grid, key=lambda w: abs(w - omega_f))
    return mp.mpf(str(nearest))


def get_teukolsky_data(l, m, s, a, omega, allow_nearest_omega=False):
    """
    Teukolsky table lookup.
    - Exact match by default (backward-compatible behavior).
    - Optional nearest-omega fallback for diagnostics.
    """
    key = (int(l), int(m), f"{float(a):.3f}", f"{float(omega):.4f}")
    if key in real_data_dict:
        return real_data_dict[key]

    legacy_key = (int(l), 2, f"{float(a):.3f}", f"{float(omega):.4f}")
    if legacy_key in real_data_dict:
        return real_data_dict[legacy_key]

    if allow_nearest_omega:
        omega_nearest = nearest_available_omega(a, omega, m=m)
        key_nearest = (int(l), int(m), f"{float(a):.3f}", f"{float(omega_nearest):.4f}")
        if key_nearest in real_data_dict:
            return real_data_dict[key_nearest]

        legacy_key_nearest = (
            int(l),
            2,
            f"{float(a):.3f}",
            f"{float(omega_nearest):.4f}",
        )
        if legacy_key_nearest in real_data_dict:
            return real_data_dict[legacy_key_nearest]

    return ZERO_RECORD


def get_teukolsky_data_mock(l, m, s, a, omega):
    """
    Real lookup interface. Name kept for compatibility with earlier stage code.
    """
    return get_teukolsky_data(l, m, s, a, omega, allow_nearest_omega=False)

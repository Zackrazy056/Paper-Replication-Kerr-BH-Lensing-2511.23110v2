import json
import os
import re

import mpmath as mp

mp.mp.dps = 80
M = mp.mpf("1.0")


class Fig4Engine:
    def __init__(self, data_file, sr_file=None):
        self.data = {}
        self.omega_by_a = {}
        self.f_sr_lookup = {}
        self.f_sr_series = {}
        self._load_data(data_file)
        if sr_file is not None:
            self._load_sr_data(sr_file)

    def _clean_str(self, s):
        s = str(s).strip().replace(" ", "").replace("*^", "e")
        s = re.sub(r"`[-+0-9.]*", "", s)
        if not s:
            return mp.mpf("0.0")
        if any(tok in s for tok in ("Infinity", "Indeterminate", "ComplexInfinity")):
            return mp.mpf("0.0")
        try:
            return mp.mpf(s)
        except Exception:
            return mp.mpf("0.0")

    def _load_data(self, path):
        paths = path if isinstance(path, list) else [path]
        for one_path in paths:
            if not os.path.exists(one_path):
                raise FileNotFoundError(f"Missing data file: {one_path}")
            raw = json.load(open(one_path, "r", encoding="utf-8"))
            for item in raw:
                if "valid" in item and not bool(item["valid"]):
                    continue
                a_key = self._a_key(item["a"])
                omega_key = self._omega_key(item["omega"])
                key = (a_key, omega_key, int(item["l"]))
                self.data[key] = {
                    "Re_C": self._clean_str(item["Re_C"]),
                    "B_ratio": mp.mpc(self._clean_str(item["B_ratio_re"]), self._clean_str(item["B_ratio_im"])),
                    "S_gamma": self._clean_str(item["S_gamma"]),
                    "S_theta_obs": self._clean_str(item["S_theta_obs"]),
                }
                if a_key not in self.omega_by_a:
                    self.omega_by_a[a_key] = set()
                self.omega_by_a[a_key].add(omega_key)

    def _load_sr_data(self, path):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing spherical SR data file: {path}")
        raw = json.load(open(path, "r", encoding="utf-8"))
        for item in raw:
            if "a" not in item or "omega" not in item:
                continue
            a_key = self._a_key(item["a"])
            omega_key = self._omega_key(item["omega"])
            f_val = mp.mpc(
                self._clean_str(item.get("f_sr_theta_re", 0)),
                self._clean_str(item.get("f_sr_theta_im", 0)),
            )
            self.f_sr_lookup[(a_key, omega_key)] = f_val
            if a_key not in self.f_sr_series:
                self.f_sr_series[a_key] = []
            self.f_sr_series[a_key].append((float(item["omega"]), f_val))
        for a_key in list(self.f_sr_series.keys()):
            self.f_sr_series[a_key].sort(key=lambda t: t[0])

    def _a_key(self, a_val):
        return f"{float(a_val):.3f}"

    def _omega_key(self, omega_val):
        return f"{float(omega_val):.12g}"

    def _omega_values(self, a_val):
        a_key = self._a_key(a_val)
        vals = sorted((float(w) for w in self.omega_by_a.get(a_key, [])))
        return vals

    def _coeff_f_map(self, a_val, omega, l_max=60):
        a_key = self._a_key(a_val)
        omega_key = self._omega_key(omega)
        omega_mp = mp.mpf(str(omega))
        pref = 2 * mp.pi / (mp.j * omega_mp)

        coeff = {}
        basis = {}
        for l in range(2, int(l_max) + 1):
            mode = self.data.get((a_key, omega_key, l))
            if mode is None:
                continue
            term1 = mp.power(-1, l + 1)
            bracket_f = term1 * mode["Re_C"] / (16 * mp.power(omega_mp, 4)) * mode["B_ratio"] - mp.mpf("1.0")
            coeff[l] = pref * mode["S_gamma"] * bracket_f
            basis[l] = mode["S_theta_obs"]
        return coeff, basis

    def _cesaro_weight(self, l, l_max, alpha):
        if alpha == 0:
            return mp.mpf("1.0")
        return mp.binomial(l_max, l) / mp.binomial(l_max + alpha, l)

    def _lookup_sr_f(self, a_val, omega):
        a_key = self._a_key(a_val)
        omega_key = self._omega_key(omega)
        exact = self.f_sr_lookup.get((a_key, omega_key), None)
        if exact is not None:
            return exact
        series = self.f_sr_series.get(a_key, [])
        if not series:
            return None
        x = [p[0] for p in series]
        yre = [float(mp.re(p[1])) for p in series]
        yim = [float(mp.im(p[1])) for p in series]
        oq = min(max(float(omega), x[0]), x[-1])
        re = float(mp.mpf(str(mp.mpf(str(__import__("numpy").interp(oq, x, yre))))))
        im = float(mp.mpf(str(mp.mpf(str(__import__("numpy").interp(oq, x, yim))))))
        return mp.mpc(re, im)

    def _series_reduce_once(self, coeff_map, m=2):
        reduced = {}
        if not coeff_map:
            return reduced
        l_values = sorted(coeff_map.keys())
        l_min = l_values[0]
        l_max = l_values[-1]
        for l in l_values:
            if l >= l_max:
                continue
            fl = coeff_map.get(l, mp.mpc("0.0"))
            fp = coeff_map.get(l + 1, mp.mpc("0.0"))
            fm = coeff_map.get(l - 1, mp.mpc("0.0"))

            term_diag = mp.mpf("1.0") - (mp.mpf("2.0") * m) / (l * (l + 1))

            a_num = (l - 1) * (l + 3) * (l - m + 1) * (l + m + 1)
            a_den = 4 * l * (l + 2) + 3
            a_l = mp.sqrt(mp.mpf(a_num) / mp.mpf(a_den)) if a_num >= 0 else mp.mpf("0.0")

            b_num = (l * l - m * m) * (l * l - 4)
            b_den = 4 * l * l - 1
            b_l = mp.sqrt(mp.mpf(b_num) / mp.mpf(b_den)) if b_num >= 0 else mp.mpf("0.0")

            # Legacy SR path (kept for backward compatibility with old Fig.4 workflow).
            reduced[l] = term_diag * fl - (a_l / (l + 1)) * fp - (b_l / l) * fm
        return reduced

    def f_abs_series_reduction(self, a_val, omega, l_max=60, k_sr=2, r_sl=100):
        coeff, basis = self._coeff_f_map(a_val, omega, l_max=l_max)
        if not coeff:
            return None

        coeff_k = dict(coeff)
        for _ in range(int(k_sr)):
            coeff_k = self._series_reduce_once(coeff_k, m=2)
            if not coeff_k:
                break
        if not coeff_k:
            return None

        theta_obs = mp.pi / 6
        denom = mp.power(1 - mp.cos(theta_obs), int(k_sr))

        sum_f = mp.mpc("0.0")
        for l, v in coeff_k.items():
            if l > l_max:
                continue
            sum_f += v * basis.get(l, mp.mpf("0.0"))
        return float(mp.fabs(sum_f / denom) / mp.mpf(str(r_sl)))

    def f_abs_cesaro(self, a_val, omega, l_max=60, alpha=2, r_sl=100):
        coeff, basis = self._coeff_f_map(a_val, omega, l_max=l_max)
        if not coeff:
            return None
        sum_f = mp.mpc("0.0")
        for l, v in coeff.items():
            if l > l_max:
                continue
            w = self._cesaro_weight(l, l_max, int(alpha))
            sum_f += w * v * basis.get(l, mp.mpf("0.0"))
        return float(mp.fabs(sum_f) / mp.mpf(str(r_sl)))

    def build_curves(self, a_val, l_max=60, k_sr=2, r_sl=100):
        omega_values = self._omega_values(a_val)
        w_out = []
        y_sr = []
        y_c2 = []
        y_c5 = []
        for omega in omega_values:
            sr = self.f_abs_series_reduction(a_val, omega, l_max=l_max, k_sr=k_sr, r_sl=r_sl)
            c2 = self.f_abs_cesaro(a_val, omega, l_max=l_max, alpha=2, r_sl=r_sl)
            c5 = self.f_abs_cesaro(a_val, omega, l_max=l_max, alpha=5, r_sl=r_sl)
            if sr is None or c2 is None or c5 is None:
                continue
            w_out.append(float(omega))
            y_sr.append(sr)
            y_c2.append(c2)
            y_c5.append(c5)
        return w_out, y_sr, y_c2, y_c5


# Backward-compatible alias for previous single-case script.
Fig4A0Engine = Fig4Engine

import json
import os
import re
from collections.abc import Sequence

import mpmath as mp

mp.mp.dps = 80
M = mp.mpf("1.0")


class Fig3Engine:
    def __init__(self, data_file):
        self.data_dict = {}
        self.f_sr_curve_lookup = {}
        self._legacy_sr_warned = set()

        paths = data_file
        if isinstance(paths, (str, os.PathLike)):
            paths = [paths]
        elif not isinstance(paths, Sequence):
            paths = [str(paths)]

        for path in paths:
            self._load_data(path)

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
        if not os.path.exists(path):
            print(f"[!] Please generate {path} first.")
            return

        with open(path, "r", encoding="utf-8") as f:
            raw = json.load(f)

        incompatible_rows = 0
        for item in raw:
            if not isinstance(item, dict):
                incompatible_rows += 1
                continue

            # Optional Appendix-C aligned spherical-basis SR rows:
            # {a, omega, k_sr, theta_deg_list, f_sr_re, f_sr_im, ...}
            if (
                "theta_deg_list" in item
                and "f_sr_re" in item
                and "f_sr_im" in item
                and "a" in item
                and "omega" in item
            ):
                theta_deg_list = item.get("theta_deg_list", [])
                f_sr_re = item.get("f_sr_re", [])
                f_sr_im = item.get("f_sr_im", [])
                if (
                    isinstance(theta_deg_list, list)
                    and isinstance(f_sr_re, list)
                    and isinstance(f_sr_im, list)
                    and len(theta_deg_list) == len(f_sr_re)
                    and len(theta_deg_list) == len(f_sr_im)
                ):
                    theta_rad = [self._clean_str(x) * mp.pi / 180 for x in theta_deg_list]
                    f_sr = [
                        mp.mpc(self._clean_str(re_x), self._clean_str(im_x))
                        for re_x, im_x in zip(f_sr_re, f_sr_im)
                    ]
                    key = (
                        f"{float(item['a']):.3f}",
                        f"{float(item['omega']):.4f}",
                    )
                    self.f_sr_curve_lookup[key] = {
                        "k_sr": int(item.get("k_sr", 2)),
                        "theta_rad": theta_rad,
                        "f_sr": f_sr,
                    }
                else:
                    incompatible_rows += 1
                continue

            if "valid" in item and not bool(item["valid"]):
                continue
            if "l" not in item or "a" not in item or "omega" not in item:
                incompatible_rows += 1
                continue

            key = (
                int(item["l"]),
                f"{float(item['a']):.3f}",
                f"{float(item['omega']):.4f}",
            )

            b_re = self._clean_str(item.get("B_ratio_re", 0))
            b_im = self._clean_str(item.get("B_ratio_im", 0))

            s_theta_raw = item.get("S_theta", [])
            s_pi_minus_theta_raw = item.get("S_pi_minus_theta", [])
            if not isinstance(s_theta_raw, list) or not isinstance(s_pi_minus_theta_raw, list):
                incompatible_rows += 1
                continue

            rec = {
                "Re_C": self._clean_str(item.get("Re_C", 0)),
                "B_ratio": mp.mpc(b_re, b_im),
                "S_0": self._clean_str(item.get("S_0", 0)),
                "S_theta": [self._clean_str(x) for x in s_theta_raw],
                "S_pi_minus_theta": [self._clean_str(x) for x in s_pi_minus_theta_raw],
            }
            prev = self.data_dict.get(key)
            if prev is None or len(rec["S_theta"]) >= len(prev.get("S_theta", [])):
                self.data_dict[key] = rec

        if incompatible_rows > 0:
            print(
                "[!] Detected incompatible rows while loading Fig3 data. "
                "Please rerun the corresponding Mathematica generators."
            )

    def _find_spherical_sr_curve(self, a_str, omega_str, theta_rad_list, k):
        row = self.f_sr_curve_lookup.get((a_str, omega_str))
        if row is None:
            return None
        if int(row.get("k_sr", 2)) != int(k):
            return None

        theta_ref = row.get("theta_rad", [])
        if len(theta_ref) != len(theta_rad_list):
            return None

        for t_ref, t_req in zip(theta_ref, theta_rad_list):
            if mp.fabs(t_ref - mp.mpf(str(t_req))) > mp.mpf("1e-12"):
                return None

        return row.get("f_sr", None)

    def _scattering_prefactor(self, omega_mp, mode="eq12"):
        if mode == "eq12":
            return 2 * mp.pi / (mp.j * omega_mp)
        if mode == "pi_over_iw":
            return mp.pi / (mp.j * omega_mp)
        if mode == "half_over_iw":
            return mp.mpf("0.5") / (mp.j * omega_mp)
        raise ValueError(f"Unknown prefactor mode: {mode}")

    def _available_l_values(self, a_str, omega_str):
        vals = []
        for (l, a_k, w_k) in self.data_dict:
            if a_k == a_str and w_k == omega_str:
                vals.append(l)
        return sorted(vals)

    def calc_angular_distribution(
        self, a_val, omega_val, num_angles, prefactor_mode="eq12", l_max=None
    ):
        """Compute d sigma / d Omega over all angles for one (a, omega)."""
        a_str = f"{float(a_val):.3f}"
        omega_str = f"{float(omega_val):.4f}"
        omega_mp = mp.mpf(str(omega_val))

        sum_f = [mp.mpc("0.0")] * num_angles
        sum_g = [mp.mpc("0.0")] * num_angles
        coeff = self._scattering_prefactor(omega_mp, mode=prefactor_mode)

        available_l = self._available_l_values(a_str, omega_str)
        if not available_l:
            return [0.0] * num_angles
        max_l_data = max(available_l)
        use_l_max = max_l_data if l_max is None else min(int(l_max), max_l_data)

        for l in range(2, use_l_max + 1):
            mode = self.data_dict.get((l, a_str, omega_str))
            if mode is None:
                continue
            if len(mode["S_theta"]) < num_angles or len(mode["S_pi_minus_theta"]) < num_angles:
                continue

            use_f_term = mp.fabs(mode["Re_C"]) >= 1e-20
            term1 = mp.power(-1, l + 1)
            if use_f_term:
                bracket_f = (
                    term1 * mode["Re_C"] / (16 * mp.power(omega_mp, 4)) * mode["B_ratio"]
                    - mp.mpf("1.0")
                )
            bracket_g = (
                term1 * 3 * mp.j * M * omega_mp / (4 * mp.power(omega_mp, 4)) * mode["B_ratio"]
            )

            for idx in range(num_angles):
                if use_f_term:
                    term_f = mode["S_0"] * mode["S_theta"][idx] * bracket_f
                    sum_f[idx] += term_f
                term_g = (
                    mp.power(-1, l)
                    * mode["S_0"]
                    * mode["S_pi_minus_theta"][idx]
                    * mp.conj(bracket_g)
                )
                sum_g[idx] += term_g

        dsdo_array = []
        for idx in range(num_angles):
            cross_section = mp.power(mp.fabs(coeff * sum_f[idx]), 2) + mp.power(
                mp.fabs(coeff * sum_g[idx]), 2
            )
            dsdo_array.append(float(cross_section))

        return dsdo_array

    def _series_reduce_once(self, coeff_map, m=2):
        reduced = {}
        if not coeff_map:
            return reduced

        l_values = sorted(coeff_map)
        l_max_current = l_values[-1]
        for l in l_values:
            # Do not use the current top mode in reduction, otherwise the missing (l+1)
            # tail is implicitly treated as zero and contaminates the reduced series.
            if l >= l_max_current:
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

            # Appendix-C Eq.(C5): reduction coefficients carry 1/(l+1) and 1/l.
            reduced[l] = term_diag * fl - (a_l / (l + 1)) * fp - (b_l / l) * fm
        return reduced

    def _cesaro_weight(self, l, l_max, alpha):
        if alpha == 0:
            return mp.mpf("1.0")
        return mp.binomial(l_max, l) / mp.binomial(l_max + alpha, l)

    def calc_angular_distribution_cesaro(
        self,
        a_val,
        omega_val,
        num_angles,
        l_max=60,
        alpha=2,
        m=2,
        prefactor_mode="eq12",
    ):
        a_str = f"{float(a_val):.3f}"
        omega_str = f"{float(omega_val):.4f}"
        omega_mp = mp.mpf(str(omega_val))
        pref = self._scattering_prefactor(omega_mp, mode=prefactor_mode)

        available_l = self._available_l_values(a_str, omega_str)
        if not available_l:
            return [0.0] * num_angles
        max_l_data = max(available_l)
        use_l_max = min(int(l_max), max_l_data)
        use_alpha = int(alpha)

        sum_f = [mp.mpc("0.0")] * num_angles
        sum_g = [mp.mpc("0.0")] * num_angles

        for l in range(2, use_l_max + 1):
            mode = self.data_dict.get((l, a_str, omega_str))
            if mode is None:
                continue
            if len(mode["S_theta"]) < num_angles or len(mode["S_pi_minus_theta"]) < num_angles:
                continue

            weight = self._cesaro_weight(l, use_l_max, use_alpha)
            use_f_term = mp.fabs(mode["Re_C"]) >= 1e-20
            term1 = mp.power(-1, l + 1)
            if use_f_term:
                bracket_f = (
                    term1 * mode["Re_C"] / (16 * mp.power(omega_mp, 4)) * mode["B_ratio"]
                    - mp.mpf("1.0")
                )
            bracket_g = term1 * 3 * mp.j * M * omega_mp / (4 * mp.power(omega_mp, 4)) * mode[
                "B_ratio"
            ]

            if use_f_term:
                coef_f = pref * mode["S_0"] * bracket_f * weight
            coef_g = pref * mp.power(-1, l) * mode["S_0"] * mp.conj(bracket_g) * weight
            for idx in range(num_angles):
                if use_f_term:
                    sum_f[idx] += coef_f * mode["S_theta"][idx]
                sum_g[idx] += coef_g * mode["S_pi_minus_theta"][idx]

        dsdo = []
        for idx in range(num_angles):
            dsdo.append(float(mp.power(mp.fabs(sum_f[idx]), 2) + mp.power(mp.fabs(sum_g[idx]), 2)))
        return dsdo

    def calc_angular_distribution_series_reduction(
        self,
        a_val,
        omega_val,
        theta_rad_list,
        l_max=60,
        k=2,
        m=2,
        prefactor_mode="eq12",
    ):
        """
        Compute d sigma / d Omega(theta) using series reduction.
        """
        a_str = f"{float(a_val):.3f}"
        omega_str = f"{float(omega_val):.4f}"
        omega_mp = mp.mpf(str(omega_val))

        pref = self._scattering_prefactor(omega_mp, mode=prefactor_mode)
        available_l = self._available_l_values(a_str, omega_str)
        if not available_l:
            return [0.0] * len(theta_rad_list)
        max_l_data = max(available_l)
        use_l_max = min(int(l_max), max_l_data)

        # Preferred path: precomputed spherical-basis SR f(theta)
        # (Appendix-C aligned), while g remains modal sum.
        f_sr_curve = self._find_spherical_sr_curve(a_str, omega_str, theta_rad_list, k=k)
        if f_sr_curve is not None:
            coeff_g = {}
            basis_g = {idx: {} for idx in range(len(theta_rad_list))}

            for l in range(2, use_l_max + 1):
                mode = self.data_dict.get((l, a_str, omega_str))
                if mode is None:
                    continue
                if len(mode["S_pi_minus_theta"]) < len(theta_rad_list):
                    continue

                term1 = mp.power(-1, l + 1)
                bracket_g = (
                    term1 * 3 * mp.j * M * omega_mp / (4 * mp.power(omega_mp, 4)) * mode["B_ratio"]
                )
                coeff_g[l] = pref * mp.power(-1, l) * mode["S_0"] * mp.conj(bracket_g)

                for idx in range(len(theta_rad_list)):
                    basis_g[idx][l] = mode["S_pi_minus_theta"][idx]

            dsdo = []
            for idx in range(len(theta_rad_list)):
                sum_g = mp.mpc("0.0")
                for l in coeff_g:
                    sum_g += coeff_g[l] * basis_g[idx].get(l, mp.mpf("0.0"))

                f_theta = f_sr_curve[idx]
                g_theta = sum_g
                dsdo.append(float(mp.power(mp.fabs(f_theta), 2) + mp.power(mp.fabs(g_theta), 2)))
            return dsdo

        warn_key = (a_str, omega_str, int(k))
        if warn_key not in self._legacy_sr_warned:
            print(
                "[!] Falling back to legacy Fig3 SR (spheroidal-index reduction). "
                "Generate data/fig3_sr_spherical.json for Appendix-C aligned SR."
            )
            self._legacy_sr_warned.add(warn_key)

        # Legacy fallback path below (not Appendix-C aligned for a*omega ~ O(1)).
        coeff_f = {}
        coeff_g = {}
        basis_f = {idx: {} for idx in range(len(theta_rad_list))}
        basis_g = {idx: {} for idx in range(len(theta_rad_list))}

        # For k-order reduction, Eq. (C5) couples to l+1 repeatedly; pad high-l input.
        build_l_max = min(max_l_data, use_l_max + int(k) + 2)
        # Avoid top-edge contamination when reduction padding is not available in data.
        safe_l_max = min(use_l_max, build_l_max - (int(k) + 2))
        if safe_l_max < 2:
            return [0.0] * len(theta_rad_list)

        for l in range(2, build_l_max + 1):
            mode = self.data_dict.get((l, a_str, omega_str))
            if mode is None:
                continue
            if len(mode["S_theta"]) < len(theta_rad_list) or len(mode["S_pi_minus_theta"]) < len(
                theta_rad_list
            ):
                continue

            use_f_term = mp.fabs(mode["Re_C"]) >= 1e-20
            term1 = mp.power(-1, l + 1)
            if use_f_term:
                bracket_f = (
                    term1 * mode["Re_C"] / (16 * mp.power(omega_mp, 4)) * mode["B_ratio"]
                    - mp.mpf("1.0")
                )
            bracket_g = (
                term1 * 3 * mp.j * M * omega_mp / (4 * mp.power(omega_mp, 4)) * mode["B_ratio"]
            )

            if use_f_term:
                coeff_f[l] = pref * mode["S_0"] * bracket_f
            coeff_g[l] = pref * mp.power(-1, l) * mode["S_0"] * mp.conj(bracket_g)

            if use_f_term:
                for idx in range(len(theta_rad_list)):
                    basis_f[idx][l] = mode["S_theta"][idx]
                    basis_g[idx][l] = mode["S_pi_minus_theta"][idx]
            else:
                for idx in range(len(theta_rad_list)):
                    basis_g[idx][l] = mode["S_pi_minus_theta"][idx]

        coeff_f_k = coeff_f
        for _ in range(int(k)):
            coeff_f_k = self._series_reduce_once(coeff_f_k, m=m)

        dsdo = []
        for idx, theta in enumerate(theta_rad_list):
            sum_f = mp.mpc("0.0")
            sum_g = mp.mpc("0.0")
            for l in coeff_f_k:
                if l > safe_l_max:
                    continue
                sum_f += coeff_f_k[l] * basis_f[idx].get(l, mp.mpf("0.0"))
            # In Appendix-C, series reduction is applied to f, while g is regular.
            for l in coeff_g:
                if l > safe_l_max:
                    continue
                sum_g += coeff_g[l] * basis_g[idx].get(l, mp.mpf("0.0"))

            denom = mp.power(mp.mpf("1.0") - mp.cos(theta), int(k))
            if denom == 0:
                f_theta = mp.mpc("0.0")
                g_theta = mp.mpc("0.0")
            else:
                f_theta = sum_f / denom
                g_theta = sum_g

            dsdo.append(float(mp.power(mp.fabs(f_theta), 2) + mp.power(mp.fabs(g_theta), 2)))

        return dsdo

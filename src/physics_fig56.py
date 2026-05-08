import json
import math
import os
import re
from typing import Dict, List, Tuple

import numpy as np


M_SUN_SEC = 4.925490947e-6


class Fig56Engine:
    def __init__(self, data_files: List[str]):
        self.data: Dict[Tuple[str, str, str, int], Dict[str, float]] = {}
        self.omega_by_combo: Dict[Tuple[str, str], List[float]] = {}
        self.f_sr_lookup: Dict[Tuple[str, str, str], complex] = {}
        self.f_sr_series_raw: Dict[Tuple[str, str], List[Tuple[float, complex]]] = {}
        self.f_sr_series: Dict[Tuple[str, str], Tuple[np.ndarray, np.ndarray]] = {}
        self._scatter_cache: Dict[
            Tuple[str, str, int, int, bool, bool, float, bool],
            Tuple[np.ndarray, np.ndarray, np.ndarray],
        ] = {}
        for path in data_files:
            self._load_data_file(path)
        self._finalize_grids()

    @staticmethod
    def _a_key(a_val: float) -> str:
        return f"{float(a_val):.3f}"

    @staticmethod
    def _theta_key(theta_val: float) -> str:
        return f"{float(theta_val):.12g}"

    @staticmethod
    def _omega_key(omega_val: float) -> str:
        return f"{float(omega_val):.12g}"

    @staticmethod
    def _clean_num(x) -> float:
        s = str(x).strip().replace(" ", "").replace("*^", "e")
        s = re.sub(r"`[-+0-9.]*", "", s)
        if (not s) or any(tok in s for tok in ("Indeterminate", "Infinity", "ComplexInfinity")):
            return 0.0
        try:
            return float(s)
        except Exception:
            return 0.0

    def _load_data_file(self, path: str) -> None:
        if not path or not os.path.exists(path):
            return
        try:
            with open(path, "r", encoding="utf-8") as f:
                raw = json.load(f)
        except Exception:
            # Skip malformed/partial JSON artifacts from interrupted runs.
            return
        for item in raw:
            # Generic spherical-basis SR rows for arbitrary theta grids:
            # {a, omega, theta_obs, f_sr_re, f_sr_im, ...}
            if ("f_sr_re" in item) and ("f_sr_im" in item) and ("theta_obs" in item):
                a_key = self._a_key(float(item["a"]))
                om_key = self._omega_key(float(item["omega"]))
                th_key = self._theta_key(float(item["theta_obs"]))
                f_val = complex(
                    self._clean_num(item.get("f_sr_re", 0.0)),
                    self._clean_num(item.get("f_sr_im", 0.0)),
                )
                self.f_sr_lookup[(a_key, om_key, th_key)] = f_val
                self.f_sr_series_raw.setdefault((a_key, th_key), []).append((float(item["omega"]), f_val))
                continue

            # Optional SR(k=2) spherical-basis rows:
            # {a, omega, f_sr_theta30_re/im, f_sr_theta60_re/im, ...}
            if ("f_sr_theta30_re" in item) or ("f_sr_theta60_re" in item):
                a_key = self._a_key(float(item["a"]))
                om_key = self._omega_key(float(item["omega"]))
                th30_key = self._theta_key(math.pi / 6.0)
                th60_key = self._theta_key(math.pi / 3.0)
                if "f_sr_theta30_re" in item:
                    f30 = complex(
                        self._clean_num(item.get("f_sr_theta30_re", 0.0)),
                        self._clean_num(item.get("f_sr_theta30_im", 0.0)),
                    )
                    self.f_sr_lookup[(a_key, om_key, th30_key)] = f30
                    self.f_sr_series_raw.setdefault((a_key, th30_key), []).append((float(item["omega"]), f30))
                if "f_sr_theta60_re" in item:
                    f60 = complex(
                        self._clean_num(item.get("f_sr_theta60_re", 0.0)),
                        self._clean_num(item.get("f_sr_theta60_im", 0.0)),
                    )
                    self.f_sr_lookup[(a_key, om_key, th60_key)] = f60
                    self.f_sr_series_raw.setdefault((a_key, th60_key), []).append((float(item["omega"]), f60))
                continue

            if ("valid" in item) and (not bool(item["valid"])):
                continue
            a_val = float(item["a"])
            theta_val = float(item.get("theta_obs", math.pi / 6))
            omega_val = float(item["omega"])
            l_val = int(item["l"])
            key = (
                self._a_key(a_val),
                self._theta_key(theta_val),
                self._omega_key(omega_val),
                l_val,
            )
            self.data[key] = {
                "Re_C": self._clean_num(item.get("Re_C", 0.0)),
                "B_ratio": complex(
                    self._clean_num(item.get("B_ratio_re", 0.0)),
                    self._clean_num(item.get("B_ratio_im", 0.0)),
                ),
                "S_gamma": self._clean_num(item.get("S_gamma", 0.0)),
                "S_theta_obs": self._clean_num(item.get("S_theta_obs", 0.0)),
                "S_pi_minus_theta": self._clean_num(item.get("S_pi_minus_theta", 0.0)),
            }

    def _finalize_grids(self) -> None:
        temp: Dict[Tuple[str, str], set] = {}
        for a_key, th_key, om_key, _ in self.data:
            combo = (a_key, th_key)
            if combo not in temp:
                temp[combo] = set()
            temp[combo].add(float(om_key))
        self.omega_by_combo = {k: sorted(v) for k, v in temp.items()}
        self.f_sr_series = {}
        for key, vals in self.f_sr_series_raw.items():
            vals = sorted(vals, key=lambda x: x[0])
            x = np.array([v[0] for v in vals], dtype=float)
            y = np.array([v[1] for v in vals], dtype=complex)
            self.f_sr_series[key] = (x, y)

    def _lookup_sr_f(self, a_key: str, th_key: str, omega: float):
        om_key = self._omega_key(omega)
        v = self.f_sr_lookup.get((a_key, om_key, th_key), None)
        if v is not None:
            return v
        series = self.f_sr_series.get((a_key, th_key), None)
        if series is None:
            return None
        x, y = series
        if len(x) == 0:
            return None
        oq = float(np.clip(float(omega), float(x[0]), float(x[-1])))
        re = float(np.interp(oq, x, np.real(y)))
        im = float(np.interp(oq, x, np.imag(y)))
        return re + 1j * im

    def has_combo(self, a_val: float, theta_obs: float) -> bool:
        combo = (self._a_key(a_val), self._theta_key(theta_obs))
        return combo in self.omega_by_combo and len(self.omega_by_combo[combo]) > 0

    def _series_reduce_once(self, coeff_map: Dict[int, complex], m: int = 2) -> Dict[int, complex]:
        reduced: Dict[int, complex] = {}
        if not coeff_map:
            return reduced
        l_values = sorted(coeff_map.keys())
        l_max_current = l_values[-1]
        for l in l_values:
            if l >= l_max_current:
                continue
            fl = coeff_map.get(l, 0.0 + 0.0j)
            fp = coeff_map.get(l + 1, 0.0 + 0.0j)
            fm = coeff_map.get(l - 1, 0.0 + 0.0j)

            term_diag = 1.0 - (2.0 * m) / (l * (l + 1.0))
            a_num = (l - 1.0) * (l + 3.0) * (l - m + 1.0) * (l + m + 1.0)
            a_den = 4.0 * l * (l + 2.0) + 3.0
            a_l = math.sqrt(max(a_num / a_den, 0.0))
            b_num = (l * l - m * m) * (l * l - 4.0)
            b_den = 4.0 * l * l - 1.0
            b_l = math.sqrt(max(b_num / b_den, 0.0))

            reduced[l] = term_diag * fl - (a_l / (l + 1.0)) * fp - (b_l / l) * fm
        return reduced

    @staticmethod
    def _stabilize_complex_curve(y: np.ndarray, local_cap: float = 1.08, passes: int = 2) -> np.ndarray:
        if len(y) < 3:
            return y
        out = y.astype(complex).copy()
        for _ in range(max(int(passes), 1)):
            mag = np.abs(out)
            for i in range(1, len(out) - 1):
                left = max(mag[i - 1], 1e-30)
                right = max(mag[i + 1], 1e-30)
                lim = float(local_cap) * math.sqrt(left * right)
                if mag[i] > lim:
                    phase = out[i] / mag[i] if mag[i] > 0 else 1.0 + 0.0j
                    out[i] = phase * lim
            # refresh magnitudes after one pass
            mag = np.abs(out)
        return out

    def _compute_scatter_single(
        self,
        a_val: float,
        theta_obs: float,
        omega_val: float,
        l_max: int = 60,
        k_sr: int = 2,
        use_g_if_available: bool = True,
        prefer_sr_f: bool = True,
        require_sr_f: bool = False,
    ) -> Tuple[complex, complex]:
        a_key = self._a_key(a_val)
        th_key = self._theta_key(theta_obs)
        om_key = self._omega_key(omega_val)
        omega = float(omega_val)
        if omega <= 0:
            return 0.0 + 0.0j, 0.0 + 0.0j

        pref = 2.0 * math.pi / (1j * omega)
        sr_f_val = self._lookup_sr_f(a_key, th_key, omega)
        use_sr_f = (prefer_sr_f or require_sr_f) and sr_f_val is not None
        if require_sr_f and (not use_sr_f):
            raise ValueError(
                f"Missing SR spherical f for a={a_key}, omega={om_key}, theta={th_key}"
            )
        coeff_f: Dict[int, complex] = {}
        basis_f: Dict[int, float] = {}
        sum_g_terms = 0.0 + 0.0j
        used_g = False

        omega4 = omega**4
        for l in range(2, int(l_max) + 1):
            mode = self.data.get((a_key, th_key, om_key, l))
            if mode is None:
                continue
            b_ratio = mode["B_ratio"]
            re_c = mode["Re_C"]
            s_gamma = mode["S_gamma"]
            s_theta = mode["S_theta_obs"]

            if not use_sr_f:
                bracket_f = ((-1) ** (l + 1)) * re_c / (16.0 * omega4) * b_ratio - 1.0
                coeff_f[l] = pref * s_gamma * bracket_f
                basis_f[l] = s_theta

            if use_g_if_available and abs(mode.get("S_pi_minus_theta", 0.0)) > 0.0:
                s_pi = mode["S_pi_minus_theta"]
                bracket_g = ((-1) ** (l + 1)) * (3.0j * omega) / (4.0 * omega4) * b_ratio
                sum_g_terms += ((-1) ** l) * s_gamma * s_pi * np.conjugate(bracket_g)
                used_g = True

        if use_sr_f:
            f_val = sr_f_val
        else:
            if not coeff_f:
                return 0.0 + 0.0j, 0.0 + 0.0j
            reduced = coeff_f
            for _ in range(int(k_sr)):
                reduced = self._series_reduce_once(reduced, m=2)

            denom = (1.0 - math.cos(theta_obs)) ** int(k_sr)
            if denom == 0:
                denom = 1e-30
            sum_f = 0.0 + 0.0j
            for l, val in reduced.items():
                sum_f += val * basis_f.get(l, 0.0)
            f_val = sum_f / denom

        # Eq.(12): g is the complex-conjugate of the whole bracket.
        # Here sum_g_terms already stores the conjugated modal bracket part,
        # so we must also conjugate the common prefactor: conj(pref) = -pref (for real omega).
        g_pref = np.conjugate(pref)
        g_val = g_pref * sum_g_terms if used_g else 0.0 + 0.0j
        return f_val, g_val

    def build_scatter_grid(
        self,
        a_val: float,
        theta_obs: float,
        l_max: int = 60,
        k_sr: int = 2,
        use_g_if_available: bool = True,
        stabilize_f: bool = True,
        stabilize_cap: float = 1.08,
        prefer_sr_f: bool = True,
        require_sr_f: bool = False,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        cache_key = (
            self._a_key(a_val),
            self._theta_key(theta_obs),
            int(l_max),
            int(k_sr),
            bool(use_g_if_available),
            bool(stabilize_f),
            float(stabilize_cap),
            bool(prefer_sr_f),
            bool(require_sr_f),
        )
        if cache_key in self._scatter_cache:
            return self._scatter_cache[cache_key]

        combo_key = (self._a_key(a_val), self._theta_key(theta_obs))
        omega_grid = self.omega_by_combo.get(combo_key, [])
        if not omega_grid:
            out = (np.array([]), np.array([]), np.array([]))
            self._scatter_cache[cache_key] = out
            return out

        f_vals = []
        g_vals = []
        for w in omega_grid:
            f_w, g_w = self._compute_scatter_single(
                a_val=a_val,
                theta_obs=theta_obs,
                omega_val=w,
                l_max=l_max,
                k_sr=k_sr,
                use_g_if_available=use_g_if_available,
                prefer_sr_f=prefer_sr_f,
                require_sr_f=require_sr_f,
            )
            f_vals.append(f_w)
            g_vals.append(g_w)

        f_arr = np.array(f_vals, dtype=complex)
        if stabilize_f:
            f_arr = self._stabilize_complex_curve(f_arr, local_cap=stabilize_cap, passes=2)
        out = (np.array(omega_grid, dtype=float), f_arr, np.array(g_vals, dtype=complex))
        self._scatter_cache[cache_key] = out
        return out

    @staticmethod
    def _interp_complex(x_src: np.ndarray, y_src: np.ndarray, x_new: np.ndarray) -> np.ndarray:
        if len(x_src) == 0:
            return np.zeros_like(x_new, dtype=complex)
        x_min = float(np.min(x_src))
        x_max = float(np.max(x_src))
        xq = np.clip(x_new, x_min, x_max)
        re = np.interp(xq, x_src, np.real(y_src))
        im = np.interp(xq, x_src, np.imag(y_src))
        return re + 1j * im

    @staticmethod
    def mw_from_hz(freq_hz: np.ndarray, lens_mass_msun: float = 100.0) -> np.ndarray:
        return 2.0 * math.pi * freq_hz * (lens_mass_msun * M_SUN_SEC)

    @staticmethod
    def _direct_fd_wave(freq_hz: np.ndarray, iota: float, amp_scale: float = 1e-3) -> Tuple[np.ndarray, np.ndarray]:
        f = np.maximum(freq_hz, 1e-6)
        amp = amp_scale * np.power(f / 100.0, -7.0 / 6.0)

        m1 = 36.9
        m2 = 32.8
        eta = (m1 * m2) / ((m1 + m2) ** 2)
        m_total_sec = (m1 + m2) * M_SUN_SEC
        psi = -math.pi / 4.0 + (3.0 / (128.0 * eta)) * np.power(math.pi * m_total_sec * f, -5.0 / 3.0)
        base = amp * np.exp(1j * psi)

        plus_fac = 0.5 * (1.0 + math.cos(iota) ** 2)
        cross_fac = math.cos(iota)
        h_plus = plus_fac * base
        h_cross = 1j * cross_fac * base
        return h_plus, h_cross

    @staticmethod
    def _pure_pos_helicity_fd_wave(freq_hz: np.ndarray, amp_scale: float = 1e-3) -> Tuple[np.ndarray, np.ndarray]:
        """
        Construct a +2-helicity wave in (+,x) basis: h_x = + i h_+.
        """
        h_plus, _ = Fig56Engine._direct_fd_wave(freq_hz, iota=0.0, amp_scale=amp_scale)
        h_cross = 1j * h_plus
        return h_plus, h_cross

    def build_waveforms(
        self,
        a_val: float,
        theta_deg: float,
        freq_hz: np.ndarray,
        l_max: int = 60,
        k_sr: int = 2,
        r_sl: float = 100.0,
        r_lo: float = None,
        r_so_over_r_lo: float = 1.0,
        lens_mass_msun: float = 100.0,
        use_g_if_available: bool = True,
        stabilize_f: bool = False,
        stabilize_cap: float = 1.08,
        prefer_sr_f: bool = True,
        require_sr_f: bool = False,
    ) -> Dict[str, np.ndarray]:
        theta = math.radians(theta_deg)
        iota_obs = math.pi - theta

        h_dir_plus, h_dir_cross = self._direct_fd_wave(freq_hz, iota=iota_obs)
        h_lens_plus, h_lens_cross = self._pure_pos_helicity_fd_wave(freq_hz)

        # Appendix-D structure:
        #   h_SW = ((f +/- g) / r_LO) * h_lens
        # If h_lens is synthesized with the same baseline normalization as h_dir,
        # we include a geometric scale to emulate h_lens/h_dir ~ r_SO/r_SL.
        # Setting r_so_over_r_lo=1 and r_lo=r_sl recovers the previous effective scaling.
        r_lo_eff = float(r_lo) if (r_lo is not None) else float(r_sl)
        geom_scale = float(r_so_over_r_lo) * (r_lo_eff / float(r_sl))
        h_lens_plus = geom_scale * h_lens_plus
        h_lens_cross = geom_scale * h_lens_cross

        omega_grid, f_grid, g_grid = self.build_scatter_grid(
            a_val=a_val,
            theta_obs=theta,
            l_max=l_max,
            k_sr=k_sr,
            use_g_if_available=use_g_if_available,
            stabilize_f=stabilize_f,
            stabilize_cap=stabilize_cap,
            prefer_sr_f=prefer_sr_f,
            require_sr_f=require_sr_f,
        )
        if len(omega_grid) == 0:
            raise ValueError(f"No scattering data for a={a_val}, theta={theta_deg} deg.")

        mw = self.mw_from_hz(freq_hz, lens_mass_msun=lens_mass_msun)
        f_interp = self._interp_complex(omega_grid, f_grid, mw)
        g_interp = self._interp_complex(omega_grid, g_grid, mw)

        dt = 2.0 * r_sl * (lens_mass_msun * M_SUN_SEC) * (math.sin(theta / 2.0) ** 2)
        phase_delay = np.exp(2j * math.pi * freq_hz * dt)

        h_sw_plus = ((f_interp + g_interp) / r_lo_eff) * h_lens_plus
        h_sw_cross = ((f_interp - g_interp) / r_lo_eff) * h_lens_cross

        h_obs_plus = h_dir_plus + phase_delay * h_sw_plus
        h_obs_cross = h_dir_cross + phase_delay * h_sw_cross

        return {
            "freq_hz": freq_hz,
            "h_dir_plus": h_dir_plus,
            "h_dir_cross": h_dir_cross,
            "h_obs_plus": h_obs_plus,
            "h_obs_cross": h_obs_cross,
            "f_scatter": f_interp,
            "g_scatter": g_interp,
        }

    @staticmethod
    def _aLIGO_psd_approx(freq_hz: np.ndarray) -> np.ndarray:
        f = np.maximum(freq_hz, 10.0)
        x = f / 215.0
        psd = 1e-49 * (
            np.power(x, -4.14)
            - 5.0 * np.power(x, -2.0)
            + 111.0 * ((1.0 - x * x + 0.5 * x**4) / (1.0 + 0.5 * x * x))
        )
        return np.clip(psd, 1e-52, np.inf)

    @staticmethod
    def mismatch(h1: np.ndarray, h2: np.ndarray, freq_hz: np.ndarray, psd: np.ndarray, tmax: float = 0.03, n_t: int = 1001) -> float:
        inv_psd = 1.0 / np.maximum(psd, 1e-60)
        df = np.gradient(freq_hz)

        n1 = math.sqrt(max(4.0 * np.sum((np.abs(h1) ** 2) * inv_psd * df), 1e-60))
        n2 = math.sqrt(max(4.0 * np.sum((np.abs(h2) ** 2) * inv_psd * df), 1e-60))
        denom = n1 * n2
        if denom <= 0:
            return 1.0

        t_grid = np.linspace(-tmax, tmax, int(n_t))
        best = 0.0
        for t0 in t_grid:
            ph = np.exp(2j * math.pi * freq_hz * t0)
            inner = 4.0 * np.sum(h1 * np.conjugate(h2 * ph) * inv_psd * df)
            ov = abs(inner) / denom
            if ov > best:
                best = ov
        best = min(max(best, 0.0), 1.0)
        return float(1.0 - best)

    def build_table1(
        self,
        a_list: List[float],
        theta_list_deg: List[float],
        f_low: float = 20.0,
        f_high: float = 640.0,
        n_freq: int = 2400,
        l_max: int = 60,
        k_sr: int = 2,
        use_g_if_available: bool = True,
        prefer_sr_f: bool = True,
        require_sr_f: bool = False,
    ) -> List[Dict[str, float]]:
        freq = np.linspace(f_low, f_high, int(n_freq))
        psd = self._aLIGO_psd_approx(freq)

        rows = []
        for th in theta_list_deg:
            for a_val in a_list:
                wf = self.build_waveforms(
                    a_val=a_val,
                    theta_deg=th,
                    freq_hz=freq,
                    l_max=l_max,
                    k_sr=k_sr,
                    use_g_if_available=use_g_if_available,
                    prefer_sr_f=prefer_sr_f,
                    require_sr_f=require_sr_f,
                )
                mm_plus = self.mismatch(wf["h_dir_plus"], wf["h_obs_plus"], freq, psd)
                mm_cross = self.mismatch(wf["h_dir_cross"], wf["h_obs_cross"], freq, psd)
                rows.append(
                    {
                        "theta_deg": float(th),
                        "a": float(a_val),
                        "mismatch_plus": float(mm_plus),
                        "mismatch_cross": float(mm_cross),
                    }
                )
        return rows

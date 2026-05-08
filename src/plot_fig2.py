import concurrent.futures
import logging
import os
import re
import shutil
import sys

import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from tqdm import tqdm

# Silence noisy loader/guards by default. Users can override with env=0.
os.environ.setdefault("KERR_GW_INTERFACE_SILENT", "1")
os.environ.setdefault("KERR_GW_PHYSICS_SILENT", "1")

try:
    from .physics import (
        M,
        calc_absorption_cross_section,
        calc_backward_scattering,
        superradiant_limit,
    )
except ImportError:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from physics import (
        M,
        calc_absorption_cross_section,
        calc_backward_scattering,
        superradiant_limit,
    )


_ABS_PARITY_MODE = os.environ.get("FIG2_ABSORPTION_PARITY", "single").strip().lower()
if _ABS_PARITY_MODE in {"sum", "both", "parity_sum"}:
    _ABSORPTION_P = None
else:
    _ABSORPTION_P = 1
_CLIP_NON_SR_NEGATIVE = os.environ.get("FIG2_CLIP_NON_SR_NEGATIVE", "1").strip() != "0"
_FORCE_NONNEGATIVE_ABSORPTION = (
    os.environ.get("FIG2_FORCE_NONNEGATIVE_ABSORPTION", "0").strip() != "0"
)
_LMAX = int(os.environ.get("FIG2_LMAX", "100"))
_OMEGA_MIN = float(os.environ.get("FIG2_OMEGA_MIN", "0.1"))
_OMEGA_MAX = float(os.environ.get("FIG2_OMEGA_MAX", "2.0"))
_OMEGA_N = int(os.environ.get("FIG2_OMEGA_N", "50"))


def worker(args):
    """
    Process worker. mpmath precision must be set per process.
    """
    a_val, omega_val = args
    mp.mp.dps = 50

    a_mp = mp.mpf(str(a_val))
    omega_mp = mp.mpf(str(omega_val))

    sigma_a = calc_absorption_cross_section(
        a_mp,
        omega_mp,
        P=_ABSORPTION_P,
        l_max=_LMAX,
        clip_non_superradiant_negative=_CLIP_NON_SR_NEGATIVE,
    )
    if _FORCE_NONNEGATIVE_ABSORPTION and sigma_a < 0:
        sigma_a = 0.0
    dsigma = calc_backward_scattering(a_mp, omega_mp, l_max=_LMAX)

    sigma_a_norm = sigma_a / (mp.pi * M**2)
    dsigma_norm = dsigma / (M**2)

    return float(sigma_a_norm), float(dsigma_norm)


def _get_results_dir():
    results_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
    os.makedirs(results_dir, exist_ok=True)
    return results_dir


def _next_run_id(results_dir):
    pattern = re.compile(r"Fig2_reproduction_(\d+)\.png$")
    max_id = 0
    for name in os.listdir(results_dir):
        match = pattern.match(name)
        if match:
            max_id = max(max_id, int(match.group(1)))
    return max_id + 1


def _setup_logger(results_dir, run_id):
    logs_dir = os.path.join(results_dir, "logs")
    os.makedirs(logs_dir, exist_ok=True)

    log_path = os.path.join(logs_dir, f"fig2_run_{run_id:02d}.log")

    logger = logging.getLogger(f"fig2_run_{run_id:02d}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False

    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")

    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    return logger, log_path


def generate_and_plot_fig2():
    results_dir = _get_results_dir()
    run_id = _next_run_id(results_dir)
    output_name = f"Fig2_reproduction_{run_id:02d}.png"
    output_path = os.path.join(results_dir, output_name)

    logger, log_path = _setup_logger(results_dir, run_id)
    logger.info("Starting Kerr GW lensing Fig.2 pipeline")
    logger.info("Run ID: %02d", run_id)
    logger.info(
        "Absorption settings: parity_mode=%s (P=%s), clip_non_superradiant_negative=%s, force_nonnegative=%s, lmax=%d",
        _ABS_PARITY_MODE,
        "sum" if _ABSORPTION_P is None else _ABSORPTION_P,
        _CLIP_NON_SR_NEGATIVE,
        _FORCE_NONNEGATIVE_ABSORPTION,
        _LMAX,
    )
    logger.info(
        "Frequency grid: Momega in [%.6g, %.6g], n=%d",
        _OMEGA_MIN,
        _OMEGA_MAX,
        _OMEGA_N,
    )

    omega_vals = np.linspace(_OMEGA_MIN, _OMEGA_MAX, _OMEGA_N)
    a_configs = [
        {"val": 0.0, "label": "a = 0", "color": "navy", "ls": "-"},
        {"val": 0.8, "label": "a = 0.8M", "color": "forestgreen", "ls": "--"},
        {"val": 0.99, "label": "a = 0.99M", "color": "crimson", "ls": "-."},
        {"val": 0.999, "label": "a = 0.999M", "color": "indigo", "ls": ":"},
        {
            "val": -0.99,
            "label": "a = -0.99M",
            "color": "goldenrod",
            "ls": (0, (3, 1, 1, 1)),
        },
    ]

    results_sigma_a = {config["val"]: [] for config in a_configs}
    results_dsigma = {config["val"]: [] for config in a_configs}

    for config in a_configs:
        a_val = config["val"]
        logger.info("Computing spin case: %s", config["label"])
        m_omega_h = float(superradiant_limit(mp.mpf(str(a_val)), m=2))
        logger.info(
            "  mOmega_H (m=2) = %.8f ; superradiant window (if any): 0 < Mω < %.8f",
            m_omega_h,
            m_omega_h,
        )

        tasks = [(a_val, w) for w in omega_vals]
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = list(
                tqdm(
                    executor.map(worker, tasks),
                    total=len(tasks),
                    desc=f"Run {run_id:02d} Processing",
                    ncols=80,
                )
            )

        for res_sig, res_dsig in results:
            results_sigma_a[a_val].append(res_sig)
            results_dsigma[a_val].append(res_dsig)

        logger.info(
            "Spin %s done: sigma[min,max]=[%.6g, %.6g], dsigma[min,max]=[%.6g, %.6g]",
            config["label"],
            min(results_sigma_a[a_val]),
            max(results_sigma_a[a_val]),
            min(results_dsigma[a_val]),
            max(results_dsigma[a_val]),
        )

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax1 = axes[0]
    for config in a_configs:
        ax1.plot(
            omega_vals,
            results_sigma_a[config["val"]],
            label=f"${config['label']}$",
            color=config["color"],
            linestyle=config["ls"],
            linewidth=2,
        )
    ax1.set_xlabel(r"$M\omega$", fontsize=14)
    ax1.set_ylabel(r"$\sigma_a / \pi M^2$", fontsize=14)
    ax1.set_xlim([_OMEGA_MIN, _OMEGA_MAX])
    ax1.grid(True, linestyle="--", alpha=0.5)
    ax1.legend(loc="lower right", fontsize=12)

    ax2 = axes[1]
    for config in a_configs:
        ax2.plot(
            omega_vals,
            results_dsigma[config["val"]],
            label=f"${config['label']}$",
            color=config["color"],
            linestyle=config["ls"],
            linewidth=2,
        )
    ax2.set_xlabel(r"$M\omega$", fontsize=14)
    ax2.set_ylabel(r"$M^{-2} d\sigma(\theta=\pi)/d\Omega$", fontsize=14)
    ax2.set_xlim([_OMEGA_MIN, _OMEGA_MAX])
    ax2.grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close(fig)

    latest_path = os.path.join(results_dir, "Fig2_reproduction.png")
    shutil.copyfile(output_path, latest_path)

    logger.info("Done. Figure saved to: %s", output_path)
    logger.info("Latest figure alias updated: %s", latest_path)
    logger.info("Run log saved to: %s", log_path)


if __name__ == "__main__":
    generate_and_plot_fig2()

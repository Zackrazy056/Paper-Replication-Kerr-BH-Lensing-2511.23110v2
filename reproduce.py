"""One-command reproduction entry point for the Kerr GW lensing figures.

The individual figure scripts remain the source of truth. This wrapper fixes
their execution order, sets conservative defaults, and provides a quick smoke
mode for environment checks.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent


STEPS = {
    "sr-audit": {
        "label": "Fig.56 spherical SR coverage audit",
        "cmd": ["scripts/audit_fig56_sr_coverage.py"],
    },
    "fig2": {
        "label": "Fig.2 absorption and backward scattering",
        "cmd": ["main.py"],
    },
    "fig3": {
        "label": "Fig.3 angular scattering distribution",
        "cmd": ["main_fig3.py"],
        "env": {"FIG3_METHOD": "sr"},
    },
    "fig3-methods": {
        "label": "Fig.3 2x2 methods comparison",
        "cmd": ["main_fig3_paper_methods_grid.py"],
    },
    "fig3-audit": {
        "label": "Fig.3 numerical data audit",
        "cmd": ["scripts/audit_fig3_numeric_status.py"],
    },
    "fig4": {
        "label": "Fig.4 SFSF frequency curves",
        "cmd": ["main_fig4.py"],
    },
    "fig56": {
        "label": "Fig.5/Fig.6/Table I waveform summary",
        "cmd": ["main_fig56.py"],
    },
    "table1-paper": {
        "label": "Paper-standard Table I dependency check",
        "cmd": ["main_table1_paper_standard.py", "--check-deps"],
    },
    "fig6-paper": {
        "label": "Paper-standard Fig.6 dependency check",
        "cmd": ["main_fig6_paper_standard.py", "--check-deps"],
    },
    "fig789-paper": {
        "label": "Paper-standard Fig.7-Fig.9 dependency check",
        "cmd": ["main_fig789.py", "--check-deps"],
    },
}


SMOKE_SNIPPETS = [
    "import main, main_fig3, main_fig4, main_fig56, main_fig6_paper_standard, main_table1_paper_standard, main_fig789; print('import smoke ok')",
    (
        "import os, math, numpy as np; "
        "from src.physics_fig3 import Fig3Engine; "
        "e=Fig3Engine([os.path.join('data','fig3_data.json'),"
        "os.path.join('data','fig3_patch_w2_a0_a099.json'),"
        "os.path.join('data','fig3_sr_spherical.json')]); "
        "y=e.calc_angular_distribution_series_reduction(0.0,2.0,"
        "np.deg2rad(np.linspace(40,180,281)),l_max=12,k=2); "
        "print('fig3 smoke points', len(y))"
    ),
    (
        "import os, math, numpy as np; "
        "from src.physics_fig56 import Fig56Engine; "
        "e=Fig56Engine([os.path.join('data','fig4_data.json'),"
        "os.path.join('data','fig56_patch_data.json'),"
        "os.path.join('data','fig56_patch_a_minus099_theta30_hi.json'),"
        "os.path.join('data','fig56_sr_spherical.json')]); "
        "wf=e.build_waveforms(0.0,30.0,np.linspace(180,190,8),l_max=12,"
        "prefer_sr_f=False,require_sr_f=False); "
        "print('fig56 smoke points', len(wf['freq_hz']))"
    ),
]


COMPILE_PATHS = [
    "reproduce.py",
    "main.py",
    "main_fig3.py",
    "main_fig3_paper_methods_grid.py",
    "main_fig3_series_test.py",
    "main_fig3_w2_a0_fast.py",
    "main_fig4.py",
    "main_fig4_a0.py",
    "main_fig56.py",
    "main_fig56_single_case.py",
    "main_fig6_paper_standard.py",
    "main_fig789.py",
    "main_table1_paper_standard.py",
    "scripts",
    "src",
]


def _env_bool(name: str, value: bool) -> str:
    return "1" if value else "0"


def build_env(args: argparse.Namespace, step_env: dict[str, str] | None = None) -> dict[str, str]:
    env = os.environ.copy()
    env.setdefault("KERR_GW_INTERFACE_SILENT", "1")
    env.setdefault("KERR_GW_PHYSICS_SILENT", "1")
    env.setdefault("MPLBACKEND", "Agg")
    env["FIG56_REQUIRE_SR_F"] = _env_bool("FIG56_REQUIRE_SR_F", args.strict_sr)

    if args.quick:
        env.setdefault("FIG2_OMEGA_N", "8")
        env.setdefault("FIG2_LMAX", "25")
        env.setdefault("FIG3_LMAX", "30")
        env.setdefault("FIG56_LMAX", "30")
        env.setdefault("FIG56_TABLE_N_FREQ", "500")

    if step_env:
        env.update(step_env)
    return env


def run_subprocess(label: str, cmd: list[str], env: dict[str, str], python: str, dry_run: bool) -> None:
    full_cmd = [python, *cmd]
    shown = " ".join(full_cmd)
    print(f"\n=== {label} ===", flush=True)
    print(f"$ {shown}", flush=True)
    if dry_run:
        return
    completed = subprocess.run(full_cmd, cwd=ROOT, env=env)
    if completed.returncode != 0:
        raise SystemExit(completed.returncode)


def run_smoke(args: argparse.Namespace) -> None:
    env = build_env(args)
    run_subprocess("Compile Python files", ["-m", "compileall", "-q", *COMPILE_PATHS], env, args.python, args.dry_run)
    for idx, snippet in enumerate(SMOKE_SNIPPETS, start=1):
        run_subprocess(f"Smoke check {idx}", ["-c", snippet], env, args.python, args.dry_run)


def expand_targets(targets: list[str]) -> list[str]:
    if not targets or "all" in targets:
        return ["fig2", "fig3", "fig4", "fig56"]
    seen: set[str] = set()
    ordered: list[str] = []
    for target in targets:
        if target == "smoke":
            ordered.append(target)
            continue
        if target not in STEPS:
            raise SystemExit(f"Unknown target: {target}")
        if target not in seen:
            ordered.append(target)
            seen.add(target)
    return ordered


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reproduce Kerr GW lensing figures and paper-standard dependency checks from cached data."
    )
    parser.add_argument(
        "targets",
        nargs="*",
        help="Targets to run: all, smoke, sr-audit, fig3-audit, fig3-methods, fig2, fig3, fig4, fig56, fig6-paper, table1-paper, fig789-paper. Default: all.",
    )
    parser.add_argument("--python", default=sys.executable, help="Python executable to use.")
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Use smaller grids for a fast wiring check. Outputs are diagnostic, not final.",
    )
    parser.add_argument(
        "--strict-sr",
        action="store_true",
        help="Require spherical-basis series-reduction data in Fig.56. Default allows fallback.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running them.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    targets = expand_targets(args.targets)

    for target in targets:
        if target == "smoke":
            run_smoke(args)
            continue
        step = STEPS[target]
        env = build_env(args, step.get("env"))
        run_subprocess(step["label"], step["cmd"], env, args.python, args.dry_run)

    if not args.dry_run:
        print("\n[+] Reproduction workflow completed.")
        print(f"[+] Outputs are under: {ROOT / 'results'}")


if __name__ == "__main__":
    main()

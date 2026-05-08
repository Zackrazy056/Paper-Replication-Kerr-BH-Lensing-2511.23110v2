# Fig.3 Cluster Runbook

Run these commands from `fig3_workspace/` after copying the workspace to a
Slurm cluster that has WolframScript, BHPToolkit packages, Python, NumPy, and
Matplotlib available.

## 1. Set BHPToolkit paths

Edit these paths for your cluster before submitting jobs:

```bash
export BHP_SW_SH_PATH="$HOME/SpinWeightedSpheroidalHarmonics/Kernel"
export BHP_KERR_GEODESICS_PATH="$HOME/KerrGeodesics/Kernel"
export BHP_TEUKOLSKY_PATH="$HOME/Teukolsky/Kernel"
```

If the executable names differ, also set:

```bash
export WOLFRAMSCRIPT=/path/to/wolframscript
export PYTHON_BIN=python3
```

## 2. Submit high-priority repair jobs

```bash
sbatch cluster/run_fig3_a099_w2_l30_70.slurm
sbatch cluster/run_fig3_a099_w06_l42_60.slurm
```

The first job fills `a=0.99, Momega=2, l=30..70` on the full 281-point angle
grid. The second fills `a=0.99, Momega=0.6, l=42..60`.

Check status:

```bash
squeue -u "$USER"
tail -f fig3_fig3-a099-w2_*.out
```

## 3. Rebuild SR cache and final Fig.3

Run this after the two patch jobs finish:

```bash
sbatch cluster/run_fig3_merge_sr_plot.slurm
```

It merges the new patch files into `data/fig3_patch_w2_a0_a099.json`, rebuilds
`data/fig3_sr_spherical.json` with `FIG3_SR_LMAX=70` and basis pad 24, reruns
the audit, redraws `results/Fig3_paper_methods_grid.png`, and writes the
single-angle convergence diagnostic:

```text
results/fig3_partial_sum_a099_w2_theta120_l70.csv
results/fig3_partial_sum_a099_w2_theta120_l70.png
```

## 4. Useful overrides

For a smaller test batch:

```bash
FIG3_PATCH_BATCH=1 FIG3_PATCH_MAX_ITERS=1 sbatch cluster/run_fig3_a099_w2_l30_70.slurm
```

For stricter radial-ratio verification:

```bash
FIG3_PATCH_CHECK_BRATIO=1 FIG3_PATCH_BRATIO_TOL='10^-8' sbatch cluster/run_fig3_a099_w2_l30_70.slurm
```

For a larger SR basis:

```bash
FIG3_SR_LMAX=70 FIG3_SR_BASIS_PAD=30 sbatch cluster/run_fig3_merge_sr_plot.slurm
```

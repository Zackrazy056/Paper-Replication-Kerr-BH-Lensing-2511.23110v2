# Fig.3 Repair Runbook

This workspace contains the Fig.3-only numerical repair path. The target is to
remove the high-frequency oscillatory artifacts by completing the missing
high-spin modes and rebuilding the spherical-basis series-reduction cache.

## Local PowerShell commands

Run from `Kerr_GW_Lensing_Fig2/fig3_workspace`.

```powershell
$env:FIG3_KERNELS = "8"
$env:FIG3_PATCH_A_LIST = "99/100"
$env:FIG3_PATCH_OMEGA_LIST = "2"
$env:FIG3_PATCH_LMIN = "30"
$env:FIG3_PATCH_LMAX = "70"
$env:FIG3_PATCH_BATCH = "2"
$env:FIG3_PATCH_OUT = "fig3_patch_w2_a099_l30_70.json"
wolframscript -file mathematica/generate_fig3_patch_w2_highspin.wls
```

```powershell
$env:FIG3_KERNELS = "8"
$env:FIG3_PATCH_A_LIST = "99/100"
$env:FIG3_PATCH_OMEGA_LIST = "6/10"
$env:FIG3_PATCH_LMIN = "42"
$env:FIG3_PATCH_LMAX = "60"
$env:FIG3_PATCH_BATCH = "4"
$env:FIG3_PATCH_OUT = "fig3_patch_w06_a099_l42_60.json"
wolframscript -file mathematica/generate_fig3_patch_w2_highspin.wls
```

After both patch files are complete:

```powershell
python scripts/merge_fig3_patch_files.py `
  --base data/fig3_patch_w2_a0_a099.json `
  --out data/fig3_patch_w2_a0_a099.json `
  data/fig3_patch_w2_a099_l30_70.json `
  data/fig3_patch_w06_a099_l42_60.json

$env:FIG3_KERNELS = "8"
$env:FIG3_SR_COMBOS = "0,0.6;0,2;99/100,0.6;99/100,2"
$env:FIG3_SR_LMAX = "70"
$env:FIG3_SR_BASIS_PAD = "24"
$env:FIG3_SR_FULL_STEPS = "800"
$env:FIG3_SR_NUM_TERMS = "80"
$env:FIG3_SR_WORKING_PREC = "80"
$env:FIG3_SR_PARALLEL = "1"
wolframscript -file mathematica/generate_fig3_sr_spherical.wls

python scripts/audit_fig3_numeric_status.py
python main_fig3_paper_methods_grid.py
python scripts/diag_fig3_partial_sum_lmax.py --a 0.99 --omega 2 --theta-deg 120 --l-max 70 --out-prefix results/fig3_partial_sum_a099_w2_theta120_l70
```

## Cluster commands

Run from `fig3_workspace`:

```bash
sbatch cluster/run_fig3_a099_w2_l30_70.slurm
sbatch cluster/run_fig3_a099_w06_l42_60.slurm
sbatch cluster/run_fig3_merge_sr_plot.slurm
```

Submit the merge/SR/plot job only after the two patch jobs finish.

## Expected checks

The audit should show that `a=0.99, Momega=2` has full 281-point modes beyond
`l=29`, ideally through `l=70`. The rebuilt `fig3_sr_spherical.json` should
report `l_mode_max=70` for the high-frequency rows and `l_basis_max` near 94
with the default basis pad.

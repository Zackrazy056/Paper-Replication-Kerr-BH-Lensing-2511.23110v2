<<<<<<< HEAD
# Kerr GW Lensing Figure Reproduction

This directory contains a cached-data reproduction workflow for the
on-axis Kerr gravitational-wave scattering/lensing figures in
`2511.23110v2.pdf`.

The current Python workflow is designed to rerun the available local results:

- Fig.2: on-axis absorption cross section and backward scattering.
- Fig.3: angular scattering distribution with series reduction.
- Fig.4: strong-field scattering factor (SFSF).
- Fig.5/Fig.6/Table I: geometry, waveform illustration, and lightweight
  mismatch table.
- Paper-standard Fig.6, Table I, and Fig.7-Fig.9: PyCBC/IMRPhenomD waveform
  and fitting-factor calculations through the WSL/Linux helpers.

## Environment

Use Python 3.11 or newer. On Windows PowerShell, from this directory:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

If the repository is checked out under a parent directory that already has a
virtual environment, use that Python executable explicitly with
`reproduce.py --python`.

In this workspace, the existing environment lives one directory above this
project and has been installed with `requirements.txt`:

```powershell
..\.venv\Scripts\python.exe reproduce.py smoke
..\.venv\Scripts\python.exe reproduce.py all
```

## One-Command Reproduction

Run all lightweight cached-data figures and the local Table I approximation:

```powershell
python reproduce.py all
```

Or with an explicit virtual-environment Python:

```powershell
.\.venv\Scripts\python.exe reproduce.py all
```

Outputs are written under `results/`.

For a fast wiring check that uses smaller numerical grids:

```powershell
python reproduce.py smoke
python reproduce.py --quick all
```

The `--quick` outputs are diagnostic only and should not be used as final
paper-comparison figures.

## Targets

Run individual targets as needed:

```powershell
python reproduce.py sr-audit
python reproduce.py fig2
python reproduce.py fig3
python reproduce.py fig4
python reproduce.py fig56
python reproduce.py fig6-paper
python reproduce.py table1-paper
python reproduce.py fig789-paper
```

`fig56` produces Fig.5, Fig.6, and `Table1_reproduction.csv/json`.
`sr-audit` writes `results/fig56_sr_coverage.json`.
`fig6-paper`, `table1-paper`, and `fig789-paper` check whether
PyCBC/LALSuite/SciPy are available in the active Python environment. On
Windows, use the WSL helpers below for the full paper-standard runs.

## Strict Series-Reduction Mode

The cached spherical-basis series-reduction data for Fig.56 is incomplete for
some spin values. By default, the one-command workflow uses the legacy modal
path for all Fig.56 spin cases, so every cached configuration reruns through the
same code path.

To prefer spherical SR rows where they exist and fall back only when missing:

```powershell
$env:FIG56_PREFER_SR_F = "1"
python reproduce.py fig56
```

To require spherical-basis SR data and fail on missing rows:

```powershell
python reproduce.py --strict-sr fig56
```

Equivalent environment variable:

```powershell
$env:FIG56_REQUIRE_SR_F = "1"
python main_fig56.py
```

## Data Generation

The Python scripts consume cached JSON files in `data/`. The Mathematica
generators under `mathematica/` are the source for those JSON files. Typical
examples:

```powershell
wolframscript -file mathematica/generate_teukolsky_data.wls
wolframscript -file mathematica/generate_fig3_sr_spherical.wls
wolframscript -file mathematica/generate_fig56_sr_spherical.wls
wolframscript -file mathematica/generate_fig789_theta_data.wls
```

The Python workflow does not require Mathematica if the JSON files are already
present.

The Fig.56 spherical SR cache can be audited with:

```powershell
python reproduce.py sr-audit
```

The Fig.3 cache and current plotted curve stability can be audited with:

```powershell
python reproduce.py fig3-audit
```

The Fig.3 paper-style 2x2 methods grid is generated with:

```powershell
python reproduce.py fig3-methods
```

If Wolfram Engine is activated, missing batches can be generated with:

```powershell
$env:FIG56_SR_A_FILTER = "0.8"
$env:FIG56_SR_BATCH = "8"
wolframscript -file mathematica/generate_fig56_sr_spherical.wls
```

Repeat with `0.99` and `-0.99` until `sr-audit` reports complete coverage.
The current cache reports `280/280 existing, 0 missing, 0 extra`.

For the Fig.3 high-frequency patch, generate batches first and then rebuild the
spherical-basis SR cache:

```powershell
$env:FIG3_PATCH_BATCH = "20"
wolframscript -file mathematica/generate_fig3_patch_w2_highspin.wls
wolframscript -file mathematica/generate_fig3_sr_spherical.wls
python reproduce.py fig3-audit
```

## Paper-Standard Table I

The paper-standard path is implemented in `main_table1_paper_standard.py`. It
is stricter than `main_fig56.py` and requires PyCBC, LALSuite/lalsimulation,
and SciPy:

```powershell
python -m pip install -r requirements-paper.txt
python main_table1_paper_standard.py --check-deps
python main_table1_paper_standard.py --strict-sr --optimize
```

On Windows, run the paper-standard path through WSL/Linux. This workspace has a
helper that downloads Linux wheels from Windows, installs them offline in a WSL
virtual environment, and runs the optimized Table I calculation:

```powershell
.\scripts\run_table1_paper_wsl.ps1
```

Useful variants:

```powershell
.\scripts\run_table1_paper_wsl.ps1 -SkipDownload
.\scripts\run_table1_paper_wsl.ps1 -SkipDownload -NoOptimize -Out results/Table1_paper_standard_unoptimized.json
```

The helper assumes WSL Python 3.12 by default and writes/uses:

- `.venv-paper-wsl/`
- `wheelhouse-paper-wsl/`
- `results/Table1_paper_standard.csv/json`

If your WSL Python is not 3.12, pass matching wheel tags, for example
`-PythonVersion 311 -Abi cp311`.

## Paper-Standard Fig.6

The paper-standard Fig.6 driver is implemented in
`main_fig6_paper_standard.py`. It uses PyCBC `IMRPhenomD` for both the direct
component and the incident lens-side waveform, then applies the cached
Teukolsky scattering factors:

```powershell
python main_fig6_paper_standard.py --check-deps
python main_fig6_paper_standard.py --strict-sr --out results/Fig6_paper_standard.png
```

On Windows, run it through WSL/Linux:

```powershell
.\scripts\run_fig6_paper_wsl.ps1 -SkipDownload
```

The helper writes:

- `results/Fig6_paper_standard.png`
- `results/Fig6_paper_standard_summary.json`

## Paper-Standard Fig.7-Fig.9

The Fig.7-Fig.9 contour driver is implemented in `main_fig789.py`. It reuses
the paper-standard PyCBC/LALSuite stack and writes JSON, CSV, and three contour
PNGs:

```powershell
python main_fig789.py --check-deps
python main_fig789.py --strict-sr --optimize --out results/Fig789_grid.json
python main_fig789.py --strict-sr --optimize --multistart --out results/Fig789_grid_multistart.json
```

On Windows, run it through WSL/Linux:

```powershell
.\scripts\run_fig789_paper_wsl.ps1 -SkipDownload
.\scripts\run_fig789_paper_wsl.ps1 -SkipDownload -MultiStart -Out results/Fig789_grid_multistart.json
```

The current verified default warm-start full-grid run writes:

- `results/Fig789_grid.json`
- `results/Fig789_grid.csv`
- `results/Fig789_grid_Fig7_mismatch_contour.png`
- `results/Fig789_grid_Fig8_snr_ratio_contour.png`
- `results/Fig789_grid_Fig9_required_snr_contour.png`

For the 9 x 5 optimized anchor-mode grid, the current ranges are:

- `mismatch_plus`: `0.0027813` to `0.3290385`
- `mismatch_cross`: `0.0035190` to `0.7127025`
- `snr_ratio_plus`: `0.9959707` to `1.5573498`
- `snr_ratio_cross`: `0.9973723` to `34.0788172`
- `required_snr_plus`: `6.0312651` to `60.0037638`
- `required_snr_cross`: `4.6689733` to `53.3546289`

The current verified multistart full-grid run writes:

- `results/Fig789_grid_multistart.json`
- `results/Fig789_grid_multistart.csv`
- `results/Fig789_grid_multistart_Fig7_mismatch_contour.png`
- `results/Fig789_grid_multistart_Fig8_snr_ratio_contour.png`
- `results/Fig789_grid_multistart_Fig9_required_snr_contour.png`

For the same 9 x 5 optimized anchor-mode grid with `--multistart`, the current
ranges are:

- `mismatch_plus`: `0.0027791` to `0.3295600`
- `mismatch_cross`: `0.0035196` to `0.4488003`
- `snr_ratio_plus`: `0.9959707` to `1.5573498`
- `snr_ratio_cross`: `0.9973723` to `34.0788172`
- `required_snr_plus`: `6.0274322` to `60.0278032`
- `required_snr_cross`: `5.3598741` to `53.3504112`

Relative to the default warm-start grid, the multistart run changes
`mismatch_plus` by `-0.0352995` to `+0.0005372` and `mismatch_cross` by
`-0.5125823` to `+0.0000050`; the largest improvement is in high-inclination
cross-polarization fits.

An exact-mode, single-point multistart check at `theta_o=30 deg, a=0` writes
`results/Fig789_exact_anchor_check_multistart.json` and currently gives
`mismatch_plus=0.0496183`, `mismatch_cross=0.0508169`. This is slightly lower
than the default warm-start full-grid value at the same point, confirming that
the main residual difference is optimizer-start sensitivity rather than a
scatter-path mismatch.

Useful quick checks:

```powershell
.\scripts\run_fig789_paper_wsl.ps1 -SkipDownload -NoOptimize -ThetaCount 3 -NFreq 400 -Out results/Fig789_smoke.json
.\scripts\run_fig789_paper_wsl.ps1 -SkipDownload -NoOptimize -ExactTheta -ThetaCount 1 -ThetaMin 30 -ThetaMax 30 -SpinValues 0 -NFreq 400 -Out results/Fig789_exact_smoke.json
.\scripts\run_fig789_paper_wsl.ps1 -SkipDownload -ExactTheta -MultiStart -ThetaCount 1 -ThetaMin 30 -ThetaMax 30 -SpinValues 0 -NFreq 1000 -Out results/Fig789_exact_anchor_check_multistart.json
```

## Useful Environment Knobs

- `FIG2_OMEGA_N`, `FIG2_LMAX`: Fig.2 frequency grid size and truncation.
- `FIG3_METHOD`: `sr`, `cesaro`, or `plain`.
- `FIG3_LMAX`: Fig.3 truncation.
- `FIG56_LMAX`, `FIG56_TABLE_N_FREQ`: Fig.56/Table I speed controls.
- `FIG56_PREFER_SR_F`: set to `1` to use spherical SR rows when available.
- `FIG56_REQUIRE_SR_F`: set to `1` for strict spherical SR data.

## Current Caveats

- `Table1_reproduction.csv/json` is the lightweight local approximation.
  `Table1_paper_standard.csv/json` is the PyCBC/IMRPhenomD plus Nelder-Mead
  fitting-factor path.
- `Fig6_reproduction.png` is the lightweight analytic-waveform plot. Use
  `Fig6_paper_standard.png` for the PyCBC/IMRPhenomD envelope comparable to
  the paper figure.
- The current paper-standard runner uses the local approximate aLIGO PSD and the
  implemented chirp-mass/q/spin optimization bounds. Exact numerical agreement
  with the paper still depends on matching any paper-specific PSD and optimizer
  settings not explicit in the cached data.
- Fig.3 is not yet paper-aligned numerically: the `Momega=0.6` panels still fall
  back to the legacy SR path, and high-spin/high-frequency angular coverage is
  incomplete above the currently densified modes. Run `python reproduce.py
  fig3-audit` before trusting a Fig.3 comparison.
- Fig.7-Fig.9 use strict SR anchor curves at `theta_o=30,60 deg` by default
  and interpolate/extrapolate the complex scattering amplitudes in `theta_o`;
  use `--theta-mode exact` only for theta values present in the scattering
  cache.
- The default Fig.7-Fig.9 optimizer uses Nelder-Mead with a baseline start and
  warm starts across the grid. Add `--multistart`/`-MultiStart` for slower but
  less order-sensitive fitting-factor checks.
=======
# Paper-Replication-Kerr-BH-Lensing-2511.23110v2
>>>>>>> origin/main

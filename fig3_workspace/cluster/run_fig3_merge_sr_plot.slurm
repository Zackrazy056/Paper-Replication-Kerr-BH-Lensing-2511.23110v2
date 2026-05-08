#!/usr/bin/env bash
#SBATCH --job-name=fig3-sr-plot
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=fig3_%x_%j.out
#SBATCH --error=fig3_%x_%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_DIR}"

mkdir -p results/logs data

if command -v module >/dev/null 2>&1; then
  module purge || true
  module load wolfram || module load Mathematica || true
  module load python || true
fi

: "${PYTHON_BIN:=python3}"
: "${WOLFRAMSCRIPT:=wolframscript}"
: "${BHP_SW_SH_PATH:=${HOME}/SpinWeightedSpheroidalHarmonics/Kernel}"

export BHP_SW_SH_PATH
export FIG3_KERNELS="${FIG3_KERNELS:-${SLURM_CPUS_PER_TASK:-4}}"
export FIG3_SR_COMBOS="${FIG3_SR_COMBOS:-0,0.6;0,2;99/100,0.6;99/100,2}"
export FIG3_SR_LMAX="${FIG3_SR_LMAX:-70}"
export FIG3_SR_BASIS_PAD="${FIG3_SR_BASIS_PAD:-24}"
export FIG3_SR_FULL_STEPS="${FIG3_SR_FULL_STEPS:-800}"
export FIG3_SR_NUM_TERMS="${FIG3_SR_NUM_TERMS:-80}"
export FIG3_SR_WORKING_PREC="${FIG3_SR_WORKING_PREC:-80}"
export FIG3_SR_PARALLEL="${FIG3_SR_PARALLEL:-1}"

echo "[*] Merging patch files into data/fig3_patch_w2_a0_a099.json"
"${PYTHON_BIN}" scripts/merge_fig3_patch_files.py \
  --base data/fig3_patch_w2_a0_a099.json \
  --out data/fig3_patch_w2_a0_a099.json \
  data/fig3_patch_w2_a099_l30_70.json \
  data/fig3_patch_w06_a099_l42_60.json

echo "[*] Rebuilding spherical SR cache"
"${WOLFRAMSCRIPT}" -file mathematica/generate_fig3_sr_spherical.wls | tee results/logs/fig3_sr_rebuild.log

echo "[*] Auditing and plotting Fig.3"
"${PYTHON_BIN}" scripts/audit_fig3_numeric_status.py | tee results/logs/fig3_numeric_audit_after_patch.log
"${PYTHON_BIN}" main_fig3_paper_methods_grid.py | tee results/logs/fig3_plot_after_patch.log
"${PYTHON_BIN}" scripts/diag_fig3_partial_sum_lmax.py \
  --a 0.99 --omega 2 --theta-deg 120 --l-max 70 \
  --out-prefix results/fig3_partial_sum_a099_w2_theta120_l70

echo "[+] Fig.3 postprocess complete."

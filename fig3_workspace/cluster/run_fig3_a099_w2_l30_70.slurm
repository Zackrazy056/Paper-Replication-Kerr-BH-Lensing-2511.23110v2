#!/usr/bin/env bash
#SBATCH --job-name=fig3-a099-w2
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
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
fi

: "${WOLFRAMSCRIPT:=wolframscript}"
: "${BHP_SW_SH_PATH:=${HOME}/SpinWeightedSpheroidalHarmonics/Kernel}"
: "${BHP_KERR_GEODESICS_PATH:=${HOME}/KerrGeodesics/Kernel}"
: "${BHP_TEUKOLSKY_PATH:=${HOME}/Teukolsky/Kernel}"

export BHP_SW_SH_PATH
export BHP_KERR_GEODESICS_PATH
export BHP_TEUKOLSKY_PATH
export FIG3_KERNELS="${FIG3_KERNELS:-${SLURM_CPUS_PER_TASK:-8}}"

export FIG3_PATCH_A_LIST="99/100"
export FIG3_PATCH_OMEGA_LIST="2"
export FIG3_PATCH_LMIN="30"
export FIG3_PATCH_LMAX="70"
export FIG3_PATCH_BATCH="${FIG3_PATCH_BATCH:-6}"
export FIG3_PATCH_OUT="fig3_patch_w2_a099_l30_70.json"
export FIG3_PATCH_BRATIO_TOL="${FIG3_PATCH_BRATIO_TOL:-10^-8}"
export FIG3_PATCH_CHECK_BRATIO="${FIG3_PATCH_CHECK_BRATIO:-0}"

max_iters="${FIG3_PATCH_MAX_ITERS:-40}"

echo "[*] Project: ${PROJECT_DIR}"
echo "[*] Wolfram: ${WOLFRAMSCRIPT}"
echo "[*] Kernels: ${FIG3_KERNELS}"
echo "[*] Output: data/${FIG3_PATCH_OUT}"

for iter in $(seq 1 "${max_iters}"); do
  log_file="results/logs/fig3_a099_w2_l30_70_iter${iter}.log"
  echo "[*] Iteration ${iter}/${max_iters}"
  "${WOLFRAMSCRIPT}" -file mathematica/generate_fig3_patch_w2_highspin.wls | tee "${log_file}"
  if grep -q "Nothing to run" "${log_file}" || grep -q "remaining=0" "${log_file}"; then
    echo "[+] Patch job complete."
    break
  fi
done

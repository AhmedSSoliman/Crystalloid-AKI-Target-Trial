#!/bin/bash
# IDEALIST LR vs NS  —  Master run script
# ==========================================
# Orchestrates: extract_cohort → full_analysis → sensitivity_analyses → figures
#
# Usage (login node, interactive):
#   bash run_idealist.sh
#
# Usage (HiPerGator sbatch):
#   sbatch run_idealist.sh
#
#SBATCH --job-name=idealist_causal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256gb
#SBATCH --time=12:00:00
#SBATCH --partition=bigmem
#SBATCH --account=prismap-ai-core
#SBATCH --output=/blue/prismap-ai-core/Ahmed/DigitalTwins/MDT/CausalAI/idealist_causal/logs/idealist_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ahmed.shokry.soliman@ufl.edu

set -euo pipefail

# ── Paths ──────────────────────────────────────────────────────────────────────
SCRIPT_DIR="/blue/prismap-ai-core/Ahmed/DigitalTwins/MDT/CausalAI/idealist_causal"
ANALYSIS="${SCRIPT_DIR}/analysis"
LOGS="${SCRIPT_DIR}/logs"
mkdir -p "${LOGS}"

echo "============================================="
echo "IDEALIST LR vs NS  —  Pipeline Start"
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "============================================="

# ── Python environment ─────────────────────────────────────────────────────────
module purge
module load conda
conda activate /blue/prismap-ai-core/Ahmed/envs/causal_env 2>/dev/null || \
  conda activate causal_env 2>/dev/null || \
  true   # fall through to system python if neither exists

PYTHON="${PYTHON:-python3}"
echo "Python: $(which ${PYTHON})"
${PYTHON} --version

# ── Step 1: Extract cohort ─────────────────────────────────────────────────────
echo ""
echo "--- STEP 1: Cohort Extraction ---"
OUT_COHORT="${ANALYSIS}/cohort_idealist.parquet"
if [[ -f "${OUT_COHORT}" ]]; then
    echo "  cohort_idealist.parquet already exists — skipping extraction"
    echo "  (delete it and rerun to force re-extraction)"
else
    ${PYTHON} "${ANALYSIS}/extract_cohort.py" --out "${OUT_COHORT}"
fi

# ── Step 2: Primary analysis (AIPW + SA1-4 + SA6a) ────────────────────────────
echo ""
echo "--- STEP 2: Primary AIPW Analysis ---"
OUT_RES="${ANALYSIS}/results_summary_idealist.csv"
OUT_BAL="${ANALYSIS}/balance_idealist.csv"
if [[ -f "${OUT_RES}" && -f "${OUT_BAL}" ]]; then
    echo "  results_summary_idealist.csv exists — skipping primary analysis"
    echo "  (delete it and rerun to force re-analysis)"
else
    ${PYTHON} "${ANALYSIS}/full_analysis.py" --parquet "${OUT_COHORT}"
fi

# ── Step 3: Sensitivity analyses ─────────────────────────────────────────────
echo ""
echo "--- STEP 3: Sensitivity Analyses ---"
# sensitivity_analyses.py reads cohort_full_ps_idealist.parquet from its directory
${PYTHON} "${ANALYSIS}/sensitivity_analyses.py"

# ── Step 4: Figures ───────────────────────────────────────────────────────────
echo ""
echo "--- STEP 4: Publication Figures ---"
${PYTHON} "${ANALYSIS}/make_figures_idealist.py"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "============================================="
echo "IDEALIST Pipeline Complete"
echo "Date: $(date)"
echo ""
echo "Key outputs:"
echo "  Cohort:   ${OUT_COHORT}"
echo "  Results:  ${OUT_RES}"
echo "  Balance:  ${OUT_BAL}"
echo "  Figures:  ${ANALYSIS}/figures/"
echo "  E-value:  ${ANALYSIS}/evalue_summary_idealist.csv"
echo "============================================="

# ── Print results summary ──────────────────────────────────────────────────────
echo ""
echo "--- Results Summary ---"
${PYTHON} - <<'PYEOF'
import pandas as pd, os, sys
csv = os.path.join(os.environ.get("ANALYSIS", "."), "results_summary_idealist.csv")
if os.path.exists(csv):
    df = pd.read_csv(csv)
    df["rd_pp"] = df["rd"] * 100
    df["lo_pp"] = df["rd_lower"] * 100
    df["hi_pp"] = df["rd_upper"] * 100
    print(df[["analysis","n","rd_pp","lo_pp","hi_pp"]].to_string(index=False))
else:
    print("results_summary_idealist.csv not found")
PYEOF


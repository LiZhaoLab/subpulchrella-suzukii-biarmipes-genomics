#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --job-name=syn_QC
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/syn_QC.%j.out
#SBATCH --error=logs/syn_QC.%j.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/synteny_4species/reanalysis_20260824

cd "${BASE}"

python3 \
  scripts/04_summarize_promer_alignment_qc.py

echo "===== Distribution summary ====="

column -t -s $'\t' \
  checks/alignment_distributions/promer_alignment_distribution_summary.tsv

echo "FINAL: PASS"

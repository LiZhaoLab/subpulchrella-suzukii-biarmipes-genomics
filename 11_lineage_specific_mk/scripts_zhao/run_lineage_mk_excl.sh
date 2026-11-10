#!/bin/bash
# Lineage-specific (polarized) McDonald-Kreitman tests after removal of fixed differences that are polymorphic in the
# other species: the steps of scripts/114_polarized_mk_excl.sbatch (job log logs/pmkx_*.out) whose
# results are in the manuscript. Written for the code folder Sato2026_github; it omits the steps of 114 that
# only collected numbers for documents (67, 71, 74, 89). The commands and settings are those of 114.
#SBATCH -J lineage_mk_excl
#SBATCH -p zhao_a
#SBATCH -A zhao_condo_bank
#SBATCH -n 5
#SBATCH -c 1
#SBATCH --mem-per-cpu=16G
#SBATCH -t 8:00:00
set -euo pipefail
unset LD_LIBRARY_PATH
P=/ru-auth/local/home/lzhao/vscode_projects/suzukii
PY=/ru-auth/local/home/lzhao/anaconda3/envs/jupyter-ml/bin/python
S=$P/results/polarized_mk_excl/summary
cd $P
mkdir -p $S logs/pmkx
# step 1: per-pool polarized counts (script 60), shared fixed differences excluded
for s in Dsub_L5_7 JP-Kan JP-Sap JP-Sap2 JP-Tok; do
  EXCLUDE_SHARED=1 srun --exact -n 1 -c 1 --mem-per-cpu=16G $PY scripts/60_polarized_mk.py $s > logs/pmkx/60_$s.log 2>&1 &
done
wait
# step 2: summaries
export PMK_VARIANT=excl
$PY scripts/64_summarize_polarized_mk.py > logs/pmkx/64.log 2>&1          # alpha, gene tests, overlap of lineages
$PY scripts/70_lineage_enrichment.py > logs/pmkx/70.log 2>&1              # lineage-specific enrichment (testis)
$PY scripts/75_per_arm_diversity_ratio.py > $S/per_arm_diversity_ratio.log   # four-fold arm ratios
$PY scripts/92_export_lineage_gene_lists.py > logs/pmkx/92.log 2>&1       # gene lists used by 124 and 130
$PY scripts/96_mk_per_site_rates.py > logs/pmkx/96.log 2>&1               # pN/pS and dN/dS per site

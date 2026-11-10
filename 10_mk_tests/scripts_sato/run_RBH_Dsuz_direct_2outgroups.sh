#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --array=0-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=0-08:00:00
#SBATCH -J RBH_Dsuz_direct
#SBATCH -o /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/RBH_Dsuz_direct.%A_%a.out
#SBATCH -e /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/RBH_Dsuz_direct.%A_%a.err

set -euo pipefail

source /ru-auth/local/home/asato/miniconda3/etc/profile.d/conda.sh
conda activate blast_env

BASE=/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817

PEPDIR=${BASE}/pep_longest
DBDIR=${BASE}/db_longest
BLASTDIR=${BASE}/blast_longest
ORTHDIR=${BASE}/orth_longest

mkdir -p \
  "${BLASTDIR}" \
  "${ORTHDIR}" \
  "${BASE}/logs" \
  "${BASE}/checks"

DSUZ=${PEPDIR}/Dsuz.IsoJpt1.0.longest.pep.fa

COMPS=("Dbia" "Dmel")
COMP=${COMPS[$SLURM_ARRAY_TASK_ID]}

case "${COMP}" in
  Dbia)
    OUTGROUP_FASTA=${PEPDIR}/Dbia.longest.pep.fa
    OUTGROUP_DB=${DBDIR}/Dbia
    ;;
  Dmel)
    OUTGROUP_FASTA=${PEPDIR}/Dmel.longest.pep.fa
    OUTGROUP_DB=${DBDIR}/Dmel
    ;;
  *)
    echo "ERROR: unsupported comparison: ${COMP}" >&2
    exit 1
    ;;
esac

DSUZ_DB=${DBDIR}/Dsuz_IsoJpt1.0

DSUZ_VS_OUT=${BLASTDIR}/Dsuz_direct_vs_${COMP}.tsv
OUT_VS_DSUZ=${BLASTDIR}/${COMP}_direct_vs_Dsuz.tsv

DSUZ_BEST=${BLASTDIR}/Dsuz.direct_best.${COMP}.tsv
OUT_BEST=${BLASTDIR}/${COMP}.direct_best.Dsuz.tsv

RBH_OUT=${ORTHDIR}/Dsuz_${COMP}.direct.RBH.longest.protein.tsv
SUMMARY=${BASE}/checks/Dsuz_${COMP}.direct_RBH.summary.tsv

for required in \
  "${DSUZ}" \
  "${OUTGROUP_FASTA}" \
  "${DSUZ_DB}.pin" \
  "${OUTGROUP_DB}.pin"
do
  if [[ ! -s "${required}" ]]; then
    echo "ERROR: required input is missing or empty: ${required}" >&2
    exit 1
  fi
done

echo "===== Direct Dsuz RBH analysis ====="
echo "Comparison: Dsuz versus ${COMP}"
echo "Dsuz FASTA: ${DSUZ}"
echo "Outgroup FASTA: ${OUTGROUP_FASTA}"
echo "Threads: ${SLURM_CPUS_PER_TASK}"
echo "BLAST E-value: 1e-5"

echo
echo "===== BLASTP: Dsuz versus ${COMP} ====="

blastp \
  -query "${DSUZ}" \
  -db "${OUTGROUP_DB}" \
  -out "${DSUZ_VS_OUT}" \
  -evalue 1e-5 \
  -outfmt 6 \
  -num_threads "${SLURM_CPUS_PER_TASK}"

echo
echo "===== BLASTP: ${COMP} versus Dsuz ====="

blastp \
  -query "${OUTGROUP_FASTA}" \
  -db "${DSUZ_DB}" \
  -out "${OUT_VS_DSUZ}" \
  -evalue 1e-5 \
  -outfmt 6 \
  -num_threads "${SLURM_CPUS_PER_TASK}"

best_hit() {
  local input_file="$1"
  local output_file="$2"

  awk 'BEGIN {OFS="\t"} {
    print $1, $2, $11, $12
  }' "${input_file}" |
  LC_ALL=C sort \
    -k1,1 \
    -k3,3g \
    -k4,4gr \
    -k2,2 |
  awk 'BEGIN {OFS="\t"} !seen[$1]++ {
    print $1, $2, $3, $4
  }' \
    > "${output_file}"
}

echo
echo "===== Select directional best hits ====="

best_hit \
  "${DSUZ_VS_OUT}" \
  "${DSUZ_BEST}"

best_hit \
  "${OUT_VS_DSUZ}" \
  "${OUT_BEST}"

extract_rbh() {
  local dsuz_best="$1"
  local outgroup_best="$2"
  local output_file="$3"

  awk '
  BEGIN {
    FS=OFS="\t"
  }

  NR==FNR {
    dsuz_to_outgroup[$1]=$2
    next
  }

  (($2 in dsuz_to_outgroup) &&
   dsuz_to_outgroup[$2] == $1) {
    print $1, $2
  }
  ' "${dsuz_best}" "${outgroup_best}" |
  LC_ALL=C sort -u \
    > "${output_file}"
}

echo
echo "===== Extract direct RBHs ====="

# Output columns:
#   outgroup protein ID
#   Dsuz protein ID

extract_rbh \
  "${DSUZ_BEST}" \
  "${OUT_BEST}" \
  "${RBH_OUT}"

RBH_COUNT=$(wc -l < "${RBH_OUT}")

DUP_OUTGROUP=$(
  cut -f1 "${RBH_OUT}" |
  sort |
  uniq -d |
  wc -l
)

DUP_DSUZ=$(
  cut -f2 "${RBH_OUT}" |
  sort |
  uniq -d |
  wc -l
)

{
  echo -e "metric\tvalue"
  echo -e "comparison\tDsuz_vs_${COMP}"
  echo -e "Dsuz_protein_records\t$(grep -c '^>' "${DSUZ}")"
  echo -e "outgroup_protein_records\t$(grep -c '^>' "${OUTGROUP_FASTA}")"
  echo -e "direct_RBH_pairs\t${RBH_COUNT}"
  echo -e "duplicated_outgroup_IDs\t${DUP_OUTGROUP}"
  echo -e "duplicated_Dsuz_IDs\t${DUP_DSUZ}"
  echo -e "RBH_output\t${RBH_OUT}"
} > "${SUMMARY}"

cat "${SUMMARY}"

if [[ "${DUP_OUTGROUP}" -ne 0 ||
      "${DUP_DSUZ}" -ne 0 ]]
then
  echo "ERROR: duplicated IDs detected in RBH table" >&2
  exit 1
fi

echo
echo "Done: direct Dsuz versus ${COMP} RBH"

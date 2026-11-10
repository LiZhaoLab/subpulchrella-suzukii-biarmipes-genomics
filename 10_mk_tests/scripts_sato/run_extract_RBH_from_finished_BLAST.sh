#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-00:10:00
#SBATCH -J extract_RBH
#SBATCH -o /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/extract_RBH.%j.out
#SBATCH -e /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/extract_RBH.%j.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817
BLASTDIR=${BASE}/blast_longest
ORTHDIR=${BASE}/orth_longest

mkdir -p "${ORTHDIR}"

extract_rbh() {
    local dsub_best="$1"
    local outgroup_best="$2"
    local output_file="$3"

    awk '
    BEGIN {
        FS=OFS="\t"
    }

    NR==FNR {
        dsub_to_outgroup[$1]=$2
        next
    }

    (($2 in dsub_to_outgroup) && dsub_to_outgroup[$2] == $1) {
        print $1, $2
    }
    ' "${dsub_best}" "${outgroup_best}" |
    LC_ALL=C sort -u \
      > "${output_file}"
}

for f in \
  "${BLASTDIR}/Dsub.best.Dsuz.tsv" \
  "${BLASTDIR}/Dsuz.best.Dsub.tsv" \
  "${BLASTDIR}/Dsub.best.Dbia.tsv" \
  "${BLASTDIR}/Dbia.best.Dsub.tsv" \
  "${BLASTDIR}/Dsub.best.Dmel.tsv" \
  "${BLASTDIR}/Dmel.best.Dsub.tsv"
do
    if [[ ! -s "${f}" ]]; then
        echo "ERROR: missing or empty best-hit file: ${f}" >&2
        exit 1
    fi
done

extract_rbh \
  "${BLASTDIR}/Dsub.best.Dsuz.tsv" \
  "${BLASTDIR}/Dsuz.best.Dsub.tsv" \
  "${ORTHDIR}/Dsub_Dsuz_IsoJpt1.0.RBH.longest.protein.tsv"

extract_rbh \
  "${BLASTDIR}/Dsub.best.Dbia.tsv" \
  "${BLASTDIR}/Dbia.best.Dsub.tsv" \
  "${ORTHDIR}/Dsub_Dbia.RBH.longest.protein.tsv"

extract_rbh \
  "${BLASTDIR}/Dsub.best.Dmel.tsv" \
  "${BLASTDIR}/Dmel.best.Dsub.tsv" \
  "${ORTHDIR}/Dsub_Dmel.RBH.longest.protein.tsv"

echo "===== RBH counts ====="

for f in \
  "${ORTHDIR}/Dsub_Dsuz_IsoJpt1.0.RBH.longest.protein.tsv" \
  "${ORTHDIR}/Dsub_Dbia.RBH.longest.protein.tsv" \
  "${ORTHDIR}/Dsub_Dmel.RBH.longest.protein.tsv"
do
    printf "%s\t" "$(basename "${f}")"
    wc -l < "${f}"
done

echo
echo "===== Duplicate-ID checks ====="

for f in \
  "${ORTHDIR}/Dsub_Dsuz_IsoJpt1.0.RBH.longest.protein.tsv" \
  "${ORTHDIR}/Dsub_Dbia.RBH.longest.protein.tsv" \
  "${ORTHDIR}/Dsub_Dmel.RBH.longest.protein.tsv"
do
    echo "$(basename "${f}")"

    printf "duplicated outgroup IDs\t"
    cut -f1 "${f}" |
      LC_ALL=C sort |
      uniq -d |
      wc -l

    printf "duplicated Dsub IDs\t"
    cut -f2 "${f}" |
      LC_ALL=C sort |
      uniq -d |
      wc -l
done

echo
echo "===== Orco check ====="

ORCO_PROTEIN=XP_037729344.1

for f in \
  "${ORTHDIR}/Dsub_Dsuz_IsoJpt1.0.RBH.longest.protein.tsv" \
  "${ORTHDIR}/Dsub_Dbia.RBH.longest.protein.tsv" \
  "${ORTHDIR}/Dsub_Dmel.RBH.longest.protein.tsv"
do
    echo "$(basename "${f}")"

    if ! grep -w "${ORCO_PROTEIN}" "${f}"; then
        echo "Orco not present in this RBH comparison"
    fi
done

echo
echo "Done."

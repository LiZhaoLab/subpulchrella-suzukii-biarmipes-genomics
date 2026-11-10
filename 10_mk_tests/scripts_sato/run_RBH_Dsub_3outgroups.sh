#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=0-08:00:00
#SBATCH -J RBH_Dsub_3out
#SBATCH -o /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/RBH_Dsub_3out.%j.out
#SBATCH -e /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/RBH_Dsub_3out.%j.err

set -euo pipefail

source /ru-auth/local/home/asato/miniconda3/etc/profile.d/conda.sh
conda activate blast_env

BASE=/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817

PEPDIR=${BASE}/pep_longest
DBDIR=${BASE}/db_longest
BLASTDIR=${BASE}/blast_longest
ORTHDIR=${BASE}/orth_longest

Dsub=${PEPDIR}/Dsub.longest_per_gene.pep.fa
Dsuz=${PEPDIR}/Dsuz.IsoJpt1.0.longest.pep.fa
Dbia=${PEPDIR}/Dbia.longest.pep.fa
Dmel=${PEPDIR}/Dmel.longest.pep.fa

mkdir -p \
  "${DBDIR}" \
  "${BLASTDIR}" \
  "${ORTHDIR}"

for f in "${Dsub}" "${Dsuz}" "${Dbia}" "${Dmel}"; do
    if [[ ! -s "${f}" ]]; then
        echo "ERROR: missing or empty protein FASTA: ${f}" >&2
        exit 1
    fi
done

echo "===== Input protein counts ====="

for f in "${Dsub}" "${Dsuz}" "${Dbia}" "${Dmel}"; do
    printf "%s\t" "$(basename "${f}")"
    grep -c '^>' "${f}"
done

echo
echo "===== Build BLAST databases ====="

makeblastdb \
  -in "${Dsub}" \
  -dbtype prot \
  -out "${DBDIR}/Dsub"

makeblastdb \
  -in "${Dsuz}" \
  -dbtype prot \
  -out "${DBDIR}/Dsuz_IsoJpt1.0"

makeblastdb \
  -in "${Dbia}" \
  -dbtype prot \
  -out "${DBDIR}/Dbia"

makeblastdb \
  -in "${Dmel}" \
  -dbtype prot \
  -out "${DBDIR}/Dmel"

echo
echo "===== Dsub versus Dsuz ====="

blastp \
  -query "${Dsub}" \
  -db "${DBDIR}/Dsuz_IsoJpt1.0" \
  -out "${BLASTDIR}/Dsub_vs_Dsuz.tsv" \
  -evalue 1e-5 \
  -outfmt 6 \
  -num_threads "${SLURM_CPUS_PER_TASK}"

blastp \
  -query "${Dsuz}" \
  -db "${DBDIR}/Dsub" \
  -out "${BLASTDIR}/Dsuz_vs_Dsub.tsv" \
  -evalue 1e-5 \
  -outfmt 6 \
  -num_threads "${SLURM_CPUS_PER_TASK}"

echo
echo "===== Dsub versus Dbia ====="

blastp \
  -query "${Dsub}" \
  -db "${DBDIR}/Dbia" \
  -out "${BLASTDIR}/Dsub_vs_Dbia.tsv" \
  -evalue 1e-5 \
  -outfmt 6 \
  -num_threads "${SLURM_CPUS_PER_TASK}"

blastp \
  -query "${Dbia}" \
  -db "${DBDIR}/Dsub" \
  -out "${BLASTDIR}/Dbia_vs_Dsub.tsv" \
  -evalue 1e-5 \
  -outfmt 6 \
  -num_threads "${SLURM_CPUS_PER_TASK}"

echo
echo "===== Dsub versus Dmel ====="

blastp \
  -query "${Dsub}" \
  -db "${DBDIR}/Dmel" \
  -out "${BLASTDIR}/Dsub_vs_Dmel.tsv" \
  -evalue 1e-5 \
  -outfmt 6 \
  -num_threads "${SLURM_CPUS_PER_TASK}"

blastp \
  -query "${Dmel}" \
  -db "${DBDIR}/Dsub" \
  -out "${BLASTDIR}/Dmel_vs_Dsub.tsv" \
  -evalue 1e-5 \
  -outfmt 6 \
  -num_threads "${SLURM_CPUS_PER_TASK}"


# ============================================================
# Select one best hit per query
#
# Same basic criterion as the previous analysis:
#   lowest BLAST E-value
#
# Deterministic tie handling:
#   highest bit score
#   lexicographically smallest subject ID
# ============================================================

best_hit() {
    local input_file="$1"
    local output_file="$2"

    awk 'BEGIN{OFS="\t"} {
        print $1, $2, $11, $12
    }' "${input_file}" |
    LC_ALL=C sort \
      -k1,1 \
      -k3,3g \
      -k4,4gr \
      -k2,2 |
    awk 'BEGIN{OFS="\t"} !seen[$1]++ {
        print $1, $2, $3, $4
    }' \
      > "${output_file}"
}


echo
echo "===== Select directional best hits ====="

best_hit \
  "${BLASTDIR}/Dsub_vs_Dsuz.tsv" \
  "${BLASTDIR}/Dsub.best.Dsuz.tsv"

best_hit \
  "${BLASTDIR}/Dsuz_vs_Dsub.tsv" \
  "${BLASTDIR}/Dsuz.best.Dsub.tsv"

best_hit \
  "${BLASTDIR}/Dsub_vs_Dbia.tsv" \
  "${BLASTDIR}/Dsub.best.Dbia.tsv"

best_hit \
  "${BLASTDIR}/Dbia_vs_Dsub.tsv" \
  "${BLASTDIR}/Dbia.best.Dsub.tsv"

best_hit \
  "${BLASTDIR}/Dsub_vs_Dmel.tsv" \
  "${BLASTDIR}/Dsub.best.Dmel.tsv"

best_hit \
  "${BLASTDIR}/Dmel_vs_Dsub.tsv" \
  "${BLASTDIR}/Dmel.best.Dsub.tsv"


# ============================================================
# Extract reciprocal best hits
#
# Output columns:
#   outgroup protein ID
#   Dsub protein ID
# ============================================================

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

echo
echo "===== Extract RBHs ====="

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

# ============================================================
# Validate RBH tables
# ============================================================

echo
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

    echo -n "  duplicated outgroup IDs: "
    cut -f1 "${f}" |
      sort |
      uniq -d |
      wc -l

    echo -n "  duplicated Dsub IDs: "
    cut -f2 "${f}" |
      sort |
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

    if grep -w "${ORCO_PROTEIN}" "${f}"; then
        :
    else
        echo "  Orco not present in this RBH comparison"
    fi
done

echo
echo "Done."

#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --array=0-11%2
#SBATCH --cpus-per-task=36
#SBATCH --mem=32G
#SBATCH --time=7-00:00:00
#SBATCH -J prank_Dsuz
#SBATCH -o /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/prank_Dsuz.%A_%a.out
#SBATCH -e /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/prank_Dsuz.%A_%a.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817

PRANK=/ru-auth/local/home/asato/bin/prank

POPS=("JP-Kan" "JP-Sap" "JP-Sap2" "JP-Tok")
COMPS=("Dsub" "Dbia" "Dmel")

POP_INDEX=$((SLURM_ARRAY_TASK_ID / 3))
COMP_INDEX=$((SLURM_ARRAY_TASK_ID % 3))

POP=${POPS[$POP_INDEX]}
COMP=${COMPS[$COMP_INDEX]}

CONDITION="${POP}.Dsuz_vs_${COMP}"

INDIR=${BASE}/combined_fasta/Dsuz/${POP}/Dsuz_vs_${COMP}
OUTDIR=${BASE}/prank_out/Dsuz/${POP}/Dsuz_vs_${COMP}
LOGDIR=${BASE}/logs/prank_gene_logs/Dsuz/${POP}/Dsuz_vs_${COMP}
CHECKDIR=${BASE}/checks/Dsuz_PRANK

mkdir -p \
    "${OUTDIR}" \
    "${LOGDIR}" \
    "${CHECKDIR}" \
    "${BASE}/logs"

if [[ ! -d "${INDIR}" ]]; then
    echo "ERROR: input directory not found: ${INDIR}" >&2
    exit 1
fi

if [[ ! -x "${PRANK}" ]]; then
    echo "ERROR: PRANK executable not found: ${PRANK}" >&2
    exit 1
fi

echo "===== Dsuz PRANK analysis ====="
echo "Population: ${POP}"
echo "Comparison: Dsuz_vs_${COMP}"
echo "Input: ${INDIR}"
echo "Output: ${OUTDIR}"
echo "PRANK: ${PRANK}"
echo "PRANK options: -codon -F"
echo "Parallel genes: ${SLURM_CPUS_PER_TASK}"

run_prank() {
    local input_fasta="$1"
    local transcript
    local output_prefix
    local best_fasta
    local gene_log
    local failed_marker

    transcript=$(
        basename "${input_fasta}" .combined.fasta
    )

    output_prefix="${OUTDIR}/${transcript}"
    best_fasta="${output_prefix}.best.fas"
    gene_log="${LOGDIR}/${transcript}.prank.log"
    failed_marker="${LOGDIR}/${transcript}.failed"

    # 完了済み遺伝子は安全にスキップ
    if [[ -s "${best_fasta}" ]]; then
        rm -f "${failed_marker}"
        return 0
    fi

    rm -f "${failed_marker}"

    if "${PRANK}" \
        -d="${input_fasta}" \
        -o="${output_prefix}" \
        -codon \
        -F \
        > "${gene_log}" 2>&1
    then
        if [[ -s "${best_fasta}" ]]; then
            rm -f "${failed_marker}"
        else
            echo \
                "PRANK returned success but best.fas is missing" \
                > "${failed_marker}"
        fi
    else
        echo \
            "PRANK command failed; see ${gene_log}" \
            > "${failed_marker}"
    fi

    # 1遺伝子の失敗で全体を止めない
    return 0
}

export PRANK
export OUTDIR
export LOGDIR
export -f run_prank

find "${INDIR}" \
    -maxdepth 1 \
    -type f \
    -name "*.combined.fasta" \
    -print0 |
sort -z |
xargs -0 \
    -P "${SLURM_CPUS_PER_TASK}" \
    -I {} \
    bash -c 'run_prank "$1"' _ {}

INPUT_LIST=${CHECKDIR}/${CONDITION}.input_IDs.txt
OUTPUT_LIST=${CHECKDIR}/${CONDITION}.output_IDs.txt
FAILED_LIST=${CHECKDIR}/${CONDITION}.failed_IDs.txt
MISSING_LIST=${CHECKDIR}/${CONDITION}.missing_IDs.txt
UNEXPECTED_LIST=${CHECKDIR}/${CONDITION}.unexpected_IDs.txt
SUMMARY=${CHECKDIR}/${CONDITION}.summary.tsv

find "${INDIR}" \
    -maxdepth 1 \
    -type f \
    -name "*.combined.fasta" \
    -printf "%f\n" |
sed 's/\.combined\.fasta$//' |
sort -u \
    > "${INPUT_LIST}"

find "${OUTDIR}" \
    -maxdepth 1 \
    -type f \
    -name "*.best.fas" \
    -size +0c \
    -printf "%f\n" |
sed 's/\.best\.fas$//' |
sort -u \
    > "${OUTPUT_LIST}"

find "${LOGDIR}" \
    -maxdepth 1 \
    -type f \
    -name "*.failed" \
    -printf "%f\n" |
sed 's/\.failed$//' |
sort -u \
    > "${FAILED_LIST}"

comm -23 \
    "${INPUT_LIST}" \
    "${OUTPUT_LIST}" \
    > "${MISSING_LIST}"

comm -13 \
    "${INPUT_LIST}" \
    "${OUTPUT_LIST}" \
    > "${UNEXPECTED_LIST}"

INPUT_COUNT=$(wc -l < "${INPUT_LIST}")
OUTPUT_COUNT=$(wc -l < "${OUTPUT_LIST}")
FAILED_COUNT=$(wc -l < "${FAILED_LIST}")
MISSING_COUNT=$(wc -l < "${MISSING_LIST}")
UNEXPECTED_COUNT=$(wc -l < "${UNEXPECTED_LIST}")

ORCO=rna-NM_001328601.1

if grep -Fxq "${ORCO}" "${OUTPUT_LIST}"; then
    ORCO_STATUS=aligned
elif grep -Fxq "${ORCO}" "${INPUT_LIST}"; then
    ORCO_STATUS=input_but_not_aligned
else
    ORCO_STATUS=not_in_input
fi

{
    printf "metric\tvalue\n"
    printf "population\t%s\n" "${POP}"
    printf "comparison\tDsuz_vs_%s\n" "${COMP}"
    printf "input_combined_FASTA\t%s\n" "${INPUT_COUNT}"
    printf "successful_PRANK_alignments\t%s\n" "${OUTPUT_COUNT}"
    printf "failed_markers\t%s\n" "${FAILED_COUNT}"
    printf "missing_alignments\t%s\n" "${MISSING_COUNT}"
    printf "unexpected_outputs\t%s\n" "${UNEXPECTED_COUNT}"
    printf "Orco_status\t%s\n" "${ORCO_STATUS}"
    printf "input_directory\t%s\n" "${INDIR}"
    printf "output_directory\t%s\n" "${OUTDIR}"
    printf "PRANK_options\t-codon -F\n"
    printf "parallel_genes\t%s\n" "${SLURM_CPUS_PER_TASK}"
} > "${SUMMARY}"

echo
echo "===== Final summary ====="
cat "${SUMMARY}"

if [[ "${MISSING_COUNT}" -eq 0 ]] \
    && [[ "${UNEXPECTED_COUNT}" -eq 0 ]]
then
    echo "FINAL: PASS"
else
    echo "FINAL: COMPLETE_WITH_FAILURES"
fi

echo "Done."

#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --array=0-2
#SBATCH --cpus-per-task=36
#SBATCH --mem=32G
#SBATCH --time=7-00:00:00
#SBATCH -J prank_Dsub_resume
#SBATCH -o /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/prank_Dsub_resume.%A_%a.out
#SBATCH -e /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/prank_Dsub_resume.%A_%a.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817
PRANK=/ru-auth/local/home/asato/bin/prank

COMPARISONS=(
    "Dsub_vs_Dsuz"
    "Dsub_vs_Dbia"
    "Dsub_vs_Dmel"
)

COMP=${COMPARISONS[$SLURM_ARRAY_TASK_ID]}

INDIR=${BASE}/combined_fasta/${COMP}

# 前回と同じ出力先
OUTDIR=${BASE}/prank_out/${COMP}

LOGDIR=${BASE}/logs/prank_gene_logs/${COMP}
CHECKDIR=${BASE}/checks/Dsub_PRANK_final

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

echo "===== Resumable Dsub PRANK ====="
echo "Comparison: ${COMP}"
echo "Input: ${INDIR}"
echo "Output: ${OUTDIR}"
echo "Parallel genes: ${SLURM_CPUS_PER_TASK}"
echo "PRANK: ${PRANK}"
echo "PRANK options: -codon -F"

run_prank() {
    local input_fasta="$1"
    local transcript
    local output_prefix
    local best_fasta
    local gene_log
    local failed_marker
    local expected_marker

    transcript=$(
        basename "${input_fasta}" .combined.fasta
    )

    output_prefix="${OUTDIR}/${transcript}"
    best_fasta="${output_prefix}.best.fas"

    gene_log="${LOGDIR}/${transcript}.prank.log"
    failed_marker="${LOGDIR}/${transcript}.failed"
    expected_marker="${LOGDIR}/${transcript}.expected_codon_exclusion"

    # 正常終了済みは再計算しない
    if [[ -s "${best_fasta}" ]]; then
        rm -f \
            "${failed_marker}" \
            "${expected_marker}"
        return 0
    fi

    # 前回すでにコドン長不適合と確定したものも再計算しない
    if [[ -s "${expected_marker}" ]]; then
        rm -f "${failed_marker}"
        return 0
    fi

    rm -f \
        "${failed_marker}" \
        "${expected_marker}" \
        "${gene_log}"

    if "${PRANK}" \
        -d="${input_fasta}" \
        -o="${output_prefix}" \
        -codon \
        -F \
        > "${gene_log}" 2>&1
    then
        if [[ -s "${best_fasta}" ]]; then
            rm -f \
                "${failed_marker}" \
                "${expected_marker}"
        else
            echo \
                "PRANK returned success but best.fas is missing" \
                > "${failed_marker}"
        fi
    else
        if grep -Fq \
            "codon sequence length is not multiple of three" \
            "${gene_log}"
        then
            {
                echo "Expected exclusion"
                echo "Reason: codon sequence length is not multiple of three"
                echo "Input: ${input_fasta}"
            } > "${expected_marker}"

            rm -f "${failed_marker}"
        else
            echo \
                "Unexpected PRANK failure; see ${gene_log}" \
                > "${failed_marker}"

            rm -f "${expected_marker}"
        fi
    fi

    # 1遺伝子の失敗で残りを止めない
    return 0
}

export PRANK
export OUTDIR
export LOGDIR
export -f run_prank

# 全入力を確認するが、成功済みは関数内で自動スキップ
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

INPUT_LIST=${CHECKDIR}/${COMP}.input_IDs.txt
SUCCESS_LIST=${CHECKDIR}/${COMP}.successful_IDs.txt
EXPECTED_LIST=${CHECKDIR}/${COMP}.expected_codon_exclusions.txt
FAILED_LIST=${CHECKDIR}/${COMP}.unexpected_failure_IDs.txt
MISSING_LIST=${CHECKDIR}/${COMP}.unaccounted_IDs.txt
UNEXPECTED_OUTPUT_LIST=${CHECKDIR}/${COMP}.unexpected_output_IDs.txt
OVERLAP_LIST=${CHECKDIR}/${COMP}.success_expected_overlap_IDs.txt
ACCOUNTED_LIST=${CHECKDIR}/${COMP}.accounted_IDs.txt
SUMMARY=${CHECKDIR}/${COMP}.summary.tsv

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
    > "${SUCCESS_LIST}"

find "${LOGDIR}" \
    -maxdepth 1 \
    -type f \
    -name "*.expected_codon_exclusion" \
    -printf "%f\n" |
sed 's/\.expected_codon_exclusion$//' |
sort -u \
    > "${EXPECTED_LIST}"

find "${LOGDIR}" \
    -maxdepth 1 \
    -type f \
    -name "*.failed" \
    -printf "%f\n" |
sed 's/\.failed$//' |
sort -u \
    > "${FAILED_LIST}"

{
    cat "${SUCCESS_LIST}"
    cat "${EXPECTED_LIST}"
} |
sort -u \
    > "${ACCOUNTED_LIST}"

comm -23 \
    "${INPUT_LIST}" \
    "${ACCOUNTED_LIST}" \
    > "${MISSING_LIST}"

comm -13 \
    "${INPUT_LIST}" \
    "${SUCCESS_LIST}" \
    > "${UNEXPECTED_OUTPUT_LIST}"

comm -12 \
    "${SUCCESS_LIST}" \
    "${EXPECTED_LIST}" \
    > "${OVERLAP_LIST}"

INPUT_COUNT=$(wc -l < "${INPUT_LIST}")
SUCCESS_COUNT=$(wc -l < "${SUCCESS_LIST}")
EXPECTED_COUNT=$(wc -l < "${EXPECTED_LIST}")
FAILED_COUNT=$(wc -l < "${FAILED_LIST}")
MISSING_COUNT=$(wc -l < "${MISSING_LIST}")
UNEXPECTED_OUTPUT_COUNT=$(wc -l < "${UNEXPECTED_OUTPUT_LIST}")
OVERLAP_COUNT=$(wc -l < "${OVERLAP_LIST}")

ACCOUNTED_COUNT=$((SUCCESS_COUNT + EXPECTED_COUNT))

ORCO=rna-XM_037873416.1
TITIN=rna-XM_037865866.1

if grep -Fxq "${ORCO}" "${SUCCESS_LIST}"; then
    ORCO_STATUS=aligned
else
    ORCO_STATUS=not_aligned
fi

if grep -Fxq "${TITIN}" "${SUCCESS_LIST}"; then
    TITIN_STATUS=aligned
else
    TITIN_STATUS=not_aligned
fi

if [[ "${FAILED_COUNT}" -eq 0 ]] \
    && [[ "${MISSING_COUNT}" -eq 0 ]] \
    && [[ "${UNEXPECTED_OUTPUT_COUNT}" -eq 0 ]] \
    && [[ "${OVERLAP_COUNT}" -eq 0 ]] \
    && [[ "${ACCOUNTED_COUNT}" -eq "${INPUT_COUNT}" ]]
then
    FINAL_STATUS=PASS
else
    FINAL_STATUS=FAIL
fi

{
    printf "metric\tvalue\n"
    printf "comparison\t%s\n" "${COMP}"
    printf "input_combined_FASTA\t%s\n" "${INPUT_COUNT}"
    printf "successful_PRANK_alignments\t%s\n" "${SUCCESS_COUNT}"
    printf "expected_codon_exclusions\t%s\n" "${EXPECTED_COUNT}"
    printf "accounted_for\t%s\n" "${ACCOUNTED_COUNT}"
    printf "unexpected_failures\t%s\n" "${FAILED_COUNT}"
    printf "unaccounted_inputs\t%s\n" "${MISSING_COUNT}"
    printf "unexpected_outputs\t%s\n" "${UNEXPECTED_OUTPUT_COUNT}"
    printf "success_expected_overlap\t%s\n" "${OVERLAP_COUNT}"
    printf "Orco_status\t%s\n" "${ORCO_STATUS}"
    printf "titin_status\t%s\n" "${TITIN_STATUS}"
    printf "PRANK_options\t-codon -F\n"
    printf "parallel_genes\t%s\n" "${SLURM_CPUS_PER_TASK}"
    printf "final_status\t%s\n" "${FINAL_STATUS}"
} > "${SUMMARY}"

echo
echo "===== Final summary ====="
cat "${SUMMARY}"

if [[ "${FINAL_STATUS}" != "PASS" ]]; then
    echo "FINAL: FAIL" >&2

    echo "Unexpected failures:"
    cat "${FAILED_LIST}"

    echo "Unaccounted inputs:"
    cat "${MISSING_LIST}"

    exit 1
fi

echo "FINAL: PASS"
echo "Done."

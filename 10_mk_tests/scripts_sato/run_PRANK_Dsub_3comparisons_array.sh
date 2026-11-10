#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --array=0-2
#SBATCH --cpus-per-task=36
#SBATCH --mem=32G
#SBATCH --time=7-00:00:00
#SBATCH -J prank_Dsub_new
#SBATCH -o /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/prank_Dsub_new.%A_%a.out
#SBATCH -e /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/prank_Dsub_new.%A_%a.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817

COMPARISONS=(
    "Dsub_vs_Dsuz"
    "Dsub_vs_Dbia"
    "Dsub_vs_Dmel"
)

COMPARISON="${COMPARISONS[$SLURM_ARRAY_TASK_ID]}"

INDIR="${BASE}/combined_fasta/${COMPARISON}"
OUTDIR="${BASE}/prank_out/${COMPARISON}"
CHECKDIR="${BASE}/checks"
TASK_LOGDIR="${BASE}/logs/prank_gene_logs/${COMPARISON}"

mkdir -p \
    "${OUTDIR}" \
    "${CHECKDIR}" \
    "${TASK_LOGDIR}"

if [[ ! -d "${INDIR}" ]]; then
    echo "ERROR: input directory not found: ${INDIR}" >&2
    exit 1
fi

if ! command -v prank >/dev/null 2>&1; then
    echo "ERROR: prank was not found in PATH" >&2
    exit 1
fi

echo "===== PRANK analysis ====="
echo "Comparison: ${COMPARISON}"
echo "Input: ${INDIR}"
echo "Output: ${OUTDIR}"
echo "Parallel processes: ${SLURM_CPUS_PER_TASK}"
echo "PRANK executable: $(command -v prank)"
echo "PRANK options: -codon -F"

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
    gene_log="${TASK_LOGDIR}/${transcript}.prank.log"
    failed_marker="${TASK_LOGDIR}/${transcript}.failed"

    # Allow safe resumption without rerunning completed genes.
    if [[ -s "${best_fasta}" ]]; then
        rm -f "${failed_marker}"
        return 0
    fi

    rm -f "${failed_marker}"

    if prank \
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

    # Return zero here so all genes are attempted.
    # The complete failure audit is performed afterward.
    return 0
}

export OUTDIR
export TASK_LOGDIR
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

INPUT_LIST="${CHECKDIR}/${COMPARISON}.PRANK_input_IDs.txt"
OUTPUT_LIST="${CHECKDIR}/${COMPARISON}.PRANK_output_IDs.txt"
MISSING_LIST="${CHECKDIR}/${COMPARISON}.PRANK_missing_IDs.txt"
UNEXPECTED_LIST="${CHECKDIR}/${COMPARISON}.PRANK_unexpected_IDs.txt"
SUMMARY="${CHECKDIR}/${COMPARISON}.PRANK_summary.tsv"

find "${INDIR}" \
    -maxdepth 1 \
    -type f \
    -name "*.combined.fasta" \
    -printf "%f\n" |
sed 's/\.combined\.fasta$//' |
LC_ALL=C sort -u \
    > "${INPUT_LIST}"

find "${OUTDIR}" \
    -maxdepth 1 \
    -type f \
    -name "*.best.fas" \
    -size +0c \
    -printf "%f\n" |
sed 's/\.best\.fas$//' |
LC_ALL=C sort -u \
    > "${OUTPUT_LIST}"

comm -23 \
    "${INPUT_LIST}" \
    "${OUTPUT_LIST}" \
    > "${MISSING_LIST}"

comm -13 \
    "${INPUT_LIST}" \
    "${OUTPUT_LIST}" \
    > "${UNEXPECTED_LIST}"

N_INPUT=$(wc -l < "${INPUT_LIST}")
N_OUTPUT=$(wc -l < "${OUTPUT_LIST}")
N_MISSING=$(wc -l < "${MISSING_LIST}")
N_UNEXPECTED=$(wc -l < "${UNEXPECTED_LIST}")

N_FAILED_MARKERS=$(
    find "${TASK_LOGDIR}" \
        -maxdepth 1 \
        -type f \
        -name "*.failed" |
    wc -l
)

ORCO="rna-XM_037873416.1"
ORCO_OUTPUT="${OUTDIR}/${ORCO}.best.fas"

if [[ -s "${ORCO_OUTPUT}" ]]; then
    ORCO_STATUS="present"
else
    ORCO_STATUS="missing"
fi

{
    echo -e \
        "comparison\tinput_combined_FASTA\tbest_fas\tmissing\tunexpected\tfailed_markers\tOrco"
    echo -e \
        "${COMPARISON}\t${N_INPUT}\t${N_OUTPUT}\t${N_MISSING}\t${N_UNEXPECTED}\t${N_FAILED_MARKERS}\t${ORCO_STATUS}"
} > "${SUMMARY}"

echo
echo "===== PRANK summary ====="
cat "${SUMMARY}"

echo
echo "Missing ID list: ${MISSING_LIST}"
echo "Unexpected ID list: ${UNEXPECTED_LIST}"
echo "Gene logs: ${TASK_LOGDIR}"

if [[ "${N_MISSING}" -ne 0 ]] \
    || [[ "${N_UNEXPECTED}" -ne 0 ]] \
    || [[ "${N_FAILED_MARKERS}" -ne 0 ]] \
    || [[ "${ORCO_STATUS}" != "present" ]]
then
    echo \
        "ERROR: PRANK validation failed for ${COMPARISON}" \
        >&2
    exit 1
fi

echo "Done: ${COMPARISON}"

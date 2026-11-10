#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --array=0-2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0-04:00:00
#SBATCH -J retry_Dsub_PRANK
#SBATCH -o /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/retry_Dsub_PRANK.%A_%a.out
#SBATCH -e /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/retry_Dsub_PRANK.%A_%a.err

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
OUTDIR=${BASE}/prank_out/${COMP}
OLD_LOGDIR=${BASE}/logs/prank_gene_logs/${COMP}
RETRY_LOGDIR=${BASE}/logs/prank_retry_logs/${COMP}
CHECKDIR=${BASE}/checks
MISSING_LIST=${CHECKDIR}/${COMP}.PRANK_missing_IDs.txt

EXPECTED_LIST=${CHECKDIR}/${COMP}.PRANK_expected_codon_exclusions.txt
UNEXPECTED_LIST=${CHECKDIR}/${COMP}.PRANK_unexpected_failures_after_retry.txt
RECOVERED_LIST=${CHECKDIR}/${COMP}.PRANK_recovered_after_retry.txt
FINAL_SUMMARY=${CHECKDIR}/${COMP}.PRANK_finalized_summary.tsv

mkdir -p \
    "${OUTDIR}" \
    "${RETRY_LOGDIR}" \
    "${CHECKDIR}"

if [[ ! -s "${MISSING_LIST}" ]]; then
    echo "No missing PRANK alignments for ${COMP}"
    : > "${EXPECTED_LIST}"
    : > "${UNEXPECTED_LIST}"
    : > "${RECOVERED_LIST}"

    {
        echo -e "metric\tvalue"
        echo -e "comparison\t${COMP}"
        echo -e "original_missing\t0"
        echo -e "recovered_after_retry\t0"
        echo -e "expected_codon_exclusions\t0"
        echo -e "unexpected_failures\t0"
        echo -e "final_status\tPASS"
    } > "${FINAL_SUMMARY}"

    echo "FINAL: PASS"
    exit 0
fi

: > "${EXPECTED_LIST}"
: > "${UNEXPECTED_LIST}"
: > "${RECOVERED_LIST}"

ORIGINAL_MISSING=$(awk 'NF{n++} END{print n+0}' "${MISSING_LIST}")

while read -r transcript; do
    [[ -z "${transcript}" ]] && continue

    input_fasta=${INDIR}/${transcript}.combined.fasta
    output_prefix=${OUTDIR}/${transcript}
    best_fasta=${output_prefix}.best.fas
    retry_log=${RETRY_LOGDIR}/${transcript}.prank.log
    failed_marker=${OLD_LOGDIR}/${transcript}.failed

    echo "===== ${transcript} ====="

    # 既に成功していれば再実行しない
    if [[ -s "${best_fasta}" ]]; then
        echo "Already completed; skip"
        echo "${transcript}" >> "${RECOVERED_LIST}"
        rm -f "${failed_marker}"
        continue
    fi

    if [[ ! -s "${input_fasta}" ]]; then
        echo "${transcript}" >> "${UNEXPECTED_LIST}"
        echo "ERROR: missing input FASTA: ${input_fasta}" >&2
        continue
    fi

    rm -f "${retry_log}"

    if "${PRANK}" \
        -d="${input_fasta}" \
        -o="${output_prefix}" \
        -codon \
        -F \
        > "${retry_log}" 2>&1
    then
        if [[ -s "${best_fasta}" ]]; then
            echo "${transcript}" >> "${RECOVERED_LIST}"
            rm -f "${failed_marker}"
            echo "Recovered successfully"
        else
            echo "${transcript}" >> "${UNEXPECTED_LIST}"
            echo "PRANK returned success but best.fas is absent" >&2
        fi
    else
        if grep -Fq \
            "codon sequence length is not multiple of three" \
            "${retry_log}"
        then
            echo "${transcript}" >> "${EXPECTED_LIST}"
            echo "Expected exclusion: codon length not divisible by three"
        else
            echo "${transcript}" >> "${UNEXPECTED_LIST}"
            echo "Unexpected PRANK failure" >&2
        fi
    fi
done < "${MISSING_LIST}"

sort -u -o "${EXPECTED_LIST}" "${EXPECTED_LIST}"
sort -u -o "${UNEXPECTED_LIST}" "${UNEXPECTED_LIST}"
sort -u -o "${RECOVERED_LIST}" "${RECOVERED_LIST}"

RECOVERED=$(
    awk 'NF{n++} END{print n+0}' "${RECOVERED_LIST}"
)

EXPECTED=$(
    awk 'NF{n++} END{print n+0}' "${EXPECTED_LIST}"
)

UNEXPECTED=$(
    awk 'NF{n++} END{print n+0}' "${UNEXPECTED_LIST}"
)

ACCOUNTED_FOR=$((RECOVERED + EXPECTED + UNEXPECTED))

if [[ "${UNEXPECTED}" -eq 0 ]] \
    && [[ "${ACCOUNTED_FOR}" -eq "${ORIGINAL_MISSING}" ]]
then
    FINAL_STATUS=PASS
else
    FINAL_STATUS=FAIL
fi

{
    echo -e "metric\tvalue"
    echo -e "comparison\t${COMP}"
    echo -e "original_missing\t${ORIGINAL_MISSING}"
    echo -e "recovered_after_retry\t${RECOVERED}"
    echo -e "expected_codon_exclusions\t${EXPECTED}"
    echo -e "unexpected_failures\t${UNEXPECTED}"
    echo -e "accounted_for\t${ACCOUNTED_FOR}"
    echo -e "final_status\t${FINAL_STATUS}"
} > "${FINAL_SUMMARY}"

cat "${FINAL_SUMMARY}"

if [[ "${FINAL_STATUS}" != "PASS" ]]; then
    echo "FINAL: FAIL" >&2
    exit 1
fi

echo "FINAL: PASS"
echo "Done."

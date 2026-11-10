#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --array=0-2
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=0-04:00:00
#SBATCH -J cmktest_Dsub_new
#SBATCH -o /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/cmktest_Dsub_new.%A_%a.out
#SBATCH -e /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/cmktest_Dsub_new.%A_%a.err

set -euo pipefail

export PATH="/ru-auth/local/home/asato/bin:${PATH}"

BASE=/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817

COMPARISONS=(
    "Dsub_vs_Dsuz"
    "Dsub_vs_Dbia"
    "Dsub_vs_Dmel"
)

OUTGROUPS=(
    "dsuz"
    "dbia"
    "dmel"
)

COMP="${COMPARISONS[$SLURM_ARRAY_TASK_ID]}"
OUTSP="${OUTGROUPS[$SLURM_ARRAY_TASK_ID]}"

WORKDIR="${BASE}/prank_out/${COMP}"
OUTDIR="${BASE}/cmktest_results/Dsub_pairwise"
TMPDIR="${OUTDIR}/tmp_${COMP}"
OUTCSV="${OUTDIR}/${COMP}_mktest.csv"
LIST="${TMPDIR}/${COMP}.bestfas.list"
FAIL_LIST="${OUTDIR}/${COMP}.cmktest_failed.txt"
SUMMARY="${OUTDIR}/${COMP}.cmktest.summary.tsv"

mkdir -p \
    "${OUTDIR}" \
    "${TMPDIR}" \
    "${BASE}/logs"

if ! command -v cmktest >/dev/null 2>&1; then
    echo "ERROR: cmktest was not found in PATH" >&2
    exit 1
fi

if [[ ! -d "${WORKDIR}" ]]; then
    echo "ERROR: PRANK output directory not found: ${WORKDIR}" >&2
    exit 1
fi

find "${WORKDIR}" \
    -maxdepth 1 \
    -type f \
    -name "*.best.fas" \
    -printf "%f\n" |
LC_ALL=C sort > "${LIST}"

INPUT_COUNT=$(awk 'NF{n++} END{print n+0}' "${LIST}")

if [[ "${INPUT_COUNT}" -eq 0 ]]; then
    echo "ERROR: no PRANK .best.fas files found in ${WORKDIR}" >&2
    exit 1
fi

echo "===== cMKtest analysis ====="
echo "Comparison: ${COMP}"
echo "Focal species: Dsub"
echo "Outgroup: ${OUTSP}"
echo "Cutoff: 0.05"
echo "Parallel processes: ${SLURM_CPUS_PER_TASK}"
echo "PRANK alignments: ${INPUT_COUNT}"
echo "Input directory: ${WORKDIR}"
echo "Output CSV: ${OUTCSV}"

run_one() {
    local filename="$1"
    local transcript
    local output
    local log
    local result

    transcript="${filename%.best.fas}"
    output="${TMPDIR}/${transcript}.csv"
    log="${TMPDIR}/${transcript}.log"

    # 再実行時には正常な既存結果をスキップする
    if [[ -s "${output}" ]]; then
        return 0
    fi

    rm -f "${log}"

    result=$(
        cd "${WORKDIR}" &&
        cmktest \
            -f "${filename}" \
            -i Dsub \
            -o "${OUTSP}" \
            -c 0.05 \
            > "${log}" 2>&1 &&
        tail -n 1 "${log}"
    ) || return 1

    if [[ -n "${result}" ]] &&
       [[ "${result}" != "Name,Ncodon,Npopulation,Pn,Dn,Ps,Ds,Alpha,Pvalue" ]]
    then
        printf "%s\n" "${result}" > "${output}"
    else
        return 1
    fi
}

export WORKDIR TMPDIR OUTSP
export -f run_one

xargs -r -I {} -P "${SLURM_CPUS_PER_TASK}" \
    bash -c 'run_one "$1"' _ {} \
    < "${LIST}" || true

: > "${FAIL_LIST}"

while read -r filename; do
    [[ -z "${filename}" ]] && continue

    transcript="${filename%.best.fas}"

    if [[ ! -s "${TMPDIR}/${transcript}.csv" ]]; then
        printf "%s\n" "${transcript}" >> "${FAIL_LIST}"
    fi
done < "${LIST}"

SUCCESS_COUNT=$(
    find "${TMPDIR}" \
        -maxdepth 1 \
        -type f \
        -name "*.csv" |
    wc -l
)

FAILED_COUNT=$(awk 'NF{n++} END{print n+0}' "${FAIL_LIST}")

{
    echo "Name,Ncodon,Npopulation,Pn,Dn,Ps,Ds,Alpha,Pvalue"

    find "${TMPDIR}" \
        -maxdepth 1 \
        -type f \
        -name "*.csv" \
        -print0 |
    LC_ALL=C sort -z |
    xargs -0 -r cat
} > "${OUTCSV}"

OUTPUT_ROWS=$(awk 'END{print NR-1}' "${OUTCSV}")

if [[ "${SUCCESS_COUNT}" -eq "${INPUT_COUNT}" ]] &&
   [[ "${FAILED_COUNT}" -eq 0 ]] &&
   [[ "${OUTPUT_ROWS}" -eq "${INPUT_COUNT}" ]]
then
    FINAL_STATUS=PASS
else
    FINAL_STATUS=FAIL
fi

{
    printf "metric\tvalue\n"
    printf "comparison\t%s\n" "${COMP}"
    printf "focal_species\tDsub\n"
    printf "outgroup\t%s\n" "${OUTSP}"
    printf "cutoff\t0.05\n"
    printf "input_PRANK_alignments\t%d\n" "${INPUT_COUNT}"
    printf "successful_cmktests\t%d\n" "${SUCCESS_COUNT}"
    printf "failed_cmktests\t%d\n" "${FAILED_COUNT}"
    printf "output_rows\t%d\n" "${OUTPUT_ROWS}"
    printf "output_CSV\t%s\n" "${OUTCSV}"
    printf "failed_list\t%s\n" "${FAIL_LIST}"
    printf "final_status\t%s\n" "${FINAL_STATUS}"
} > "${SUMMARY}"

cat "${SUMMARY}"

if [[ "${FINAL_STATUS}" != "PASS" ]]; then
    echo "FINAL: FAIL — inspect ${FAIL_LIST}" >&2
    exit 1
fi

echo "FINAL: PASS"
echo "Done."

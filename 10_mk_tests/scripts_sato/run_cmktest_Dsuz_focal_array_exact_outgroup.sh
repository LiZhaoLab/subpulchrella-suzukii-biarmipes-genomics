#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --array=0-11
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=0-04:00:00
#SBATCH -J cmktest_Dsuz_fix
#SBATCH -o /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/cmktest_Dsuz_fix.%A_%a.out
#SBATCH -e /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/cmktest_Dsuz_fix.%A_%a.err

set -euo pipefail

export PATH="/ru-auth/local/home/asato/bin:${PATH}"

BASE=/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817

TASKS=(
  "JP-Kan Dsuz_vs_Dsub dsub"
  "JP-Kan Dsuz_vs_Dbia dbia"
  "JP-Kan Dsuz_vs_Dmel dmel"
  "JP-Sap Dsuz_vs_Dsub dsub"
  "JP-Sap Dsuz_vs_Dbia dbia"
  "JP-Sap Dsuz_vs_Dmel dmel"
  "JP-Sap2 Dsuz_vs_Dsub dsub"
  "JP-Sap2 Dsuz_vs_Dbia dbia"
  "JP-Sap2 Dsuz_vs_Dmel dmel"
  "JP-Tok Dsuz_vs_Dsub dsub"
  "JP-Tok Dsuz_vs_Dbia dbia"
  "JP-Tok Dsuz_vs_Dmel dmel"
)

read -r POP COMP OUTSP <<< "${TASKS[$SLURM_ARRAY_TASK_ID]}"

WORKDIR=${BASE}/prank_out/Dsuz/${POP}/${COMP}
OUTDIR=${BASE}/cmktest_results/Dsuz_pairwise/${POP}

# This new directory prevents reuse of the incorrect pre-fix per-gene CSVs.
TMPDIR=${OUTDIR}/tmp_${COMP}_exact_outgroup_header
INPUTDIR=${TMPDIR}/normalized_inputs

OUTCSV=${OUTDIR}/${COMP}_mktest.csv
LIST=${TMPDIR}/${COMP}.bestfas.list
FAILED=${OUTDIR}/${COMP}.cmktest_failed.txt
SUMMARY=${OUTDIR}/${COMP}.cmktest.summary.tsv

mkdir -p "${OUTDIR}" "${TMPDIR}" "${INPUTDIR}" "${BASE}/logs"

if ! command -v cmktest >/dev/null 2>&1; then
    echo "ERROR: cmktest not found in PATH" >&2
    exit 1
fi

if [[ ! -d "${WORKDIR}" ]]; then
    echo "ERROR: missing WORKDIR: ${WORKDIR}" >&2
    exit 1
fi

echo "===== Corrected Dsuz focal cMKtest ====="
echo "POP=${POP}"
echo "COMP=${COMP}"
echo "FOCAL=${POP}"
echo "OUTGROUP=${OUTSP}"
echo "WORKDIR=${WORKDIR}"
echo "TMPDIR=${TMPDIR}"
echo "OUTCSV=${OUTCSV}"
echo "cutoff=0.05"
echo "header_fix=>${OUTSP}|accession converted to >${OUTSP} in temporary input only"

find "${WORKDIR}" \
  -maxdepth 1 \
  -type f \
  -name "*.best.fas" \
  -printf "%f\n" |
sort > "${LIST}"

input_count=$(awk 'NF{n++} END{print n+0}' "${LIST}")

if [[ "${input_count}" -eq 0 ]]; then
    echo "ERROR: no PRANK best.fas inputs found in ${WORKDIR}" >&2
    exit 1
fi

run_one() {
    local f="$1"
    local base
    local source_fasta
    local normalized_fasta
    local out
    local failed_marker
    local original_outgroup_count
    local normalized_outgroup_count
    local result

    base="${f%.best.fas}"
    source_fasta=${WORKDIR}/${f}
    normalized_fasta=${INPUTDIR}/${f}
    out=${TMPDIR}/${base}.csv
    failed_marker=${TMPDIR}/${base}.failed

    # Resume safely only from outputs created by this corrected run.
    if [[ -s "${out}" ]]; then
        return 0
    fi

    rm -f "${failed_marker}" "${normalized_fasta}"

    if [[ ! -s "${source_fasta}" ]]; then
        echo "missing_input" > "${failed_marker}"
        return 0
    fi

    # Accept either the exact header or the accession-bearing header.
    original_outgroup_count=$(
      awk -v out="${OUTSP}" '
        $0 == ">" out || index($0, ">" out "|") == 1 {n++}
        END {print n+0}
      ' "${source_fasta}"
    )

    if [[ "${original_outgroup_count}" -ne 1 ]]; then
        echo "original_outgroup_header_count=${original_outgroup_count}" \
          > "${failed_marker}"
        return 0
    fi

    # Preserve all sequences and alignment columns; change only the one
    # outgroup FASTA header to the exact value passed with cmktest -o.
    awk -v out="${OUTSP}" '
      $0 == ">" out || index($0, ">" out "|") == 1 {
          print ">" out
          next
      }
      {print}
    ' "${source_fasta}" > "${normalized_fasta}"

    normalized_outgroup_count=$(
      awk -v out="${OUTSP}" '$0 == ">" out {n++} END {print n+0}' \
        "${normalized_fasta}"
    )

    if [[ "${normalized_outgroup_count}" -ne 1 ]]; then
        echo "normalized_outgroup_header_count=${normalized_outgroup_count}" \
          > "${failed_marker}"
        rm -f "${normalized_fasta}"
        return 0
    fi

    result=$(
      cmktest \
        -f "${normalized_fasta}" \
        -i "${POP}" \
        -o "${OUTSP}" \
        -c 0.05 |
      tail -n 1 || true
    )

    if [[ -n "${result}" ]] \
       && [[ "${result}" != "Name,Ncodon,Npopulation,Pn,Dn,Ps,Ds,Alpha,Pvalue" ]] \
       && [[ "$(awk -F, '{print NF}' <<< "${result}")" -eq 9 ]]
    then
        echo "${result}" > "${out}"
    else
        echo "invalid_or_missing_cmktest_result" > "${failed_marker}"
    fi

    rm -f "${normalized_fasta}"
}

export POP OUTSP WORKDIR TMPDIR INPUTDIR
export -f run_one

awk 'NF' "${LIST}" |
xargs -I {} -P "${SLURM_CPUS_PER_TASK}" bash -c 'run_one "$@"' _ {}

find "${TMPDIR}" \
  -maxdepth 1 \
  -type f \
  -name "*.failed" \
  -printf "%f\n" |
sed 's/\.failed$//' |
sort > "${FAILED}"

echo "Name,Ncodon,Npopulation,Pn,Dn,Ps,Ds,Alpha,Pvalue" > "${OUTCSV}"

while read -r f; do
    [[ -z "${f}" ]] && continue
    base=${f%.best.fas}
    per_gene_csv=${TMPDIR}/${base}.csv

    if [[ -s "${per_gene_csv}" ]]; then
        cat "${per_gene_csv}" >> "${OUTCSV}"
    fi
done < "${LIST}"

output_rows=$(awk 'END{print NR-1}' "${OUTCSV}")
failed_count=$(awk 'NF{n++} END{print n+0}' "${FAILED}")

{
    echo -e "metric\tvalue"
    echo -e "population\t${POP}"
    echo -e "comparison\t${COMP}"
    echo -e "focal_species\t${POP}"
    echo -e "outgroup\t${OUTSP}"
    echo -e "cutoff\t0.05"
    echo -e "outgroup_header_normalization\t>${OUTSP}|accession to >${OUTSP}"
    echo -e "input_PRANK_alignments\t${input_count}"
    echo -e "successful_cmktests\t${output_rows}"
    echo -e "failed_cmktests\t${failed_count}"
    echo -e "output_rows\t${output_rows}"
    echo -e "output_CSV\t${OUTCSV}"
    echo -e "failed_list\t${FAILED}"
} > "${SUMMARY}"

if [[ "${failed_count}" -eq 0 ]] \
   && [[ "${output_rows}" -eq "${input_count}" ]]
then
    echo -e "final_status\tPASS" >> "${SUMMARY}"
    cat "${SUMMARY}"
    echo "FINAL: PASS"
    echo "Done."
else
    echo -e "final_status\tFAIL" >> "${SUMMARY}"
    cat "${SUMMARY}"
    echo "FINAL: FAIL" >&2
    exit 1
fi

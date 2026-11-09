#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --job-name=syn_promer4
#SBATCH --array=0-2%3
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/syn_promer4.%A_%a.out
#SBATCH --error=logs/syn_promer4.%A_%a.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/synteny_4species/reanalysis_20260824
MUMMER_ENV=/ru-auth/local/home/asato/scratch/conda_envs/mummer4_env

PROMER=${BASE}/software_wrappers/promer_mummer4_4.0.1
DELTA_FILTER=${MUMMER_ENV}/bin/delta-filter
SHOW_COORDS=${MUMMER_ENV}/bin/show-coords

COMPARISONS=(
    Dsuz_vs_Dsub
    Dsub_vs_Dbia
    Dbia_vs_Dmel
)

REFERENCES=(
    "${BASE}/repeatmasker/Dsuz/Dsuz.GCF_043229965.1.fna.masked"
    "${BASE}/repeatmasker/Dsub/Dsub.GCF_014743375.2.fna.masked"
    "${BASE}/repeatmasker/Dbia/Dbia.GCF_025231255.1.fna.masked"
)

QUERIES=(
    "${BASE}/repeatmasker/Dsub/Dsub.GCF_014743375.2.fna.masked"
    "${BASE}/repeatmasker/Dbia/Dbia.GCF_025231255.1.fna.masked"
    "${BASE}/repeatmasker/Dmel/Dmel.GCF_000001215.4.fna.masked"
)

task=${SLURM_ARRAY_TASK_ID}

comparison=${COMPARISONS[$task]}
reference=${REFERENCES[$task]}
query=${QUERIES[$task]}

outdir=${BASE}/alignments/${comparison}
checkdir=${BASE}/checks/alignment_qc

prefix=${outdir}/${comparison}.promer
raw_delta=${prefix}.delta

many_delta=${outdir}/${comparison}.promer.many_to_many.L100.delta
one_delta=${outdir}/${comparison}.promer.one_to_one.L100.delta

many_coords=${outdir}/${comparison}.promer.many_to_many.L100.coords.tsv
one_coords=${outdir}/${comparison}.promer.one_to_one.L100.coords.tsv

summary=${checkdir}/${comparison}.promer.summary.tsv

mkdir -p \
    "${outdir}" \
    "${checkdir}"

for required in \
    "${PROMER}" \
    "${DELTA_FILTER}" \
    "${SHOW_COORDS}" \
    "${reference}" \
    "${query}"
do
    if [[ ! -s "${required}" ]]; then
        echo "ERROR: missing or empty input: ${required}" >&2
        exit 1
    fi
done

echo "===== MUMmer4 PROmer synteny ====="
echo "comparison=${comparison}"
echo "reference=${reference}"
echo "query=${query}"
echo "output_directory=${outdir}"
echo "MUMmer_package=4.0.1"
echo "promer_anchor_mode=mumreference"
echo "minimum_filtered_alignment_length_bp=100"
echo "many_to_many_filter=-m"
echo "one_to_one_filter=-1"
echo "cpus=${SLURM_CPUS_PER_TASK}"
date

echo "===== Software ====="

"${PROMER}" --version 2>&1 | head -5 || true
"${MUMMER_ENV}/bin/nucmer" --version 2>&1 | head -5 || true

echo "===== PROmer alignment ====="

rm -f \
    "${raw_delta}" \
    "${many_delta}" \
    "${one_delta}" \
    "${many_coords}" \
    "${one_coords}" \
    "${summary}"

"${PROMER}" \
    --mumreference \
    -p "${prefix}" \
    "${reference}" \
    "${query}"

if [[ ! -s "${raw_delta}" ]]; then
    echo "ERROR: PROmer produced no delta output" >&2
    exit 1
fi

echo "===== Many-to-many filtering ====="

"${DELTA_FILTER}" \
    -m \
    -l 100 \
    "${raw_delta}" \
    > "${many_delta}"

echo "===== One-to-one filtering ====="

"${DELTA_FILTER}" \
    -1 \
    -l 100 \
    "${raw_delta}" \
    > "${one_delta}"

for delta in \
    "${many_delta}" \
    "${one_delta}"
do
    if [[ ! -s "${delta}" ]]; then
        echo "ERROR: empty filtered delta file: ${delta}" >&2
        exit 1
    fi
done

echo "===== Coordinate tables ====="

"${SHOW_COORDS}" \
    -T \
    -H \
    -r \
    -c \
    -l \
    "${many_delta}" \
    > "${many_coords}"

"${SHOW_COORDS}" \
    -T \
    -H \
    -r \
    -c \
    -l \
    "${one_delta}" \
    > "${one_coords}"

for coords in \
    "${many_coords}" \
    "${one_coords}"
do
    if [[ ! -s "${coords}" ]]; then
        echo "ERROR: empty coordinate table: ${coords}" >&2
        exit 1
    fi
done

many_alignments=$(wc -l < "${many_coords}")
one_alignments=$(wc -l < "${one_coords}")

{
    printf "metric\tvalue\n"
    printf "comparison\t%s\n" "${comparison}"
    printf "reference\t%s\n" "${reference}"
    printf "query\t%s\n" "${query}"
    printf "mummer_package_version\t4.0.1\n"
    printf "promer_internal_version\t3.07\n"
    printf "anchor_mode\tmumreference\n"
    printf "minimum_alignment_length_bp\t100\n"
    printf "many_to_many_alignments\t%s\n" "${many_alignments}"
    printf "one_to_one_alignments\t%s\n" "${one_alignments}"
    printf "raw_delta\t%s\n" "${raw_delta}"
    printf "many_to_many_delta\t%s\n" "${many_delta}"
    printf "one_to_one_delta\t%s\n" "${one_delta}"
    printf "many_to_many_coords\t%s\n" "${many_coords}"
    printf "one_to_one_coords\t%s\n" "${one_coords}"
    printf "status\tPASS\n"
} > "${summary}"

echo "===== Output summary ====="
cat "${summary}"

echo "===== Output sizes ====="

ls -lh \
    "${raw_delta}" \
    "${many_delta}" \
    "${one_delta}" \
    "${many_coords}" \
    "${one_coords}" \
    "${summary}"

echo "FINAL: PASS"
date

#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --array=0-9%4
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=0-04:00:00
#SBATCH -J poly_regionBAM
#SBATCH -o /ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821/logs/poly_regionBAM.%A_%a.out
#SBATCH -e /ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821/logs/poly_regionBAM.%A_%a.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821

SAMPLETAB=${BASE}/config/samples.tsv
REGIONDIR=${BASE}/config/regions
INDIR=${BASE}/bam_rmdup
OUTDIR=${BASE}/bam_filtered
CHECKDIR=${BASE}/checks/bam_filtered

SAMTOOLS=/ru-auth/local/home/asato/scratch/conda_envs/augustus_env/bin/samtools

SAMPLES=(
    Dsub_L5_7
    JP-Kan
    JP-Sap
    JP-Sap2
    JP-Tok
)

CLASSES=(
    autosome
    X
)

sample_index=$((SLURM_ARRAY_TASK_ID / 2))
class_index=$((SLURM_ARRAY_TASK_ID % 2))

SAMPLE=${SAMPLES[$sample_index]}
CLASS=${CLASSES[$class_index]}

SPECIES=$(
    awk -F'\t' -v sample="${SAMPLE}" '
        NR > 1 && $1 == sample {
            print $2
            exit
        }
    ' "${SAMPLETAB}"
)

[[ -n "${SPECIES}" ]] || {
    echo "ERROR: species not found for ${SAMPLE}" >&2
    exit 1
}

INPUT_BAM=${INDIR}/${SAMPLE}.rmdup.bam
REGION_BED=${REGIONDIR}/${SPECIES}.${CLASS}.bed

OUTPUT_BAM=${OUTDIR}/${SAMPLE}.${CLASS}.filtered.bam
TMP_BAM=${OUTPUT_BAM}.tmp.${SLURM_JOB_ID}.${SLURM_ARRAY_TASK_ID}
FLAGSTAT=${CHECKDIR}/${SAMPLE}.${CLASS}.filtered.flagstat.txt

mkdir -p "${OUTDIR}" "${CHECKDIR}"

cleanup() {
    rm -f "${TMP_BAM}" "${TMP_BAM}.bai"
}
trap cleanup EXIT

for required in \
    "${INPUT_BAM}" \
    "${INPUT_BAM}.bai" \
    "${REGION_BED}"
do
    [[ -s "${required}" ]] || {
        echo "ERROR: missing input: ${required}" >&2
        exit 1
    }
done

echo "sample=${SAMPLE}"
echo "species=${SPECIES}"
echo "analysis_class=${CLASS}"
echo "input_bam=${INPUT_BAM}"
echo "region_bed=${REGION_BED}"
echo "output_bam=${OUTPUT_BAM}"
echo "required_flags=2"
echo "excluded_flags=3844"
echo "minimum_MAPQ=20"

"${SAMTOOLS}" view \
    -@ "${SLURM_CPUS_PER_TASK}" \
    -bh \
    -f 2 \
    -F 3844 \
    -q 20 \
    -L "${REGION_BED}" \
    -o "${TMP_BAM}" \
    "${INPUT_BAM}"

"${SAMTOOLS}" quickcheck -v "${TMP_BAM}"

mv "${TMP_BAM}" "${OUTPUT_BAM}"

"${SAMTOOLS}" index \
    -@ "${SLURM_CPUS_PER_TASK}" \
    "${OUTPUT_BAM}"

"${SAMTOOLS}" flagstat \
    -@ "${SLURM_CPUS_PER_TASK}" \
    "${OUTPUT_BAM}" \
    > "${FLAGSTAT}"

echo "output_bam=${OUTPUT_BAM}"
echo "flagstat=${FLAGSTAT}"
echo "FINAL: PASS"

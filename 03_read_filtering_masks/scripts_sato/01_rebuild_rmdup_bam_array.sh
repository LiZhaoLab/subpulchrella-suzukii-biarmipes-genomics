#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --array=0-4%2
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH -J poly_rmdup
#SBATCH -o /ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821/logs/poly_rmdup.%A_%a.out
#SBATCH -e /ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821/logs/poly_rmdup.%A_%a.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821
SAMPLE_TABLE=${BASE}/config/samples.tsv

SAMTOOLS=/ru-auth/local/home/asato/scratch/conda_envs/augustus_env/bin/samtools
THREADS=${SLURM_CPUS_PER_TASK}

TASK_LINE=$((SLURM_ARRAY_TASK_ID + 2))

IFS=$'\t' read -r \
    SAMPLE \
    SPECIES \
    INPUT_BAM \
    REFERENCE_FASTA \
    AUTOSOMAL_POOL_SIZE \
    X_POOL_SIZE \
    < <(sed -n "${TASK_LINE}p" "${SAMPLE_TABLE}")

if [[ -z "${SAMPLE:-}" ]]; then
    echo "ERROR: no sample found for array task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

if [[ ! -s "${INPUT_BAM}" ]]; then
    echo "ERROR: missing input BAM: ${INPUT_BAM}" >&2
    exit 1
fi

OUTDIR=${BASE}/bam_rmdup
TMPDIR=${BASE}/tmp_bam/${SAMPLE}
CHECKDIR=${BASE}/checks/rmdup/${SAMPLE}

NAME_BAM=${TMPDIR}/${SAMPLE}.name_sorted.bam
FIXMATE_BAM=${TMPDIR}/${SAMPLE}.fixmate.bam
POS_BAM=${TMPDIR}/${SAMPLE}.position_sorted.bam

OUT_BAM=${OUTDIR}/${SAMPLE}.rmdup.bam
OUT_BAI=${OUT_BAM}.bai

mkdir -p \
    "${OUTDIR}" \
    "${TMPDIR}" \
    "${CHECKDIR}"

echo "===== Sample ====="
echo "sample=${SAMPLE}"
echo "species=${SPECIES}"
echo "input_bam=${INPUT_BAM}"
echo "output_bam=${OUT_BAM}"
echo "samtools=${SAMTOOLS}"
echo "threads=${THREADS}"

"${SAMTOOLS}" --version |
head -3

# Validな完成済み出力があれば再計算しない
if [[ -s "${OUT_BAM}" && -s "${OUT_BAI}" ]] &&
   "${SAMTOOLS}" quickcheck -q "${OUT_BAM}"
then
    echo "Existing valid rmdup BAM found; skip rebuilding."
    echo "FINAL: PASS"
    exit 0
fi

# 部分的な出力だけを削除
rm -f \
    "${OUT_BAM}" \
    "${OUT_BAI}" \
    "${NAME_BAM}" \
    "${FIXMATE_BAM}" \
    "${POS_BAM}"

echo "===== Input BAM check ====="

"${SAMTOOLS}" quickcheck -v "${INPUT_BAM}"

"${SAMTOOLS}" view -H "${INPUT_BAM}" \
    > "${CHECKDIR}/${SAMPLE}.input.header.sam"

"${SAMTOOLS}" flagstat \
    -@ "${THREADS}" \
    "${INPUT_BAM}" \
    > "${CHECKDIR}/${SAMPLE}.input.flagstat.txt"

echo "===== Name sort ====="

"${SAMTOOLS}" sort \
    -n \
    -@ "${THREADS}" \
    -m 2G \
    -T "${TMPDIR}/${SAMPLE}.namesort.tmp" \
    -o "${NAME_BAM}" \
    "${INPUT_BAM}"

echo "===== Fix mate information ====="

"${SAMTOOLS}" fixmate \
    -m \
    -@ "${THREADS}" \
    "${NAME_BAM}" \
    "${FIXMATE_BAM}"

echo "===== Coordinate sort ====="

"${SAMTOOLS}" sort \
    -@ "${THREADS}" \
    -m 2G \
    -T "${TMPDIR}/${SAMPLE}.coordsort.tmp" \
    -o "${POS_BAM}" \
    "${FIXMATE_BAM}"

echo "===== Remove duplicates ====="

"${SAMTOOLS}" markdup \
    -r \
    -s \
    -@ "${THREADS}" \
    "${POS_BAM}" \
    "${OUT_BAM}" \
    2> "${CHECKDIR}/${SAMPLE}.markdup.stats.txt"

echo "===== Index and validate ====="

"${SAMTOOLS}" index \
    -@ "${THREADS}" \
    "${OUT_BAM}"

"${SAMTOOLS}" quickcheck -v "${OUT_BAM}"

"${SAMTOOLS}" view -H "${OUT_BAM}" \
    > "${CHECKDIR}/${SAMPLE}.rmdup.header.sam"

"${SAMTOOLS}" flagstat \
    -@ "${THREADS}" \
    "${OUT_BAM}" \
    > "${CHECKDIR}/${SAMPLE}.rmdup.flagstat.txt"

"${SAMTOOLS}" idxstats \
    "${OUT_BAM}" \
    > "${CHECKDIR}/${SAMPLE}.rmdup.idxstats.tsv"

# 完成後、一時BAMだけ削除
rm -f \
    "${NAME_BAM}" \
    "${FIXMATE_BAM}" \
    "${POS_BAM}"

rmdir "${TMPDIR}" 2>/dev/null || true

echo "output_bam=${OUT_BAM}"
echo "output_index=${OUT_BAI}"
echo "FINAL: PASS"
echo "Done ${SAMPLE}"

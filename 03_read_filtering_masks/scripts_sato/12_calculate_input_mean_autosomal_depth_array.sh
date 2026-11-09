#!/usr/bin/env bash
#SBATCH --job-name=inputMeanDepth
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=0-4%2
#SBATCH --output=logs/inputMeanDepth.%A_%a.out
#SBATCH --error=logs/inputMeanDepth.%A_%a.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821
SAMTOOLS=/ru-auth/local/home/asato/scratch/conda_envs/augustus_env/bin/samtools
SAMPLES=${BASE}/config/samples.tsv
REGION_DIR=${BASE}/config/regions
OUTDIR=${BASE}/checks/input_mean_depth

mkdir -p "${OUTDIR}"

mapfile -t sample_rows < <(tail -n +2 "${SAMPLES}")

row=${sample_rows[${SLURM_ARRAY_TASK_ID}]}

IFS=$'\t' read -r \
  sample \
  species \
  input_bam \
  reference_fasta \
  autosomal_pool_size \
  X_pool_size \
  <<< "${row}"

region_bed=${REGION_DIR}/${species}.autosome.bed
output=${OUTDIR}/${sample}.mean_autosomal_depth.tsv

[[ -s "${input_bam}" ]] || {
  echo "ERROR: BAM missing: ${input_bam}" >&2
  exit 1
}

[[ -s "${region_bed}" ]] || {
  echo "ERROR: BED missing: ${region_bed}" >&2
  exit 1
}

echo "sample=${sample}"
echo "species=${species}"
echo "input_bam=${input_bam}"
echo "region_bed=${region_bed}"

mean_depth=$(
  "${SAMTOOLS}" bedcov \
    "${region_bed}" \
    "${input_bam}" |
  awk '
    {
      interval_length = $3 - $2
      total_length += interval_length
      total_depth += $NF
    }
    END {
      if (total_length > 0) {
        printf "%.3f", total_depth / total_length
      } else {
        print "NA"
      }
    }
  '
)

printf "sample\tspecies\tmean_autosomal_depth\n" > "${output}"
printf "%s\t%s\t%s\n" \
  "${sample}" \
  "${species}" \
  "${mean_depth}" \
  >> "${output}"

cat "${output}"

echo "FINAL: PASS"

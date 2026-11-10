#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --array=0-3
#SBATCH -J vcf_Dsuz4_Iso
#SBATCH -o logs/vcf_Dsuz4_Iso.%A_%a.out
#SBATCH -e logs/vcf_Dsuz4_Iso.%A_%a.err

set -euo pipefail

BCFTOOLS=/programs/x86_64-linux/bcftools/1.20/bin/bcftools
SAMTOOLS=/programs/x86_64-linux/samtools/1.22/bin/samtools

BASE=/ru-auth/local/home/asato/scratch/MKtest/Dsuz4pop_MKtest
REF=${BASE}/ref/Dsuz_IsoJpt1.0.fa
BAMDIR=${BASE}/bam_iso
VCFDIR=${BASE}/vcf_IsoJpt1

mkdir -p "${VCFDIR}" logs

POPS=(JP-Kan JP-Sap JP-Sap2 JP-Tok)
POP=${POPS[$SLURM_ARRAY_TASK_ID]}

BAM=${BAMDIR}/${POP}.sorted.bam

RAW=${VCFDIR}/${POP}.IsoJpt1.raw.vcf
BIALLELIC=${VCFDIR}/${POP}.IsoJpt1.biallelic.vcf
ALT=${VCFDIR}/${POP}.IsoJpt1.altfiltered.vcf
ALT_GZ=${VCFDIR}/${POP}.IsoJpt1.altfiltered.vcf.gz

if [ ! -f "${REF}.fai" ]; then
  ${SAMTOOLS} faidx "${REF}"
fi

echo "Processing ${POP}"
echo "BAM: ${BAM}"

${BCFTOOLS} mpileup -Ou -f "${REF}" "${BAM}" | \
  ${BCFTOOLS} call -mv -Ov -o "${RAW}"

${BCFTOOLS} view -v snps -m2 -M2 "${RAW}" -Ov -o "${BIALLELIC}"

${BCFTOOLS} filter \
  -i 'DP>=20 && (INFO/DP4[2] + INFO/DP4[3]) >= 3 && INFO/DP4[2] > 0 && INFO/DP4[3] > 0 && (INFO/DP4[2] + INFO/DP4[3]) / (INFO/DP4[0] + INFO/DP4[1] + INFO/DP4[2] + INFO/DP4[3]) > 0.05' \
  "${BIALLELIC}" \
  -Ov \
  -o "${ALT}"

${BCFTOOLS} view "${ALT}" -Oz -o "${ALT_GZ}"
${BCFTOOLS} index -t "${ALT_GZ}"

echo "Done ${POP}"

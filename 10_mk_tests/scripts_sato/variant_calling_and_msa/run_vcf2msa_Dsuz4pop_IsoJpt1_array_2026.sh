#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00
#SBATCH --array=0-3
#SBATCH -J vcf2msa26
#SBATCH -o /ru-auth/local/home/asato/scratch/MKtest/Dsuz4pop_MKtest/scripts/logs/vcf2msa26.%A_%a.out
#SBATCH -e /ru-auth/local/home/asato/scratch/MKtest/Dsuz4pop_MKtest/scripts/logs/vcf2msa26.%A_%a.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/MKtest/Dsuz4pop_MKtest

POPS=("JP-Kan" "JP-Sap" "JP-Sap2" "JP-Tok")
POP="${POPS[$SLURM_ARRAY_TASK_ID]}"

VCF2MSA=/ru-auth/local/home/jpeng/bin/python/vcf2msa_2026.py
VCFDIR=${BASE}/vcf_IsoJpt1
CDS=${BASE}/cds/Dsuz_IsoJpt1.0.cds.coords.fa
VCF=${VCFDIR}/${POP}.IsoJpt1.altfiltered.vcf.gz

OUTDIR=${BASE}/vcf2msa_2026/${POP}_msa
mkdir -p "${OUTDIR}"

echo "POP: ${POP}"
echo "VCF2MSA: ${VCF2MSA}"
echo "VCF: ${VCF}"
echo "CDS: ${CDS}"
echo "OUTDIR: ${OUTDIR}"

python "${VCF2MSA}" \
  -vcf "${VCF}" \
  -fasta "${CDS}" \
  -out "${OUTDIR}" \
  -species "${POP}"

echo -n "msa fasta files: "
find "${OUTDIR}" -maxdepth 1 -type f -name "msa_*.fasta" | wc -l

echo "Done: ${POP}"

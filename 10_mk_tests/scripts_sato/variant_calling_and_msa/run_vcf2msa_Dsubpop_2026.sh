#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0-04:00:00
#SBATCH -J vcf2msa_Dsub
#SBATCH -o /ru-auth/local/home/asato/scratch/MKtest/Dsubpop/scripts/logs/vcf2msa_Dsub.%j.out
#SBATCH -e /ru-auth/local/home/asato/scratch/MKtest/Dsubpop/scripts/logs/vcf2msa_Dsub.%j.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/MKtest/Dsubpop

VCF2MSA=/ru-auth/local/home/jpeng/bin/python/vcf2msa_2026.py
VCF=${BASE}/vcf/Dsub.altfiltered.vcf
CDS=${BASE}/cds/Dsub.cds.coords.fa
OUTDIR=${BASE}/vcf2msa_2026/Dsub_msa

mkdir -p "${OUTDIR}"

echo "Focal species: Dsub"
echo "VCF2MSA: ${VCF2MSA}"
echo "VCF: ${VCF}"
echo "CDS: ${CDS}"
echo "OUTDIR: ${OUTDIR}"

python "${VCF2MSA}" \
  -vcf "${VCF}" \
  -fasta "${CDS}" \
  -out "${OUTDIR}" \
  -species Dsub

echo -n "msa fasta files: "
find "${OUTDIR}" -maxdepth 1 -type f -name "msa_*.fasta" | wc -l

echo -n "cds_snp_index rows: "
tail -n +2 "${OUTDIR}/cds_snp_index.csv" | wc -l

echo "SNP counts:"
tail -n +2 "${OUTDIR}/cds_snp_index.csv" \
  | awk -F',' '{count[$3]++} END{for (s in count) print "SNP="s, count[s]}'

echo "Done."

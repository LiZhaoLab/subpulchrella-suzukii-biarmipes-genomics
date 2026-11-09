#!/bin/bash
#SBATCH -p zhao_a
#SBATCH -A zhao_condo_bank
#SBATCH -J suzQC_chr3mm2
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lzhao@rockefeller.edu
#SBATCH -o /ru-auth/local/home/lzhao/vscode_projects/suzukii/logs/121_chr3_mm2.%j.out
#SBATCH -e /ru-auth/local/home/lzhao/vscode_projects/suzukii/logs/121_chr3_mm2.%j.err
# Script 121: is the gene-poor, TE-rich centre of D. suzukii IsoJpt1.0 chromosome 3 (NC_092082.1,
# 27-62 Mb) present in the lab's earlier D. suzukii assembly (Dsuz_RU_1.0) and in D. subpulchrella (RU_Dsub_v1.1)?
# Unmasked nucleotide alignment (PROmer ran on hard-masked genomes, so a TE-rich region cannot align there).
# minimap2 2.30 (Airi's popgen env), preset asm20 for both targets, primary alignments only. Coverage: script 122.
set -euo pipefail
G=/ru-auth/local/home/asato/scratch/genome
MM=/ru-auth/local/home/asato/scratch/conda_envs/popgen/bin/minimap2
OUT=/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/synteny_chr3
awk '/^>/{p=($1==">NC_092082.1")} p' $G/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.fna > $OUT/Dsuz_chr3.fa
grep -c '>' $OUT/Dsuz_chr3.fa
$MM --version
for tgt in RU1:$G/Dsuz_RU_1.0/GCF_037355615.1_Dsuz_RU_1.0_genomic.fna Dsub:$G/Dsub_v1.1/GCF_014743375.2_RU_Dsub_v1.1_genomic.fna; do
  name=${tgt%%:*}; fa=${tgt#*:}
  $MM -x asm20 -t 8 --secondary=no $fa $OUT/Dsuz_chr3.fa > $OUT/chr3_vs_${name}.paf
  wc -l $OUT/chr3_vs_${name}.paf
done
echo FINAL: PASS

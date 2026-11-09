#!/bin/bash
# Purity of the D. suzukii pools from mitochondrial reads (analysis 4.2, third method).
# 1. D. subpulchrella mtDNA: the D. subpulchrella assembly has no annotated mitochondrial genome.
#    Its BAM has no unmapped mtDNA reads (0 of 5,452,866 unmapped reads align to the
#    D. suzukii mtDNA), so mtDNA reads must be mapped to mtDNA-like sequence in the assembly.
#    Instead, the raw read subsample of 63_competitive_mapping.sh (2 million pairs of lane L5, not
#    filtered against any reference) is mapped (single-end, bwa mem) to the D. suzukii mitochondrial
#    genome NC_060762.1 (part of the IsoJpt1.0 assembly).
# 2. Per-position base counts (samtools mpileup, MAPQ >= 20, base quality >= 20, no BAQ, no depth cap)
#    for the D. subpulchrella reads and for the mitochondrial reads in each D. suzukii pool BAM.
# 3. Diagnostic sites and frequencies: scripts/79_summarize_purity.py.
# Tools: project env envs/te (bwa 0.7.18, samtools 1.21). Output: $QC/mito/
set -euo pipefail
unset LD_LIBRARY_PATH
P=/ru-auth/local/home/lzhao/vscode_projects/suzukii
export PATH=$P/envs/te/bin:/usr/bin:/bin
QC=/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/mito
BAM=/ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821/bam_rmdup
FNA=/ru-auth/local/home/asato/scratch/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.fna
MITO=NC_060762.1
mkdir -p $QC && cd $QC
samtools faidx $FNA $MITO > dsuz_mito.fa
samtools faidx dsuz_mito.fa
bwa index dsuz_mito.fa 2> /dev/null

# D. subpulchrella pool: raw lane-L5 read subsample -> single-end mapping to D. suzukii mtDNA
CM=/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/compmap
cat $CM/L5_1.fq $CM/L5_2.fq | \
  bwa mem -t 2 dsuz_mito.fa - 2> dsub_mito.bwa.log | samtools view -b -F 0x904 -q 20 - | \
  samtools sort -o Dsub_L5_7.mito.bam -
samtools index Dsub_L5_7.mito.bam
samtools mpileup -B -Q 20 -q 20 -d 0 -f dsuz_mito.fa Dsub_L5_7.mito.bam 2> /dev/null | cut -f 1-5 > Dsub_L5_7.mito.pileup

[ -n "${SKIP_DSUZ:-}" ] && exit 0
for s in JP-Kan JP-Sap JP-Sap2 JP-Tok; do
  samtools mpileup -B -Q 20 -q 20 -d 0 --ff UNMAP,SECONDARY,QCFAIL,DUP -r $MITO -f dsuz_mito.fa \
    $BAM/$s.rmdup.bam 2> /dev/null | cut -f 1-5 > $s.mito.pileup
done
wc -l *.pileup

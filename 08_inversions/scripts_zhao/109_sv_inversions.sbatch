#!/bin/bash
# Inversion breakpoints in the D. subpulchrella pool, with the reference strain as control (
# inversion analysis, approach 2).
# 1. Reference-strain Illumina reads (SRR12476582; strain Dsub33 F10#4, inbred, mixed adult females;
#    the strain of the RU_Dsub_v1.1 assembly): 35% random subsample of pairs (seqtk, seed 20261030),
#    bwa mem 0.7.18 to the same reference as the pool BAM (te/ref/Dsub.fa; sequence names and
#    lengths identical to GCF_014743375.2), fixmate, sort, duplicates marked.
# 2. Delly 1.7.3 inversion discovery (-t INV, PE mapping quality >= 20) in the pool BAM (Airi's
#    duplicate-removed BAM) and in the reference-strain BAM; sites merged (delly merge) and
#    genotyped in both BAMs, so every candidate has read counts in the pool and in the reference strain.
# A breakpoint supported in the reference strain is an assembly or mapping artefact, not a polymorphism.
#SBATCH -J suzQC_svinv
#SBATCH -p zhao_a
#SBATCH -A zhao_condo_bank
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o /ru-auth/local/home/lzhao/vscode_projects/suzukii/logs/svinv_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lzhao@rockefeller.edu
set -euo pipefail
unset LD_LIBRARY_PATH
P=/ru-auth/local/home/lzhao/vscode_projects/suzukii
TEB=$P/envs/te/bin
SV=$P/envs/sv/bin
REF=/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/te/ref/Dsub.fa
POOL=/ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821/bam_rmdup/Dsub_L5_7.rmdup.bam
RAW=/ru-auth/local/home/lzhao/Data_scratch/zhao/Genome/SRR12476582_D_sub_genome_Kmer
O=/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/sv
mkdir -p $O/tmp
cd $O
echo "$(date) start on $(hostname); delly $($SV/delly --version 2>&1 | head -1)"
# pool discovery in the background
( $SV/delly call -t INV -q 20 -h 4 -g $REF -o $O/pool.inv.bcf $POOL > $O/delly_pool.log 2>&1
  echo "$(date) pool discovery done" ) &
# reference-strain reads: subsample, map, mark duplicates
$TEB/bwa mem -t 12 -R "@RG\tID:refstrain\tSM:refstrain" $REF \
    <($TEB/seqtk sample -s 20261030 $RAW/SRR12476582_1.fastq.gz 0.35) \
    <($TEB/seqtk sample -s 20261030 $RAW/SRR12476582_2.fastq.gz 0.35) 2> $O/bwa_refstrain.log | \
  $TEB/samtools fixmate -m -u - - | \
  $TEB/samtools sort -@ 4 -m 2G -T $O/tmp/s -u - | \
  $TEB/samtools markdup -@ 4 -T $O/tmp/m - $O/refstrain.bam
$TEB/samtools index $O/refstrain.bam
$TEB/samtools flagstat $O/refstrain.bam > $O/refstrain.flagstat
echo "$(date) reference-strain BAM done"
$SV/delly call -t INV -q 20 -h 4 -g $REF -o $O/refstrain.inv.bcf $O/refstrain.bam > $O/delly_refstrain.log 2>&1
wait
echo "$(date) discovery done"
# merge filters off: delly merge defaults (max size 1 Mb, min ALT fraction 0.15, min quality 200) remove
# large and low-frequency inversions, which matter in pool data. The first run used the
# defaults (608 sites; kept as *_defaultmerge.*); merge and genotyping were repeated with these options.
$SV/delly merge -y 0 -a 0 -v 1 -m 0 -n 500000000 -b 1000 -r 0.8 -o $O/sites.inv.bcf $O/pool.inv.bcf $O/refstrain.inv.bcf > $O/delly_merge.log 2>&1
$SV/delly call -t INV -q 20 -h 4 -g $REF -v $O/sites.inv.bcf -o $O/geno.inv.bcf $POOL $O/refstrain.bam > $O/delly_geno.log 2>&1
$SV/bcftools index -f $O/geno.inv.bcf
$SV/bcftools query -H -f '%CHROM\t%POS\t%INFO/END\t%ID\t%FILTER\t%INFO/CT\t%INFO/PE\t%INFO/SR\t%INFO/PRECISE\t%QUAL[\t%DR\t%DV\t%RR\t%RV]\n' \
    $O/geno.inv.bcf > $O/geno.inv.tsv
rm -rf $O/tmp
echo "$(date) done"; wc -l $O/geno.inv.tsv

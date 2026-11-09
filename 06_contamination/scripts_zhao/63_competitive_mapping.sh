#!/bin/bash
#SBATCH --job-name=suzQC_compmap
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=4:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lzhao@rockefeller.edu
#SBATCH --output=/ru-auth/local/home/lzhao/vscode_projects/suzukii/logs/compmap.%j.out
#SBATCH --error=/ru-auth/local/home/lzhao/vscode_projects/suzukii/logs/compmap.%j.err
# Competitive mapping of D. subpulchrella pool reads to a combined four-species reference
# (D. subpulchrella, D. suzukii, D. biarmipes, D. melanogaster), with simulated-read controls
# for the misassignment rate. bwa 0.7.19 (conda env; the linuxbrew bwa needs LD_LIBRARY_PATH), samtools 1.9. Seeds 20261021 (Dsub), 20261022 (Dsuz).
set -euo pipefail
unset LD_LIBRARY_PATH
A=/ru-auth/local/home/asato/scratch
W=/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/compmap
R=/ru-auth/local/home/lzhao/Data_scratch/zhao/Dsub_pool_genome
P=/ru-auth/local/home/lzhao/vscode_projects/suzukii
PY=/ru-auth/local/home/lzhao/anaconda3/envs/jupyter-ml/bin/python
BWA=/ru-auth/local/home/asato/scratch/conda_envs/popgen/bin/bwa
NPAIRS=2000000
mkdir -p $W && cd $W

declare -A G=(
  [Dsub]=$A/genome/Dsub_v1.1/GCF_014743375.2_RU_Dsub_v1.1_genomic.fna
  [Dsuz]=$A/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.fna
  [Dbia]=$A/genome/Dbia/GCF_025231255.1_RU_DBia_V1.1_genomic.fna
  [Dmel]=$A/synteny_4species/reanalysis_20260824/input_genomes/Dmel.GCF_000001215.4.fna
)
if [[ ! -s combined.fa.bwt ]]; then
  : > combined.fa
  for sp in Dsub Dsuz Dbia Dmel; do
    awk -v sp=$sp '/^>/{print ">" sp "|" substr($1,2); next} {print}' ${G[$sp]} >> combined.fa
  done
  samtools faidx combined.fa
  $BWA index -a bwtsw combined.fa
fi

# inputs: first NPAIRS read pairs of each lane, and simulated controls
L5=$R/Dsub_pool_CKDN250010106-1A_22WVMFLT4_L5
L7=$R/Dsub_pool_CKDN250010106-1A_22WVNJLT4_L7
# zcat receives SIGPIPE when head stops reading (exit 141); tolerate it, then check line counts.
for lane in L5 L7; do
  base=${!lane}
  for m in 1 2; do
    out=${lane}_${m}.fq
    if [[ ! -s $out || $(wc -l < $out) -ne $((4*NPAIRS)) ]]; then
      { zcat ${base}_${m}.fq.gz.trim.gz || true; } | head -n $((4*NPAIRS)) > $out.tmp
      [[ $(wc -l < $out.tmp) -eq $((4*NPAIRS)) ]] || { echo "short read file $out.tmp" >&2; exit 1; }
      mv $out.tmp $out
    fi
  done
done
[[ -s simDsub_1.fq.gz ]] || $PY $P/scripts/62_simulate_reads.py ${G[Dsub]} 500000 20261021 simDsub
[[ -s simDsuz_1.fq.gz ]] || $PY $P/scripts/62_simulate_reads.py ${G[Dsuz]} 500000 20261022 simDsuz

# map and tally: primary alignments (exclude unmapped 4, secondary 256, supplementary 2048)
for s in L5 L7 simDsub simDsuz; do
  r1=${s}_1.fq; r2=${s}_2.fq
  [[ -s $r1 ]] || { r1=${s}_1.fq.gz; r2=${s}_2.fq.gz; }
  $BWA mem -t 8 combined.fa $r1 $r2 2> $s.bwa.log | samtools view -F 0x900 - | \
    awk -v s=$s 'BEGIN{OFS="\t"}
      { n++; if (and($2,4)) {u++; next}
        split($3,a,"|"); sp=a[1];
        if ($5>=20) q20[sp]++; else low[sp]++ }
      END{ print s, "reads", n; print s, "unmapped", u+0;
           for (k in q20) print s, "MAPQ>=20", k, q20[k];
           for (k in low) print s, "MAPQ<20", k, low[k] }' > $s.tally.tsv
done
cat L5.tally.tsv L7.tally.tsv simDsub.tally.tsv simDsuz.tally.tsv > $P/results/compmap_tally.tsv
echo done

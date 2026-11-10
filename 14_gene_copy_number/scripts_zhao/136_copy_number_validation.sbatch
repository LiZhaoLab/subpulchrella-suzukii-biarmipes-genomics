#!/bin/bash
#SBATCH -J c2_validate
#SBATCH -p zhao_a
#SBATCH -A zhao_condo_bank
#SBATCH -c 28
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH -o /ru-auth/local/home/lzhao/vscode_projects/suzukii/logs/%x_%j.out
#SBATCH -e /ru-auth/local/home/lzhao/vscode_projects/suzukii/logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lzhao@rockefeller.edu
# C2 validation (scripts/136_copy_number_validation.py).
# Part 1: samtools bedcov over the CDS of every annotated gene, 6 BAMs x MAPQ >= 0 and >= 20
#         (12 single-threaded processes in the background).
# Part 2: miniprot (16 threads) of D. melanogaster, D. biarmipes, D. subpulchrella and D. suzukii
#         IsoJpt1.0 longest proteins to the four assemblies.
set -euo pipefail
S=/ru-auth/local/home/lzhao/vscode_projects/suzukii/envs/sv/bin/samtools
MP=/lustre/fs8/zhao_lab/scratch/lzhao/envs/miniprot/bin/miniprot
O=/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/gene_family_cn
B=/ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821/bam_rmdup
G=/ru-auth/local/home/asato/scratch/genome
mkdir -p $O/depth $O/mp
$S --version | head -1
echo "miniprot $($MP --version)"
date

pids=()
run_bedcov() {  # sample species bam mapq
  $S bedcov -Q $4 -j -c --max-depth 1000000 $O/$2_cds.bed $3 > $O/depth/$1.Q$4.bedcov.tmp
  mv $O/depth/$1.Q$4.bedcov.tmp $O/depth/$1.Q$4.bedcov
}
for q in 0 20; do
  run_bedcov Dsub_L5_7 Dsub $B/Dsub_L5_7.rmdup.bam $q & pids+=($!)
  run_bedcov refstrain Dsub /ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/sv/refstrain.bam $q & pids+=($!)
  for p in JP-Kan JP-Sap JP-Sap2 JP-Tok; do
    run_bedcov $p DsuzI $B/$p.rmdup.bam $q & pids+=($!)
  done
done

declare -A FNA=(
  [Dbia]=$G/Dbia/GCF_025231255.1_RU_DBia_V1.1_genomic.fna
  [Dsub]=$G/Dsub_v1.1/GCF_014743375.2_RU_Dsub_v1.1_genomic.fna
  [DsuzI]=$G/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.fna
  [DsuzR]=$G/Dsuz_RU_1.0/GCF_037355615.1_Dsuz_RU_1.0_genomic.fna
)
for g in Dbia Dsub DsuzI DsuzR; do
  $MP -t 16 -d $O/mp/$g.mpi ${FNA[$g]} 2> $O/mp/$g.index.log
  $MP -t 16 --gff -P $g -p 0.4 -N 50 --outn=50 --outs=0.4 --outc=0.3 $O/mp/$g.mpi $O/mp_queries.faa \
      > $O/mp/$g.gff.tmp 2> $O/mp/$g.map.log
  mv $O/mp/$g.gff.tmp $O/mp/$g.gff
  echo "miniprot $g done $(date)"
done

fail=0
for pid in "${pids[@]}"; do wait $pid || fail=1; done
ls -la $O/depth
[ $fail -eq 0 ] || { echo "a bedcov process failed"; exit 1; }
echo "all done $(date)"

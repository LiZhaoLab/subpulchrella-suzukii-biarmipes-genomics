#!/bin/bash
#SBATCH -J c2_locidepth
#SBATCH -p zhao_a
#SBATCH -A zhao_condo_bank
#SBATCH -c 6
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH -o /ru-auth/local/home/lzhao/vscode_projects/suzukii/logs/%x_%j.out
#SBATCH -e /ru-auth/local/home/lzhao/vscode_projects/suzukii/logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lzhao@rockefeller.edu
# Read depth over the exons of every miniprot locus (BEDs written by scripts/137_copy_number_analysis.py), MAPQ >= 0.
set -euo pipefail
S=/ru-auth/local/home/lzhao/vscode_projects/suzukii/envs/sv/bin/samtools
O=/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/gene_family_cn
B=/ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821/bam_rmdup
$S --version | head -1
run() {  # sample genome bam
  $S bedcov -Q 0 -j -c --max-depth 1000000 $O/mp/$2_loci_exons.bed $3 > $O/depth/loci_$1.Q0.bedcov.tmp
  mv $O/depth/loci_$1.Q0.bedcov.tmp $O/depth/loci_$1.Q0.bedcov
}
pids=()
run Dsub_L5_7 Dsub $B/Dsub_L5_7.rmdup.bam & pids+=($!)
run refstrain Dsub /ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/sv/refstrain.bam & pids+=($!)
for p in JP-Kan JP-Sap JP-Sap2 JP-Tok; do run $p DsuzI $B/$p.rmdup.bam & pids+=($!); done
fail=0; for pid in "${pids[@]}"; do wait $pid || fail=1; done
ls -la $O/depth/loci_*
[ $fail -eq 0 ] || { echo "a bedcov process failed"; exit 1; }
echo "all done $(date)"

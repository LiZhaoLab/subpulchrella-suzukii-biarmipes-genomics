#!/bin/bash
#SBATCH -J c2_diamond
#SBATCH -p zhao_a
#SBATCH -A zhao_condo_bank
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH -o /ru-auth/local/home/lzhao/vscode_projects/suzukii/logs/%x_%j.out
#SBATCH -e /ru-auth/local/home/lzhao/vscode_projects/suzukii/logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lzhao@rockefeller.edu
# C2 first pass (scripts/135_gene_family_copy_number.py): DIAMOND 2.1.15 blastp of longest proteins
# against D. melanogaster (family assignment) and against the same genome (near-identical paralogs).
set -euo pipefail
D=/ru-auth/local/home/asato/scratch/conda_envs/mcscan_env/bin/diamond
O=/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/gene_family_cn
F="6 qseqid sseqid pident length qlen slen evalue bitscore"
$D version
$D makedb --in $O/Dmel_longest.faa -d $O/Dmel_longest -p 16
for sp in Dmel Dbia Dsub DsuzI DsuzR; do
  $D blastp -q $O/${sp}_longest.faa -d $O/Dmel_longest --more-sensitive -e 1e-5 -k 25 -p 16 \
     --outfmt $F -o $O/${sp}_vs_Dmel.tsv --quiet
  $D makedb --in $O/${sp}_longest.faa -d $O/${sp}_self -p 16 --quiet
  $D blastp -q $O/${sp}_longest.faa -d $O/${sp}_self -e 1e-10 -k 50 --id 90 -p 16 \
     --outfmt $F -o $O/${sp}_self.tsv --quiet
  echo "$sp done $(date)"
done

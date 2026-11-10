#!/bin/bash
# Analysis F: SweepFinder2 composite likelihood ratio (CLR) on one 1-Mb core of one
# chromosome. SNPs from the core plus 100-kb flanks; grid every 2 kb over the core only (-lu);
# background folded spectrum of the pool and chromosome class (scripts/98_sweep_input.sh).
# Usage: 99_sweep_chunk.sh "POOL CHROM START END CLS"
set -euo pipefail
unset LD_LIBRARY_PATH
set -- $*
POOL=$1; CHROM=$2; S=$3; E=$4; CLS=$5
O=/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/sweeps/$POOL
SF2=/ru-auth/local/home/lzhao/vscode_projects/suzukii/envs/popgen/bin/SweepFinder2
mkdir -p $O/clr $O/tmp
T=$O/tmp/${CHROM}_$S
awk -v a=$((S - 100000)) -v b=$((E + 100000)) 'NR == 1 || ($1 >= a && $1 <= b)' $O/$CHROM.sf2.txt > $T.in
seq $S 2000 $E > $T.grid
if [ $(wc -l < $T.in) -lt 50 ]; then echo "too few SNPs"; rm -f $T.in $T.grid; exit 0; fi
$SF2 -lu $T.grid $T.in $O/$CLS.spect.txt $T.out > /dev/null
mv $T.out $O/clr/${CHROM}_$S.txt
rm -f $T.in $T.grid

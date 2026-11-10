#!/bin/bash
# Spanning-read coverage (MAPQ >= 20, primary; read must extend >= 20 bp past the point on each side) at
# (a) six points around each retrocopy (500, 1000, 1500 bp outside each end) and
# (b) points every 1 Mb on IsoJpt1.0 chr 2L, 2R, 3 (autosomes), X and Y, in each D. suzukii pool. No random numbers.
set -euo pipefail
W=/lustre/fs8/home/lzhao/vscode_projects/suzukii/results/retrocopies_six_genes
S=/lustre/fs8/home/lzhao/vscode_projects/suzukii/envs/sv/bin/samtools
FAI=/ru-auth/local/home/asato/scratch/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.fna.fai
B=/ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821/bam_rmdup
: > $W/points.tsv
while IFS=$'\t' read n c s e; do
  for d in 500 1000 1500; do printf "%s\t%s\t%d\n" "retro_$n" $c $((s-d)) >> $W/points.tsv; printf "%s\t%s\t%d\n" "retro_$n" $c $((e+d)) >> $W/points.tsv; done
done < $W/retro.tsv
awk -F'\t' 'BEGIN{cls["NC_092080.1"]="A";cls["NC_092081.1"]="A";cls["NC_092082.1"]="A";cls["NC_092084.1"]="X";cls["NC_092085.1"]="Y"}
  ($1 in cls){for(p=500000;p<$2-500000;p+=1000000) printf "%s\t%s\t%d\n", cls[$1], $1, p}' $FAI >> $W/points.tsv
for p in JP-Kan JP-Sap JP-Sap2 JP-Tok; do
  while IFS=$'\t' read cl c x; do
    k=$($S view -q 20 -F 0x904 $B/$p.rmdup.bam "$c:$x-$x" | awk -v x=$x '
      { a = $4; len = 0; y = $6
        while (match(y, /^[0-9]+[MIDNSHP=X]/)) { m = substr(y, RSTART, RLENGTH); if (m ~ /[MDN=X]$/) len += m + 0; y = substr(y, RLENGTH + 1) }
        if (a <= x - 20 && a + len - 1 >= x + 20) k++ } END { print k + 0 }')
    printf "%s\t%s\t%s\t%d\t%d\n" $p "$cl" $c $x $k
  done < $W/points.tsv
done > $W/point_cov.tsv
wc -l $W/point_cov.tsv

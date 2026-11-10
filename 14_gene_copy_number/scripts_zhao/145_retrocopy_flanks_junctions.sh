#!/bin/bash
# Retrocopy checks for RpL36, RpS14, CCT1 (T-cp1) and VhaAC39-1 in D. suzukii.
# 1) Flank mapping: retrocopy + 2 kb flanks from IsoJpt1.0 -> RU_1.0 and Dsub_v1.1 (minimap2 2.30, asm20).
# 2) Junction-spanning reads (MAPQ >= 20, primary) in the four D. suzukii pools on IsoJpt1.0, relative to reads
#    spanning a flank point 1 kb outside each junction.
set -euo pipefail
W=/lustre/fs8/home/lzhao/vscode_projects/suzukii/results/retrocopies_six_genes
mkdir -p $W
S=/lustre/fs8/home/lzhao/vscode_projects/suzukii/envs/sv/bin/samtools
MM=/ru-auth/local/home/asato/scratch/conda_envs/popgen/bin/minimap2
GI=/ru-auth/local/home/asato/scratch/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.fna
GR=/ru-auth/local/home/asato/scratch/genome/Dsuz_RU_1.0/GCF_037355615.1_Dsuz_RU_1.0_genomic.fna
GS=/ru-auth/local/home/asato/scratch/genome/Dsub_v1.1/Dsub_v1.1_clean.fa
B=/ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821/bam_rmdup
FL=2000
# name seq start end (1-based, inclusive; miniprot alignment span)
cat > $W/retro.tsv <<EOF
RpL36	NW_027255939.1	1918069	1918413
RpS14	NW_027255927.1	2662680	2663100
CCT1	NW_027255899.1	1547553	1549223
VhaAC39-1	NW_027255896.1	24702032	24703088
EOF
: > $W/query.fa
while IFS=$'\t' read n c s e; do
  $S faidx $GI "$c:$((s-FL))-$((s-1))" | sed "1s/.*/>${n}_Lflank/" >> $W/query.fa
  $S faidx $GI "$c:$s-$e" | sed "1s/.*/>${n}_insert/" >> $W/query.fa
  $S faidx $GI "$c:$((e+1))-$((e+FL))" | sed "1s/.*/>${n}_Rflank/" >> $W/query.fa
  $S faidx $GI "$c:$((s-FL))-$((e+FL))" | sed "1s/.*/>${n}_whole/" >> $W/query.fa
done < $W/retro.tsv
for t in RU:$GR Dsub:$GS; do
  $MM -t 1 -x asm20 -c --secondary=yes -N 5 ${t#*:} $W/query.fa 2> $W/mm_${t%%:*}.log > $W/map_${t%%:*}.paf
done
echo "== junction-spanning reads (MAPQ>=20, primary): span_L span_R flankL flankR"
while IFS=$'\t' read n c s e; do
  for p in JP-Kan JP-Sap JP-Sap2 JP-Tok; do
    $S view -q 20 -F 0x904 $B/$p.rmdup.bam "$c:$((s-1500))-$((e+1500))" | \
    awk -v n=$n -v p=$p -v s=$s -v e=$e '
      function rend(pos, cig,   x, m, len) { len = 0; x = cig
        while (match(x, /^[0-9]+[MIDNSHP=X]/)) { m = substr(x, RSTART, RLENGTH); if (m ~ /[MDN=X]$/) len += m + 0; x = substr(x, RLENGTH + 1) }
        return pos + len - 1 }
      { a = $4; b = rend($4, $6)
        if (a <= s - 20 && b >= s + 20) jl++
        if (a <= e - 20 && b >= e + 20) jr++
        if (a <= s - 1000 - 20 && b >= s - 1000 + 20) fl++
        if (a <= e + 1000 - 20 && b >= e + 1000 + 20) fr++ }
      END { printf "%s\t%s\t%d\t%d\t%d\t%d\n", n, p, jl, jr, fl, fr }'
  done
done < $W/retro.tsv

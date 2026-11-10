#!/bin/bash
# Script 146f: sequence inputs of the retrocopy checks (scripts/146b, 146e) and the BLAST searches for
# copies of the parental coding sequences (CDS) of the seven D. suzukii retrogenes of Chiu et al. 2013 (parents
# RpS14a, Sce, T-cp1 = CCT1, betaTub97EF, RpL36, VhaAC39-1) in IsoJpt1.0 and in two assemblies of the WT3 strain:
# GCA_000472105.1 (Dsuzukii.v01, Chiu et al. 2013, short reads) and GCA_013340165.1 (LBDM_Dsuz_2.1.pri, Paris et al.
# 2020, long reads). These commands were first run interactively; this script collects them.
# Outputs (results/retrocopies_six_genes/): cds6.fa, cds6_vs_IsoJpt.tsv, parent_cds.fa,
#   rps14_cds_vs_IsoJpt.tsv, parent_cds_vs_<assembly>.tsv, retro_regions.fa, retro_regions_3asm.fa.
# BLAST+ 2.11.0 (blastn -task dc-megablast, E <= 1e-10); samtools 1.24 (envs/sv). No random numbers.
set -euo pipefail
P=/ru-auth/local/home/lzhao/vscode_projects/suzukii
R=$P/results/retrocopies_six_genes
D=/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC/wt3_assemblies      # downloaded assemblies (untrusted input)
S=$P/envs/sv/bin/samtools
G=/ru-auth/local/home/asato/scratch/genome
GI=$G/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.fna          # has a BLAST database
CDS=$G/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_cds_from_genomic.fna
FMT="6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore"
mkdir -p "$R" "$D"

# 1. Parental CDS from the IsoJpt1.0 RefSeq annotation (first CDS record per gene). The RpS14 parent is annotated as
#    LOC108018733 in IsoJpt1.0, so it is extracted separately.
awk 'BEGIN{RS=">"; ORS=""} NR>1{h=substr($0,1,index($0,"\n")-1); if(h~/\[gene=(Sce|betaTub97EF|RpL36|RpS14a|CCT1|VhaAC39-1)\]/){match(h,/\[gene=[^]]+\]/); g=substr(h,RSTART+6,RLENGTH-7); if(!(g in seen)){seen[g]=1; print ">"g"\n"substr($0,index($0,"\n")+1)}}}' \
  "$CDS" > "$R/cds6.fa"
awk 'BEGIN{RS=">"; ORS=""} NR>1{h=substr($0,1,index($0,"\n")-1); if(h~/\[gene=LOC108018733\]/ && !s){s=1; print ">RpS14\n"substr($0,index($0,"\n")+1)}}' \
  "$CDS" > "$R/parent_cds.fa"
cat "$R/cds6.fa" >> "$R/parent_cds.fa"

# 2. Copies in IsoJpt1.0
blastn -task dc-megablast -query "$R/cds6.fa" -db "$GI" -evalue 1e-10 -outfmt "$FMT" -max_target_seqs 50 -max_hsps 50 \
  -num_threads 1 > "$R/cds6_vs_IsoJpt.tsv"
awk 'BEGIN{RS=">";ORS=""} NR>1 && /^RpS14\n/{print ">"$0}' "$R/parent_cds.fa" > "$D/rps14_cds.fa"
blastn -task dc-megablast -query "$D/rps14_cds.fa" -db "$GI" -evalue 1e-10 -outfmt "$FMT" -num_threads 1 \
  > "$R/rps14_cds_vs_IsoJpt.tsv"

# 3. WT3 assemblies (NCBI FTP), BLAST databases and copies of the parental CDS
for a in GCA_000472105.1_Dsuzukii.v01 GCA_013340165.1_LBDM_Dsuz_2.1.pri; do
  acc=${a%%_[A-Z]*}; n=${acc#GCA_}; n=${n%.*}
  url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/${n:0:3}/${n:3:3}/${n:6:3}/$a/${a}_genomic.fna.gz"
  [ -s "$D/$a.fna" ] || { curl -s --max-time 900 -o "$D/${a}_genomic.fna.gz" "$url"; gzip -t "$D/${a}_genomic.fna.gz";
                          gunzip -c "$D/${a}_genomic.fna.gz" > "$D/$a.fna"; }
  [ -s "$D/$a.nsq" ] || makeblastdb -in "$D/$a.fna" -dbtype nucl -out "$D/$a" > /dev/null
  [ -s "$D/$a.fna.fai" ] || $S faidx "$D/$a.fna"
  blastn -task dc-megablast -query "$R/parent_cds.fa" -db "$D/$a" -evalue 1e-10 -outfmt "$FMT" -max_target_seqs 50 \
    -max_hsps 50 -num_threads 1 > "$R/parent_cds_vs_$a.tsv"
done

# 4. Retrocopy regions +- 300 nt. Coordinates: IsoJpt1.0 loci of scripts/137 (VhaAC39-1: its alignment, see
#    scripts/145_hidden_nested_copies.py); WT3 copies = full-length single-block hits of step 3.
add() { $S faidx "$1" "$2:$(( $3 - 300 ))-$(( $4 + 300 ))" | sed "1s/.*/>$5/"; }
{ add "$GI" NW_027255939.1 1918069 1918413 RpL36; add "$GI" NW_027255927.1 2662680 2663100 RpS14
  add "$GI" NW_027255899.1 1547553 1549223 CCT1;  add "$GI" NW_027255896.1 24702032 24703088 VhaAC39-1; } > "$R/retro_regions.fa"
V=$D/GCA_000472105.1_Dsuzukii.v01.fna; L=$D/GCA_013340165.1_LBDM_Dsuz_2.1.pri.fna
{ add "$GI" NW_027255939.1 1918069 1918413 RpL36_IsoJpt; add "$GI" NW_027255927.1 2662680 2663100 RpS14_IsoJpt
  add "$GI" NW_027255899.1 1547553 1549223 CCT1_IsoJpt;  add "$GI" NW_027255896.1 24702032 24703088 VhaAC39-1_IsoJpt
  add "$V" KI420146.1 7977 8324 RpL36_WT3v01;  add "$V" KI419779.1 121568 121991 RpS14_WT3v01
  add "$V" KI419368.1 189207 190880 CCT1_WT3v01; add "$V" KI419582.1 44463 45509 VhaAC39-1_WT3v01
  add "$L" WWNF01000046.1 557307 557654 RpL36_WT3LBDM;  add "$L" WWNF01000053.1 534970 535393 RpS14_WT3LBDM
  add "$L" WWNF01000023.1 1763947 1765620 CCT1_WT3LBDM; add "$L" WWNF01000140.1 231265 232324 VhaAC39-1_WT3LBDM
} > "$R/retro_regions_3asm.fa"
echo "done"

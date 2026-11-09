#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH -J dfam39_p1
#SBATCH -o logs/dfam39_p1.%j.out
#SBATCH -e logs/dfam39_p1.%j.err

set -euo pipefail

OUTDIR=/ru-auth/local/home/asato/scratch/synteny_4species/reanalysis_20260824/databases/Dfam_3.9

URL=https://www.dfam.org/releases/Dfam_3.9/families/FamDB/dfam39_full.1.h5.gz
MD5_URL=https://www.dfam.org/releases/Dfam_3.9/families/FamDB/dfam39_full.1.h5.gz.md5

mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "===== Download Dfam 3.9 partition 1: Brachycera ====="
date
hostname

/usr/bin/wget \
  -c \
  --progress=dot:giga \
  "$URL"

/usr/bin/wget \
  -O dfam39_full.1.h5.gz.md5 \
  "$MD5_URL"

echo "===== Validate checksum ====="

md5sum -c dfam39_full.1.h5.gz.md5

ls -lh \
  dfam39_full.1.h5.gz \
  dfam39_full.1.h5.gz.md5

echo "FINAL: PASS"
date

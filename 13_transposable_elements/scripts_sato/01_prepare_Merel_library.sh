#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --job-name=TE_MerelPrep
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only/logs/TE_MerelPrep.%j.out
#SBATCH --error=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only/logs/TE_MerelPrep.%j.err

set -euo pipefail

OLD=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20260824
NEW=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only

SOURCE_FA="${OLD}/resources/DSu-TE/TE_db/all_TEs.fa"
SOURCE_TXT="${OLD}/resources/DSu-TE/TE_db/all_TEs.txt"

FULL="${NEW}/shared_library/Merel.all_TEs.original.fa"
ANNOTATION="${NEW}/shared_library/Merel.all_TEs.original.txt"
MIN500="${NEW}/shared_library/Merel.min500.fa"
SUMMARY="${NEW}/checks/01_Merel_library.summary.tsv"

for INPUT in "$SOURCE_FA" "$SOURCE_TXT"; do
    [[ -s "$INPUT" ]] || {
        echo "ERROR: missing or empty input: $INPUT" >&2
        exit 1
    }
done

for OUTPUT in "$FULL" "$ANNOTATION" "$MIN500" "$SUMMARY"; do
    [[ ! -e "$OUTPUT" ]] || {
        echo "ERROR: output already exists; will not overwrite: $OUTPUT" >&2
        exit 1
    }
done

mkdir -p "${NEW}/shared_library" "${NEW}/checks"

TMP_FULL="${FULL}.tmp.${SLURM_JOB_ID}"
TMP_TXT="${ANNOTATION}.tmp.${SLURM_JOB_ID}"
TMP_MIN="${MIN500}.tmp.${SLURM_JOB_ID}"
TMP_COUNTS="${SUMMARY}.counts.tmp.${SLURM_JOB_ID}"
TMP_SUMMARY="${SUMMARY}.tmp.${SLURM_JOB_ID}"

cleanup() {
    rm -f -- "$TMP_FULL" "$TMP_TXT" "$TMP_MIN" "$TMP_COUNTS" "$TMP_SUMMARY"
}
trap cleanup EXIT

cp -- "$SOURCE_FA" "$TMP_FULL"
cp -- "$SOURCE_TXT" "$TMP_TXT"
cmp -s "$SOURCE_FA" "$TMP_FULL"
cmp -s "$SOURCE_TXT" "$TMP_TXT"

# 元のヘッダーと配列行を保持したまま、長さ500 bp以上だけ抽出。
awk -v counts="$TMP_COUNTS" '
function finish_record() {
    if (!have_record) return
    total++
    if (length_bp == 0) empty++
    if (length_bp >= 500) {
        printf "%s", record
        kept++
    } else {
        excluded++
    }
}
/^>/ {
    finish_record()
    have_record = 1
    record = $0 ORS
    length_bp = 0
    next
}
{
    if (!have_record) {
        malformed++
        next
    }
    sequence_line = $0
    gsub(/[[:space:]]/, "", sequence_line)
    length_bp += length(sequence_line)
    record = record $0 ORS
}
END {
    finish_record()
    printf "%d\t%d\t%d\t%d\t%d\n",
        total, kept, excluded, empty, malformed > counts
    if (total == 0 || empty > 0 || malformed > 0) exit 1
}
' "$TMP_FULL" > "$TMP_MIN"

read -r TOTAL KEPT EXCLUDED EMPTY MALFORMED < "$TMP_COUNTS"

# 前回確認したMérel原本の件数と照合。異なれば停止する。
[[ "$TOTAL" == 3975 && "$KEPT" == 2960 && "$EXCLUDED" == 1015 ]] || {
    echo "ERROR: unexpected sequence counts: total=$TOTAL kept=$KEPT excluded=$EXCLUDED" >&2
    exit 1
}

[[ -s "$TMP_MIN" ]] || {
    echo "ERROR: empty min500 FASTA" >&2
    exit 1
}

{
    printf 'metric\tvalue\n'
    printf 'source_fasta\t%s\n' "$SOURCE_FA"
    printf 'source_annotation\t%s\n' "$SOURCE_TXT"
    printf 'original_sequences\t%s\n' "$TOTAL"
    printf 'min500_sequences\t%s\n' "$KEPT"
    printf 'shorter_than500_sequences\t%s\n' "$EXCLUDED"
    printf 'empty_sequences\t%s\n' "$EMPTY"
    printf 'malformed_lines\t%s\n' "$MALFORMED"
    printf 'original_fasta_sha256\t%s\n' "$(sha256sum "$TMP_FULL" | awk '{print $1}')"
    printf 'min500_fasta_sha256\t%s\n' "$(sha256sum "$TMP_MIN" | awk '{print $1}')"
    printf 'status\tPASS\n'
} > "$TMP_SUMMARY"

mv -- "$TMP_FULL" "$FULL"
mv -- "$TMP_TXT" "$ANNOTATION"
mv -- "$TMP_MIN" "$MIN500"
mv -- "$TMP_SUMMARY" "$SUMMARY"

cat "$SUMMARY"
echo "FINAL: PASS"

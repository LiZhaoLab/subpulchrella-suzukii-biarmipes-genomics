#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --job-name=TE_RM_Merel
#SBATCH --array=0-1%2
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --output=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only/logs/TE_RM_Merel.%A_%a.out
#SBATCH --error=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only/logs/TE_RM_Merel.%A_%a.err

set -euo pipefail

NEW=/ru-auth/local/home/asato/scratch/TE_reanalysis/reanalysis_20261022_Merel_only
RM=/ru-auth/local/home/asato/scratch/conda_envs/earlgrey_env/bin/RepeatMasker
MASK_LIBRARY="${NEW}/shared_library/Merel.all_TEs.original.fa"
CALL_LIBRARY="${NEW}/shared_library/Merel.min500.fa"

SPECIES=(Dsub Dsuz)
GENOMES=(
    /ru-auth/local/home/asato/scratch/genome/Dsub_v1.1/GCF_014743375.2_RU_Dsub_v1.1_genomic.fna
    /ru-auth/local/home/asato/scratch/genome/Dsuz_IsoJpt1.0/GCF_043229965.1_CBGP_Dsuzu_IsoJpt1.0_genomic.fna
)

INDEX=${SLURM_ARRAY_TASK_ID}
[[ "$INDEX" =~ ^[01]$ ]] || {
    echo "ERROR: unexpected array task ID: $INDEX" >&2
    exit 1
}

SPECIES_NAME=${SPECIES[$INDEX]}
GENOME=${GENOMES[$INDEX]}
GENOME_NAME=$(basename "$GENOME")

OUTDIR="${NEW}/shared_library_annotation_fullmask/${SPECIES_NAME}"
CHECKDIR="${NEW}/checks/repeatmasker_fullmask"
SUMMARY="${CHECKDIR}/${SPECIES_NAME}.summary.tsv"
MASKED="${OUTDIR}/${GENOME_NAME}.masked"

for FILE in "$RM" "$MASK_LIBRARY" "$CALL_LIBRARY" "$GENOME"; do
    [[ -s "$FILE" ]] || {
        echo "ERROR: missing or empty input: $FILE" >&2
        exit 1
    }
done

MASK_RECORDS=$(grep -c '^>' "$MASK_LIBRARY")
CALL_RECORDS=$(grep -c '^>' "$CALL_LIBRARY")
[[ "$MASK_RECORDS" -eq 3975 && "$CALL_RECORDS" -eq 2960 ]] || {
    echo "ERROR: unexpected TE counts: mask=$MASK_RECORDS call=$CALL_RECORDS" >&2
    exit 1
}

MASK_SHA256=$(sha256sum "$MASK_LIBRARY" | awk '{print $1}')
CALL_SHA256=$(sha256sum "$CALL_LIBRARY" | awk '{print $1}')
[[ "$MASK_SHA256" == e6257f1d260c047525733e091adfdab6242cf580a7d01ebc69cfe6fadc606685 ]] || {
    echo "ERROR: full library checksum changed" >&2
    exit 1
}
[[ "$CALL_SHA256" == 16ed7549009f80fc5a0197e98af1633ec397a29b64d7ba8fd1357aee48e421a3 ]] || {
    echo "ERROR: min500 library checksum changed" >&2
    exit 1
}

[[ ! -e "$OUTDIR" && ! -e "$SUMMARY" ]] || {
    echo "ERROR: output already exists; will not overwrite" >&2
    exit 1
}
mkdir -p "$OUTDIR" "$CHECKDIR"

# RMBlastでは -pa 8 が概ね32 CPUに対応する。
RM_PARALLEL=$((SLURM_CPUS_PER_TASK / 4))

echo "species=$SPECIES_NAME"
echo "genome=$GENOME"
echo "mask_library=$MASK_LIBRARY"
echo "mask_library_records=$MASK_RECORDS"
echo "later_mapping_TE_library=$CALL_LIBRARY"
echo "later_mapping_TE_records=$CALL_RECORDS"
echo "output_directory=$OUTDIR"
echo "options=-e rmblast -no_is -nolow -norna -s -gff"
echo "rm_parallel=$RM_PARALLEL"
"$RM" -v 2>&1 | head -n 1
date

"$RM" \
    -e rmblast \
    -lib "$MASK_LIBRARY" \
    -no_is \
    -nolow \
    -norna \
    -s \
    -gff \
    -pa "$RM_PARALLEL" \
    -dir "$OUTDIR" \
    "$GENOME"

for SUFFIX in masked out out.gff; do
    FILE="${OUTDIR}/${GENOME_NAME}.${SUFFIX}"
    [[ -s "$FILE" ]] || {
        echo "ERROR: expected output missing or empty: $FILE" >&2
        exit 1
    }
done

fasta_lengths() {
    awk '
        /^>/ {
            if (id != "") print id "\t" length_bp
            id = substr($1, 2)
            length_bp = 0
            next
        }
        {
            gsub(/[[:space:]]/, "")
            length_bp += length($0)
        }
        END {
            if (id != "") print id "\t" length_bp
        }
    ' "$1"
}

if ! diff -q <(fasta_lengths "$GENOME") <(fasta_lengths "$MASKED"); then
    echo "ERROR: genome sequence names, lengths, or order changed" >&2
    exit 1
fi

{
    printf 'metric\tvalue\n'
    printf 'species\t%s\n' "$SPECIES_NAME"
    printf 'input_genome\t%s\n' "$GENOME"
    printf 'mask_library\t%s\n' "$MASK_LIBRARY"
    printf 'mask_library_sequences\t%s\n' "$MASK_RECORDS"
    printf 'mask_library_sha256\t%s\n' "$MASK_SHA256"
    printf 'later_mapping_TE_library\t%s\n' "$CALL_LIBRARY"
    printf 'later_mapping_TE_sequences\t%s\n' "$CALL_RECORDS"
    printf 'masked_genome\t%s\n' "$MASKED"
    printf 'repeatmasker_options\t-e rmblast -no_is -nolow -norna -s -gff\n'
    printf 'genome_names_lengths_order\tMATCH\n'
    printf 'status\tPASS\n'
} > "$SUMMARY"

cat "$SUMMARY"
echo "FINAL: PASS"
date

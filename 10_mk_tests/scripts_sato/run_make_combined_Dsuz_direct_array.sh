#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --array=0-11%4
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=0-02:00:00
#SBATCH -J combine_Dsuz
#SBATCH -o /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/combine_Dsuz.%A_%a.out
#SBATCH -e /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/combine_Dsuz.%A_%a.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817

POPS=("JP-Kan" "JP-Sap" "JP-Sap2" "JP-Tok")
COMPS=("Dsub" "Dbia" "Dmel")

POP_INDEX=$((SLURM_ARRAY_TASK_ID / 3))
COMP_INDEX=$((SLURM_ARRAY_TASK_ID % 3))

POP=${POPS[$POP_INDEX]}
COMP=${COMPS[$COMP_INDEX]}

INPUT_MSA=${BASE}/vcf2msa_longest_RBH/Dsuz/${POP}/Dsuz_vs_${COMP}
OUTPUT_DIR=${BASE}/combined_fasta/Dsuz/${POP}/Dsuz_vs_${COMP}
CHECK_DIR=${BASE}/checks/Dsuz_direct_combined

PROTEIN_RBH=${BASE}/orth_longest/Dsuz_${COMP}.direct.RBH.longest.protein.tsv
MRNA_RBH=${BASE}/orth_longest/Dsuz_${COMP}.direct.RBH.longest.mRNA.tsv

case "${COMP}" in
    Dsub)
        CDS_FASTA=/ru-auth/local/home/asato/scratch/genome/Dsub_v1.1/GCF_014743375.2_RU_Dsub_v1.1_cds_from_genomic.fna
        ;;
    Dbia)
        CDS_FASTA=/ru-auth/local/home/asato/scratch/genome/Dbia/GCF_025231255.1_RU_DBia_V1.1_cds_from_genomic.fna
        ;;
    Dmel)
        CDS_FASTA=/ru-auth/local/home/asato/scratch/genome/Dmel/GCF_000001215.4_Release_6_plus_ISO1_MT_cds_from_genomic.fna
        ;;
    *)
        echo "ERROR: unsupported comparison: ${COMP}" >&2
        exit 1
        ;;
esac

mkdir -p \
    "${OUTPUT_DIR}" \
    "${CHECK_DIR}" \
    "${BASE}/logs"

export POP
export COMP
export INPUT_MSA
export OUTPUT_DIR
export CHECK_DIR
export PROTEIN_RBH
export MRNA_RBH
export CDS_FASTA

python3 - <<'PY'
from pathlib import Path
import os
import re

pop = os.environ["POP"]
comp = os.environ["COMP"]

input_dir = Path(os.environ["INPUT_MSA"])
output_dir = Path(os.environ["OUTPUT_DIR"])
check_dir = Path(os.environ["CHECK_DIR"])

protein_rbh = Path(os.environ["PROTEIN_RBH"])
mrna_rbh = Path(os.environ["MRNA_RBH"])
cds_fasta = Path(os.environ["CDS_FASTA"])

required_paths = (
    input_dir,
    protein_rbh,
    mrna_rbh,
    cds_fasta,
)

for path in required_paths:
    if not path.exists():
        raise FileNotFoundError(path)

output_dir.mkdir(parents=True, exist_ok=True)
check_dir.mkdir(parents=True, exist_ok=True)

# この新規出力ディレクトリ内にある既存のcombined FASTAのみ削除
for old_file in output_dir.glob("*.combined.fasta"):
    old_file.unlink()


def read_pairs(path):
    pairs = []

    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.rstrip("\n")

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) < 2:
                raise RuntimeError(
                    f"Malformed row in {path}, line {line_number}"
                )

            pairs.append((fields[0], fields[1]))

    return pairs


protein_pairs = read_pairs(protein_rbh)
mrna_pairs = read_pairs(mrna_rbh)

if len(protein_pairs) != len(mrna_pairs):
    raise RuntimeError(
        "Protein and mRNA RBH row counts differ: "
        f"{len(protein_pairs)} versus {len(mrna_pairs)}"
    )

# 両RBH表の列順：
# 1列目 = outgroup
# 2列目 = focal Dsuz
dsuz_mrna_to_outgroup_protein = {}
rbh_mapping_conflicts = []

for protein_pair, mrna_pair in zip(protein_pairs, mrna_pairs):
    outgroup_protein, dsuz_protein = protein_pair
    outgroup_mrna, dsuz_mrna = mrna_pair

    previous = dsuz_mrna_to_outgroup_protein.get(dsuz_mrna)

    if previous is not None and previous != outgroup_protein:
        rbh_mapping_conflicts.append(
            (
                dsuz_mrna,
                previous,
                outgroup_protein,
            )
        )
        continue

    dsuz_mrna_to_outgroup_protein[dsuz_mrna] = (
        outgroup_protein
    )

if rbh_mapping_conflicts:
    raise RuntimeError(
        "Conflicting RBH mappings were found: "
        f"{len(rbh_mapping_conflicts)}"
    )


def read_fasta(path):
    header = None
    sequence = []

    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)

                header = line[1:]
                sequence = []
            else:
                sequence.append(line.strip())

        if header is not None:
            yield header, "".join(sequence)


def protein_id_from_header(header):
    if comp in {"Dsub", "Dbia"}:
        match = re.search(
            r"\[protein_id=([^\]]+)\]",
            header,
        )
    elif comp == "Dmel":
        match = re.search(
            r"FLYBASE:(FBpp[0-9]+)",
            header,
        )
    else:
        match = None

    if match:
        return match.group(1)

    return None


needed_proteins = set(
    dsuz_mrna_to_outgroup_protein.values()
)

protein_to_cds = {}
conflicting_cds = set()

for header, sequence in read_fasta(cds_fasta):
    protein_id = protein_id_from_header(header)

    if protein_id is None:
        continue

    if protein_id not in needed_proteins:
        continue

    sequence = sequence.upper()

    if protein_id in protein_to_cds:
        if protein_to_cds[protein_id] != sequence:
            conflicting_cds.add(protein_id)
    else:
        protein_to_cds[protein_id] = sequence

if conflicting_cds:
    raise RuntimeError(
        "Protein IDs associated with different CDS sequences: "
        + ",".join(sorted(conflicting_cds)[:20])
    )

manifest = []
missing_cds = []
invalid = []
not_in_rbh = []

input_files = sorted(input_dir.glob("*.fasta"))

for msa_file in input_files:
    dsuz_mrna = msa_file.stem

    outgroup_protein = (
        dsuz_mrna_to_outgroup_protein.get(dsuz_mrna)
    )

    if outgroup_protein is None:
        not_in_rbh.append(
            (
                dsuz_mrna,
                str(msa_file),
            )
        )
        continue

    outgroup_cds = protein_to_cds.get(outgroup_protein)

    if outgroup_cds is None:
        missing_cds.append(
            (
                dsuz_mrna,
                outgroup_protein,
                "protein_not_found_in_complete_CDS_FASTA",
            )
        )
        continue

    if len(outgroup_cds) == 0:
        invalid.append(
            (
                dsuz_mrna,
                outgroup_protein,
                str(len(outgroup_cds)),
                "empty_outgroup_CDS",
            )
        )
        continue

    if len(outgroup_cds) % 3 != 0:
        invalid.append(
            (
                dsuz_mrna,
                outgroup_protein,
                str(len(outgroup_cds)),
                "outgroup_CDS_length_not_multiple_of_three",
            )
        )
        continue

    focal_records = list(read_fasta(msa_file))

    if not focal_records:
        invalid.append(
            (
                dsuz_mrna,
                outgroup_protein,
                str(len(outgroup_cds)),
                "empty_focal_MSA",
            )
        )
        continue

    focal_lengths = {
        len(sequence)
        for header, sequence in focal_records
    }

    if len(focal_lengths) != 1:
        invalid.append(
            (
                dsuz_mrna,
                outgroup_protein,
                str(len(outgroup_cds)),
                "unequal_focal_MSA_lengths",
            )
        )
        continue

    focal_length = next(iter(focal_lengths))

    if focal_length == 0:
        invalid.append(
            (
                dsuz_mrna,
                outgroup_protein,
                str(len(outgroup_cds)),
                "empty_focal_MSA",
            )
        )
        continue

    if focal_length % 3 != 0:
        invalid.append(
            (
                dsuz_mrna,
                outgroup_protein,
                str(len(outgroup_cds)),
                "focal_MSA_length_not_multiple_of_three",
            )
        )
        continue

    destination = (
        output_dir
        / f"{dsuz_mrna}.combined.fasta"
    )

    with destination.open("w") as output:
        for header, sequence in focal_records:
            output.write(f">{header}\n")

            for i in range(0, len(sequence), 80):
                output.write(
                    sequence[i:i + 80] + "\n"
                )

        output.write(
            f">{comp.lower()}|{outgroup_protein}\n"
        )

        for i in range(0, len(outgroup_cds), 80):
            output.write(
                outgroup_cds[i:i + 80] + "\n"
            )

    manifest.append(
        (
            dsuz_mrna,
            outgroup_protein,
            str(len(focal_records)),
            str(focal_length),
            str(len(outgroup_cds)),
            str(msa_file),
            str(destination),
        )
    )

prefix = f"{pop}.Dsuz_vs_{comp}"

manifest_file = (
    check_dir
    / f"{prefix}.combined_manifest.tsv"
)

missing_file = (
    check_dir
    / f"{prefix}.missing_CDS.tsv"
)

invalid_file = (
    check_dir
    / f"{prefix}.invalid_combined.tsv"
)

not_rbh_file = (
    check_dir
    / f"{prefix}.MSA_not_in_RBH.tsv"
)

summary_file = (
    check_dir
    / f"{prefix}.summary.tsv"
)

with manifest_file.open("w") as output:
    output.write(
        "Dsuz_mRNA\toutgroup_protein\t"
        "focal_sequence_count\t"
        "focal_alignment_length\t"
        "outgroup_CDS_length\t"
        "input_MSA\tcombined_FASTA\n"
    )

    for row in manifest:
        output.write("\t".join(row) + "\n")

with missing_file.open("w") as output:
    output.write(
        "Dsuz_mRNA\toutgroup_protein\treason\n"
    )

    for row in missing_cds:
        output.write("\t".join(row) + "\n")

with invalid_file.open("w") as output:
    output.write(
        "Dsuz_mRNA\toutgroup_protein\t"
        "outgroup_CDS_length\treason\n"
    )

    for row in invalid:
        output.write("\t".join(row) + "\n")

with not_rbh_file.open("w") as output:
    output.write(
        "Dsuz_mRNA\tinput_MSA\n"
    )

    for row in not_in_rbh:
        output.write("\t".join(row) + "\n")

# 3の倍数でない配列はPRANK -codonに入力できないため、
# 追跡可能なexpected exclusionとして記録する。
expected_invalid_reasons = {
    "outgroup_CDS_length_not_multiple_of_three",
    "focal_MSA_length_not_multiple_of_three",
}

unexpected_invalid = [
    row
    for row in invalid
    if row[3] not in expected_invalid_reasons
]

expected_codon_exclusions = (
    len(invalid) - len(unexpected_invalid)
)

fatal_problems = (
    len(missing_cds)
    + len(not_in_rbh)
    + len(unexpected_invalid)
)

orco = "rna-NM_001328601.1"

orco_combined = (
    output_dir
    / f"{orco}.combined.fasta"
).exists()

summary = [
    (
        "population",
        pop,
    ),
    (
        "comparison",
        f"Dsuz_vs_{comp}",
    ),
    (
        "input_direct_RBH_pairs",
        len(protein_pairs),
    ),
    (
        "available_focal_MSA_files",
        len(input_files),
    ),
    (
        "selected_outgroup_CDS",
        len(protein_to_cds),
    ),
    (
        "combined_FASTA_written",
        len(manifest),
    ),
    (
        "missing_outgroup_CDS",
        len(missing_cds),
    ),
    (
        "expected_codon_incompatible_exclusions",
        expected_codon_exclusions,
    ),
    (
        "unexpected_invalid_combined",
        len(unexpected_invalid),
    ),
    (
        "MSA_not_in_RBH",
        len(not_in_rbh),
    ),
    (
        "fatal_problems",
        fatal_problems,
    ),
    (
        "Orco_combined",
        "present" if orco_combined else "not_present",
    ),
    (
        "output_directory",
        str(output_dir),
    ),
    (
        "manifest",
        str(manifest_file),
    ),
    (
        "expected_exclusion_table",
        str(invalid_file),
    ),
]

with summary_file.open("w") as output:
    output.write("metric\tvalue\n")

    for key, value in summary:
        output.write(f"{key}\t{value}\n")

print(
    f"===== {pop}: Dsuz versus {comp} ====="
)

for key, value in summary:
    print(f"{key}\t{value}")

if fatal_problems != 0:
    raise RuntimeError(
        f"Found {fatal_problems} fatal problems; "
        f"inspect {check_dir}"
    )

if (
    len(manifest)
    + expected_codon_exclusions
    != len(input_files)
):
    raise RuntimeError(
        "Combined FASTA plus expected exclusions "
        "does not equal focal MSA count"
    )

print("Done.")
PY

#!/bin/bash
#SBATCH --partition=zhao_a
#SBATCH --account=zhao_condo_bank
#SBATCH --array=0-11
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0-04:00:00
#SBATCH -J ext_Dsuz_direct
#SBATCH -o /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/extract_Dsuz_direct.%A_%a.out
#SBATCH -e /ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817/logs/extract_Dsuz_direct.%A_%a.err

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/MKtest/reanalysis_longest_fix_20260817
DSUZBASE=/ru-auth/local/home/asato/scratch/MKtest/Dsuz4pop_MKtest
MKBASE=/ru-auth/local/home/asato/scratch/MKtest

POPS=(
  "JP-Kan" "JP-Kan" "JP-Kan"
  "JP-Sap" "JP-Sap" "JP-Sap"
  "JP-Sap2" "JP-Sap2" "JP-Sap2"
  "JP-Tok" "JP-Tok" "JP-Tok"
)

COMPS=(
  "Dsub" "Dbia" "Dmel"
  "Dsub" "Dbia" "Dmel"
  "Dsub" "Dbia" "Dmel"
  "Dsub" "Dbia" "Dmel"
)

POP=${POPS[$SLURM_ARRAY_TASK_ID]}
COMP=${COMPS[$SLURM_ARRAY_TASK_ID]}

RAW_MSA=${DSUZBASE}/vcf2msa_2026/${POP}_msa
INDEXCSV=${RAW_MSA}/cds_snp_index.csv

CDS_MAPPING=${DSUZBASE}/cds/Dsuz_IsoJpt1.0.cds_header_mapping.tsv
DSUZ_MRNA_MAP=${MKBASE}/orth_longest/Dsuz.IsoJpt1.0.protein_to_mrna.tsv

PROTEIN_RBH=${BASE}/orth_longest/Dsuz_${COMP}.direct.RBH.longest.protein.tsv
MRNA_RBH=${BASE}/orth_longest/Dsuz_${COMP}.direct.RBH.longest.mRNA.tsv

OUTDIR=${BASE}/vcf2msa_longest_RBH/Dsuz/${POP}/Dsuz_vs_${COMP}
CHECKDIR=${BASE}/checks/Dsuz_direct_VCF2MSA

mkdir -p \
  "${OUTDIR}" \
  "${CHECKDIR}" \
  "${BASE}/logs"

echo "===== Dsuz direct-RBH VCF2MSA extraction ====="
echo "Population: ${POP}"
echo "Comparison: Dsuz versus ${COMP}"
echo "Raw MSA: ${RAW_MSA}"
echo "Protein RBH: ${PROTEIN_RBH}"
echo "mRNA RBH: ${MRNA_RBH}"
echo "Output: ${OUTDIR}"

python3 - <<PY
from pathlib import Path
from collections import defaultdict, Counter
import csv
import shutil

pop = "${POP}"
comp = "${COMP}"

raw_msa = Path("${RAW_MSA}")
index_file = Path("${INDEXCSV}")
mapping_file = Path("${CDS_MAPPING}")
dsuz_mrna_map_file = Path("${DSUZ_MRNA_MAP}")

protein_rbh_file = Path("${PROTEIN_RBH}")
mrna_rbh_file = Path("${MRNA_RBH}")

outdir = Path("${OUTDIR}")
checkdir = Path("${CHECKDIR}")

required = [
    raw_msa,
    index_file,
    mapping_file,
    dsuz_mrna_map_file,
    protein_rbh_file,
    mrna_rbh_file,
]

for path in required:
    if not path.exists():
        raise RuntimeError(
            f"Required input does not exist: {path}"
        )

    if path.is_file() and path.stat().st_size == 0:
        raise RuntimeError(
            f"Required input is empty: {path}"
        )

outdir.mkdir(parents=True, exist_ok=True)
checkdir.mkdir(parents=True, exist_ok=True)

# Remove only previous FASTA outputs from this exact
# population/comparison directory.
for old_fasta in outdir.glob("*.fasta"):
    old_fasta.unlink()


def normalize_mrna(value):
    value = value.strip()

    if value.startswith("rna-"):
        return value

    if value.startswith(("NM_", "XM_", "NR_", "XR_")):
        return "rna-" + value

    return value


# ----------------------------------------------------------
# Read Dsuz protein-to-mRNA map
# ----------------------------------------------------------

protein_to_mrna = {}
map_conflicts = []

with dsuz_mrna_map_file.open() as fh:
    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in reader:
        protein = row["protein_id"].strip()
        mrna = normalize_mrna(
            row["mrna_id"]
        )

        previous = protein_to_mrna.get(protein)

        if previous is not None and previous != mrna:
            map_conflicts.append(
                (protein, previous, mrna)
            )
            continue

        protein_to_mrna[protein] = mrna

if map_conflicts:
    raise RuntimeError(
        "Conflicting Dsuz protein-to-mRNA mappings: "
        f"{len(map_conflicts)}"
    )


# ----------------------------------------------------------
# Read direct protein RBH
#
# Columns:
#   outgroup protein
#   focal Dsuz protein
# ----------------------------------------------------------

protein_pairs = []

with protein_rbh_file.open() as fh:
    for line in fh:
        fields = line.rstrip("\n").split("\t")

        if len(fields) < 2:
            continue

        outgroup_protein = fields[0].strip()
        dsuz_protein = fields[1].strip()

        protein_pairs.append(
            (outgroup_protein, dsuz_protein)
        )


# ----------------------------------------------------------
# Read direct mRNA RBH
#
# Columns:
#   outgroup mRNA
#   focal Dsuz mRNA
# ----------------------------------------------------------

mrna_pair_by_dsuz = {}

with mrna_rbh_file.open() as fh:
    for line in fh:
        fields = line.rstrip("\n").split("\t")

        if len(fields) < 2:
            continue

        outgroup_mrna = normalize_mrna(fields[0])
        dsuz_mrna = normalize_mrna(fields[1])

        if dsuz_mrna in mrna_pair_by_dsuz:
            raise RuntimeError(
                f"Duplicated Dsuz mRNA in RBH: {dsuz_mrna}"
            )

        mrna_pair_by_dsuz[dsuz_mrna] = outgroup_mrna


# ----------------------------------------------------------
# Read CDS mapping, keyed by exact Dsuz protein ID
# ----------------------------------------------------------

mapping_by_protein = defaultdict(list)

with mapping_file.open() as fh:
    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in reader:
        protein = row["protein_id"].strip()

        if protein:
            mapping_by_protein[protein].append(row)


# ----------------------------------------------------------
# Read all VCF2MSA indices for each coordinate
# ----------------------------------------------------------

indices_by_header = defaultdict(list)

with index_file.open() as fh:
    reader = csv.DictReader(fh)

    for row in reader:
        coord = row["header"]
        index = row["index"]
        snp = row["SNP"]

        indices_by_header[coord].append(
            (index, snp)
        )

for coord in indices_by_header:
    indices_by_header[coord].sort(
        key=lambda x: int(x[0])
    )


# ----------------------------------------------------------
# Extract exactly one population MSA for every eligible
# direct-RBH Dsuz protein
# ----------------------------------------------------------

manifest = []
missing = []
fatal = []
seen_dsuz_mrna = set()

for outgroup_protein, dsuz_protein in protein_pairs:
    dsuz_mrna = protein_to_mrna.get(dsuz_protein)

    if dsuz_mrna is None:
        fatal.append(
            (
                outgroup_protein,
                dsuz_protein,
                "",
                "Dsuz_protein_not_in_mRNA_map",
                "",
                "",
            )
        )
        continue

    expected_outgroup_mrna = mrna_pair_by_dsuz.get(
        dsuz_mrna
    )

    if expected_outgroup_mrna is None:
        fatal.append(
            (
                outgroup_protein,
                dsuz_protein,
                dsuz_mrna,
                "Dsuz_mRNA_not_in_mRNA_RBH",
                "",
                "",
            )
        )
        continue

    mapping_rows = mapping_by_protein.get(
        dsuz_protein,
        []
    )

    if not mapping_rows:
        fatal.append(
            (
                outgroup_protein,
                dsuz_protein,
                dsuz_mrna,
                "Dsuz_protein_not_in_CDS_mapping",
                "",
                "",
            )
        )
        continue

    candidates = []

    for mapping_row in mapping_rows:
        coord = mapping_row["coord_header"]

        for index, snp in indices_by_header.get(
            coord,
            []
        ):
            source = raw_msa / f"msa_{index}.fasta"

            candidates.append(
                {
                    "coord": coord,
                    "index": index,
                    "snp": snp,
                    "source": source,
                }
            )

    retained = [
        row
        for row in candidates
        if row["snp"] == "1"
        and row["source"].is_file()
        and row["source"].stat().st_size > 0
    ]

    snp1_missing_file = [
        row
        for row in candidates
        if row["snp"] == "1"
        and not row["source"].is_file()
    ]

    if not retained:
        if snp1_missing_file:
            fatal.append(
                (
                    outgroup_protein,
                    dsuz_protein,
                    dsuz_mrna,
                    "SNP1_index_but_MSA_file_missing",
                    ";".join(
                        row["index"]
                        for row in snp1_missing_file
                    ),
                    ";".join(
                        row["coord"]
                        for row in snp1_missing_file
                    ),
                )
            )
        else:
            missing.append(
                (
                    outgroup_protein,
                    dsuz_protein,
                    dsuz_mrna,
                    "no_retained_SNP",
                    ";".join(
                        row["index"]
                        for row in candidates
                    ),
                    ";".join(
                        sorted(set(
                            row["coord"]
                            for row in candidates
                        ))
                    ),
                )
            )

        continue

    retained.sort(
        key=lambda row: (
            int(row["index"]),
            row["coord"],
        )
    )

    selected = retained[0]

    if dsuz_mrna in seen_dsuz_mrna:
        fatal.append(
            (
                outgroup_protein,
                dsuz_protein,
                dsuz_mrna,
                "duplicated_output_Dsuz_mRNA",
                selected["index"],
                selected["coord"],
            )
        )
        continue

    seen_dsuz_mrna.add(dsuz_mrna)

    destination = outdir / f"{dsuz_mrna}.fasta"

    shutil.copy2(
        selected["source"],
        destination,
    )

    manifest.append(
        (
            outgroup_protein,
            expected_outgroup_mrna,
            dsuz_protein,
            dsuz_mrna,
            selected["coord"],
            ";".join(
                row["index"]
                for row in retained
            ),
            selected["index"],
            selected["snp"],
            str(selected["source"]),
            str(destination),
        )
    )


# ----------------------------------------------------------
# Write audit files
# ----------------------------------------------------------

prefix = (
    f"{pop}.Dsuz_vs_{comp}.direct_RBH"
)

manifest_file = (
    checkdir
    / f"{prefix}.extraction_manifest.tsv"
)

missing_file = (
    checkdir
    / f"{prefix}.no_retained_SNP.tsv"
)

fatal_file = (
    checkdir
    / f"{prefix}.fatal_problems.tsv"
)

summary_file = (
    checkdir
    / f"{prefix}.summary.tsv"
)

with manifest_file.open("w") as out:
    out.write(
        "outgroup_protein\toutgroup_mRNA\t"
        "Dsuz_protein\tDsuz_mRNA\t"
        "CDS_coordinate\teligible_indices\t"
        "selected_index\tSNP\tsource_MSA\toutput_MSA\n"
    )

    for row in manifest:
        out.write("\t".join(row) + "\n")

problem_header = (
    "outgroup_protein\tDsuz_protein\tDsuz_mRNA\t"
    "reason\tcandidate_indices\tCDS_coordinates\n"
)

with missing_file.open("w") as out:
    out.write(problem_header)

    for row in missing:
        out.write("\t".join(row) + "\n")

with fatal_file.open("w") as out:
    out.write(problem_header)

    for row in fatal:
        out.write("\t".join(row) + "\n")

output_ids = {
    path.stem
    for path in outdir.glob("*.fasta")
}

manifest_ids = {
    row[3]
    for row in manifest
}

if output_ids != manifest_ids:
    raise RuntimeError(
        "Manifest and output FASTA IDs do not match"
    )

if len(output_ids) != len(manifest):
    raise RuntimeError(
        "Duplicated IDs detected in extraction manifest"
    )

orco_dsuz = "rna-NM_001328601.1"

if comp == "Dsub":
    if orco_dsuz in output_ids:
        orco_status = "extracted"
    elif any(
        row[2] == orco_dsuz
        for row in missing
    ):
        orco_status = "no_retained_SNP"
    elif any(
        row[2] == orco_dsuz
        for row in fatal
    ):
        orco_status = "fatal_problem"
    else:
        orco_status = "not_in_direct_RBH"
else:
    orco_status = "not_applicable"

summary_rows = [
    ("population", pop),
    ("comparison", f"Dsuz_vs_{comp}"),
    ("input_direct_RBH_pairs", len(protein_pairs)),
    ("extracted_MSA_files", len(manifest)),
    ("expected_no_retained_SNP", len(missing)),
    ("fatal_problems", len(fatal)),
    ("unique_output_Dsuz_mRNA", len(output_ids)),
    ("Orco_status", orco_status),
    ("output_directory", str(outdir)),
    ("manifest", str(manifest_file)),
    ("no_retained_SNP_table", str(missing_file)),
    ("fatal_problem_table", str(fatal_file)),
]

with summary_file.open("w") as out:
    out.write("metric\tvalue\n")

    for key, value in summary_rows:
        out.write(f"{key}\t{value}\n")

print(summary_file.read_text())

if fatal:
    raise RuntimeError(
        f"{pop} Dsuz_vs_{comp}: "
        f"{len(fatal)} fatal extraction problems"
    )

print("Done.")
PY
